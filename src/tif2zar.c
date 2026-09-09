/* tif2zar - sequential TIFF CT volume -> OME-Zarr v0.4 (Zarr v2) converter
 *
 *   tif2zar tifDir out.zarr {--chunk N} {--pixel P} {--min A --max B}
 *                           {--levels L} {--codec C} {--clevel N}
 *
 *     tifDir    directory holding one numbered grayscale TIFF series
 *               (uint16 or float32, e.g. rec00000.tif ...); all slices the
 *               same size and depth
 *     out.zarr  output directory (created; must not already exist)
 *     --chunk   isotropic chunk edge, default 128 (must be even)
 *     --pixel   pixel size override [um]; default: ImageDescription, else 1.0
 *     --min/max quantisation range override for float32 input
 *     --levels  pyramid levels; default: auto, halve until the largest
 *               dimension is <= 256
 *
 * uint16 input is copied through.  float32 input is linearly quantised to
 * uint16: q = 65535 (v - min) / (max - min), clipped; min/max come from a
 * metadata-only prepass over the per-slice ImageDescription min/max fields
 * (global min of mins / max of maxs), or from --min/--max (CLI wins).
 *
 * The resolution pyramid is built during the same single streaming pass:
 * every level keeps one chunk-deep band; a full band is written out as
 * chunks and simultaneously 2x2x2-averaged into the band of the next
 * level (odd edges average over the voxels that exist).
 *
 * Output: OME-Zarr v0.4 on Zarr v2, axes t,c,z,y,x (t = c = 1), dtype <u2,
 * dimension_separator "/", blosc-compressed (zstd by default) and
 * OpenMP-parallel.  Usage, output layout and measured performance:
 * 20260909_tif2zar.md.
 * Open in Fiji: Plugins > BigDataViewer > OME ZARR, or the N5 importer
 * with a file:/// URI (plain C:\ paths trip a known n5 dialog bug).
 *
 * --- ct-rec TIFF metadata convention (fixed 2026-09, from the writers) ---
 * Everything is TAB-separated text in ImageDescription (TIFF tag 270);
 * there are no custom binary tags.  Field layouts by producing program:
 *   hp_tg / ct_rec (rec*.tif, 32-bit float), 6 fields:
 *     pixel[um] \t axis_pos \t Nproj \t angle_offset \t min \t max
 *   tif_f2i (8/16-bit, normalised from the above), 8 fields:
 *     ...same 6... \t norm_min \t norm_max
 *     (norm_min/max are the physical values mapped to 0 and 65535/255)
 *   ct_prj_f (p*.tif, 32-bit float), 2 fields:  min \t max
 * pixel size = field 1 (um) when 6 or 8 fields are present.
 * Slice min/max = fields 5-6 (6/8 fields) or 1-2 (2 fields).
 * The quantisation range in use is written to .zattrs as
 * "physical_min"/"physical_max" (for uint16 input: fields 7-8 if present)
 * so raw values can be mapped back to physical ones.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include <errno.h>
#include "tiffio.h"
#include "blosc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WINDOWS
#include "msdirent.h"
#include <direct.h>
#define MKDIR(p)	_mkdir(p)
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(p)	mkdir(p, 0755)
#endif

#define MAX_SLICES	65536
#define MAX_LEVELS	12
#define NAME_LEN	512
#define DESC_FIELDS	16

static char	prefix[NAME_LEN];	/* series name parts: prefix + %0*d + .tif */
static int	numwidth;
static int	num0;			/* first slice number */
static int	nZ;			/* slices */
static char	tifdir[NAME_LEN], zardir[NAME_LEN];

static uint32_t	NX, NY;
static int	chunk = 128;
static int	is_float = 0;		/* input: 0 = uint16, 1 = float32 */
static int	nlev = 0;		/* pyramid levels; 0 = auto */

static double	pixel_um = 1.0;		/* voxel pitch, isotropic */
static int	pixel_from = 0;		/* 0=default 1=tag 2=CLI */
static double	phys_min = 0.0, phys_max = 0.0;
static int	have_phys = 0;
static double	cli_min = 0.0, cli_max = 0.0;
static int	have_cli_range = 0;
static char	desc0[2048];		/* first slice's ImageDescription */

static long	n_clip_lo = 0, n_clip_hi = 0, n_nan = 0;
static long	n_chunkfiles = 0;
static uint16_t	g_lo = 65535, g_hi = 0;	/* written data range, for omero */

static char	codec[16] = "zstd";	/* blosc cname; "none" = uncompressed */
static int	clevel = 3;
static int	nthr = 1;		/* OpenMP threads in use */
static long long n_rawbytes = 0, n_cmpbytes = 0;

/* one pyramid level: dimensions, chunk grid, one chunk-deep slice band */
typedef struct {
	int		X, Y, Z;	/* level dimensions */
	int		nyc, nxc;	/* chunk grid in y, x */
	uint16_t	*band;		/* chunk * Y * X voxels */
	int		fill;		/* slices currently in the band */
	int		zc;		/* next band (z-chunk) index */
	int		zout;		/* slices produced so far (global) */
} Level;
static Level	lv[MAX_LEVELS];
static uint16_t	*cbufs;			/* nthr chunk buffers */
static char	*zbufs;			/* nthr compression buffers */
static size_t	zbytes;			/* size of one compression buffer */
static uint16_t	*dsbuf;			/* one downsampled slice, shared */
static float	*frows;			/* nthr float rows (float32 input) */

/*----------------------------------------------------------------------*/
static void die(const char *fmt, const char *arg)
{
	fprintf(stderr, "tif2zar: ");
	fprintf(stderr, fmt, arg);
	fputc('\n', stderr);
	exit(1);
}

/*----------------------------------------------------------------------*/
/* mkdir that tolerates "already exists" */
static void mkdir_ok(const char *path)
{
	if (MKDIR(path) != 0 && errno != EEXIST)
		die("cannot create directory %s", path);
}

/*----------------------------------------------------------------------*/
/* scan tifDir for one numbered .tif series; fill prefix/num0/nZ.        */
/* Numbering gaps are fatal and every missing file name is printed.      */
static void scan_series(void)
{
	DIR		*d;
	struct dirent	*e;
	static char	pfx[MAX_SLICES][64];	/* prefix per file (short) */
	static int	num[MAX_SLICES], wid[MAX_SLICES];
	int		n = 0, i, j, k, best, bestcnt, miss;
	char		bp[64];

	if ((d = opendir(tifdir)) == NULL) die("cannot open directory %s", tifdir);
	while ((e = readdir(d)) != NULL) {
		const char	*s = e->d_name;
		int		L = (int)strlen(s), dg;

		if (L < 5 || L >= 64) continue;
		if (strcmp(s + L - 4, ".tif") != 0 && strcmp(s + L - 4, ".TIF") != 0)
			continue;
		for (dg = L - 4; dg > 0 && s[dg-1] >= '0' && s[dg-1] <= '9'; --dg) ;
		if (dg == L - 4) continue;		/* no trailing digits */
		if (n >= MAX_SLICES) die("more than %s slices", "65536");
		memcpy(pfx[n], s, dg); pfx[n][dg] = '\0';
		wid[n] = L - 4 - dg;
		num[n] = atoi(s + dg);
		++n;
	}
	closedir(d);
	if (n == 0) die("no numbered .tif files in %s", tifdir);

	/* dominant prefix (there is normally exactly one) */
	best = 0; bestcnt = 0;
	for (i = 0; i < n; ++i) {
		for (k = 0, j = 0; j < n; ++j)
			if (strcmp(pfx[j], pfx[i]) == 0 && wid[j] == wid[i]) ++k;
		if (k > bestcnt) { bestcnt = k; best = i; }
	}
	strcpy(bp, pfx[best]);
	numwidth = wid[best];
	if (bestcnt != n)
		fprintf(stderr, "tif2zar: warning: ignoring %d file(s) not matching "
		    "'%s%%0%dd.tif'\n", n - bestcnt, bp, numwidth);

	/* min/max slice number of the dominant group */
	num0 = 0x7fffffff; k = -1;
	for (i = 0; i < n; ++i) {
		if (strcmp(pfx[i], bp) != 0 || wid[i] != numwidth) continue;
		if (num[i] < num0) num0 = num[i];
		if (num[i] > k)    k    = num[i];
	}
	nZ = k - num0 + 1;
	if (nZ != bestcnt) {			/* gaps: list every missing file */
		static char	seen[MAX_SLICES];
		memset(seen, 0, sizeof(seen));
		for (i = 0; i < n; ++i)
			if (strcmp(pfx[i], bp) == 0 && wid[i] == numwidth)
				seen[num[i] - num0] = 1;
		fprintf(stderr, "tif2zar: the series %s%%0%dd.tif has %d gap(s):\n",
		    bp, numwidth, nZ - bestcnt);
		for (miss = 0, i = 0; i < nZ; ++i)
			if (!seen[i]) {
				fprintf(stderr, "  missing: %s/%s%0*d.tif\n",
				    tifdir, bp, numwidth, num0 + i);
				if (++miss >= 50) {
					fprintf(stderr, "  ... (%d more)\n",
					    nZ - bestcnt - miss);
					break;
				}
			}
		exit(1);
	}
	snprintf(prefix, sizeof(prefix), "%s", bp);
	fprintf(stderr, "tif2zar: series %s%0*d.tif .. %s%0*d.tif (%d slices)\n",
	    prefix, numwidth, num0, prefix, numwidth, num0 + nZ - 1, nZ);
}

/*----------------------------------------------------------------------*/
static void slice_path(char *buf, size_t len, int iz)
{
	snprintf(buf, len, "%s/%s%0*d.tif", tifdir, prefix, numwidth, num0 + iz);
}

/*----------------------------------------------------------------------*/
/* split the TAB-separated ImageDescription into numeric fields */
static int desc_fields(const char *desc, double *f, int maxf)
{
	int	nf = 0;
	char	tmp[2048], *tok, *end;

	if (desc == NULL) return 0;
	snprintf(tmp, sizeof(tmp), "%s", desc);
	for (tok = strtok(tmp, "\t"); tok != NULL && nf < maxf;
	     tok = strtok(NULL, "\t")) {
		f[nf] = strtod(tok, &end);
		if (end == tok) return 0;	/* not our numeric format */
		++nf;
	}
	return nf;
}

/* per-slice min/max by field count (see header comment); 0 on success */
static int desc_minmax(const char *desc, double *mn, double *mx)
{
	double	f[DESC_FIELDS];
	int	nf = desc_fields(desc, f, DESC_FIELDS);

	if (nf >= 6)      { *mn = f[4]; *mx = f[5]; return 0; }
	else if (nf == 2) { *mn = f[0]; *mx = f[1]; return 0; }
	return -1;
}

/* first slice: pixel size and (uint16) physical range */
static void parse_desc0(const char *desc)
{
	double	f[DESC_FIELDS];
	int	nf;

	snprintf(desc0, sizeof(desc0), "%s", desc == NULL ? "" : desc);
	nf = desc_fields(desc, f, DESC_FIELDS);
	if (nf >= 6 && pixel_from == 0 && f[0] > 0.0) {
		pixel_um = f[0];
		pixel_from = 1;
	}
	if (nf >= 8 && !is_float) {		/* tif_f2i: normalisation range */
		phys_min = f[6];
		phys_max = f[7];
		have_phys = 1;
	}
}

/*----------------------------------------------------------------------*/
/* open one slice, verify geometry, return open TIFF (caller closes) */
static TIFF *open_slice(int iz, int first)
{
	char		path[NAME_LEN + 64];
	TIFF		*tif;
	uint32_t	w, h;
	uint16_t	bps = 0, spp = 1, fmt = SAMPLEFORMAT_UINT;
	char		*desc = NULL;

	slice_path(path, sizeof(path), iz);
	if ((tif = TIFFOpen(path, "r")) == NULL) die("cannot open %s", path);
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,  &w);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &fmt);

	if (spp != 1) die("not single-channel grayscale: %s", path);
	if (first) {
		NX = w; NY = h;
		if (bps == 32 && fmt == SAMPLEFORMAT_IEEEFP)      is_float = 1;
		else if (bps == 16 && fmt != SAMPLEFORMAT_IEEEFP) is_float = 0;
		else die("only uint16 or float32 grayscale input: %s", path);
		if (TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &desc) == 1)
			parse_desc0(desc);
	} else {
		if (w != NX || h != NY)
			die("slice size differs from the first slice: %s", path);
		if (( is_float && !(bps == 32 && fmt == SAMPLEFORMAT_IEEEFP)) ||
		    (!is_float && !(bps == 16 && fmt != SAMPLEFORMAT_IEEEFP)))
			die("pixel type differs from the first slice: %s", path);
	}
	return tif;
}

/*----------------------------------------------------------------------*/
/* float32 prepass: ImageDescription min/max of every slice, no pixels.  */
static void prepass_range(void)
{
	char	path[NAME_LEN + 64];
	int	iz, bad = 0;
	double	mn, mx;

	if (have_cli_range) {
		phys_min = cli_min;
		phys_max = cli_max;
		have_phys = 1;
		fprintf(stderr, "tif2zar: quantisation range %g .. %g (from --min/--max)\n",
		    phys_min, phys_max);
		return;
	}
	phys_min =  1.0e300;
	phys_max = -1.0e300;
	for (iz = 0; iz < nZ; ++iz) {
		TIFF	*tif;
		char	*desc = NULL;

		slice_path(path, sizeof(path), iz);
		if ((tif = TIFFOpen(path, "r")) == NULL) die("cannot open %s", path);
		if (TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &desc) != 1 ||
		    desc_minmax(desc, &mn, &mx) != 0) {
			if (++bad <= 10)
				fprintf(stderr, "tif2zar: no min/max in the "
				    "ImageDescription of %s\n", path);
		} else {
			if (mn < phys_min) phys_min = mn;
			if (mx > phys_max) phys_max = mx;
		}
		TIFFClose(tif);
	}
	if (bad > 0)
		die("float32 input needs min/max in every slice's description; "
		    "give --min and --max to continue%s", "");
	if (!(phys_max > phys_min))
		die("degenerate range: all slices have min == max; "
		    "give --min and --max%s", "");
	have_phys = 1;
	fprintf(stderr, "tif2zar: quantisation range %g .. %g (from %d descriptions)\n",
	    phys_min, phys_max, nZ);
}

/*----------------------------------------------------------------------*/
/* JSON string escape into out (tabs etc.) */
static void jesc(const char *s, char *out, size_t len)
{
	size_t	o = 0;

	for (; *s != '\0' && o + 8 < len; ++s) {
		unsigned char c = (unsigned char)*s;
		if (c == '"' || c == '\\') { out[o++] = '\\'; out[o++] = (char)c; }
		else if (c == '\t') { out[o++] = '\\'; out[o++] = 't'; }
		else if (c == '\n') { out[o++] = '\\'; out[o++] = 'n'; }
		else if (c == '\r') { out[o++] = '\\'; out[o++] = 'r'; }
		else if (c < 0x20) { o += (size_t)sprintf(out + o, "\\u%04x", c); }
		else out[o++] = (char)c;
	}
	out[o] = '\0';
}

/*----------------------------------------------------------------------*/
static void write_root_json(void)
{
	char	path[NAME_LEN + 32], esc[4096];
	FILE	*f;
	int	l;

	snprintf(path, sizeof(path), "%s/.zgroup", zardir);
	if ((f = fopen(path, "wb")) == NULL) die("cannot write %s", path);
	fprintf(f, "{\n    \"zarr_format\": 2\n}");
	fclose(f);

	jesc(desc0, esc, sizeof(esc));
	snprintf(path, sizeof(path), "%s/.zattrs", zardir);
	if ((f = fopen(path, "wb")) == NULL) die("cannot write %s", path);
	fprintf(f,
	    "{\n"
	    "    \"multiscales\": [{\n"
	    "        \"version\": \"0.4\",\n"
	    "        \"name\": \"ct volume\",\n"
	    "        \"axes\": [\n"
	    "            {\"name\": \"t\", \"type\": \"time\"},\n"
	    "            {\"name\": \"c\", \"type\": \"channel\"},\n"
	    "            {\"name\": \"z\", \"type\": \"space\", \"unit\": \"micrometer\"},\n"
	    "            {\"name\": \"y\", \"type\": \"space\", \"unit\": \"micrometer\"},\n"
	    "            {\"name\": \"x\", \"type\": \"space\", \"unit\": \"micrometer\"}\n"
	    "        ],\n"
	    "        \"datasets\": [");
	for (l = 0; l < nlev; ++l) {
		double	s = pixel_um * (double)(1 << l);
		fprintf(f,
		    "%s{\n"
		    "            \"path\": \"%d\",\n"
		    "            \"coordinateTransformations\": [{\n"
		    "                \"type\": \"scale\",\n"
		    "                \"scale\": [1.0, 1.0, %.9g, %.9g, %.9g]\n"
		    "            }]\n"
		    "        }", l == 0 ? "" : ", ", l, s, s, s);
	}
	fprintf(f,
	    "]\n"
	    "    }],\n"
	    "    \"omero\": {\n"
	    "        \"id\": 1,\n"
	    "        \"channels\": [{\n"
	    "            \"color\": \"FFFFFF\",\n"
	    "            \"label\": \"CT\",\n"
	    "            \"active\": true,\n"
	    "            \"window\": {\"min\": 0.0, \"max\": 65535.0, "
	                             "\"start\": %.1f, \"end\": %.1f}\n"
	    "        }]\n"
	    "    },\n"
	    "    \"ct_rec\": {\n"
	    "        \"source\": \"tif2zar\",\n"
	    "        \"input_type\": \"%s\",\n"
	    "        \"pixel_size_um\": %.9g,\n"
	    "        \"pixel_size_from\": \"%s\",\n",
	    (double)g_lo, (double)(g_hi > g_lo ? g_hi : g_lo + 1),
	    is_float ? "float32" : "uint16",
	    pixel_um, pixel_from == 2 ? "cli" : pixel_from == 1 ? "tag" : "default");
	if (have_phys)
		fprintf(f,
		    "        \"physical_min\": %.9g,\n"
		    "        \"physical_max\": %.9g,\n", phys_min, phys_max);
	fprintf(f,
	    "        \"source_description\": \"%s\"\n"
	    "    }\n"
	    "}", esc);
	fclose(f);
}

/*----------------------------------------------------------------------*/
static void write_zarray(int l)
{
	char	path[NAME_LEN + 32];
	FILE	*f;

	snprintf(path, sizeof(path), "%s/%d/.zarray", zardir, l);
	if ((f = fopen(path, "wb")) == NULL) die("cannot write %s", path);
	fprintf(f,
	    "{\n"
	    "    \"zarr_format\": 2,\n"
	    "    \"shape\": [1, 1, %d, %d, %d],\n"
	    "    \"chunks\": [1, 1, %d, %d, %d],\n"
	    "    \"dtype\": \"<u2\",\n",
	    lv[l].Z, lv[l].Y, lv[l].X, chunk, chunk, chunk);
	if (strcmp(codec, "none") == 0)
		fprintf(f, "    \"compressor\": null,\n");
	else
		fprintf(f,
		    "    \"compressor\": {\"id\": \"blosc\", \"cname\": \"%s\", "
		    "\"clevel\": %d, \"shuffle\": %d, \"blocksize\": 0},\n",
		    codec, clevel, BLOSC_SHUFFLE);
	fprintf(f,
	    "    \"fill_value\": 0,\n"
	    "    \"order\": \"C\",\n"
	    "    \"filters\": null,\n"
	    "    \"dimension_separator\": \"/\"\n"
	    "}");
	fclose(f);
}

/*----------------------------------------------------------------------*/
/* read one slice into dst (quantising float32 on the fly) */
static void read_slice(int iz, uint16_t *dst, float *frow)
{
	TIFF	*tif = open_slice(iz, 0);
	char	path[NAME_LEN + 64];
	int	y, x;
	long	lnan = 0, llo = 0, lhi = 0;

	for (y = 0; y < (int)NY; ++y) {
		if (!is_float) {
			if (TIFFReadScanline(tif, dst + (size_t)y * NX, y, 0) < 0)
				goto readerr;
		} else {
			double	sc = 65535.0 / (phys_max - phys_min);

			if (TIFFReadScanline(tif, frow, y, 0) < 0)
				goto readerr;
			for (x = 0; x < (int)NX; ++x) {
				double	v = frow[x];
				double	q;

				if (v != v) { ++lnan; q = 0.0; }   /* NaN */
				else {
					q = (v - phys_min) * sc + 0.5;
					/* count only real out-of-range pixels
					   (>= 1/2 LSB outside), not the text
					   rounding of the min/max fields */
					if (q < 0.0)       ++llo;
					if (q >= 65536.0)  ++lhi;
					if (q < 0.0)     q = 0.0;
					if (q > 65535.0) q = 65535.0;
				}
				dst[(size_t)y * NX + x] = (uint16_t)q;
			}
		}
	}
	TIFFClose(tif);
	if (lnan + llo + lhi > 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
		{ n_nan += lnan; n_clip_lo += llo; n_clip_hi += lhi; }
	}
	return;
readerr:
	slice_path(path, sizeof(path), iz);
	die("error reading %s", path);
}

/*----------------------------------------------------------------------*/
static void push_slice(int l, const uint16_t *src);

/* write the current band of level l as chunk files, then feed the       */
/* 2x2x2-averaged slices into level l+1, and reset the band              */
static void flush_band(int l)
{
	Level	*L = &lv[l];
	int	nz = L->fill, yc, s;
	char	path[NAME_LEN + 96];

	if (nz == 0) return;

	snprintf(path, sizeof(path), "%s/%d/0/0/%d", zardir, l, L->zc);
	mkdir_ok(path);
	for (yc = 0; yc < L->nyc; ++yc) {
		snprintf(path, sizeof(path), "%s/%d/0/0/%d/%d", zardir, l, L->zc, yc);
		mkdir_ok(path);
	}

	/* chunks of this band: copy, (blosc) compress and write in parallel */
	{
		int		job, njobs = L->nyc * L->nxc;
		long long	raw = 0, cmp = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(nthr) reduction(+: raw, cmp)
#endif
		{
			int		tid = 0, s2, y2, x2;
			uint16_t	lo = 65535, hi = 0;
			uint16_t	*cb;
			char		*zb, cpath[NAME_LEN + 96];
			FILE		*cf;
			size_t		cn = (size_t)chunk * chunk * chunk;

#ifdef _OPENMP
			tid = omp_get_thread_num();
#endif
			cb = cbufs + (size_t)tid * cn;
			zb = zbufs + (size_t)tid * zbytes;
#ifdef _OPENMP
#pragma omp for
#endif
			for (job = 0; job < njobs; ++job) {
				int	jyc = job / L->nxc, jxc = job % L->nxc;
				int	ny = (jyc + 1) * chunk <= L->Y ? chunk : L->Y - jyc * chunk;
				int	nx = (jxc + 1) * chunk <= L->X ? chunk : L->X - jxc * chunk;
				const void	*out;
				size_t		outn;

				memset(cb, 0, cn * sizeof(uint16_t));
				for (s2 = 0; s2 < nz; ++s2)
					for (y2 = 0; y2 < ny; ++y2)
						memcpy(cb + ((size_t)s2 * chunk + y2) * chunk,
						    L->band + ((size_t)s2 * L->Y + jyc * chunk + y2) * L->X
						            + (size_t)jxc * chunk,
						    (size_t)nx * sizeof(uint16_t));
				if (l == 0)	/* data range (frame region only) */
					for (s2 = 0; s2 < nz; ++s2)
						for (y2 = 0; y2 < ny; ++y2) {
							const uint16_t *r =
							    cb + ((size_t)s2 * chunk + y2) * chunk;
							for (x2 = 0; x2 < nx; ++x2) {
								if (r[x2] < lo) lo = r[x2];
								if (r[x2] > hi) hi = r[x2];
							}
						}
				if (strcmp(codec, "none") != 0) {
					int	cs = blosc_compress_ctx(clevel,
						    BLOSC_SHUFFLE, sizeof(uint16_t),
						    cn * sizeof(uint16_t), cb, zb, zbytes,
						    codec, 0, 1);
					if (cs <= 0)
						die("blosc compression failed (%s)", codec);
					out = zb; outn = (size_t)cs;
				} else {
					out = cb; outn = cn * sizeof(uint16_t);
				}
				raw += (long long)(cn * sizeof(uint16_t));
				cmp += (long long)outn;
				snprintf(cpath, sizeof(cpath), "%s/%d/0/0/%d/%d/%d",
				    zardir, l, L->zc, jyc, jxc);
				if ((cf = fopen(cpath, "wb")) == NULL)
					die("cannot write %s", cpath);
				if (fwrite(out, 1, outn, cf) != outn)
					die("short write to %s", cpath);
				fclose(cf);
#ifdef _OPENMP
#pragma omp atomic
#endif
				++n_chunkfiles;
			}
			if (l == 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					if (lo < g_lo) g_lo = lo;
					if (hi > g_hi) g_hi = hi;
				}
			}
		}
		n_rawbytes += raw;
		n_cmpbytes += cmp;
	}

	/* cascade: average slice pairs into the next level (odd edges use
	   the voxels that exist; the chunk depth is even, so pairs never
	   straddle two bands) */
	if (l + 1 < nlev) {
		Level	*M = &lv[l + 1];

		for (s = 0; s < nz; s += 2) {
			int	z2 = (s + 1 < nz) ? 2 : 1;	/* volume's last
				slice may be unpaired (only in the last band) */
			int	xo, yo;

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthr) private(xo)
#endif
			for (yo = 0; yo < M->Y; ++yo) {
				int	y1 = 2 * yo, yn = (y1 + 1 < L->Y) ? 2 : 1;

				for (xo = 0; xo < M->X; ++xo) {
					int	x1 = 2 * xo, xn = (x1 + 1 < L->X) ? 2 : 1;
					int	dz, dy, dx, cnt = 0;
					unsigned long	sum = 0;

					for (dz = 0; dz < z2; ++dz)
					for (dy = 0; dy < yn; ++dy)
					for (dx = 0; dx < xn; ++dx) {
						sum += L->band[((size_t)(s + dz) * L->Y
						    + y1 + dy) * L->X + x1 + dx];
						++cnt;
					}
					dsbuf[(size_t)yo * M->X + xo] =
					    (uint16_t)((double)sum / cnt + 0.5);
				}
			}
			push_slice(l + 1, dsbuf);
		}
	}
	L->fill = 0;
	L->zc += 1;
}

/*----------------------------------------------------------------------*/
static void push_slice(int l, const uint16_t *src)
{
	Level	*L = &lv[l];

	memcpy(L->band + (size_t)L->fill * L->Y * L->X, src,
	    (size_t)L->Y * L->X * sizeof(uint16_t));
	L->fill += 1;
	L->zout += 1;
	if (L->fill == chunk) flush_band(l);
}

/*----------------------------------------------------------------------*/
int main(int argc, char **argv)
{
	int	i, a, l, mx;
	char	path[NAME_LEN + 96];
	FILE	*f;

	setlocale(LC_NUMERIC, "C");	/* fixed decimal point in JSON */

	if (argc < 3) {
		fprintf(stderr,
		    "usage : %s tifDir out.zarr {--chunk N} {--pixel P} {--min A --max B} {--levels L}\n"
		    "  numbered grayscale TIFF series (uint16 / float32) -> OME-Zarr v0.4\n"
		    "  --chunk : isotropic chunk edge (default 128, must be even)\n"
		    "  --pixel : pixel size [um], overrides the ImageDescription\n"
		    "  --min/--max : float32 quantisation range, overrides the tags\n"
		    "  --levels : pyramid levels (default: halve until <= 256 px)\n"
		    "  --codec : zstd (default) / zlib / lz4 / none (blosc, shuffle on)\n"
		    "  --clevel : blosc compression level 1-9 (default 3)\n",
		    argv[0]);
		return 1;
	}
	snprintf(tifdir, sizeof(tifdir), "%s", argv[1]);
	snprintf(zardir, sizeof(zardir), "%s", argv[2]);
	for (a = 3; a < argc; ++a) {
		if      (strcmp(argv[a], "--chunk") == 0 && a + 1 < argc) {
			chunk = atoi(argv[++a]);
			if (chunk < 8 || chunk > 1024 || (chunk & 1))
				die("bad chunk size %s (8..1024, even)", argv[a]);
		}
		else if (strcmp(argv[a], "--pixel") == 0 && a + 1 < argc) {
			pixel_um = atof(argv[++a]);
			if (!(pixel_um > 0.0)) die("bad pixel size %s", argv[a]);
			pixel_from = 2;
		}
		else if (strcmp(argv[a], "--min") == 0 && a + 1 < argc) {
			cli_min = atof(argv[++a]); have_cli_range |= 1;
		}
		else if (strcmp(argv[a], "--max") == 0 && a + 1 < argc) {
			cli_max = atof(argv[++a]); have_cli_range |= 2;
		}
		else if (strcmp(argv[a], "--levels") == 0 && a + 1 < argc) {
			nlev = atoi(argv[++a]);
			if (nlev < 1 || nlev > MAX_LEVELS)
				die("bad level count %s", argv[a]);
		}
		else if (strcmp(argv[a], "--codec") == 0 && a + 1 < argc) {
			++a;
			if (strcmp(argv[a], "zstd") && strcmp(argv[a], "zlib") &&
			    strcmp(argv[a], "lz4") && strcmp(argv[a], "none"))
				die("bad codec %s (zstd / zlib / lz4 / none)", argv[a]);
			snprintf(codec, sizeof(codec), "%s", argv[a]);
		}
		else if (strcmp(argv[a], "--clevel") == 0 && a + 1 < argc) {
			clevel = atoi(argv[++a]);
			if (clevel < 1 || clevel > 9) die("bad clevel %s", argv[a]);
		}
		else die("unknown or incomplete option %s", argv[a]);
	}
	if (have_cli_range != 0 && have_cli_range != 3)
		die("give both --min and --max%s", "");
	if (have_cli_range == 3) {
		have_cli_range = 1;
		if (!(cli_max > cli_min)) die("--max must be greater than --min%s", "");
	}

	{	/* Fiji / MoBIE / napari pick their reader by the suffix */
		size_t	zl = strlen(zardir);

		if (zl < 5 ||
		    (strcmp(zardir + zl - 5, ".zarr") != 0 &&
		     strcmp(zardir + zl - 5, ".ZARR") != 0))
			fprintf(stderr,
			    "tif2zar: warning: the output name '%s' does not end in "
			    "'.zarr'; Fiji/MoBIE detect the format by that suffix and "
			    "will try to open the folder as something else\n", zardir);
	}

	scan_series();

	{	/* first slice: geometry, pixel type, metadata */
		TIFF *tif = open_slice(0, 1);
		TIFFClose(tif);
	}
	if (pixel_from == 0)
		fprintf(stderr, "tif2zar: warning: no pixel size in the "
		    "ImageDescription; using 1.0 um (override with --pixel)\n");
	if (is_float)
		prepass_range();	/* metadata-only pass; sets phys_min/max */

	/* pyramid geometry */
	lv[0].X = (int)NX; lv[0].Y = (int)NY; lv[0].Z = nZ;
	if (nlev == 0) {		/* auto: halve until <= 256 px */
		nlev = 1;
		for (;;) {
			Level *L = &lv[nlev - 1];
			mx = L->X > L->Y ? L->X : L->Y;
			if (L->Z > mx) mx = L->Z;
			if (mx <= 256 || nlev == MAX_LEVELS) break;
			lv[nlev].X = (lv[nlev-1].X + 1) / 2;
			lv[nlev].Y = (lv[nlev-1].Y + 1) / 2;
			lv[nlev].Z = (lv[nlev-1].Z + 1) / 2;
			++nlev;
		}
	} else {
		for (l = 1; l < nlev; ++l) {
			lv[l].X = (lv[l-1].X + 1) / 2;
			lv[l].Y = (lv[l-1].Y + 1) / 2;
			lv[l].Z = (lv[l-1].Z + 1) / 2;
		}
	}
	fprintf(stderr, "tif2zar: %u x %u x %d %s, chunk %d, pixel %.6g um, "
	    "%d level%s (smallest %d x %d x %d)\n",
	    NX, NY, nZ, is_float ? "float32 -> uint16" : "uint16", chunk, pixel_um,
	    nlev, nlev == 1 ? "" : "s",
	    lv[nlev-1].X, lv[nlev-1].Y, lv[nlev-1].Z);

	if (MKDIR(zardir) != 0)
		die("cannot create %s (must not already exist)", zardir);
	/* .zattrs is written at the end: the omero display window needs the
	   actual data range, known only after the pass */
	for (l = 0; l < nlev; ++l) {
		Level	*L = &lv[l];

		L->nyc = (L->Y + chunk - 1) / chunk;
		L->nxc = (L->X + chunk - 1) / chunk;
		L->band = (uint16_t *)malloc((size_t)chunk * L->X * L->Y * sizeof(uint16_t));
		if (L->band == NULL) die("no memory for the level %s band", "pyramid");
		snprintf(path, sizeof(path), "%s/%d", zardir, l);     mkdir_ok(path);
		snprintf(path, sizeof(path), "%s/%d/0", zardir, l);   mkdir_ok(path);
		snprintf(path, sizeof(path), "%s/%d/0/0", zardir, l); mkdir_ok(path);
		write_zarray(l);
	}
#ifdef _OPENMP
	{
		const char	*e = getenv("TIF2ZAR_THREADS");

		nthr = (e != NULL && atoi(e) > 0) ? atoi(e) : omp_get_max_threads();
		if (nthr < 1) nthr = 1;
	}
#endif
	zbytes = (size_t)chunk * chunk * chunk * sizeof(uint16_t) + BLOSC_MAX_OVERHEAD;
	cbufs = (uint16_t *)malloc((size_t)nthr * chunk * chunk * chunk * sizeof(uint16_t));
	zbufs = (char *)malloc((size_t)nthr * zbytes);
	dsbuf = (uint16_t *)malloc(nlev > 1 ?
	    (size_t)lv[1].X * lv[1].Y * sizeof(uint16_t) : sizeof(uint16_t));
	if (is_float) frows = (float *)malloc((size_t)nthr * NX * sizeof(float));
	if (cbufs == NULL || zbufs == NULL || dsbuf == NULL ||
	    (is_float && frows == NULL))
		die("no memory for the %s buffers", "work");

	{	/* stream: read one chunk-deep band (slices in parallel), then
		   write + compress its chunks and cascade into the pyramid */
		int	zc0, nzc0 = (nZ + chunk - 1) / chunk, s0;

		for (zc0 = 0; zc0 < nzc0; ++zc0) {
			int	nz0 = (zc0 + 1) * chunk <= nZ ? chunk : nZ - zc0 * chunk;

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthr)
#endif
			for (s0 = 0; s0 < nz0; ++s0) {
				int	t = 0;
#ifdef _OPENMP
				t = omp_get_thread_num();
#endif
				read_slice(zc0 * chunk + s0,
				    lv[0].band + (size_t)s0 * NY * NX,
				    is_float ? frows + (size_t)t * NX : NULL);
			}
			lv[0].fill = nz0;
			lv[0].zout += nz0;
			flush_band(0);
			fprintf(stderr, "band %d / %d\r", zc0 + 1, nzc0);
		}
	}
	for (l = 1; l < nlev; ++l)	/* cascade the partial bands */
		flush_band(l);
	write_root_json();		/* omero window = actual data range */

	fprintf(stderr, "\ntif2zar: done (%ld chunk files in %d levels, %d thread%s)\n",
	    n_chunkfiles, nlev, nthr, nthr == 1 ? "" : "s");
	if (strcmp(codec, "none") != 0 && n_cmpbytes > 0)
		fprintf(stderr, "tif2zar: blosc-%s clevel %d: %.2f x compression "
		    "(%.1f -> %.1f MB)\n", codec, clevel,
		    (double)n_rawbytes / (double)n_cmpbytes,
		    (double)n_rawbytes / 1048576.0, (double)n_cmpbytes / 1048576.0);
	if (is_float && (n_clip_lo + n_clip_hi + n_nan) > 0)
		fprintf(stderr, "tif2zar: quantisation clipped %ld low / %ld high, "
		    "%ld NaN -> 0\n", n_clip_lo, n_clip_hi, n_nan);
	for (l = 0; l < nlev; ++l) free(lv[l].band);
	free(cbufs); free(zbufs); free(dsbuf); free(frows);

	/* command history, ct-rec convention */
	if ((f = fopen("cmd-hst.log", "a")) != NULL) {
		for (i = 0; i < argc; ++i) fprintf(f, "%s ", argv[i]);
		fprintf(f, "\n");
		fclose(f);
	}
	return 0;
}
