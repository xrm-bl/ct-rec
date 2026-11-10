/*
 * ct_gmem - GPU memory check for a CT data set before reconstructing it.
 *
 *   cd HiPic ; ct_gmem {layer} {center}
 *
 *   Run it inside the data directory (like his_rec / his_info): output.log,
 *            dark.img and q????.img, or dark.tif and q????.tif are read from
 *            the current directory (the format is taken from the dark frame,
 *            as the readers of hp_tg / ofct_srec do; RHP_O / RHP_D / RHP_Q
 *            override the file names as usual).
 *   layer    the layer whose rows are read for the axis estimate of a
 *            360-degree offset scan (default: the middle row).  Not used
 *            for a 180-degree scan.
 *   center   the rotation axis position of an offset scan (fractional
 *            allowed; "-" or omitted = estimate it from the layer with the
 *            one-layer search of his_rec / his_sg, see axis1.h).  Ignored for
 *            a 180-degree scan, whose GPU demand does not depend on the axis.
 *
 * The scan type is read from the angle column of output.log (sample views
 * reaching 350 degrees = 360-degree offset scan, the test ofct_rec and the
 * his_* programs apply).  The reconstruction size follows ct_rec / ofct_rec:
 *   180 degrees : N = Nx (detector width),      M = number of sample views
 *   offset      : N = 2 * (integer part of center), M = sample views / 2
 * and the GPU demand of one slice is what cbp.cu allocates in InitCBP
 * (G + PQ + F + SC + cuFFT work areas, EstimateCBPMemory) plus the three
 * buffers of the GPU ring removal (sort_filter_g.cu, unless KERNEL_SIZE
 * turns it off).  That demand is the same for every GPU program built on
 * cbp.cu (ct_rec / ofct_rec / his_rec, hp_tg / ofct_srec / his_srec, act_tg /
 * ofact_srec), which process one slice at a time; the whole-volume programs
 * only need more host memory (see the three-argument form of ofct_srec).
 * The CUDA context itself (a few hundred MB) is not included.  The result is
 * compared with the free memory of the device CUDA_GPU selects.
 *
 * PAD_THRESH (truncation pad, default on) enlarges the FFT length and is
 * honoured; so is KERNEL_SIZE.  Built without -DUSE_GPU the program reports
 * the buffer sizes from the same formula but without the cuFFT work areas
 * and without a device query.
 *
 * Build (makefileGPU.bat / MakefileGPU): ct_gmem.c error.c rhp_c.c
 * sort_filter_g.obj cbp.obj (any filter) + libtiff + cufft + cudart, with
 * -DONLY_CT_VIEWS -DFOM=float -DFloat=float like hp_tg.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "rhp.h"
#include "cbp.h"
#include "axis1.h"
#ifdef	USE_GPU
#include "sort_filter_g.h"
#endif

#define OFFSET_ANGLE	350.0		/* as ofct_rec / his_main.c */

extern void	Error(char *msg);

/* ---- row reader for axis1.h: -log transmittance of view t at one layer ---- */
typedef struct {
	HiPic		*hp;
	int		layer;
	FOM		*row;
	unsigned short	*raw;
} GmemCtx;

static void	GmemRow(void *vctx,int view,double *out)
{
	GmemCtx	*c=(GmemCtx *)vctx;
	FOM	*dst[1];
	int	x;

	dst[0]=c->row;
	ReadHiPicBand(c->hp,view,c->layer,c->layer,dst,c->raw);
	for (x=0; x<c->hp->Nx; x++) out[x]=(c->row[x]>0.0)?-log((double)c->row[x]):0.0;
}

static const char	*MiB(size_t b,char *buf)
{
	(void)sprintf(buf,"%.0f MiB",(double)b/1048576.0);
	return buf;
}

int	main(int argc,char **argv)
{
	HiPic		hp;
	struct stat	st;
	char		*env,b1[32],b2[32],b3[32],b4[32];
	const char	*fmt;
	double		amax=-1.0e300,center=0.0;
	int		t,offset,layer=-1,N,M,have_center=0;
	size_t		cbp=0,ring=0,need,hostp,hostf;

	if (argc>3 || (argc>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")))) {
	    (void)fprintf(stderr,
		"usage : ct_gmem {layer} {center}     (run inside the HiPic directory)\n"
		"        GPU memory one slice of ct_rec / ofct_rec (and every other GPU program\n"
		"        built on cbp.cu) needs for the data set in the current directory\n"
		"        (output.log, dark.img + q????.img or dark.tif + q????.tif).  180-degree\n"
		"        or 360-degree offset scan is read from the angles in output.log; for an\n"
		"        offset scan the rotation axis is taken from center or estimated from\n"
		"        layer (default: the middle row) as his_rec does when center is omitted.\n");
	    return 1;
	}
	if (argc>1 && strcmp(argv[1],"-")!=0) layer=atoi(argv[1]);
	if (argc>2 && strcmp(argv[2],"-")!=0) { center=atof(argv[2]); have_center=1; }

	/* which format the reader will pick (same rule as rhp_c.c) */
	if ((env=getenv("RHP_D"))!=NULL && *env!='\0') {
	    size_t	l=strlen(env);
	    fmt=(l>=4 && !strcmp(env+l-4,".tif")) || (l>=5 && !strcmp(env+l-5,".tiff"))?"TIFF (RHP_D)":"HiPic .img (RHP_D)";
	}
	else if (stat("dark.img",&st)==0) fmt="HiPic .img (dark.img)";
	else if (stat("dark.tif",&st)==0) fmt="TIFF (dark.tif)";
	else fmt="(no dark.img / dark.tif in the current directory)";

	InitReadHiPic(".",&hp);		/* reads output.log, dark and the I0 frames of the current directory */

	for (t=0; t<hp.Nt; t++) if (hp.OL[t].a>amax) amax=hp.OL[t].a;
	offset=(amax>=OFFSET_ANGLE);

	/* ---- the two numbers that matter: scan type / axis, and the GPU demand ---- */
	if (offset) {
	    if (hp.Nt%2) (void)fprintf(stderr,"warning: odd number of sample views (%d); the last one is not paired.\n",hp.Nt);
	    M=hp.Nt/2;
	    if (!have_center) {
		GmemCtx	c;
		if (layer<0) layer=hp.Ny/2;
		if (layer>=hp.Ny) { (void)fprintf(stderr,"layer %d is outside the frame (0 .. %d)\n",layer,hp.Ny-1); return 1; }
		c.hp=&hp; c.layer=layer;
		if ((c.row=(FOM *)malloc(sizeof(FOM)*(size_t)hp.Nx))==NULL ||
		    (c.raw=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)hp.Nx*(size_t)hp.Ny))==NULL)
		    Error("memory allocation error.");
		Axis1Quiet=1;			/* the summary goes into the details below */
		center=EstimateOffsetAxis1(hp.Nx,M,layer,GmemRow,&c);
		free(c.raw); free(c.row);
	    }
	    if (center<=0.0 || center>=hp.Nx) { (void)fprintf(stderr,"center %.3f is outside the detector (0 .. %d)\n",center,hp.Nx); return 1; }
	    N=2*(int)floor(center);		/* NN = 2*idcnt of ofct_rec */
	    (void)printf("rotation angle: 360 degree\ncenter of rotation: %.3f pixel\n",center);
	}
	else {
	    N=hp.Nx; M=hp.Nt;
	    (void)printf("rotation angle: 180 degree\ncenter of rotation: - (does not change the GPU demand)\n");
	}

	hostp=sizeof(Float)*(size_t)N*(size_t)M;
	hostf=sizeof(Float)*(size_t)N*(size_t)N;
#ifdef	USE_GPU
	{
	    char	name[256];
	    size_t	freeB,totalB;
	    int		ks=get_kernel_size_from_env();

	    cbp=EstimateCBPMemory(N,M);
	    ring=sort_filter_gpu_bytes((unsigned)N,(unsigned)M);
	    need=cbp+ring;
	    (void)printf("GPU memory: %s\n",MiB(need,b1));
	    (void)printf("\ndetails\n");
	    if (CBPDeviceInfo(name,sizeof name,&freeB,&totalB)==0)
		(void)printf("  device     : %s (CUDA_GPU=%s), %s free of %s -> %s\n",name,((env=getenv("CUDA_GPU"))!=NULL)?env:"0",
			     MiB(freeB,b2),MiB(totalB,b3),
			     (need<=freeB)?"OK":"NOT ENOUGH (Windows: the driver may spill to system RAM, very slow; Linux: allocation fails)");
	    else
		(void)printf("  device     : no CUDA device could be queried\n");
	    (void)printf("  GPU        : per slice %s (%.2f GiB) = CBP buffers + cuFFT %s + ring removal %s (KERNEL_SIZE %d); the CUDA context (a few hundred MB) is not included\n",
			 MiB(need,b4),(double)need/1073741824.0,MiB(cbp,b1),MiB(ring,b2),ks);
	}
#else
	{	/* the formula of cbp.cu without the cuFFT work areas */
	    int	l,nl=((env=getenv("PAD_THRESH"))!=NULL?atof(env):PAD_THRESH_DEFAULT)>0.0?N+N/2:N;
	    int	ks=((env=getenv("KERNEL_SIZE"))!=NULL && atoi(env)>-1)?atoi(env):5;
	    for (l=1; l<nl; l*=2) ;
	    cbp=8*(size_t)(l+1)*(size_t)(M+1)+4*(size_t)N*(size_t)N+8*(size_t)M;
	    ring=(ks<2)?0:12*(size_t)N*(size_t)M;
	    need=cbp+ring;
	    (void)printf("GPU memory: %s\n",MiB(need,b1));
	    (void)printf("\ndetails\n");
	    (void)printf("  device     : CPU build, no device query\n");
	    (void)printf("  GPU        : per slice %s (%.2f GiB) = CBP buffers %s + ring removal %s (KERNEL_SIZE %d); cuFFT work areas and the CUDA context are not included\n",
			 MiB(need,b4),(double)need/1073741824.0,MiB(cbp,b1),MiB(ring,b2),ks);
	}
#endif
	(void)printf("  slice      : %d x %d (N x M)%s, PAD_THRESH %s\n",N,M,
		     offset?" = stitched width 2 * (integer part of the axis) x half the sample views":" = detector width x sample views",
		     ((env=getenv("PAD_THRESH"))!=NULL)?env:"default");
	(void)printf("  host       : pinned memory per slice: projections %s, image %s\n",MiB(hostp,b3),MiB(hostf,b4));
	(void)printf("  data       : %s, %d x %d pixels, %d sample views, %d I0 frames, angles up to %.1f degrees\n",
		     fmt,hp.Nx,hp.Ny,hp.Nt,hp.Ni,amax);
	if (offset && !have_center)
	    (void)printf("  axis search: layer %d, %d view pairs, shift %d%+.3f, Rc %d%s (one-layer search of his_rec; give center to skip it)\n",
			 layer,Axis1Pairs,Axis1Shift,Axis1Dx,Axis1Rc,Axis1Edge?" - AT THE EDGE OF THE SEARCH RANGE, unreliable":"");
	else if (offset)
	    (void)printf("  axis       : center given on the command line\n");
	(void)printf("  programs   : the same demand for every GPU program built on cbp.cu (%s)\n",
		     offset?"ofct_rec / ofct_srec / ofact_srec, his_rec / his_srec on a 360-degree his":"ct_rec / hp_tg / act_tg, his_rec / his_srec on a 180-degree his");

	TermReadHiPic(&hp);
	return 0;
}
