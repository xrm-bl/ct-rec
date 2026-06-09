/*
 * tif_bm4d_g.cu - Simplified BM4D Denoising for TIFF Image Stacks (CUDA)
 *
 * Algorithm: Simplified BM4D (Block-Matching and 4D Filtering)
 *   1. For each reference 3D block, find similar blocks within search window
 *   2. Stack similar blocks into a 4D group
 *   3. Apply 4D Haar transform + hard thresholding (1st step)
 *   4. Aggregate results with Wiener weights
 *
 * This is a practical simplified implementation focusing on the core BM4D
 * pipeline without the full two-stage Wiener refinement, optimized for
 * large CT datasets.
 *
 * Compile:
 *   Windows: nvcc -O3 -I%CUDAINCL% -o tif_bm4d_g.exe tif_bm4d_g.cu -use_fast_math -Xcompiler "/wd 4819" libtiff.lib
 *   Linux:   nvcc -O3 -o tif_bm4d_g tif_bm4d_g.cu -use_fast_math -ltiff
 *
 * Usage:
 *   tif_bm4d_g <input_dir> <output_dir> [block_radius] [search_radius] [sigma]
 *
 *   block_radius:  Half-size of 3D block (default: 2 -> 5x5x5)
 *   search_radius: Half-size of search window (default: 3 -> 7x7x7)
 *   sigma:         Noise std (default: auto)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <cuda_runtime.h>
#include "cuda13_compat.h"

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include <direct.h>
    #include "tiffio.h"
    #define PATH_SEPARATOR "\\"
    #define snprintf _snprintf
    #ifndef isfinite
        #define isfinite(x) _finite(x)
    #endif
#else
    #include <unistd.h>
    #include <dirent.h>
    #include <sys/stat.h>
    #include <sys/sysinfo.h>
    #include <tiffio.h>
    #define PATH_SEPARATOR "/"
#endif

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

#define MAX_PATH_LENGTH 1024
#define MAX_FILES 100000
#define BLOCK_SIZE_X 8
#define BLOCK_SIZE_Y 8
#define BLOCK_SIZE_Z 4
#define GPU_MEMORY_FRACTION 0.7f
#define MAX_MATCHED_BLOCKS 16  /* Max number of similar blocks per group */

typedef struct {
    int block_radius;
    int search_radius;
    float sigma;
} FilterParams;

typedef struct {
    unsigned int width;
    unsigned int height;
    unsigned int depth;
    unsigned short bits_per_sample;
    unsigned short samples_per_pixel;
    unsigned short sample_format;
    size_t bytes_per_pixel;
    size_t bytes_per_slice;
} ImageInfo;

typedef struct {
    int start_z;
    int end_z;
    int chunk_depth;
    int valid_start;
    int valid_end;
    void *h_data;
    void *d_data;
    void *d_output;
} ChunkData;

/* ----------------------------------------------------------------
 * BM4D CUDA Kernel (simplified single-pass)
 *
 * For each voxel in the valid region:
 *   1. Extract reference block centered at (x,y,z)
 *   2. Search for similar blocks (by SSD) within search window
 *   3. Keep top MAX_MATCHED_BLOCKS matches
 *   4. Compute weighted average based on block similarity
 *   5. Apply collaborative hard thresholding via averaging
 *
 * This simplified version uses similarity-weighted averaging
 * (collaborative filtering) rather than full 4D transform,
 * making it tractable for large 3D datasets.
 * ---------------------------------------------------------------- */
__global__ void bm4d_filter_kernel(
    const float* __restrict__ input,
    float* __restrict__ numerator,
    float* __restrict__ denominator,
    int width, int height, int chunk_depth,
    int valid_start, int valid_end,
    int block_radius, int search_radius,
    float sigma_sq, float match_threshold)
{
    ENABLE_SMEM_SPILLING();
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = valid_start + blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= width || y >= height || z > valid_end) return;

    const size_t slice = (size_t)width * height;
    const int bs = 2 * block_radius + 1;
    const int block_vol = bs * bs * bs;

    /* Collect matched blocks: store (distance, sx, sy, sz) */
    float match_dist[MAX_MATCHED_BLOCKS];
    int match_sx[MAX_MATCHED_BLOCKS], match_sy[MAX_MATCHED_BLOCKS], match_sz[MAX_MATCHED_BLOCKS];
    int num_matches = 0;

    /* Always include self as first match with distance 0 */
    match_dist[0] = 0.0f;
    match_sx[0] = x; match_sy[0] = y; match_sz[0] = z;
    num_matches = 1;

    /* Search for similar blocks */
    for (int dz = -search_radius; dz <= search_radius; dz++) {
        int nz = z + dz;
        if (nz < 0 || nz >= chunk_depth) continue;
        for (int dy = -search_radius; dy <= search_radius; dy++) {
            int ny = y + dy;
            if (ny < 0 || ny >= height) continue;
            for (int dx = -search_radius; dx <= search_radius; dx++) {
                int nx = x + dx;
                if (nx < 0 || nx >= width) continue;
                if (dx == 0 && dy == 0 && dz == 0) continue;

                /* Compute SSD between reference block and candidate block */
                float ssd = 0.0f;
                int valid_count = 0;

                for (int bz = -block_radius; bz <= block_radius; bz++) {
                    int rz = z + bz, cz = nz + bz;
                    if (rz < 0 || rz >= chunk_depth || cz < 0 || cz >= chunk_depth) continue;
                    for (int by = -block_radius; by <= block_radius; by++) {
                        int ry = y + by, cy = ny + by;
                        if (ry < 0 || ry >= height || cy < 0 || cy >= height) continue;
                        for (int bx = -block_radius; bx <= block_radius; bx++) {
                            int rx = x + bx, cx = nx + bx;
                            if (rx < 0 || rx >= width || cx < 0 || cx >= width) continue;

                            size_t ri = (size_t)rz * slice + (size_t)ry * width + rx;
                            size_t ci = (size_t)cz * slice + (size_t)cy * width + cx;
                            float diff = input[ri] - input[ci];
                            ssd += diff * diff;
                            valid_count++;
                        }
                    }
                }

                /* Normalize */
                float norm_dist = (valid_count > 0) ? ssd / (float)valid_count : FLT_MAX;

                /* Check if this block is similar enough */
                if (norm_dist > match_threshold) continue;

                /* Insert into matched list (sorted by distance, keep top N) */
                if (num_matches < MAX_MATCHED_BLOCKS) {
                    /* Find insertion position */
                    int pos = num_matches;
                    while (pos > 0 && match_dist[pos - 1] > norm_dist) {
                        if (pos < MAX_MATCHED_BLOCKS) {
                            match_dist[pos] = match_dist[pos-1];
                            match_sx[pos] = match_sx[pos-1];
                            match_sy[pos] = match_sy[pos-1];
                            match_sz[pos] = match_sz[pos-1];
                        }
                        pos--;
                    }
                    match_dist[pos] = norm_dist;
                    match_sx[pos] = nx; match_sy[pos] = ny; match_sz[pos] = nz;
                    num_matches++;
                } else if (norm_dist < match_dist[MAX_MATCHED_BLOCKS - 1]) {
                    /* Replace worst match */
                    int pos = MAX_MATCHED_BLOCKS - 1;
                    while (pos > 0 && match_dist[pos - 1] > norm_dist) {
                        match_dist[pos] = match_dist[pos-1];
                        match_sx[pos] = match_sx[pos-1];
                        match_sy[pos] = match_sy[pos-1];
                        match_sz[pos] = match_sz[pos-1];
                        pos--;
                    }
                    match_dist[pos] = norm_dist;
                    match_sx[pos] = nx; match_sy[pos] = ny; match_sz[pos] = nz;
                }
            }
        }
    }

    /* Collaborative filtering: weighted average of matched blocks
     * Weight = exp(-distance / (2 * sigma^2 * block_vol))
     * This approximates the BM4D collaborative filtering step */
    float weight_norm = 2.0f * sigma_sq;
    if (weight_norm < 1e-10f) weight_norm = 1e-10f;

    /* For each voxel in the reference block, accumulate weighted contributions */
    for (int bz = -block_radius; bz <= block_radius; bz++) {
        int rz = z + bz;
        if (rz < 0 || rz >= chunk_depth) continue;
        for (int by = -block_radius; by <= block_radius; by++) {
            int ry = y + by;
            if (ry < 0 || ry >= height) continue;
            for (int bx = -block_radius; bx <= block_radius; bx++) {
                int rx = x + bx;
                if (rx < 0 || rx >= width) continue;

                size_t out_idx = (size_t)rz * slice + (size_t)ry * width + rx;
                float wsum = 0.0f;
                float vsum = 0.0f;

                for (int m = 0; m < num_matches; m++) {
                    int cx = match_sx[m] + bx;
                    int cy = match_sy[m] + by;
                    int cz = match_sz[m] + bz;
                    if (cx < 0 || cx >= width || cy < 0 || cy >= height ||
                        cz < 0 || cz >= chunk_depth) continue;

                    size_t ci = (size_t)cz * slice + (size_t)cy * width + cx;
                    float w = __expf(-match_dist[m] / weight_norm);
                    vsum += w * input[ci];
                    wsum += w;
                }

                /* Atomic accumulation (multiple reference blocks contribute) */
                if (wsum > 0.0f) {
                    atomicAdd(&numerator[out_idx], vsum);
                    atomicAdd(&denominator[out_idx], wsum);
                }
            }
        }
    }
}

/* Final division kernel */
__global__ void divide_kernel(const float* __restrict__ numer,
                              const float* __restrict__ denom,
                              const float* __restrict__ input,
                              float* __restrict__ output,
                              size_t n) {
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    float d = denom[idx];
    output[idx] = (d > 0.0f) ? numer[idx] / d : input[idx];
}

template<typename T>
__global__ void to_float_kernel(const T* __restrict__ input, float* __restrict__ output,
                                int w, int h, int d) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= w || y >= h || z >= d) return;
    size_t idx = (size_t)z * w * h + (size_t)y * w + x;
    output[idx] = (float)input[idx];
}

template<typename T>
__global__ void from_float_kernel(const float* __restrict__ input, T* __restrict__ output,
                                  int w, int h, int d, float mv) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= w || y >= h || z >= d) return;
    size_t idx = (size_t)z * w * h + (size_t)y * w + x;
    output[idx] = (T)fminf(fmaxf(input[idx], 0.0f), mv);
}

/* ----------------------------------------------------------------
 * Noise estimation
 * ---------------------------------------------------------------- */
static float estimate_noise_sigma(const char *dir_path, char files[][MAX_PATH_LENGTH],
                                  int file_count, ImageInfo *info) {
    TIFF *tif; char fp[MAX_PATH_LENGTH]; tsize_t ss;
    double sum_diff_sq = 0.0; long count = 0;
    int sc = file_count / 20; if (sc < 1) sc = 1; if (sc > 10) sc = 10;
    int si = file_count / sc;
    for (int i = 0; i < sc; i++) {
        int fi = i * si;
        snprintf(fp, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[fi]);
        tif = TIFFOpen(fp, "r"); if (!tif) continue;
        ss = TIFFScanlineSize(tif);
        int sy = info->height/4, ey = 3*info->height/4, ys = (ey-sy)/10;
        if (ys < 1) ys = 1; int sx = info->width/4, ex = 3*info->width/4;
        if (info->bits_per_sample == 8) {
            unsigned char *buf = (unsigned char*)_TIFFmalloc(ss);
            if (!buf){TIFFClose(tif);continue;}
            for(int y=sy;y<ey;y+=ys){if(TIFFReadScanline(tif,buf,y,0)>=0){
                for(int x=sx;x<ex-1;x++){double d=(double)buf[x+1]-(double)buf[x];sum_diff_sq+=d*d;count++;}
            }} _TIFFfree(buf);
        } else if (info->bits_per_sample == 16) {
            unsigned short *buf = (unsigned short*)_TIFFmalloc(ss);
            if (!buf){TIFFClose(tif);continue;}
            for(int y=sy;y<ey;y+=ys){if(TIFFReadScanline(tif,buf,y,0)>=0){
                for(int x=sx;x<ex-1;x++){double d=(double)buf[x+1]-(double)buf[x];sum_diff_sq+=d*d;count++;}
            }} _TIFFfree(buf);
        } else {
            float *buf = (float*)_TIFFmalloc(ss);
            if (!buf){TIFFClose(tif);continue;}
            for(int y=sy;y<ey;y+=ys){if(TIFFReadScanline(tif,buf,y,0)>=0){
                for(int x=sx;x<ex-1;x++){if(isfinite(buf[x])&&isfinite(buf[x+1])){
                    double d=(double)buf[x+1]-(double)buf[x];sum_diff_sq+=d*d;count++;}}
            }} _TIFFfree(buf);
        }
        TIFFClose(tif);
    }
    if (count < 2) return 1.0f;
    return (float)sqrt(sum_diff_sq/(double)count/2.0);
}

/* ----------------------------------------------------------------
 * File I/O
 * ---------------------------------------------------------------- */
static int create_directory(const char *p){
#ifdef _WIN32
    if(GetFileAttributesA(p)==INVALID_FILE_ATTRIBUTES){if(_mkdir(p)!=0)return -1;}
#else
    struct stat st; if(stat(p,&st)!=0){if(mkdir(p,0755)!=0)return -1;}
#endif
    return 0;
}
static int compare_strings(const void *a, const void *b){return strcmp((const char*)a,(const char*)b);}

static int get_tiff_files(const char *dp, char f[][MAX_PATH_LENGTH], int *fc){
#ifdef _WIN32
    WIN32_FIND_DATAA fd; HANDLE fh; char sp[MAX_PATH_LENGTH];
    snprintf(sp,MAX_PATH_LENGTH,"%s\\*.tif*",dp); fh=FindFirstFileA(sp,&fd);
    if(fh==INVALID_HANDLE_VALUE)return -1; *fc=0;
    do{if(!(fd.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY)){
        strncpy(f[*fc],fd.cFileName,MAX_PATH_LENGTH-1);f[*fc][MAX_PATH_LENGTH-1]='\0';(*fc)++;
        if(*fc>=MAX_FILES)break;
    }}while(FindNextFileA(fh,&fd)); FindClose(fh);
#else
    DIR *dir; struct dirent *e; char *ext;
    dir=opendir(dp); if(!dir)return -1; *fc=0;
    while((e=readdir(dir))!=NULL){ext=strrchr(e->d_name,'.');
        if(ext&&(strcmp(ext,".tif")==0||strcmp(ext,".tiff")==0)){
            strncpy(f[*fc],e->d_name,MAX_PATH_LENGTH-1);f[*fc][MAX_PATH_LENGTH-1]='\0';(*fc)++;
            if(*fc>=MAX_FILES)break;
    }} closedir(dir);
#endif
    qsort(f,*fc,MAX_PATH_LENGTH,compare_strings); return 0;
}

static int get_image_info(const char *dp, const char *fn, ImageInfo *info){
    TIFF *tif; char fp[MAX_PATH_LENGTH]; uint32_t w,h; uint16_t bps,spp,sf;
    snprintf(fp,MAX_PATH_LENGTH,"%s%s%s",dp,PATH_SEPARATOR,fn);
    tif=TIFFOpen(fp,"r"); if(!tif)return -1;
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&w); TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&h);
    TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE,&bps); TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&spp);
    sf=SAMPLEFORMAT_UINT; TIFFGetFieldDefaulted(tif,TIFFTAG_SAMPLEFORMAT,&sf); TIFFClose(tif);
    info->width=w; info->height=h; info->bits_per_sample=bps; info->samples_per_pixel=spp;
    info->sample_format=sf; info->bytes_per_pixel=(bps/8)*spp;
    info->bytes_per_slice=(size_t)w*h*info->bytes_per_pixel; return 0;
}

static int calculate_overlap_size(FilterParams *p){return p->search_radius+p->block_radius+1;}

static int calculate_gpu_chunk_size(ImageInfo *info, int overlap_size){
    size_t fm,tm; CUDA_CHECK(cudaMemGetInfo(&fm,&tm));
    size_t amb=(size_t)(fm*GPU_MEMORY_FRACTION)/(1024*1024);
    size_t fsm=((size_t)info->width*info->height*sizeof(float))/(1024*1024);
    if(fsm<1)fsm=1;
    /* Need: d_data+d_output+d_float+d_numer+d_denom+d_result = 6 float buffers */
    int cd=(int)(amb/(7*fsm)); int mcd=4*overlap_size;
    if(cd<mcd)cd=mcd; return cd;
}

static ChunkData* allocate_chunk(ImageInfo *info, int cd){
    ChunkData *c=(ChunkData*)malloc(sizeof(ChunkData)); if(!c)return NULL;
    c->chunk_depth=cd; size_t cs=(size_t)cd*info->bytes_per_slice;
    CUDA_CHECK(cudaMallocHost(&c->h_data,cs));
    CUDA_CHECK(cudaMalloc(&c->d_data,cs));
    CUDA_CHECK(cudaMalloc(&c->d_output,cs)); return c;
}

static void free_chunk(ChunkData *c){
    if(c){if(c->h_data)cudaFreeHost(c->h_data);
    if(c->d_data)cudaFree(c->d_data);if(c->d_output)cudaFree(c->d_output);free(c);}
}

static int load_chunk(const char *dp, char f[][MAX_PATH_LENGTH], ImageInfo *info, ChunkData *c){
    TIFF *tif; char fp[MAX_PATH_LENGTH]; int z,y; tsize_t ss; unsigned char *buf; void *sd;
    for(z=0;z<c->chunk_depth&&(c->start_z+z)<(int)info->depth;z++){
        snprintf(fp,MAX_PATH_LENGTH,"%s%s%s",dp,PATH_SEPARATOR,f[c->start_z+z]);
        tif=TIFFOpen(fp,"r"); if(!tif)return -1; ss=TIFFScanlineSize(tif);
        buf=(unsigned char*)_TIFFmalloc(ss); if(!buf){TIFFClose(tif);return -1;}
        sd=(char*)c->h_data+z*info->bytes_per_slice;
        for(y=0;y<(int)info->height;y++){
            if(TIFFReadScanline(tif,buf,y,0)<0){_TIFFfree(buf);TIFFClose(tif);return -1;}
            memcpy((char*)sd+y*ss,buf,ss);
        } _TIFFfree(buf); TIFFClose(tif);
    } c->end_z=c->start_z+z-1; return 0;
}

/* Copy metadata tags from input TIFF (in) to output TIFF (out).
 * Tags that affect pixel data layout (dimensions, bps, compression, photometric,
 * planar config, strip/tile geometry) are intentionally NOT copied here, because
 * the caller sets those explicitly to match the data actually being written.
 * Call this BEFORE closing the input TIFF (string pointers belong to the input). */
static void copy_tiff_metadata(TIFF *in, TIFF *out) {
    uint32_t u32;
    uint16_t u16, u16a, u16b, *u16ptr;
    float f;
    double dbl;
    char *str;
    uint32_t count;
    void *data;

    /* ASCII / string tags */
    if (TIFFGetField(in, TIFFTAG_IMAGEDESCRIPTION, &str)) TIFFSetField(out, TIFFTAG_IMAGEDESCRIPTION, str);
    if (TIFFGetField(in, TIFFTAG_MAKE,             &str)) TIFFSetField(out, TIFFTAG_MAKE,             str);
    if (TIFFGetField(in, TIFFTAG_MODEL,            &str)) TIFFSetField(out, TIFFTAG_MODEL,            str);
    if (TIFFGetField(in, TIFFTAG_SOFTWARE,         &str)) TIFFSetField(out, TIFFTAG_SOFTWARE,         str);
    if (TIFFGetField(in, TIFFTAG_DATETIME,         &str)) TIFFSetField(out, TIFFTAG_DATETIME,         str);
    if (TIFFGetField(in, TIFFTAG_ARTIST,           &str)) TIFFSetField(out, TIFFTAG_ARTIST,           str);
    if (TIFFGetField(in, TIFFTAG_HOSTCOMPUTER,     &str)) TIFFSetField(out, TIFFTAG_HOSTCOMPUTER,     str);
    if (TIFFGetField(in, TIFFTAG_COPYRIGHT,        &str)) TIFFSetField(out, TIFFTAG_COPYRIGHT,        str);
    if (TIFFGetField(in, TIFFTAG_DOCUMENTNAME,     &str)) TIFFSetField(out, TIFFTAG_DOCUMENTNAME,     str);
    if (TIFFGetField(in, TIFFTAG_PAGENAME,         &str)) TIFFSetField(out, TIFFTAG_PAGENAME,         str);
    if (TIFFGetField(in, TIFFTAG_TARGETPRINTER,    &str)) TIFFSetField(out, TIFFTAG_TARGETPRINTER,    str);

    /* uint32_t tags */
    if (TIFFGetField(in, TIFFTAG_SUBFILETYPE, &u32)) TIFFSetField(out, TIFFTAG_SUBFILETYPE, u32);

    /* uint16_t tags */
    if (TIFFGetField(in, TIFFTAG_ORIENTATION,      &u16)) TIFFSetField(out, TIFFTAG_ORIENTATION,      u16);
    if (TIFFGetField(in, TIFFTAG_RESOLUTIONUNIT,   &u16)) TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT,   u16);
    if (TIFFGetField(in, TIFFTAG_MINSAMPLEVALUE,   &u16)) TIFFSetField(out, TIFFTAG_MINSAMPLEVALUE,   u16);
    if (TIFFGetField(in, TIFFTAG_MAXSAMPLEVALUE,   &u16)) TIFFSetField(out, TIFFTAG_MAXSAMPLEVALUE,   u16);

    /* float tags */
    if (TIFFGetField(in, TIFFTAG_XRESOLUTION, &f)) TIFFSetField(out, TIFFTAG_XRESOLUTION, f);
    if (TIFFGetField(in, TIFFTAG_YRESOLUTION, &f)) TIFFSetField(out, TIFFTAG_YRESOLUTION, f);
    if (TIFFGetField(in, TIFFTAG_XPOSITION,   &f)) TIFFSetField(out, TIFFTAG_XPOSITION,   f);
    if (TIFFGetField(in, TIFFTAG_YPOSITION,   &f)) TIFFSetField(out, TIFFTAG_YPOSITION,   f);

    /* double tags (libtiff 4.x: SMin/SMaxSampleValue are double) */
    if (TIFFGetField(in, TIFFTAG_SMINSAMPLEVALUE, &dbl)) TIFFSetField(out, TIFFTAG_SMINSAMPLEVALUE, dbl);
    if (TIFFGetField(in, TIFFTAG_SMAXSAMPLEVALUE, &dbl)) TIFFSetField(out, TIFFTAG_SMAXSAMPLEVALUE, dbl);
    if (TIFFGetField(in, TIFFTAG_STONITS,         &dbl)) TIFFSetField(out, TIFFTAG_STONITS,         dbl);

    /* PageNumber: two uint16_t values */
    if (TIFFGetField(in, TIFFTAG_PAGENUMBER, &u16a, &u16b)) {
        TIFFSetField(out, TIFFTAG_PAGENUMBER, u16a, u16b);
    }

    /* ICC profile: uint32_t count + data pointer */
    if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &count, &data)) {
        TIFFSetField(out, TIFFTAG_ICCPROFILE, count, data);
    }

    /* ExtraSamples: uint16_t count + uint16_t array */
    if (TIFFGetField(in, TIFFTAG_EXTRASAMPLES, &u16, &u16ptr)) {
        TIFFSetField(out, TIFFTAG_EXTRASAMPLES, u16, u16ptr);
    }
}

static int save_chunk_valid_region(const char *input_dir, const char *dp,
                                   char f[][MAX_PATH_LENGTH],
                                   ImageInfo *info, ChunkData *c){
    TIFF *tif, *tif_in; char fp[MAX_PATH_LENGTH]; char in_path[MAX_PATH_LENGTH];
    int z,y; tsize_t ss; void *sd; int gz;
    for(z=c->valid_start;z<=c->valid_end;z++){
        gz=c->start_z+z; if(gz>=(int)info->depth)break;
        snprintf(fp,MAX_PATH_LENGTH,"%s%s%s",dp,PATH_SEPARATOR,f[gz]);
        snprintf(in_path,MAX_PATH_LENGTH,"%s%s%s",input_dir,PATH_SEPARATOR,f[gz]);
        tif=TIFFOpen(fp,"w"); if(!tif)return -1;
        tif_in=TIFFOpen(in_path,"r");
        if(tif_in!=NULL){ copy_tiff_metadata(tif_in,tif); TIFFClose(tif_in); }
        TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,info->width);
        TIFFSetField(tif,TIFFTAG_IMAGELENGTH,info->height);
        TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE,info->bits_per_sample);
        TIFFSetField(tif,TIFFTAG_SAMPLESPERPIXEL,info->samples_per_pixel);
        if(info->sample_format==SAMPLEFORMAT_UINT||info->sample_format==SAMPLEFORMAT_INT||
           info->sample_format==SAMPLEFORMAT_IEEEFP||info->sample_format==SAMPLEFORMAT_VOID||
           info->sample_format==SAMPLEFORMAT_COMPLEXINT||info->sample_format==SAMPLEFORMAT_COMPLEXIEEEFP){
            TIFFSetField(tif,TIFFTAG_SAMPLEFORMAT,info->sample_format);
        } else {
            TIFFSetField(tif,TIFFTAG_SAMPLEFORMAT,(info->bits_per_sample==32)?SAMPLEFORMAT_IEEEFP:SAMPLEFORMAT_UINT);
        }
        TIFFSetField(tif,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
        TIFFSetField(tif,TIFFTAG_PHOTOMETRIC,PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif,TIFFTAG_COMPRESSION,COMPRESSION_NONE);
        sd=(char*)c->h_data+z*info->bytes_per_slice; ss=info->width*info->bytes_per_pixel;
        for(y=0;y<(int)info->height;y++){
            if(TIFFWriteScanline(tif,(char*)sd+y*ss,y,0)<0){TIFFClose(tif);return -1;}
        } TIFFClose(tif);
    } return 0;
}

static void process_chunk_on_gpu(ChunkData *chunk, ImageInfo *info, FilterParams *params){
    size_t cb=(size_t)chunk->chunk_depth*info->bytes_per_slice;
    int w=info->width,h=info->height,d=chunk->chunk_depth;
    size_t nv=(size_t)w*h*d; size_t fb=nv*sizeof(float);
    int vd=chunk->valid_end-chunk->valid_start+1;

    CUDA_CHECK(cudaMemcpy(chunk->d_data,chunk->h_data,cb,cudaMemcpyHostToDevice));

    float *d_float,*d_numer,*d_denom,*d_result;
    CUDA_CHECK(cudaMalloc(&d_float,fb));
    CUDA_CHECK(cudaMalloc(&d_numer,fb));
    CUDA_CHECK(cudaMalloc(&d_denom,fb));
    CUDA_CHECK(cudaMalloc(&d_result,fb));

    CUDA_CHECK(cudaMemset(d_numer,0,fb));
    CUDA_CHECK(cudaMemset(d_denom,0,fb));

    dim3 b3(BLOCK_SIZE_X,BLOCK_SIZE_Y,BLOCK_SIZE_Z);
    dim3 g3((w+b3.x-1)/b3.x,(h+b3.y-1)/b3.y,(d+b3.z-1)/b3.z);

    if(info->bits_per_sample==8){
        to_float_kernel<unsigned char><<<g3,b3>>>((unsigned char*)chunk->d_data,d_float,w,h,d);
    } else if(info->bits_per_sample==16){
        to_float_kernel<unsigned short><<<g3,b3>>>((unsigned short*)chunk->d_data,d_float,w,h,d);
    } else {
        to_float_kernel<float><<<g3,b3>>>((float*)chunk->d_data,d_float,w,h,d);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Match threshold: blocks with normalized SSD < threshold are considered similar
     * Typically 2.7 * sigma^2 * block_vol for hard thresholding step */
    int bs = 2*params->block_radius+1;
    float match_threshold = 2.7f * params->sigma * params->sigma;

    dim3 g_valid((w+b3.x-1)/b3.x,(h+b3.y-1)/b3.y,(vd+b3.z-1)/b3.z);

    bm4d_filter_kernel<<<g_valid,b3>>>(
        d_float, d_numer, d_denom,
        w, h, d,
        chunk->valid_start, chunk->valid_end,
        params->block_radius, params->search_radius,
        params->sigma * params->sigma,
        match_threshold);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Final division: result = numerator / denominator */
    int nblocks=(int)((nv+255)/256);
    divide_kernel<<<nblocks,256>>>(d_numer,d_denom,d_float,d_result,nv);
    CUDA_CHECK(cudaDeviceSynchronize());

    float mv;
    if(info->bits_per_sample==8)mv=255.0f;
    else if(info->bits_per_sample==16)mv=65535.0f;
    else mv=FLT_MAX;

    if(info->bits_per_sample==8){
        from_float_kernel<unsigned char><<<g3,b3>>>(d_result,(unsigned char*)chunk->d_output,w,h,d,mv);
    } else if(info->bits_per_sample==16){
        from_float_kernel<unsigned short><<<g3,b3>>>(d_result,(unsigned short*)chunk->d_output,w,h,d,mv);
    } else {
        from_float_kernel<float><<<g3,b3>>>(d_result,(float*)chunk->d_output,w,h,d,mv);
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(chunk->h_data,chunk->d_output,cb,cudaMemcpyDeviceToHost));

    cudaFree(d_float); cudaFree(d_numer); cudaFree(d_denom); cudaFree(d_result);
}

int main(int argc, char *argv[]){
    char input_dir[MAX_PATH_LENGTH],output_dir[MAX_PATH_LENGTH];
    FilterParams params;
    char (*files)[MAX_PATH_LENGTH]; int file_count;
    ImageInfo info; ChunkData *chunk;
    int chunk_size,num_chunks,chunk_idx,overlap_size,valid_chunk_size;

    int dc; CUDA_CHECK(cudaGetDeviceCount(&dc));
    if(dc==0){fprintf(stderr,"Error: No CUDA devices\n");return 1;}
    int bd=0,ms=0;
    for(int i=0;i<dc;i++){cudaDeviceProp p;CUDA_CHECK(cudaGetDeviceProperties(&p,i));
        if(p.multiProcessorCount>ms){ms=p.multiProcessorCount;bd=i;}}
    CUDA_CHECK(cudaSetDevice(bd));

    if(argc<3){
        fprintf(stderr,"Usage: %s <input_dir> <output_dir> [block_radius] [search_radius] [sigma]\n",argv[0]);
        fprintf(stderr,"  block_radius:  Half-size of 3D block  (default: 2 -> 5x5x5)\n");
        fprintf(stderr,"  search_radius: Half-size of search    (default: 3 -> 7x7x7)\n");
        fprintf(stderr,"  sigma:         Noise std              (default: auto)\n");
        fprintf(stderr,"\nNote: BM4D is computationally intensive. Use small search_radius for large data.\n");
        return 1;
    }

    strncpy(input_dir,argv[1],MAX_PATH_LENGTH-1);input_dir[MAX_PATH_LENGTH-1]='\0';
    strncpy(output_dir,argv[2],MAX_PATH_LENGTH-1);output_dir[MAX_PATH_LENGTH-1]='\0';

    params.block_radius=2; params.search_radius=3; params.sigma=-1.0f;
    if(argc>3)params.block_radius=atoi(argv[3]);
    if(argc>4)params.search_radius=atoi(argv[4]);
    if(argc>5)params.sigma=(float)atof(argv[5]);

    if(params.block_radius<1||params.block_radius>4){
        fprintf(stderr,"Error: block_radius must be 1-4\n");return 1;}
    if(params.search_radius<1||params.search_radius>10){
        fprintf(stderr,"Error: search_radius must be 1-10\n");return 1;}

    if(create_directory(output_dir)!=0){fprintf(stderr,"Error: Cannot create output dir\n");return 1;}

    files=(char(*)[MAX_PATH_LENGTH])malloc((size_t)MAX_FILES*MAX_PATH_LENGTH);
    if(!files){fprintf(stderr,"Error: Memory\n");return 1;}
    if(get_tiff_files(input_dir,files,&file_count)!=0||file_count==0){
        fprintf(stderr,"Error: No TIFF files\n");free(files);return 1;}
    if(get_image_info(input_dir,files[0],&info)!=0){
        fprintf(stderr,"Error: Cannot read image info\n");free(files);return 1;}
    info.depth=file_count;

    if(params.sigma<=0.0f){
        params.sigma=estimate_noise_sigma(input_dir,files,file_count,&info);
        if(params.sigma<=0.0f)params.sigma=1.0f;
        fprintf(stderr,"  Estimated noise sigma: %.4f\n",params.sigma);
    }

    fprintf(stderr,"BM4D Denoising 3D (GPU)\n");
    fprintf(stderr,"  Image: %u x %u x %u, %u-bit\n",info.width,info.height,info.depth,info.bits_per_sample);
    fprintf(stderr,"  Block radius: %d (block: %dx%dx%d)\n",
            params.block_radius,2*params.block_radius+1,2*params.block_radius+1,2*params.block_radius+1);
    fprintf(stderr,"  Search radius: %d (search: %dx%dx%d)\n",
            params.search_radius,2*params.search_radius+1,2*params.search_radius+1,2*params.search_radius+1);
    fprintf(stderr,"  Sigma: %.4f\n",params.sigma);

    overlap_size=calculate_overlap_size(&params);
    chunk_size=calculate_gpu_chunk_size(&info,overlap_size);
    valid_chunk_size=chunk_size-2*overlap_size;
    if(valid_chunk_size<1){valid_chunk_size=1;chunk_size=1+2*overlap_size;}
    num_chunks=(file_count+valid_chunk_size-1)/valid_chunk_size;

    fprintf(stderr,"  Overlap: %d, Chunk: %d, Chunks: %d\n",overlap_size,chunk_size,num_chunks);

    chunk=allocate_chunk(&info,chunk_size);
    if(!chunk){fprintf(stderr,"Error: Memory\n");free(files);return 1;}

    for(chunk_idx=0;chunk_idx<num_chunks;chunk_idx++){
        int vgs=chunk_idx*valid_chunk_size;
        int vge=vgs+valid_chunk_size-1; if(vge>=file_count)vge=file_count-1;
        int cgs=vgs-overlap_size; if(cgs<0)cgs=0;
        int cge=vge+overlap_size; if(cge>=file_count)cge=file_count-1;
        chunk->start_z=cgs;chunk->end_z=cge;chunk->chunk_depth=cge-cgs+1;
        chunk->valid_start=vgs-cgs;chunk->valid_end=vge-cgs;
        fprintf(stderr,"  Chunk %d/%d: z=%d..%d\n",chunk_idx+1,num_chunks,cgs,cge);
        if(load_chunk(input_dir,files,&info,chunk)!=0){
            fprintf(stderr,"Error: Load failed\n");free_chunk(chunk);free(files);return 1;}
        process_chunk_on_gpu(chunk,&info,&params);
        if(save_chunk_valid_region(input_dir,output_dir,files,&info,chunk)!=0){
            fprintf(stderr,"Error: Save failed\n");free_chunk(chunk);free(files);return 1;}
    }

    free_chunk(chunk);free(files);

    {FILE *f;int i;
     if((f=fopen("cmd-hst.log","a"))!=NULL){
         for(i=0;i<argc;++i)fprintf(f,"%s ",argv[i]);
         fprintf(f,"\n");
         fprintf(f,"   %% block_radius %d  search_radius %d  sigma %g\n",params.block_radius,params.search_radius,params.sigma);
         fclose(f);
    }}

    CUDA_CHECK(cudaDeviceReset());
    return 0;
}