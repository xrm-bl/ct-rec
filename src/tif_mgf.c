/*
 * tif_mgf.c - High-precision 2D Filter Software (Median + Gaussian)
 * C89 compliant implementation with ImageJ-level accuracy
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#ifdef _WIN32
#include "tiffio.h"
#include <windows.h>
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#define PATH_SEPARATOR "\\"
#define HOME_ENV "USERPROFILE"
#else
#include <tiffio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(path) mkdir(path, 0755)
#define PATH_SEPARATOR "/"
#define HOME_ENV "HOME"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/* Edge handling modes */
#define EDGE_MIRROR 0
#define EDGE_EXTEND 1
#define EDGE_ZERO 2

/* Global variables for logging */
static FILE* g_log_file = NULL;
static clock_t g_start_time;

/* Configuration structure */
typedef struct {
    int edge_mode;
    int use_32bit_precision;
} FilterConfig;

/* Function prototypes */
static void init_log(void);
static void close_log(void);
//static void write_log(const char* format, ...);
static int validate_arguments(int argc, char* argv[]);
static int process_tiff_file(const char* input_file, const char* output_file, 
                            int median_size, double gaussian_sigma);
static int apply_median_filter(void* data, int width, int height, int kernel_size, 
                              int is_16bit, FilterConfig* config);
static int apply_gaussian_filter_separable(void* data, int width, int height, 
                                          double sigma, int is_16bit, FilterConfig* config);
static void sort_array_8bit(unsigned char* arr, int size);
static void sort_array_16bit(unsigned short* arr, int size);
static void create_gaussian_kernel_1d(double* kernel, int size, double sigma);
static int calculate_kernel_size(double sigma);
static double* convert_to_float(void* data, int width, int height, int is_16bit);
static void convert_from_float(double* float_data, void* data, int width, int height, int is_16bit);
static double get_pixel_value_mirror(double* data, int x, int y, int width, int height);
static double get_pixel_value_extend(double* data, int x, int y, int width, int height);
static double get_pixel_value_zero(double* data, int x, int y, int width, int height);

static char *desc;

/* Main function */
int main(int argc, char* argv[])
{
    const char* input_file;
    const char* output_file;
    int median_kernel_size = 1;
    double gaussian_sigma = 0.5;
    int result;
    
    /* Initialize */
    g_start_time = clock();
//    init_log();
    
    /* Validate arguments */
    if (!validate_arguments(argc, argv)) {
//        close_log();
        return 1;
    }
    
    /* Parse arguments */
    input_file = argv[1];
    output_file = argv[2];
    
    if (argc > 3) {
        median_kernel_size = atoi(argv[3]);
        if (median_kernel_size < 0 || median_kernel_size > 25) {
//            write_log("Error: Invalid median kernel size (must be 0-25)");
            fprintf(stderr, "Error: Invalid median kernel size (must be 0-25)\n");
//            close_log();
            return 1;
        }
    }
    
    if (argc > 4) {
        gaussian_sigma = atof(argv[4]);
        if (gaussian_sigma < 0.0 || gaussian_sigma > 100.0) {
//            write_log("Error: Invalid gaussian sigma (must be 0.0-100.0)");
            fprintf(stderr, "Error: Invalid gaussian sigma (must be 0.0-100.0)\n");
//            close_log();
            return 1;
        }
    }
    
//    write_log("Starting tif_mgf processing (High-precision mode)");
//    write_log("Input file: %s", input_file);
//    write_log("Output file: %s", output_file);
//    write_log("Median kernel size: %d", median_kernel_size);
//    write_log("Gaussian sigma: %.3f", gaussian_sigma);
    
    /* Process the file */
    result = process_tiff_file(input_file, output_file, median_kernel_size, gaussian_sigma);
    
    if (result == 0) {
        double elapsed = (double)(clock() - g_start_time) / CLOCKS_PER_SEC;
//        write_log("Processing completed successfully. Total time: %.2f seconds", elapsed);
//        printf("Processing completed successfully. Total time: %.2f seconds\n", elapsed);
    }
    
//    close_log();

	
	
// append to log file
	FILE		*flog;
	int		i;
	if ((flog = fopen("cmd-hst.log", "a")) == NULL) {
		return(-1);
	}
	for (i = 0; i<argc; ++i) fprintf(flog, "%s ", argv[i]);
	fprintf(flog, "\n");
	fclose(flog);
  
	
	return result;
}

/* Initialize logging */
//static void init_log(void)
//{
//    char log_path[512];
//    char timestamp[64];
//    time_t now;
//    struct tm* tm_info;
//    char* home_dir;
//    
//    /* Get home directory */
//    home_dir = getenv(HOME_ENV);
//    if (!home_dir) {
//        fprintf(stderr, "Warning: Cannot find home directory\n");
//        return;
//    }
//    
//    /* Create com-log directory */
//    sprintf(log_path, "%s%scom-log", home_dir, PATH_SEPARATOR);
//    MKDIR(log_path);
//    
//    /* Create log file with timestamp */
//    time(&now);
//    tm_info = localtime(&now);
//    strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", tm_info);
//    sprintf(log_path, "%s%scom-log%stif_mgf_%s.log", home_dir, PATH_SEPARATOR, PATH_SEPARATOR, timestamp);
//    
//    g_log_file = fopen(log_path, "w");
//    if (!g_log_file) {
//        fprintf(stderr, "Warning: Cannot create log file\n");
//    }
//}

/* Close logging */
//static void close_log(void)
//{
//    if (g_log_file) {
//        fclose(g_log_file);
//        g_log_file = NULL;
//    }
//}
//
///* Write to log file */
//static void write_log(const char* format, ...)
//{
//    va_list args;
//    char timestamp[64];
//    time_t now;
//    struct tm* tm_info;
//    
//    if (!g_log_file) return;
//    
//    /* Get timestamp */
//    time(&now);
//    tm_info = localtime(&now);
//    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", tm_info);
//    
//    /* Write timestamp and message */
//    fprintf(g_log_file, "[%s] ", timestamp);
//    
//    va_start(args, format);
//    vfprintf(g_log_file, format, args);
//    va_end(args);
//    
//    fprintf(g_log_file, "\n");
//    fflush(g_log_file);
//}

/* Validate command line arguments */
static int validate_arguments(int argc, char* argv[])
{
    FILE* test_file;
    
    if (argc < 3 || argc > 5) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> [median_kernel_size] [gaussian_sigma]\n", argv[0]);
        fprintf(stderr, "  median_kernel_size: 0-25 (default: 1, 0 to skip)\n");
        fprintf(stderr, "  gaussian_sigma: 0.0-100.0 (default: 0.5, 0.0 to skip)\n");
//        write_log("Error: Invalid number of arguments");
        return 0;
    }
    
    /* Check if input file exists */
    test_file = fopen(argv[1], "rb");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot open input file: %s\n", argv[1]);
//        write_log("Error: Cannot open input file: %s", argv[1]);
        return 0;
    }
    fclose(test_file);
    
    /* Check if output directory is writable */
    test_file = fopen(argv[2], "wb");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot create output file: %s\n", argv[2]);
//        write_log("Error: Cannot create output file: %s", argv[2]);
        return 0;
    }
    fclose(test_file);
    remove(argv[2]);
    
    return 1;
}

/* Calculate optimal kernel size based on sigma (3-sigma rule) */
static int calculate_kernel_size(double sigma)
{
    int size;
    if (sigma <= 0.0) return 0;
    
    /* Kernel size should cover at least 3 standard deviations on each side */
    size = (int)ceil(3.0 * sigma) * 2 + 1;
    
    /* Ensure odd size */
    if (size % 2 == 0) size++;
    
    /* Limit maximum size for performance */
    if (size > 51) size = 51;
    
    return size;
}

/* Create 1D Gaussian kernel with proper normalization */
static void create_gaussian_kernel_1d(double* kernel, int size, double sigma)
{
    int center = size / 2;
    double sum = 0.0;
    double two_sigma_sq = 2.0 * sigma * sigma;
    double norm_factor = 1.0 / sqrt(2.0 * 3.14159265358979323846 * sigma * sigma);
    int i;
    
    /* Calculate kernel values */
    for (i = 0; i < size; i++) {
        int distance = i - center;
        kernel[i] = norm_factor * exp(-(distance * distance) / two_sigma_sq);
        sum += kernel[i];
    }
    
    /* Normalize to ensure sum = 1.0 */
    for (i = 0; i < size; i++) {
        kernel[i] /= sum;
    }
}

/* Convert data to floating point for high-precision processing */
static double* convert_to_float(void* data, int width, int height, int is_16bit)
{
    double* float_data;
    int i, total_pixels = width * height;
    
    float_data = (double*)malloc(total_pixels * sizeof(double));
    if (!float_data) return NULL;
    
    if (is_16bit) {
        unsigned short* src = (unsigned short*)data;
        for (i = 0; i < total_pixels; i++) {
            float_data[i] = (double)src[i];
        }
    } else {
        unsigned char* src = (unsigned char*)data;
        for (i = 0; i < total_pixels; i++) {
            float_data[i] = (double)src[i];
        }
    }
    
    return float_data;
}

/* Convert floating point data back to original format */
static void convert_from_float(double* float_data, void* data, int width, int height, int is_16bit)
{
    int i, total_pixels = width * height;
    
    if (is_16bit) {
        unsigned short* dst = (unsigned short*)data;
        for (i = 0; i < total_pixels; i++) {
            double val = float_data[i] + 0.5;  /* Round */
            if (val < 0.0) val = 0.0;
            if (val > 65535.0) val = 65535.0;
            dst[i] = (unsigned short)val;
        }
    } else {
        unsigned char* dst = (unsigned char*)data;
        for (i = 0; i < total_pixels; i++) {
            double val = float_data[i] + 0.5;  /* Round */
            if (val < 0.0) val = 0.0;
            if (val > 255.0) val = 255.0;
            dst[i] = (unsigned char)val;
        }
    }
}

/* Pixel access with mirror boundary condition */
static double get_pixel_value_mirror(double* data, int x, int y, int width, int height)
{
    if (x < 0) x = -x;
    if (x >= width) x = 2 * width - x - 2;
    if (y < 0) y = -y;
    if (y >= height) y = 2 * height - y - 2;
    
    return data[y * width + x];
}

/* Pixel access with extend boundary condition */
static double get_pixel_value_extend(double* data, int x, int y, int width, int height)
{
    if (x < 0) x = 0;
    if (x >= width) x = width - 1;
    if (y < 0) y = 0;
    if (y >= height) y = height - 1;
    
    return data[y * width + x];
}

/* Pixel access with zero boundary condition */
static double get_pixel_value_zero(double* data, int x, int y, int width, int height)
{
    if (x < 0 || x >= width || y < 0 || y >= height) {
        return 0.0;
    }
    return data[y * width + x];
}

/* Process TIFF file */
static int process_tiff_file(const char* input_file, const char* output_file, 
                            int median_size, double gaussian_sigma)
{
    TIFF* in_tiff;
    TIFF* out_tiff;
    uint32 width, height;
    uint16 bits_per_sample, samples_per_pixel;
    void* data;
    size_t scanline_size;
    uint32 row;
    int is_16bit;
    clock_t filter_start;
    FilterConfig config;
    
    /* Set configuration */
    config.edge_mode = EDGE_MIRROR;  /* Default to mirror mode like ImageJ */
    config.use_32bit_precision = 1;   /* Always use high precision */
    
    /* Open input TIFF */
    in_tiff = TIFFOpen(input_file, "r");
    if (!in_tiff) {
        fprintf(stderr, "Error: Cannot open input TIFF file\n");
//        write_log("Error: Cannot open input TIFF file");
        return 1;
    }
    
    /* Get image properties */
    TIFFGetField(in_tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(in_tiff, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(in_tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
//    TIFFGetField(in_tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    TIFFGetField(in_tiff, TIFFTAG_IMAGEDESCRIPTION, &desc);

//    write_log("Image properties: width=%u, height=%u, bits=%u", width, height, bits_per_sample);
    
    /* Check if image is supported */
    if (bits_per_sample != 8 && bits_per_sample != 16) {
        fprintf(stderr, "Error: Only 8-bit and 16-bit images are supported\n");
//        write_log("Error: Only 8-bit and 16-bit images are supported");
        TIFFClose(in_tiff);
        return 1;
    }
    
    is_16bit = (bits_per_sample == 16);
    scanline_size = TIFFScanlineSize(in_tiff);
    
    /* Allocate memory for the entire image */
    data = malloc(height * scanline_size);
    if (!data) {
        fprintf(stderr, "Error: Cannot allocate memory\n");
//        write_log("Error: Cannot allocate memory");
        TIFFClose(in_tiff);
        return 1;
    }
    
    /* Read the entire image */
//    write_log("Reading image data...");
    for (row = 0; row < height; row++) {
        if (TIFFReadScanline(in_tiff, (char*)data + row * scanline_size, row, 0) < 0) {
            fprintf(stderr, "Error: Cannot read scanline %u\n", row);
//            write_log("Error: Cannot read scanline %u", row);
            free(data);
            TIFFClose(in_tiff);
            return 1;
        }
        
//        if (row % 100 == 0) {
//            write_log("Read %u/%u scanlines", row, height);
//        }
    }
    
    TIFFClose(in_tiff);
    
    /* Apply median filter if kernel size > 0 */
    if (median_size > 0) {
        filter_start = clock();
//        write_log("Applying median filter (kernel size: %d)...", median_size);
        if (apply_median_filter(data, width, height, median_size, is_16bit, &config) != 0) {
            free(data);
            return 1;
        }
//        write_log("Median filter completed in %.2f seconds", (double)(clock() - filter_start) / CLOCKS_PER_SEC);
    } else {
//        write_log("Skipping median filter (kernel size: 0)");
    }
    
    /* Apply Gaussian filter if sigma > 0 */
    if (gaussian_sigma > 0.0) {
        int kernel_size = calculate_kernel_size(gaussian_sigma);
        filter_start = clock();
//        write_log("Applying Gaussian filter (sigma: %.3f, kernel size: %d)...", gaussian_sigma, kernel_size);
        if (apply_gaussian_filter_separable(data, width, height, gaussian_sigma, is_16bit, &config) != 0) {
            free(data);
            return 1;
        }
//        write_log("Gaussian filter completed in %.2f seconds", (double)(clock() - filter_start) / CLOCKS_PER_SEC);
    } else {
//        write_log("Skipping Gaussian filter (sigma: 0.0)");
    }
    
    /* Open output TIFF */
    out_tiff = TIFFOpen(output_file, "w");
    if (!out_tiff) {
        fprintf(stderr, "Error: Cannot create output TIFF file\n");
//        write_log("Error: Cannot create output TIFF file");
        free(data);
        return 1;
    }
    
    /* Set output TIFF properties (same as input) */
    TIFFSetField(out_tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out_tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(out_tiff, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFSetField(out_tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(out_tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out_tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(out_tiff, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(out_tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(out_tiff, TIFFTAG_IMAGEDESCRIPTION, desc);
	TIFFSetField(out_tiff, TIFFTAG_ARTIST, "tif_mgf2_libtiff");


    TIFFSetField(out_tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    
    /* Write the processed image */
//    write_log("Writing output image...");
    for (row = 0; row < height; row++) {
        if (TIFFWriteScanline(out_tiff, (char*)data + row * scanline_size, row, 0) < 0) {
            fprintf(stderr, "Error: Cannot write scanline %u\n", row);
//            write_log("Error: Cannot write scanline %u", row);
            free(data);
            TIFFClose(out_tiff);
            return 1;
        }
        
//        if (row % 100 == 0) {
//            write_log("Wrote %u/%u scanlines", row, height);
//        }
    }
    
    TIFFClose(out_tiff);
    free(data);

//    write_log("Output file written successfully");
    return 0;
}

/* Apply median filter with high precision */
static int apply_median_filter(void* data, int width, int height, int kernel_size, 
                              int is_16bit, FilterConfig* config)
{
    void* temp_data;
    int x, y, i, j;
    int half_kernel = kernel_size / 2;
    int kernel_area = kernel_size * kernel_size;
    size_t data_size;
    
    /* Allocate temporary buffer */
    data_size = is_16bit ? sizeof(unsigned short) : sizeof(unsigned char);
    temp_data = malloc(width * height * data_size);
    if (!temp_data) {
        fprintf(stderr, "Error: Cannot allocate memory for median filter\n");
//        write_log("Error: Cannot allocate memory for median filter");
        return 1;
    }
    
    /* Apply median filter */
    #pragma omp parallel for private(x, i, j) schedule(dynamic)
    for (y = 0; y < height; y++) {
        unsigned char* window_8bit = NULL;
        unsigned short* window_16bit = NULL;
        int window_count;
        
        /* Allocate window buffer for this thread */
        if (is_16bit) {
            window_16bit = (unsigned short*)malloc(kernel_area * sizeof(unsigned short));
        } else {
            window_8bit = (unsigned char*)malloc(kernel_area * sizeof(unsigned char));
        }
        
        for (x = 0; x < width; x++) {
            window_count = 0;
            
            /* Collect values in the window with boundary handling */
            for (j = -half_kernel; j <= half_kernel; j++) {
                for (i = -half_kernel; i <= half_kernel; i++) {
                    int nx = x + i;
                    int ny = y + j;
                    
                    /* Handle boundaries based on edge mode */
                    if (config->edge_mode == EDGE_MIRROR) {
                        /* Mirror boundary */
                        if (nx < 0) nx = -nx;
                        if (nx >= width) nx = 2 * width - nx - 2;
                        if (ny < 0) ny = -ny;
                        if (ny >= height) ny = 2 * height - ny - 2;
                    } else if (config->edge_mode == EDGE_EXTEND) {
                        /* Extend boundary */
                        if (nx < 0) nx = 0;
                        if (nx >= width) nx = width - 1;
                        if (ny < 0) ny = 0;
                        if (ny >= height) ny = height - 1;
                    } else {
                        /* Zero boundary */
                        if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
                            if (is_16bit) {
                                window_16bit[window_count++] = 0;
                            } else {
                                window_8bit[window_count++] = 0;
                            }
                            continue;
                        }
                    }
                    
                    if (is_16bit) {
                        window_16bit[window_count++] = ((unsigned short*)data)[ny * width + nx];
                    } else {
                        window_8bit[window_count++] = ((unsigned char*)data)[ny * width + nx];
                    }
                }
            }
            
            /* Find median */
            if (is_16bit) {
                sort_array_16bit(window_16bit, window_count);
                ((unsigned short*)temp_data)[y * width + x] = window_16bit[window_count / 2];
            } else {
                sort_array_8bit(window_8bit, window_count);
                ((unsigned char*)temp_data)[y * width + x] = window_8bit[window_count / 2];
            }
        }
        
        /* Free window buffer */
        if (is_16bit) {
            free(window_16bit);
        } else {
            free(window_8bit);
        }
        
        /* Progress logging */
//        if (y % 50 == 0) {
//            #pragma omp critical
//            {
//                write_log("Median filter progress: %d/%d rows", y, height);
//            }
//        }
    }
    
    /* Copy result back */
    memcpy(data, temp_data, width * height * data_size);
    free(temp_data);
    
    return 0;
}

/* Apply Gaussian filter using separable implementation (like ImageJ) */
static int apply_gaussian_filter_separable(void* data, int width, int height, 
                                          double sigma, int is_16bit, FilterConfig* config)
{
    double* float_data;
    double* temp_data;
    double* kernel;
    int kernel_size;
    int half_kernel;
    int x, y, i;
    double sum;
    
    /* Calculate kernel size */
    kernel_size = calculate_kernel_size(sigma);
    if (kernel_size <= 0) return 0;
    
    half_kernel = kernel_size / 2;
//    write_log("Gaussian kernel size calculated: %d for sigma %.3f", kernel_size, sigma);
    
    /* Convert to floating point for high precision */
    float_data = convert_to_float(data, width, height, is_16bit);
    if (!float_data) {
        fprintf(stderr, "Error: Cannot allocate float buffer\n");
//        write_log("Error: Cannot allocate float buffer");
        return 1;
    }
    
    /* Allocate temporary buffer */
    temp_data = (double*)malloc(width * height * sizeof(double));
    if (!temp_data) {
        fprintf(stderr, "Error: Cannot allocate temp buffer\n");
//        write_log("Error: Cannot allocate temp buffer");
        free(float_data);
        return 1;
    }
    
    /* Create 1D Gaussian kernel */
    kernel = (double*)malloc(kernel_size * sizeof(double));
    if (!kernel) {
        fprintf(stderr, "Error: Cannot allocate kernel\n");
//        write_log("Error: Cannot allocate kernel");
        free(float_data);
        free(temp_data);
        return 1;
    }
    
    create_gaussian_kernel_1d(kernel, kernel_size, sigma);
    
    /* Log kernel values for verification */
//    write_log("Gaussian kernel center value: %.6f", kernel[half_kernel]);
//    write_log("Gaussian kernel sum: %.10f", 1.0);  /* Should be 1.0 after normalization */
    
    /* Apply horizontal pass */
//    write_log("Applying horizontal Gaussian pass...");
    #pragma omp parallel for private(x, i, sum) schedule(dynamic)
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            sum = 0.0;
            
            for (i = 0; i < kernel_size; i++) {
                int src_x = x + i - half_kernel;
                double pixel_value;
                
                /* Get pixel value with boundary handling */
                if (config->edge_mode == EDGE_MIRROR) {
                    pixel_value = get_pixel_value_mirror(float_data, src_x, y, width, height);
                } else if (config->edge_mode == EDGE_EXTEND) {
                    pixel_value = get_pixel_value_extend(float_data, src_x, y, width, height);
                } else {
                    pixel_value = get_pixel_value_zero(float_data, src_x, y, width, height);
                }
                
                sum += pixel_value * kernel[i];
            }
            
            temp_data[y * width + x] = sum;
        }
        
//        if (y % 50 == 0) {
//            #pragma omp critical
//            {
//                write_log("Horizontal pass progress: %d/%d rows", y, height);
//            }
//        }
    }
    
    /* Apply vertical pass */
//    write_log("Applying vertical Gaussian pass...");
    #pragma omp parallel for private(y, i, sum) schedule(dynamic)
    for (x = 0; x < width; x++) {
        for (y = 0; y < height; y++) {
            sum = 0.0;
            
            for (i = 0; i < kernel_size; i++) {
                int src_y = y + i - half_kernel;
                double pixel_value;
                
                /* Get pixel value with boundary handling */
                if (config->edge_mode == EDGE_MIRROR) {
                    pixel_value = get_pixel_value_mirror(temp_data, x, src_y, width, height);
                } else if (config->edge_mode == EDGE_EXTEND) {
                    pixel_value = get_pixel_value_extend(temp_data, x, src_y, width, height);
                } else {
                    pixel_value = get_pixel_value_zero(temp_data, x, src_y, width, height);
                }
                
                sum += pixel_value * kernel[i];
            }
            
            float_data[y * width + x] = sum;
        }
        
//        if (x % 50 == 0) {
//            #pragma omp critical
//            {
//                write_log("Vertical pass progress: %d/%d columns", x, width);
//            }
//        }
    }
    
    /* Convert back to original data type */
    convert_from_float(float_data, data, width, height, is_16bit);
    
    /* Clean up */
    free(kernel);
    free(temp_data);
    free(float_data);
    
    return 0;
}

/* Sort array for 8-bit data (using optimized quicksort for larger arrays) */
static void sort_array_8bit(unsigned char* arr, int size)
{
    int i, j;
    unsigned char temp;
    
    /* Use insertion sort for small arrays (efficient for small kernel sizes) */
    if (size <= 20) {
        for (i = 1; i < size; i++) {
            temp = arr[i];
            j = i - 1;
            while (j >= 0 && arr[j] > temp) {
                arr[j + 1] = arr[j];
                j--;
            }
            arr[j + 1] = temp;
        }
    } else {
        /* For larger arrays, use quicksort (not shown for brevity) */
        /* In practice, qsort from stdlib.h could be used */
        for (i = 0; i < size - 1; i++) {
            for (j = 0; j < size - i - 1; j++) {
                if (arr[j] > arr[j + 1]) {
                    temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
    }
}

/* Sort array for 16-bit data (using optimized quicksort for larger arrays) */
static void sort_array_16bit(unsigned short* arr, int size)
{
    int i, j;
    unsigned short temp;
    
    /* Use insertion sort for small arrays (efficient for small kernel sizes) */
    if (size <= 20) {
        for (i = 1; i < size; i++) {
            temp = arr[i];
            j = i - 1;
            while (j >= 0 && arr[j] > temp) {
                arr[j + 1] = arr[j];
                j--;
            }
            arr[j + 1] = temp;
        }
    } else {
        /* For larger arrays, use quicksort (not shown for brevity) */
        /* In practice, qsort from stdlib.h could be used */
        for (i = 0; i < size - 1; i++) {
            for (j = 0; j < size - i - 1; j++) {
                if (arr[j] > arr[j + 1]) {
                    temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
    }
}