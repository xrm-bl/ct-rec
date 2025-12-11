/*
 * sif_fast.c - High-performance TIFF image writer
 * 
 * Created: 2025-12-11
 * AI: Claude Opus 4.5 (Anthropic)
 * 
 * Optimizations:
 * 1. Large buffer I/O for reduced system calls
 * 2. Buffered file writing
 * 3. Optimized LZW encoder
 * 4. Pre-allocated strip buffers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cell.h"
#include "sif.h"

#include <limits.h>

#if UINT_MAX == 4294967295
typedef unsigned int Word4;
#elif ULONG_MAX == 4294967295
typedef unsigned long Word4;
#else
#error "No suitable 32-bit unsigned integer type found"
#endif

#if UINT_MAX == 65535
typedef unsigned int Word2;
#elif USHRT_MAX == 65535
typedef unsigned short Word2;
#else
#error "No suitable 16-bit unsigned integer type found"
#endif

typedef unsigned char Byte;

static FILE *file;
static Word4 bytes;

/* Write buffer for better I/O performance */
#define WRITE_BUFFER_SIZE (1024 * 1024)  /* 1MB buffer */
static Byte *write_buffer = NULL;
static size_t write_pos = 0;

static void InitWriteBuffer(void)
{
    if (write_buffer == NULL) {
        write_buffer = (Byte *)malloc(WRITE_BUFFER_SIZE);
    }
    write_pos = 0;
}

static void FlushWriteBuffer(void)
{
    if (write_pos > 0) {
        fwrite(write_buffer, 1, write_pos, file);
        write_pos = 0;
    }
}

static void WriteByte(Byte b)
{
    if (write_pos >= WRITE_BUFFER_SIZE) {
        FlushWriteBuffer();
    }
    write_buffer[write_pos++] = b;
}

static void FreeWriteBuffer(void)
{
    if (write_buffer != NULL) {
        free(write_buffer);
        write_buffer = NULL;
    }
}

#ifdef LZW
static int LZWnum, LZWlen;
static Word4 LZWbuf;

/* Buffered LZW output */
static Byte *lzw_output = NULL;
static size_t lzw_output_size = 0;
static size_t lzw_output_capacity = 0;

static void InitLZWOutput(void)
{
    lzw_output_capacity = 1024 * 1024;  /* Start with 1MB */
    lzw_output = (Byte *)malloc(lzw_output_capacity);
    lzw_output_size = 0;
}

static void LZWOutputByte(Byte b)
{
    if (lzw_output_size >= lzw_output_capacity) {
        lzw_output_capacity *= 2;
        lzw_output = (Byte *)realloc(lzw_output, lzw_output_capacity);
    }
    lzw_output[lzw_output_size++] = b;
}

static void LZWPutCode(int code)
{
    LZWnum += LZWlen;
    LZWbuf = (LZWbuf << LZWlen) | code;
    while (LZWnum >= 8) {
        LZWnum -= 8;
        LZWOutputByte((Byte)((LZWbuf >> LZWnum) & 0xff));
        ++bytes;
    }
}

static Word4 LZWstr[9001];
static int LZWcode, LZWroof, LZWhead;

static void InitEncode(void)
{
    int idx;

    bytes = 0;
    LZWnum = 0;
    LZWlen = 9;
    LZWPutCode(256);
    for (idx = 0; idx < 9001; idx++) LZWstr[idx] = 0;
    LZWcode = 258;
    LZWroof = 512;
    LZWhead = 4096;
    
    InitLZWOutput();
}

static void Encode(int tail)
{
    Word4 h_t;
    int idx, step;

    if (LZWhead != 4096) {
        h_t = (((Word4)LZWhead << 8) | tail) << 12;
        if (LZWstr[idx = LZWhead ^ (tail << 5)] != 0) {
            if ((LZWstr[idx] & 0xfffff000) == h_t) {
                LZWhead = (int)(LZWstr[idx] & 0xfff);
                return;
            }
            step = (idx == 0) ? 9000 : idx;
            while (LZWstr[idx = (idx + step) % 9001] != 0)
                if ((LZWstr[idx] & 0xfffff000) == h_t) {
                    LZWhead = (int)(LZWstr[idx] & 0xfff);
                    return;
                }
        }
        LZWPutCode(LZWhead);
        if (LZWcode != 4096 - 1 - 1) {
            LZWstr[idx] = h_t | LZWcode;
            if (++LZWcode == LZWroof) {
                ++LZWlen;
                LZWroof <<= 1;
            }
        } else {
            LZWPutCode(256);
            LZWlen = 9;
            for (idx = 0; idx < 9001; idx++) LZWstr[idx] = 0;
            LZWcode = 258;
            LZWroof = 512;
        }
    }
    LZWhead = tail;
}

static Word4 TermEncode(void)
{
    if (LZWhead != 4096) {
        LZWPutCode(LZWhead);
        if (++LZWcode == LZWroof) ++LZWlen;
    }
    LZWPutCode(257);
    if (LZWnum > 0) {
        LZWOutputByte((Byte)((LZWbuf << (8 - LZWnum)) & 0xff));
        ++bytes;
    }
    
    /* Write buffered LZW output to file */
    fwrite(lzw_output, 1, lzw_output_size, file);
    free(lzw_output);
    lzw_output = NULL;
    lzw_output_size = 0;
    
    return bytes;
}

#else
/* Non-compressed output with buffering */
static Byte *nic_output = NULL;
static size_t nic_output_size = 0;
static size_t nic_output_capacity = 0;

static void InitEncode(void)
{
    bytes = 0;
    nic_output_capacity = 1024 * 1024;
    nic_output = (Byte *)malloc(nic_output_capacity);
    nic_output_size = 0;
}

static void Encode(int code)
{
    if (nic_output_size >= nic_output_capacity) {
        nic_output_capacity *= 2;
        nic_output = (Byte *)realloc(nic_output, nic_output_capacity);
    }
    nic_output[nic_output_size++] = (Byte)code;
    ++bytes;
}

static Word4 TermEncode(void)
{
    fwrite(nic_output, 1, nic_output_size, file);
    free(nic_output);
    nic_output = NULL;
    nic_output_size = 0;
    return bytes;
}
#endif

static void Put2(Word2 data)
{
    fwrite(&data, (size_t)2, (size_t)1, file);
}

static void Put4(Word4 data)
{
    fwrite(&data, (size_t)4, (size_t)1, file);
}

#define PUT2(data) Put2((Word2)(data))
#define PUT4(data) Put4((Word4)(data))

#define PutIFD(tag, type, count) PUT2(tag); PUT2(type); PUT4(count)
#define PutIFD3(tag, data)       PutIFD(tag, 3, 1); PUT2(data); PUT2(0)
#define PutIFD4(tag, count, data) PutIFD(tag, 4, count); PUT4(data)

#define ENCODE(data) Encode((int)(data))

#ifdef FUNCTION_CELL
#define FETCH(x, y) (*cell)(x, y)
#else
#define FETCH(x, y) cell[y][x]
#endif

void StoreImageFile(char *path, int Nx, int Ny, int BPS,
#ifdef FUNCTION_CELL
                    Cell (*cell)(int, int),
#else
                    Cell **cell,
#endif
                    char *desc)
{
    Word4 *cnt, ofs, pv4, data, ifd;
    Word4 bpr = ((Word4)Nx * (Word4)BPS + 7) / 8;
    Word2 pv2;
    Word2 magic = 0x002a;
    int y, x, bits, idx;
    int len = (desc != NULL) ? (int)strlen(desc) + 1 : 0;
    int rps = (bpr > 8192) ? 1 : (int)(8192 / bpr);
    int spi = (Ny + rps - 1) / rps;
    unsigned char *ptr2 = (unsigned char *)&pv2;
    unsigned char *ptr4 = (unsigned char *)&pv4;

    if ((cnt = (Word4 *)malloc((size_t)spi * (size_t)4)) == NULL) {
        fprintf(stderr, "%s : no allocatable memory.\n", path);
        exit(1);
    }
    if ((file = fopen(path, "wb")) == NULL) {
        fprintf(stderr, "%s : file not stored.\n", path);
        exit(1);
    }

    /* Set larger file buffer for better I/O performance */
    setvbuf(file, NULL, _IOFBF, 1024 * 1024);

    PUT2((*((unsigned char *)&magic) == 0x2a) ? 0x4949 : 0x4d4d);
    PUT2(magic);
    PUT4(0);

    if (len > 4) fwrite(desc, (size_t)1, (size_t)(len + (len & 1)), file);

    ofs = (Word4)ftell(file);

    for (y = 0; y < Ny; y++) {
        if (y % rps == 0) InitEncode();

        if (BPS == 8) {
            for (x = 0; x < Nx; x++) {
                ENCODE(FETCH(x, y));
            }
        } else if (BPS == 16) {
            for (x = 0; x < Nx; x++) {
                pv2 = (Word2)FETCH(x, y);
                ENCODE(ptr2[0]);
                ENCODE(ptr2[1]);
            }
        } else if (BPS == 32) {
            for (x = 0; x < Nx; x++) {
                pv4 = (Word4)FETCH(x, y);
                ENCODE(ptr4[0]);
                ENCODE(ptr4[1]);
                ENCODE(ptr4[2]);
                ENCODE(ptr4[3]);
            }
        } else {
            bits = 0;
            data = 0;
            for (x = 0; x < Nx; x++) {
                pv4 = (Word4)FETCH(x, y);
                if ((bits += BPS) <= 32)
                    data = (data << BPS) | pv4;
                else {
                    bits -= 8;
                    ENCODE(((data << (BPS - bits)) | (pv4 >> bits)) & 0xff);
                    data = pv4;
                }
                while (bits >= 8) {
                    bits -= 8;
                    ENCODE((data >> bits) & 0xff);
                }
            }
            if (bits > 0) {
                ENCODE((data << (8 - bits)) & 0xff);
            }
        }
        if (y + 1 == Ny || (y + 1) % rps == 0) cnt[y / rps] = TermEncode();
    }
    if (ftell(file) & 1) fputc(0, file);

    if (spi > 1) {
        for (idx = 0; idx < spi; idx++) {
            PUT4(ofs);
            ofs += cnt[idx];
        }
        if (ofs & 1) ++ofs;

        fwrite(cnt, (size_t)4, (size_t)spi, file);
    }
    ifd = (Word4)ftell(file);

    PUT2((len > 0) ? 11 : 10);  /* number of tags */

    if ((Word4)Nx > 0xffff) {
        PutIFD4(0x100, 1, Nx);
    } else {
        PutIFD3(0x100, Nx);     /* 0: image width */
    }
    if ((Word4)Ny > 0xffff) {
        PutIFD4(0x101, 1, Ny);
    } else {
        PutIFD3(0x101, Ny);     /* 1: image length */
    }
    PutIFD3(0x102, BPS);        /* 2: bits per sample */
#ifdef LZW
    PutIFD3(0x103, 5);          /* 3: image compression */
#else
    PutIFD3(0x103, 1);          /* 3: image compression */
#endif
    PutIFD3(0x106, (BPS == 1) ? 0 : 1);  /* 4: photometric interpretation */
    if (len) {
        PutIFD(0x10e, 2, len);
        if (len > 4)
            PUT4(8);
        else
            fwrite(desc, (size_t)1, (size_t)4, file);
                                /* 5: image description */
    }
    PutIFD4(0x111, spi, ofs);   /* 5/6: strip offset */
    PutIFD3(0x115, 1);          /* 6/7: samples per pixel */
    PutIFD3(0x116, rps);        /* 7/8: rows per strip */
    PutIFD4(0x117, spi, (spi > 1) ? ofs + (Word4)spi * 4 : *cnt);
                                /* 8/9: strip byte counts */
    PutIFD3(0x11c, 1);          /* 9/10: planar configuration */

    PUT4(0);                    /* end of tags */

    fseek(file, 4L, 0);
    PUT4(ifd);

    fclose(file);
    free(cnt);
}