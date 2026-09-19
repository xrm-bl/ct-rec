#include <stdio.h>
#include <stdlib.h>

/*
	macros originally written in "*_SDK/C/common/inc/cutil.h"
*/

/*
	extra diagnostics printed when the GPU runs out of memory
*/

static void	CudaOOMInfo(void)
{
	size_t	freeB=0,totalB=0;

	if (cudaMemGetInfo(&freeB,&totalB)==cudaSuccess)
	    (void)fprintf(stderr,
			  "  GPU memory: free %.1f MiB / total %.1f MiB\n",
			  (double)freeB/1048576.0,(double)totalB/1048576.0);
	(void)fputs("  hint: close other GPU programs, reduce the data size, "
		    "or use the CPU variant.\n",stderr);
}

/*
	cufft.h may not be included by every user of this header,
	so take the error code as a plain int
*/

static const char	*CufftErrorName(int err)
{
	switch (err) {
	case  1: return "CUFFT_INVALID_PLAN";
	case  2: return "CUFFT_ALLOC_FAILED";
	case  3: return "CUFFT_INVALID_TYPE";
	case  4: return "CUFFT_INVALID_VALUE";
	case  5: return "CUFFT_INTERNAL_ERROR";
	case  6: return "CUFFT_EXEC_FAILED";
	case  7: return "CUFFT_SETUP_FAILED";
	case  8: return "CUFFT_INVALID_SIZE";
	case  9: return "CUFFT_UNALIGNED_DATA";
	case 10: return "CUFFT_INCOMPLETE_PARAMETER_LIST";
	case 11: return "CUFFT_INVALID_DEVICE";
	case 12: return "CUFFT_PARSE_ERROR";
	case 13: return "CUFFT_NO_WORKSPACE";
	case 14: return "CUFFT_NOT_IMPLEMENTED";
	case 15: return "CUFFT_LICENSE_ERROR";
	case 16: return "CUFFT_NOT_SUPPORTED";
	}
	return "unknown cufftResult";
}

#define	CUDA_SAFE_CALL_NO_SYNC(call)	\
{\
	cudaError	err=call;\
\
	if (cudaSuccess!=err) {\
	    (void)fprintf(stderr,\
			  "Cuda error in file '%s' in line %i : %s.\n",\
			  __FILE__,__LINE__,cudaGetErrorString(err));\
	    if (cudaErrorMemoryAllocation==err) CudaOOMInfo();\
	    exit(EXIT_FAILURE);\
	}\
}

#define	CUDA_SAFE_CALL(call)	CUDA_SAFE_CALL_NO_SYNC(call)

#define	CUFFT_SAFE_CALL(call)	\
{\
	cufftResult	err=call;\
\
	if (CUFFT_SUCCESS!=err) {\
	    (void)fprintf(stderr,\
			  "CUFFT error in file '%s' in line %i : %s (%d).\n",\
			  __FILE__,__LINE__,\
			  CufftErrorName((int)err),(int)err);\
	    /* CUDA 13 reports a workspace allocation failure as\
	       CUFFT_INTERNAL_ERROR, so show the memory state for both */\
	    if (CUFFT_ALLOC_FAILED==err ||\
		CUFFT_INTERNAL_ERROR==err) CudaOOMInfo();\
	    exit(EXIT_FAILURE);\
	}\
}

#ifdef	_DEBUG
#if	CUDART_VERSION>=4000
#define	CUT_DEVICE_SYNCHRONIZE()	cudaDeviceSynchronize()
#else
#define	CUT_DEVICE_SYNCHRONIZE()	cudaThreadSynchronize()
#endif

#define	CUT_CHECK_ERROR(msg)	\
{\
	cudaError_t	err=cudaGetLastError();\
\
	if (cudaSuccess!=err) {\
	    (void)fprintf(stderr,\
			  "Cuda error: %s in file '%s' in line %i : %s.\n",\
			  msg,__FILE__,__LINE__,cudaGetErrorString(err));\
	    if (cudaErrorMemoryAllocation==err) CudaOOMInfo();\
	    exit(EXIT_FAILURE);\
	}\
	err=CUT_DEVICE_SYNCHRONIZE();\
	if (cudaSuccess!=err) {\
	    (void)fprintf(stderr,\
			  "Cuda error: %s in file '%s' in line %i : %s.\n",\
			  msg,__FILE__,__LINE__,cudaGetErrorString(err));\
	    if (cudaErrorMemoryAllocation==err) CudaOOMInfo();\
	    exit(EXIT_FAILURE);\
	}\
}
#else
#define	CUT_CHECK_ERROR(msg)	\
{\
	cudaError_t	err=cudaGetLastError();\
\
	if (cudaSuccess!=err) {\
	    (void)fprintf(stderr,\
			  "Cuda error: %s in file '%s' in line %i : %s.\n",\
			  msg,__FILE__,__LINE__,cudaGetErrorString(err));\
	    if (cudaErrorMemoryAllocation==err) CudaOOMInfo();\
	    exit(EXIT_FAILURE);\
	}\
}
#endif

/*
	set up device assigned by the environmental variable "CUDA_GPU"
*/

#ifdef	__DEVICE_EMULATION__
#define SETUP_CUDA_GPU()	/* nothing */
#else
#ifndef	CUDA_GPU
#define	CUDA_GPU	0
#endif

#define	SETUP_CUDA_GPU()	\
{\
	char		*ev=getenv("CUDA_GPU");\
	const char	*nd="no device assigned by CUDA_GPU.\n";\
	int		cg,dc;\
        cudaDeviceProp	dp;\
\
	if (ev==NULL)\
	    cg=CUDA_GPU;\
	else\
	    if (sscanf(ev,"%d",&cg)!=1) {\
		(void)fputs(nd,stderr); exit(EXIT_FAILURE);\
	    }\
\
	if (cg<0) {\
	    (void)fputs(nd,stderr); exit(EXIT_FAILURE);\
	}\
\
	CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&dc));\
	if (cg>=dc) {\
	    (void)fputs(nd,stderr); exit(EXIT_FAILURE);\
	}\
\
	CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&dp,cg));\
	if (dp.major<1) {\
	    (void)fputs("device not supported by CUDA.\n",stderr);\
	    exit(EXIT_FAILURE);\
	}\
\
	CUDA_SAFE_CALL(cudaSetDevice(cg));\
}
#endif
