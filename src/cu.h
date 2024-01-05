
#include <stdio.h>
#include <stdlib.h>

/*
	macros originally written in "*_SDK/C/common/inc/cutil.h"
*/

#define	CUDA_SAFE_CALL_NO_SYNC(call)	\
{\
	cudaError	err=call;\
\
	if (cudaSuccess!=err) {\
	    (void)fprintf(stderr,\
			  "Cuda error in file '%s' in line %i : %s.\n",\
			  __FILE__,__LINE__,cudaGetErrorString(err));\
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
			  "CUFFT error in file '%s' in line %i.\n",\
			  __FILE__,__LINE__);\
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
	    exit(EXIT_FAILURE);\
	}\
	err=CUT_DEVICE_SYNCHRONIZE();\
	if (cudaSuccess!=err) {\
	    (void)fprintf(stderr,\
			  "Cuda error: %s in file '%s' in line %i : %s.\n",\
			  msg,__FILE__,__LINE__,cudaGetErrorString(err));\
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
	char		*ev=getenv("CUDA_GPU"),\
			*nd="no device assigned by CUDA_GPU.\n";\
	int		cg,dc;\
        cudaDeviceProp	dp;\
\
	if (ev==NULL)\
	    cg=CUDA_GPU;\
	else\
	    if (sscanf(ev,"%d",&cg)!=1) (void)fputs(nd,stderr);\
\
	if (cg<0) (void)fputs(nd,stderr);\
\
	CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&dc));\
	if (cg>=dc) (void)fputs(nd,stderr);\
\
	CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&dp,cg));\
	if (dp.major<1) (void)fputs("device not supported by CUDA.\n",stderr);\
\
	CUDA_SAFE_CALL(cudaSetDevice(cg));\
}
#endif
