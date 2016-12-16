/* Copyright 2013 The MathWorks, Inc. */
#ifndef _IPPGEOTRANS_H_
#define _IPPGEOTRANS_H_

#ifndef EXTERN_C
#  ifdef __cplusplus
#    define EXTERN_C extern "C"
#  else
#    define EXTERN_C extern
#  endif
#endif

#ifndef LIBMWIPPGEOTRANS_API
#    define LIBMWIPPGEOTRANS_API
#endif

#ifdef MATLAB_MEX_FILE
#include "tmwtypes.h"
#else
#include "rtwtypes.h"
#endif

EXTERN_C LIBMWIPPGEOTRANS_API void ippgeotransCaller(real32_T *pDst, real64_T *dstSizeDouble, const real64_T ndims,
                                                     const real32_T *pSrc,  real64_T *srcSize, const real64_T numelSrc,
                                                     const real64_T *tformPtr,  int8_T interpMethodEnum, const real32_T *fillVal);

#endif /* _IPPGEOTRANS_H_ */
