/* Copyright 2014 The MathWorks, Inc. */

#ifdef SUPPORTS_PRAGMA_ONCE
#pragma once
#endif

#ifndef SL_EXCEL_READER_C_API_H
#define SL_EXCEL_READER_C_API_H

#include "tmwtypes.h"

#ifndef DLL_EXPORT_SYM
#ifdef SL_INTERNAL
#include "package.h"
#else
#define DLL_EXPORT_SYM
#endif
#endif

#ifdef __cplusplus
#define SL_EXCELREADER_EXPORT_EXTERN_C extern "C" DLL_EXPORT_SYM
#else
#define SL_EXCELREADER_EXPORT_EXTERN_C extern DLL_EXPORT_SYM
#endif

SL_EXCELREADER_EXPORT_EXTERN_C CHAR16_T * rtwExcelLoaderGetUnicodeStringFromChars(const char * str);

SL_EXCELREADER_EXPORT_EXTERN_C void rtwExcelLoaderFreeLabel(CHAR16_T * str);

SL_EXCELREADER_EXPORT_EXTERN_C const char *rtwExcelLoaderCreateInstance(const CHAR16_T * fileName,
                                                                        const CHAR16_T * sheetName,
                                                                        const int    extrapolationBeforeFirstDataPointInt,
                                                                        const int    interpolationWithinTimeRangeInt,
                                                                        const int    extrapolationAfterLastDataPointInt,
                                                                        const unsigned char *ground,
                                                                        const int    iZeroCrossingSwitch,
                                                                        const int    signalNumber,
                                                                        const char * spreadsheetIOImpl,
                                                                        void        **outExcelLoader);
SL_EXCELREADER_EXPORT_EXTERN_C const char *terminate(void * readerObj);

SL_EXCELREADER_EXPORT_EXTERN_C const char *GetZeroCrossingSignal(
                                                                 void         *pExcelFileLoader,
                                                                 const double  t,
                                                                 const int     iMajorTimeStep,
                                                                 void         *outZeroCrossingSignal);

SL_EXCELREADER_EXPORT_EXTERN_C const char *getOutput(void ** outOutputValue,
                                                     void * pExcelLoader,
                                                     const double t,
                                                     const int isMajorTimeStep);

#endif