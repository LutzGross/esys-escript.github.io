
/*****************************************************************************
*
* Copyright (c) 2010-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************/

/*    Functions for C error handling  (and timing)*/

/************************************************************************************/


#ifndef INC_ESYS_ERROR
#define INC_ESYS_ERROR

#include "system_dep.h"
#include "types.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>	/* For FILENAME_MAX */
#define LenString_MAX FILENAME_MAX*2
#define LenErrorMsg_MAX LenString_MAX

/************************************************************************************/

typedef enum {
  NO_ERROR,
  WARNING,
  DIVERGED,
  VALUE_ERROR,
  TYPE_ERROR,
  MEMORY_ERROR,
  IO_ERROR,
  ZERO_DIVISION_ERROR,
  EOF_ERROR,
  FLOATING_POINT_ERROR,
  INDEX_ERROR,
  OS_ERROR,
  OVERFLOW_ERROR,
  SYSTEM_ERROR,
  ESYS_MPI_ERROR,
  NO_PROGRESS_ERROR
} Esys_ErrorCodeType;

/* interfaces */


ESYSUTILS_DLL_API
double Esys_timer(void);

ESYSUTILS_DLL_API
bool Esys_checkPtr(void*);

ESYSUTILS_DLL_API
void Esys_resetError(void);

ESYSUTILS_DLL_API
void Esys_setError(Esys_ErrorCodeType err,__const char* msg);

ESYSUTILS_DLL_API
bool Esys_noError(void);

ESYSUTILS_DLL_API
Esys_ErrorCodeType Esys_getErrorType(void);

ESYSUTILS_DLL_API
char* Esys_getErrorMessage(void);

#ifndef _OPENMP
int serial_get_max_threads(void);
int serial_get_thread_num(void);

/* Nasty hack to get 3.2 out */
#define omp_get_max_threads serial_get_max_threads
#define omp_get_thread_num serial_get_thread_num
#endif


#endif /* #ifndef INC_ESYS_ERROR */
