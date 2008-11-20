
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    Paso finite element solver */

/**************************************************************/

/*  Copyrights by ACcESS Australia, 2003,2004,2005 */
/*  Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO
#define INC_PASO

#include "Common.h"

/**************************************************************/

enum Paso_ErrorCodeType {
  NO_ERROR,
  WARNING,
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
  PASO_MPI_ERROR 
};

typedef enum Paso_ErrorCodeType Paso_ErrorCodeType;

/* interfaces */


PASO_DLL_API
double Paso_timer(void);

PASO_DLL_API
bool_t Paso_checkPtr(void*);

PASO_DLL_API
void Paso_resetError(void);

PASO_DLL_API
void Paso_setError(Paso_ErrorCodeType err,__const char* msg);

PASO_DLL_API
bool_t Paso_noError(void);

PASO_DLL_API
Paso_ErrorCodeType Paso_getErrorType(void);

PASO_DLL_API
char* Paso_getErrorMessage(void);

#ifndef _OPENMP 
int omp_get_max_threads(void);
#endif

#endif /* #ifndef INC_PASO */
