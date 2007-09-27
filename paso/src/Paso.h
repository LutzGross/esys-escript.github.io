/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

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

double Paso_timer(void);
bool_t Paso_checkPtr(void*);
void Paso_resetError(void);
void Paso_setError(Paso_ErrorCodeType err,char* msg);
bool_t Paso_noError(void);
Paso_ErrorCodeType Paso_getErrorType(void);
char* Paso_getErrorMessage(void);

#endif /* #ifndef INC_PASO */
