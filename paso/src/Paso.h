/* $Id$ */


/*
********************************************************************************
*               Copyright © 2006 by ACcESS MNRF                                *
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
#include "Options.h"
#include "SystemMatrix.h"

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
  SYSTEM_ERROR
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
void Paso_solve(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options);
void Paso_solve_free(Paso_SystemMatrix* in);

#endif /* #ifndef INC_PASO */

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
