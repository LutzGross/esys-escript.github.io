/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrix: interface to SGI SCSL sparse solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_SCSL
#define INC_SCSL

#include "System.h"

void Finley_SCSL_iterative_free(Finley_SystemMatrix* A);
void Finley_SCSL_iterative(Finley_SystemMatrix* A, double* out,double* in,Finley_SolverOptions* options);
void Finley_SCSL_solve_free(Finley_SystemMatrix* A);
void Finley_SCSL_solve(Finley_SystemMatrix* A, double* out, double* in, Finley_SolverOptions* options);

#endif
/*
* $Log$
* Revision 1.3  2004/12/15 03:48:47  jgs
* *** empty log message ***
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1  2004/07/02 04:21:14  gross
* Finley C code has been included
*
*
*/
