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

void Finley_SCSL_free(Finley_SystemMatrix* A);
void Finley_SCSL(Finley_SystemMatrix* A, double* out, double* in, Finley_SolverOptions* options);
void Finley_SCSL_iterative_free(Finley_SystemMatrix* A);
void Finley_SCSL_iterative(Finley_SystemMatrix* A, double* out,double* in,Finley_SolverOptions* options);
void Finley_SCSL_direct_free(Finley_SystemMatrix* A);
void Finley_SCSL_direct(Finley_SystemMatrix* A, double* out, double* in, Finley_SolverOptions* options);

#endif
