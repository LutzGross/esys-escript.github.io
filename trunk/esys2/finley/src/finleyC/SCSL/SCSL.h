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
/*
* $Log$
* Revision 1.2  2004/12/14 05:39:31  jgs
* *** empty log message ***
*
* Revision 1.1.1.1.2.1  2004/11/12 06:58:20  gross
* a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1  2004/07/02 04:21:14  gross
* Finley C code has been included
*
*
*/
