/* $Id: SCSL.h 150 2005-09-15 03:44:45Z jgs $ */

/**************************************************************/

/* Paso: SystemMatrix: interface to SGI SCSL sparse solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SCSL
#define INC_PASO_SCSL

#include "SystemMatrix.h"

void Paso_SCSL_free(Paso_SystemMatrix* A);
void Paso_SCSL(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options);
void Paso_SCSL_iterative_free(Paso_SystemMatrix* A);
void Paso_SCSL_iterative(Paso_SystemMatrix* A, double* out,double* in,Paso_Options* options);
void Paso_SCSL_direct_free(Paso_SystemMatrix* A);
void Paso_SCSL_direct(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options);

#endif
