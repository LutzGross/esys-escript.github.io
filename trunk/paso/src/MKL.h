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

/* Paso: SystemMatrix: interface to intel MKL sparse solver */

/**************************************************************/

/* Copyrights by ACcESS Australia 2006 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_MKL
#define INC_PASO_MKL

#include "SystemMatrix.h"
#include "performance.h"

# if defined(_WIN32) || defined(_WIN64)
#define PARDISO pardiso
#else
#define PARDISO pardiso_
#endif

#ifdef MKL
#include "mkl_solver.h"
#endif


#define MKL_ERROR_NO 0
#define MKL_MTYPE_SYM -2
#define MKL_MTYPE_UNSYM 11

#define MKL_REORDERING_MINIMUM_DEGREE 0
#define MKL_REORDERING_NESTED_DISSECTION 2
#define MKL_PHASE_SYMBOLIC_FACTORIZATION 11
#define MKL_PHASE_FACTORIZATION 22
#define MKL_PHASE_SOLVE 33
#define MKL_PHASE_RELEASE_MEMORY -1

/* extern int PARDISO
#         (void *, int *, int *, int *, int *, int *,
#         double *, int *, int *, int *, int *, int *,
#         int *, double *, double *, int *);
*/


void Paso_MKL_free(Paso_SystemMatrix* A);
void Paso_MKL(Paso_SystemMatrix* A, double* out, double* in, Paso_Options* options,Paso_Performance* pp);
#endif
