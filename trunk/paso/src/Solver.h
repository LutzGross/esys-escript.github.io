
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#ifndef INC_SOLVER
#define INC_SOLVER

#include "SystemMatrix.h"
#include "performance.h"
#include "Functions.h"

#define PASO_TRACE
/* error codes used in the solver */
#define SOLVER_NO_ERROR 0
#define SOLVER_MAXITER_REACHED 1
#define SOLVER_INPUT_ERROR -1
#define SOLVER_MEMORY_ERROR -9
#define SOLVER_BREAKDOWN -10
#define SOLVER_NEGATIVE_NORM_ERROR -11
#define SOLVER_DIVERGENCE -12

#define TOLERANCE_FOR_SCALARS (double)(0.)

PASO_DLL_API
void Paso_Solver(Paso_SystemMatrix*,double*,double*,Paso_Options*,Paso_Performance* pp);

PASO_DLL_API
void Paso_Solver_free(Paso_SystemMatrix*);

err_t Paso_Solver_BiCGStab( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_PCG( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_TFQMR( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_MINRES( Paso_SystemMatrix * A, double* B, double * X, dim_t *iter, double * tolerance, Paso_Performance* pp);
err_t Paso_Solver_GMRES(Paso_SystemMatrix * A, double * r, double * x, dim_t *num_iter, double * tolerance,dim_t length_of_recursion,dim_t restart, Paso_Performance* pp);
err_t Paso_Solver_GMRES2(Paso_Function * F, const double* f0, const double* x0, double * x, dim_t *iter, double* tolerance, Paso_Performance* pp);
err_t Paso_Solver_NewtonGMRES(Paso_Function *F, double *x, Paso_Options* options, Paso_Performance* pp);

Paso_Function * Paso_Function_LinearSystem_alloc(Paso_SystemMatrix* A, double* b, Paso_Options* options);
err_t Paso_Function_LinearSystem_call(Paso_Function * F,double* value, const double* arg, Paso_Performance *pp);
void Paso_Function_LinearSystem_free(Paso_Function * F);
err_t Paso_Function_LinearSystem_setInitialGuess(Paso_SystemMatrix* A, double* x, Paso_Performance *pp);

#endif /* #ifndef INC_SOLVER */

