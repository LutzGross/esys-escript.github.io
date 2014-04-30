
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#ifndef __PASO_SOLVER_H__
#define __PASO_SOLVER_H__

#include "SystemMatrix.h"
#include "performance.h"
#include "Functions.h"

namespace paso {

// error codes used in the solver
#define SOLVER_NO_ERROR 0
#define SOLVER_MAXITER_REACHED 1
#define SOLVER_INPUT_ERROR -1
#define SOLVER_MEMORY_ERROR -9
#define SOLVER_BREAKDOWN -10
#define SOLVER_NEGATIVE_NORM_ERROR -11
#define SOLVER_DIVERGENCE -12

#define TOLERANCE_FOR_SCALARS (double)(0.)

void solve(SystemMatrix_ptr A, double* out, double* in, Options* options);

void solve_free(SystemMatrix* A);

PASO_DLL_API
void Solver(SystemMatrix_ptr, double*, double*, Options*, Performance*);

PASO_DLL_API
void Solver_free(SystemMatrix*);

err_t Solver_BiCGStab(SystemMatrix_ptr A, double* B, double* X, dim_t* iter,
                      double* tolerance, Performance* pp);

err_t Solver_PCG(SystemMatrix_ptr A, double* B, double* X, dim_t* iter,
                 double* tolerance, Performance* pp);

err_t Solver_TFQMR(SystemMatrix_ptr A, double* B, double* X, dim_t* iter,
                   double* tolerance, Performance* pp);

err_t Solver_MINRES(SystemMatrix_ptr A, double* B, double* X, dim_t* iter,
                    double* tolerance, Performance* pp);

err_t Solver_GMRES(SystemMatrix_ptr A, double* r, double* x, dim_t* num_iter,
                   double* tolerance, dim_t length_of_recursion, dim_t restart,
                   Performance* pp);

err_t Solver_GMRES2(Function* F, const double* f0, const double* x0, double* x,
                    dim_t* iter, double* tolerance, Performance* pp);

err_t Solver_NewtonGMRES(Function* F, double* x, Options* options,
                         Performance* pp);

} // namespace paso

#endif // __PASO_SOLVER_H__

