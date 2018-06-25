
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#ifndef __PASO_SOLVER_H__
#define __PASO_SOLVER_H__

#include "Paso.h"
#include "Functions.h"
#include "performance.h"
#include "SystemMatrix.h"

namespace paso {

#define TOLERANCE_FOR_SCALARS (double)(0.)

void solve_free(SystemMatrix* A);

SolverResult Solver(SystemMatrix_ptr, double*, double*, Options*, Performance*);

void Solver_free(SystemMatrix*);

SolverResult Solver_BiCGStab(SystemMatrix_ptr A, double* B, double* X,
                             dim_t* iter, double* tolerance, Performance* pp);

SolverResult Solver_PCG(SystemMatrix_ptr A, double* B, double* X, dim_t* iter,
                        double* tolerance, Performance* pp);

SolverResult Solver_TFQMR(SystemMatrix_ptr A, double* B, double* X, dim_t* iter,
                          double* tolerance, Performance* pp);

SolverResult Solver_MINRES(SystemMatrix_ptr A, double* B, double* X,
                           dim_t* iter, double* tolerance, Performance* pp);

SolverResult Solver_GMRES(SystemMatrix_ptr A, double* r, double* x,
                          dim_t* num_iter, double* tolerance,
                          dim_t length_of_recursion, dim_t restart,
                          Performance* pp);

SolverResult Solver_GMRES2(Function* F, const double* f0, const double* x0,
                           double* x, dim_t* iter, double* tolerance,
                           Performance* pp);

SolverResult Solver_NewtonGMRES(Function* F, double* x, Options* options,
                                Performance* pp);

} // namespace paso

#endif // __PASO_SOLVER_H__

