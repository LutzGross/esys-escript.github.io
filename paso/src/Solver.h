
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#ifndef __PASO_SOLVER_H__
#define __PASO_SOLVER_H__

#include "Paso.h"
#include "Functions.h"
#include "performance.h"
#include "SystemMatrix.h"

namespace paso {

#define TOLERANCE_FOR_SCALARS (double)(0.)

struct Function;

template <typename T>
void solve_free(SystemMatrix<T>* A);

SolverResult Solver(SystemMatrix_ptr<double>, double*, double*, Options*, Performance*);
void PASO_DLL_API Solver_free(SystemMatrix<double>*);

SolverResult Solver(SystemMatrix_ptr<cplx_t>, cplx_t*, cplx_t*, Options*, Performance*);
void PASO_DLL_API Solver_free(SystemMatrix<cplx_t>*);

SolverResult Solver_BiCGStab(SystemMatrix_ptr<double> A, double* B, double* X,
                             dim_t* iter, double* tolerance, Performance* pp);

SolverResult Solver_PCG(SystemMatrix_ptr<double> A, double* B, double* X, dim_t* iter,
                        double* tolerance, Performance* pp);

SolverResult Solver_TFQMR(SystemMatrix_ptr<double> A, double* B, double* X, dim_t* iter,
                          double* tolerance, Performance* pp);

SolverResult Solver_MINRES(SystemMatrix_ptr<double> A, double* B, double* X,
                           dim_t* iter, double* tolerance, Performance* pp);

SolverResult Solver_GMRES(SystemMatrix_ptr<double> A, double* r, double* x,
                          dim_t* num_iter, double* tolerance,
                          dim_t length_of_recursion, dim_t restart,
                          Performance* pp);

SolverResult Solver_GMRES2(Function* F, const double* f0, const double* x0,
                           double* x, dim_t* iter, double* tolerance,
                           Performance* pp);

SolverResult Solver_NewtonGMRES(Function* F, double* x, Options* options,
                                Performance* pp);

} // namespace paso

#include "Preconditioner.h"
#include "MKL.h"
#include "UMFPACK.h"
#include "MUMPS.h"

namespace paso {

struct Preconditioner_Smoother;
void PASO_DLL_API Preconditioner_Smoother_free(Preconditioner_Smoother * in);

template <typename T>
void solve_free(SystemMatrix<T>* in)
{
    if (!in) return;

    switch(in->solver_package) {
        case PASO_PASO:
            Solver_free(in);
            break;

        case PASO_SMOOTHER:
            Preconditioner_Smoother_free((Preconditioner_Smoother*) in->solver_p);
            break;

        case PASO_MKL:
            MKL_free(in->mainBlock.get());
            break;

        case PASO_UMFPACK:
            UMFPACK_free(in->mainBlock.get());
            break;

        case PASO_MUMPS:
            MUMPS_free(in->mainBlock.get());
            break;
   }
}

} // namespace paso

#endif // __PASO_SOLVER_H__

