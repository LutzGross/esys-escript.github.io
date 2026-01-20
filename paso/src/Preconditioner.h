
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __PASO_PRECONDITIONER_H__
#define __PASO_PRECONDITIONER_H__

#include "Paso.h"
#include "SystemMatrix.h"

namespace paso {

struct MergedSolver;
struct Preconditioner;
typedef boost::shared_ptr<Preconditioner> Preconditioner_ptr;
typedef boost::shared_ptr<const Preconditioner> const_Preconditioner_ptr;

struct Preconditioner_Smoother;
struct Solver_ILU;
struct Solver_RILU;

// general preconditioner interface
struct Preconditioner
{
    dim_t type;
    dim_t sweeps;
    /// Jacobi preconditioner
    Preconditioner_Smoother* jacobi;
    /// Gauss-Seidel preconditioner
    Preconditioner_Smoother* gs;
    /// ILU preconditioner
    Solver_ILU* ilu;
    /// RILU preconditioner
    Solver_RILU* rilu;
};

void Preconditioner_free(Preconditioner*);
Preconditioner* Preconditioner_alloc(SystemMatrix_ptr<double> A, Options* options);
void Preconditioner_solve(Preconditioner* prec, SystemMatrix_ptr<double> A, double*, double*);


// GAUSS SEIDEL & Jacobi
struct Preconditioner_LocalSmoother
{
    bool Jacobi;
    double* diag;
    double* buffer;
    index_t* pivot;
};

struct Preconditioner_Smoother
{
    Preconditioner_LocalSmoother* localSmoother;
    bool is_local;
};

void PASO_DLL_API Preconditioner_Smoother_free(Preconditioner_Smoother * in);
void PASO_DLL_API Preconditioner_LocalSmoother_free(Preconditioner_LocalSmoother * in);

Preconditioner_Smoother* Preconditioner_Smoother_alloc(
        SystemMatrix_ptr<double> A, bool jacobi, bool is_local, bool verbose);

Preconditioner_LocalSmoother* Preconditioner_LocalSmoother_alloc(
        SparseMatrix_ptr<double> A, bool jacobi, bool verbose);

void Preconditioner_Smoother_solve(SystemMatrix_ptr<double> A,
        Preconditioner_Smoother* gs, double* x, const double* b,
        dim_t sweeps, bool x_is_initial);

void Preconditioner_LocalSmoother_solve(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* gs, double* x, const double* b,
        dim_t sweeps, bool x_is_initial);

SolverResult Preconditioner_Smoother_solve_byTolerance(SystemMatrix_ptr<double> A,
                    Preconditioner_Smoother* gs, double* x, const double* b,
                    double atol, dim_t* sweeps, bool x_is_initial);

void Preconditioner_LocalSmoother_Sweep(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* gs, double* x);

void Preconditioner_LocalSmoother_Sweep_sequential(
        SparseMatrix_ptr<double> A, Preconditioner_LocalSmoother* gs,
        double* x);

void Preconditioner_LocalSmoother_Sweep_tiled(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* gs, double* x);

void Preconditioner_LocalSmoother_Sweep_colored(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* gs, double* x);

/// ILU preconditioner
struct Solver_ILU
{
    double* factors;
};

/// RILU preconditioner
struct Solver_RILU
{
    dim_t n;
    dim_t n_block;
    dim_t n_F;
    dim_t n_C;
    double* inv_A_FF;
    index_t* A_FF_pivot;
    SparseMatrix_ptr<double> A_FC;
    SparseMatrix_ptr<double> A_CF;
    index_t* rows_in_F;
    index_t* rows_in_C;
    index_t* mask_F;
    index_t* mask_C;
    double* x_F;
    double* b_F;
    double* x_C;
    double* b_C;
    Solver_RILU* RILU_of_Schur;
};

void Solver_ILU_free(Solver_ILU * in);
Solver_ILU* Solver_getILU(SparseMatrix_ptr<double> A, bool verbose);
void Solver_solveILU(SparseMatrix_ptr<double> A, Solver_ILU* ilu, double* x, const double* b);

void Solver_RILU_free(Solver_RILU* in);
Solver_RILU* Solver_getRILU(SparseMatrix_ptr<double> A, bool verbose);
void Solver_solveRILU(Solver_RILU* rilu, double* x, double* b);

void Solver_updateIncompleteSchurComplement(SparseMatrix_ptr<double> A_CC,
        SparseMatrix_ptr<double> A_CF, double* invA_FF, index_t* A_FF_pivot,
        SparseMatrix_ptr<double> A_FC);

} // namespace paso

#endif // __PASO_PRECONDITIONER_H__

