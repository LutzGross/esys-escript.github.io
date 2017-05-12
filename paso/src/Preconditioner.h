
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#ifndef __PASO_PRECONDITIONER_H__
#define __PASO_PRECONDITIONER_H__

#include "Paso.h"
#include "BOOMERAMG.h"
#include "SystemMatrix.h"

namespace paso {

template <class T>
struct MergedSolver;

template <class T>
struct Preconditioner;

template <class T>
using Preconditioner_ptr = boost::shared_ptr<Preconditioner<T> >;

template <class T>
using const_Preconditioner_ptr = boost::shared_ptr<const Preconditioner<T> >;

struct Preconditioner_Smoother;

template<class T>
struct Preconditioner_AMG_Root;
struct Solver_ILU;
struct Solver_RILU;

// general preconditioner interface
template <class T>
struct Preconditioner
{
    dim_t type;
    dim_t sweeps;
    /// Jacobi preconditioner
    Preconditioner_Smoother* jacobi;
    /// Gauss-Seidel preconditioner
    Preconditioner_Smoother* gs;
    /// AMG preconditioner
    Preconditioner_AMG_Root<T>* amg;
    /// ILU preconditioner
    Solver_ILU* ilu;
    /// RILU preconditioner
    Solver_RILU* rilu;
};

template <class T>
void Preconditioner_free(Preconditioner<T>*);

template<class T>
Preconditioner<T>* Preconditioner_alloc(SystemMatrix_ptr<T> A, Options* options);

template<class T>
void Preconditioner_solve(Preconditioner<T>* prec, SystemMatrix_ptr<T> A, double*, double*);


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

void Preconditioner_Smoother_free(Preconditioner_Smoother * in);
void Preconditioner_LocalSmoother_free(Preconditioner_LocalSmoother * in);

template <class T>
Preconditioner_Smoother* Preconditioner_Smoother_alloc(
        SystemMatrix_ptr<T> A, bool jacobi, bool is_local, bool verbose);

Preconditioner_LocalSmoother* Preconditioner_LocalSmoother_alloc(
        SparseMatrix_ptr A, bool jacobi, bool verbose);

template <class T>
void Preconditioner_Smoother_solve(SystemMatrix_ptr<T> A,
        Preconditioner_Smoother* gs, double* x, const double* b,
        dim_t sweeps, bool x_is_initial);

void Preconditioner_LocalSmoother_solve(SparseMatrix_ptr A,
        Preconditioner_LocalSmoother* gs, double* x, const double* b,
        dim_t sweeps, bool x_is_initial);

template <class T>
SolverResult Preconditioner_Smoother_solve_byTolerance(SystemMatrix_ptr<T> A,
                    Preconditioner_Smoother* gs, double* x, const double* b,
                    double atol, dim_t* sweeps, bool x_is_initial);

void Preconditioner_LocalSmoother_Sweep(SparseMatrix_ptr A,
        Preconditioner_LocalSmoother* gs, double* x);

void Preconditioner_LocalSmoother_Sweep_sequential(
        SparseMatrix_ptr A, Preconditioner_LocalSmoother* gs,
        double* x);

void Preconditioner_LocalSmoother_Sweep_tiled(SparseMatrix_ptr A,
        Preconditioner_LocalSmoother* gs, double* x);

void Preconditioner_LocalSmoother_Sweep_colored(SparseMatrix_ptr A,
        Preconditioner_LocalSmoother* gs, double* x);


typedef enum
{
    PASO_AMG_UNDECIDED=-1,
    PASO_AMG_IN_F=0,
    PASO_AMG_IN_C=1
} AMGBlockSelect;

/// Local preconditioner
template <class T>
struct Preconditioner_AMG
{
    int level;
    /// coarse level matrix
    SystemMatrix_ptr<T> A_C;
    /// prolongation n x n_C
    SystemMatrix_ptr<T> P;
    /// restriction  n_C x n
    SystemMatrix_ptr<T> R;

    Preconditioner_Smoother* Smoother;
    int post_sweeps;
    int pre_sweeps;
    /// used in direct solver
    dim_t options_smoother;
    /// used in direct solver
    bool verbose;
    /// applied reordering in direct solver
    index_t reordering;
    /// number of refinements in direct solver (typically =0)
    int refinements;
    /// buffer for residual
    double* r;
    /// solution of coarse level system
    double* x_C;
    /// right hand side of coarse level system
    double* b_C;
    /// used on the coarsest level
    MergedSolver<T>* merged_solver;
    Preconditioner_AMG* AMG_C;
};

template<class T>
Preconditioner_AMG<T>* Preconditioner_AMG_alloc(SystemMatrix_ptr<T> A, int level,
                                             Options* options);

template<class T>
void Preconditioner_AMG_free(Preconditioner_AMG<T>* in);

template<class T>
void Preconditioner_AMG_solve(SystemMatrix_ptr<T> A, Preconditioner_AMG<T>* amg,
                              double* x, double* b);

template<class T>
void Preconditioner_AMG_setStrongConnections(SystemMatrix_ptr<T> A,
                        dim_t* degree_S, index_t* offset_S, index_t* S,
                        double theta, double tau);

template<class T>
void Preconditioner_AMG_setStrongConnections_Block(SystemMatrix_ptr<T> A,
                        dim_t* degree_S, index_t* offset_S, index_t* S,
                        double theta, double tau);

template<class T>
SystemMatrix_ptr<T> Preconditioner_AMG_getProlongation(SystemMatrix_ptr<T> A,
                        const index_t* offset_S, const dim_t* degree_S,
                        const index_t* S, dim_t n_C, index_t* counter_C,
                        index_t interpolation_method);

template<class T>
void Preconditioner_AMG_setClassicProlongation(SystemMatrix_ptr<T> P,
                        SystemMatrix_ptr<T> A, const index_t* offset_S,
                        const dim_t* degree_S, const index_t* S,
                        const index_t* counter_C);

template<class T>
void Preconditioner_AMG_setClassicProlongation_Block(SystemMatrix_ptr<T> P,
                        SystemMatrix_ptr<T> A, const index_t* offset_S,
                        const dim_t* degree_S, const index_t* S,
                        const index_t* counter_C);

template<class T>
void Preconditioner_AMG_setDirectProlongation(SystemMatrix_ptr<T> P,
                        SystemMatrix_ptr<T> A, const index_t* offset_S,
                        const dim_t* degree_S, const index_t* S,
                        const index_t* counter_C);

template<class T>
void Preconditioner_AMG_setDirectProlongation_Block(SystemMatrix_ptr<T> P,
                        SystemMatrix_ptr<T> A, const index_t* offset_S,
                        const dim_t* degree_S, const index_t* S,
                        const index_t* counter_C);

template<class T>
double Preconditioner_AMG_getCoarseLevelSparsity(const Preconditioner_AMG<T>* in);

template<class T>
dim_t Preconditioner_AMG_getNumCoarseUnknowns(const Preconditioner_AMG<T>* in);

template<class T>
int Preconditioner_AMG_getMaxLevel(const Preconditioner_AMG<T>* in);

void Preconditioner_AMG_transposeStrongConnections(dim_t n,
                        const dim_t* degree_S, const index_t* offset_S,
                        const index_t* S, dim_t nT, dim_t* degree_ST,
                        index_t* offset_ST, index_t* ST);

void Preconditioner_AMG_CIJPCoarsening(dim_t n, dim_t my_n,
                        AMGBlockSelect* split_marker, const dim_t* degree_S,
                        const index_t* offset_S, const index_t* S,
                        const dim_t* degree_ST, const index_t* offset_ST,
                        const index_t* ST, const_Connector_ptr col_connector,
                        escript::const_Distribution_ptr col_dist);

template<class T>
SystemMatrix_ptr<T> Preconditioner_AMG_getRestriction(SystemMatrix_ptr<T> P);

template<class T>
SystemMatrix_ptr<T> Preconditioner_AMG_buildInterpolationOperator(
        SystemMatrix_ptr<T> A, SystemMatrix_ptr<T> P, SystemMatrix_ptr<T> R);

template<class T>
SystemMatrix_ptr<T> Preconditioner_AMG_buildInterpolationOperatorBlock(
        SystemMatrix_ptr<T> A, SystemMatrix_ptr<T> P, SystemMatrix_ptr<T> R);

template<class T>
SparseMatrix_ptr Preconditioner_AMG_mergeSystemMatrix(SystemMatrix_ptr<T> A);

template<class T>
void Preconditioner_AMG_mergeSolve(Preconditioner_AMG<T>* amg);

/// Local AMG preconditioner
struct Preconditioner_LocalAMG
{
    dim_t level;
    SparseMatrix_ptr A_C;  // coarse level matrix
    SparseMatrix_ptr P;    // prolongation n x n_C
    SparseMatrix_ptr R;    // restriction  n_C x n

    Preconditioner_LocalSmoother* Smoother;
    int post_sweeps;
    int pre_sweeps;
    index_t reordering; // applied reordering in direct solver
    int refinements;    // number of refinements in direct solver (typically=0)
    double* r;          // buffer for residual
    double* x_C;        // solution of coarse level system
    double* b_C;        // right hand side of coarse level system
    struct Preconditioner_LocalAMG* AMG_C;
};

Preconditioner_LocalAMG* Preconditioner_LocalAMG_alloc(SparseMatrix_ptr A,
                                             int level, Options* options);
void Preconditioner_LocalAMG_free(Preconditioner_LocalAMG* in);
void Preconditioner_LocalAMG_solve(SparseMatrix_ptr A,
                         Preconditioner_LocalAMG* amg, double* x, double* b);

void Preconditioner_LocalAMG_RungeStuebenSearch(dim_t n, const index_t* offset,
                         const dim_t* degree, const index_t* S,
                         AMGBlockSelect* split_marker, bool usePanel);

void Preconditioner_LocalAMG_setStrongConnections_Block(SparseMatrix_ptr A,
                         dim_t* degree, index_t* S, double theta, double tau);

void Preconditioner_LocalAMG_setStrongConnections(SparseMatrix_ptr A,
                         dim_t* degree, index_t* S, double theta, double tau);

SparseMatrix_ptr Preconditioner_LocalAMG_getProlongation(SparseMatrix_ptr A,
                         const index_t* offset_S, const dim_t* degree_S,
                         const index_t* S, dim_t n_C, const index_t* counter_C,
                         index_t interpolation_method);

void Preconditioner_LocalAMG_setDirectProlongation_Block(SparseMatrix_ptr P,
                         const_SparseMatrix_ptr A, const index_t* counter_C);

void Preconditioner_LocalAMG_setDirectProlongation(SparseMatrix_ptr P,
                         const_SparseMatrix_ptr A, const index_t* counter_C);

void Preconditioner_LocalAMG_setClassicProlongation(SparseMatrix_ptr P,
                         SparseMatrix_ptr A, const index_t* offset_S,
                         const dim_t* degree_S, const index_t* S,
                         const index_t* counter_C);

void Preconditioner_LocalAMG_setClassicProlongation_Block(SparseMatrix_ptr P,
                         SparseMatrix_ptr A, const index_t* offset_S,
                         const dim_t* degree_S, const index_t* S,
                         const index_t* counter_C);

int Preconditioner_LocalAMG_getMaxLevel(const Preconditioner_LocalAMG* in);
double Preconditioner_LocalAMG_getCoarseLevelSparsity(const Preconditioner_LocalAMG* in);
dim_t Preconditioner_LocalAMG_getNumCoarseUnknowns(const Preconditioner_LocalAMG* in);
void Preconditioner_LocalAMG_enforceFFConnectivity(dim_t n,
                         const index_t* offset_S, const dim_t* degree_S,
                         const index_t* S, AMGBlockSelect* split_marker);


template <class T>
struct Preconditioner_AMG_Root
{
    bool is_local;
    Preconditioner_AMG<T>* amg;
    Preconditioner_LocalAMG* localamg;
    Preconditioner_BoomerAMG* boomeramg;
    int sweeps;
    Preconditioner_Smoother* amgsubstitute;
};

template <class T>
Preconditioner_AMG_Root<T>* Preconditioner_AMG_Root_alloc(SystemMatrix_ptr<T> A,
                                                       Options* options);

template <class T>
void Preconditioner_AMG_Root_free(Preconditioner_AMG_Root<T>* in);

template <class T>
void Preconditioner_AMG_Root_solve(SystemMatrix_ptr<T> A,
                         Preconditioner_AMG_Root<T>* amg, double* x, double* b);

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
    SparseMatrix_ptr A_FC;
    SparseMatrix_ptr A_CF;
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
Solver_ILU* Solver_getILU(SparseMatrix_ptr A, bool verbose);
void Solver_solveILU(SparseMatrix_ptr A, Solver_ILU* ilu, double* x, const double* b);

void Solver_RILU_free(Solver_RILU* in);
Solver_RILU* Solver_getRILU(SparseMatrix_ptr A, bool verbose);
void Solver_solveRILU(Solver_RILU* rilu, double* x, double* b);

void Solver_updateIncompleteSchurComplement(SparseMatrix_ptr A_CC,
        SparseMatrix_ptr A_CF, double* invA_FF, index_t* A_FF_pivot,
        SparseMatrix_ptr A_FC);

} // namespace paso

#endif // __PASO_PRECONDITIONER_H__

