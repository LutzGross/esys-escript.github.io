
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

#ifndef __PASO_PRECONDITIONER_H__
#define __PASO_PRECONDITIONER_H__

#include "SystemMatrix.h"
#include "BOOMERAMG.h"
#include "MergedSolver.h"

#define PRECONDITIONER_NO_ERROR 0
#define PRECONDITIONER_MAXITER_REACHED 1
#define PRECONDITIONER_INPUT_ERROR -1
#define PRECONDITIONER_MEMORY_ERROR -9
#define PRECONDITIONER_BREAKDOWN -10
#define PRECONDITIONER_NEGATIVE_NORM_ERROR -11
#define PRECONDITIONER_DIVERGENCE -12

namespace paso {

struct Preconditioner;
typedef boost::shared_ptr<Preconditioner> Preconditioner_ptr;
typedef boost::shared_ptr<const Preconditioner> const_Preconditioner_ptr;

struct Preconditioner_Smoother;
struct Preconditioner_AMG_Root;
struct Solver_ILU;
struct Solver_RILU;

// general preconditioner interface
struct Preconditioner
{
    dim_t type;
    dim_t sweeps;
    /* jacobi preconditioner */
    Preconditioner_Smoother* jacobi;
    /* Gauss-Seidel preconditioner */
    Preconditioner_Smoother* gs;  
    /* amg preconditioner */
    Preconditioner_AMG_Root *amg;

    /* ilu preconditioner */
    Solver_ILU* ilu;
    /* rilu preconditioner */
    Solver_RILU* rilu;
};

void Preconditioner_free(Preconditioner*);
Preconditioner* Preconditioner_alloc(SystemMatrix_ptr A, Options* options);
void Preconditioner_solve(Preconditioner* prec, SystemMatrix_ptr A, double*, double*);


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

void Paso_solve(SystemMatrix_ptr A, double* out, double* in, Options* options);

void Paso_solve_free(SystemMatrix* A);

Preconditioner_Smoother* Preconditioner_Smoother_alloc(
        SystemMatrix_ptr A, bool jacobi, bool is_local, bool verbose);

Preconditioner_LocalSmoother* Preconditioner_LocalSmoother_alloc(
        SparseMatrix_ptr A, bool jacobi, bool verbose);

void Preconditioner_Smoother_solve(SystemMatrix_ptr A,
        Preconditioner_Smoother* gs, double* x, const double* b,
        dim_t sweeps, bool x_is_initial);

void Preconditioner_LocalSmoother_solve(SparseMatrix_ptr A,
        Preconditioner_LocalSmoother* gs, double* x, const double* b,
        dim_t sweeps, bool x_is_initial);

err_t Preconditioner_Smoother_solve_byTolerance(SystemMatrix_ptr A,
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


/* Local preconditioner */
struct Preconditioner_AMG
{
    dim_t level;
    SystemMatrix_ptr A_C; /* coarse level matrix */
    SystemMatrix_ptr P;   /* prolongation n x n_C*/ 
    SystemMatrix_ptr R;   /* restriction  n_C x n */
   
    Preconditioner_Smoother* Smoother;
    dim_t post_sweeps;
    dim_t pre_sweeps;
    dim_t options_smoother;  /* used in direct solver */
    bool verbose;            /* used in direct solver */
    index_t reordering;  /* applied reordering in direct solver */
    dim_t refinements;  /* number of refinements in direct solver (typically =0) */
    double* r;         /* buffer for residual */
    double* x_C;       /* solution of coarse level system */
    double* b_C;       /* right hand side of coarse level system */
    Paso_MergedSolver* merged_solver; /* used on the coarsest level */
    Preconditioner_AMG* AMG_C;
};

void Preconditioner_AMG_free(Preconditioner_AMG * in);
Preconditioner_AMG* Preconditioner_AMG_alloc(SystemMatrix_ptr A, dim_t level, Options* options);
void Preconditioner_AMG_solve(SystemMatrix_ptr A, Preconditioner_AMG * amg, double * x, double * b);
void Preconditioner_AMG_setStrongConnections(SystemMatrix_ptr A,  dim_t *degree_S, index_t* offset_S, index_t *S, const double theta, const double tau);
void Preconditioner_AMG_setStrongConnections_Block(SystemMatrix_ptr A, dim_t *degree_S, index_t* offset_S, index_t *S, const double theta, const double tau);
SystemMatrix_ptr Preconditioner_AMG_getProlongation(SystemMatrix_ptr A_p, const index_t* offset_S, const dim_t* degree_S, const index_t* S, const dim_t n_C, index_t* counter_C, const index_t interpolation_method);
void Preconditioner_AMG_setClassicProlongation(SystemMatrix_ptr P, SystemMatrix_ptr A, const index_t* offset_S, const dim_t* degree_S, const index_t* S,const index_t *counter_C);
void Preconditioner_AMG_setClassicProlongation_Block(SystemMatrix_ptr P, SystemMatrix_ptr A, const index_t* offset_S, const dim_t* degree_S, const index_t* S,const index_t *counter_C);
void Preconditioner_AMG_setDirectProlongation(SystemMatrix_ptr P, SystemMatrix_ptr A, const index_t* offset_S, const dim_t* degree_S, const index_t* S,const index_t *counter_C);
void Preconditioner_AMG_setDirectProlongation_Block(SystemMatrix_ptr P, SystemMatrix_ptr A, const index_t* offset_S, const dim_t* degree_S, const index_t* S,const index_t *counter_C);
double Preconditioner_AMG_getCoarseLevelSparsity(const Preconditioner_AMG * in);
dim_t Preconditioner_AMG_getNumCoarseUnknowns(const Preconditioner_AMG * in);
index_t Preconditioner_AMG_getMaxLevel(const Preconditioner_AMG * in);
void Preconditioner_AMG_transposeStrongConnections(const dim_t n, const dim_t* degree_S, const index_t* offset_S, const index_t* S, const dim_t nT, dim_t *degree_ST, index_t* offset_ST,index_t* ST);
void Preconditioner_AMG_CIJPCoarsening(const dim_t n, const dim_t my_n, AMGBlockSelect*split_marker,
                                            const dim_t* degree_S, const index_t* offset_S, const index_t* S,
                                            const dim_t* degree_ST, const index_t* offset_ST, const index_t* ST,
                                            Connector_ptr col_connector, const_Distribution_ptr col_dist);
SystemMatrix_ptr Preconditioner_AMG_getRestriction(SystemMatrix_ptr P);
SystemMatrix_ptr Preconditioner_AMG_buildInterpolationOperator(
        SystemMatrix_ptr A, SystemMatrix_ptr P, SystemMatrix_ptr R);

SystemMatrix_ptr Preconditioner_AMG_buildInterpolationOperatorBlock(
        SystemMatrix_ptr A, SystemMatrix_ptr P, SystemMatrix_ptr R);

SparseMatrix_ptr Preconditioner_AMG_mergeSystemMatrix(SystemMatrix_ptr A);

void Preconditioner_AMG_mergeSolve(Preconditioner_AMG* amg);

/* Local AMG preconditioner */
struct Preconditioner_LocalAMG {
   dim_t level;
   SparseMatrix_ptr A_C;  /* coarse level matrix */
   SparseMatrix_ptr P;    /* prolongation n x n_C*/ 
   SparseMatrix_ptr R;    /* restriction  n_C x n */

   Preconditioner_LocalSmoother* Smoother;
   dim_t post_sweeps;
   dim_t pre_sweeps;
   index_t reordering;  /* applied reordering in direct solver */
   dim_t refinements;  /* number of refinements in direct solver (typically =0) */
   double* r;         /* buffer for residual */
   double* x_C;       /* solution of coarse level system */
   double* b_C;       /* right hand side of coarse level system */
   struct Preconditioner_LocalAMG * AMG_C;
};

void Preconditioner_LocalAMG_free(Preconditioner_LocalAMG * in);
Preconditioner_LocalAMG* Preconditioner_LocalAMG_alloc(SparseMatrix_ptr A, dim_t level, Options* options);
void Preconditioner_LocalAMG_solve(SparseMatrix_ptr A, Preconditioner_LocalAMG * amg, double * x, double * b);

void Preconditioner_LocalAMG_RungeStuebenSearch(const dim_t n, const index_t* offset, const dim_t* degree, const index_t* S, AMGBlockSelect*split_marker, const bool usePanel);
void Preconditioner_LocalAMG_setStrongConnections_Block(SparseMatrix_ptr A, dim_t *degree, index_t *S, const double theta, const double tau);
void Preconditioner_LocalAMG_setStrongConnections(SparseMatrix_ptr A, dim_t *degree, index_t *S, const double theta, const double tau);

SparseMatrix_ptr Preconditioner_LocalAMG_getProlongation(
        SparseMatrix_ptr A_p, const index_t* offset_S,
        const dim_t* degree_S, const index_t* S, dim_t n_C,
        const index_t* counter_C, index_t interpolation_method);

void Preconditioner_LocalAMG_setDirectProlongation_Block(SparseMatrix_ptr P_p, const_SparseMatrix_ptr A_p, const index_t *counter_C);

void Preconditioner_LocalAMG_setDirectProlongation(SparseMatrix_ptr P_p, const_SparseMatrix_ptr A_p, const index_t *counter_C);
void Preconditioner_LocalAMG_setClassicProlongation(SparseMatrix_ptr P_p, SparseMatrix_ptr A_p, const index_t* offset_S, const dim_t* degree_S, const index_t* S, const index_t *counter_C);
void Preconditioner_LocalAMG_setClassicProlongation_Block(SparseMatrix_ptr P_p, SparseMatrix_ptr A_p, const index_t* offset_S, const dim_t* degree_S, const index_t* S, const index_t *counter_C);
index_t Preconditioner_LocalAMG_getMaxLevel(const Preconditioner_LocalAMG * in);
double Preconditioner_LocalAMG_getCoarseLevelSparsity(const Preconditioner_LocalAMG * in);
dim_t Preconditioner_LocalAMG_getNumCoarseUnknowns(const Preconditioner_LocalAMG * in);
void Preconditioner_LocalAMG_enforceFFConnectivity(const dim_t n, const index_t* offset_S, const dim_t* degree_S, const index_t* S, AMGBlockSelect*split_marker);


struct Preconditioner_BoomerAMG
{
    Paso_BOOMERAMG_Handler* pt;
};

void Preconditioner_BoomerAMG_free(Preconditioner_BoomerAMG * in);
Preconditioner_BoomerAMG* Preconditioner_BoomerAMG_alloc(SystemMatrix_ptr A, Options* options);
void Preconditioner_BoomerAMG_solve(SystemMatrix_ptr A, Preconditioner_BoomerAMG * amg, double * x, double * b);


struct Preconditioner_AMG_Root 
{
  bool is_local;
  Preconditioner_AMG* amg;
  Preconditioner_LocalAMG* localamg;
  Preconditioner_BoomerAMG* boomeramg;
  dim_t sweeps;
  Preconditioner_Smoother* amgsubstitute;
};

Preconditioner_AMG_Root* Preconditioner_AMG_Root_alloc(SystemMatrix_ptr A, Options* options);
void Preconditioner_AMG_Root_free(Preconditioner_AMG_Root * in);
void Preconditioner_AMG_Root_solve(SystemMatrix_ptr A, Preconditioner_AMG_Root * amg, double * x, double * b);

/*===============================================*/
/* ILU preconditioner */
struct Solver_ILU
{
    double* factors;
};

/* RILU preconditioner */
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


/*******************************************/
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

