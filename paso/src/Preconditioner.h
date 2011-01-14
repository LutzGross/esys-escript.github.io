
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#ifndef INC_PRECS
#define INC_PRECS

#include "SystemMatrix.h"
#include "performance.h"


#define PASO_AMG_UNDECIDED -1
#define PASO_AMG_IN_F 0
#define PASO_AMG_IN_C 1

/* GAUSS SEIDEL & Jacobi */
typedef struct Paso_Preconditioner_LocalSmoother {
   bool_t Jacobi;
   double* diag;
   double* buffer;
   index_t* pivot;
} Paso_Preconditioner_LocalSmoother;

typedef struct Paso_Preconditioner_Smoother {
   Paso_Preconditioner_LocalSmoother* localSmoother;
   bool_t is_local;
} Paso_Preconditioner_Smoother;

void Paso_Preconditioner_Smoother_free(Paso_Preconditioner_Smoother * in);
void Paso_Preconditioner_LocalSmoother_free(Paso_Preconditioner_LocalSmoother * in);

Paso_Preconditioner_Smoother* Paso_Preconditioner_Smoother_alloc(Paso_SystemMatrix * A_p, const bool_t jacobi, const bool_t is_local, const bool_t verbose);
Paso_Preconditioner_LocalSmoother* Paso_Preconditioner_LocalSmoother_alloc(Paso_SparseMatrix * A_p,const bool_t jacobi, const bool_t verbose);

void Paso_Preconditioner_Smoother_solve(Paso_SystemMatrix* A, Paso_Preconditioner_Smoother * gs, double * x, const double * b, const dim_t sweeps, const bool_t x_is_initial);
void Paso_Preconditioner_LocalSmoother_solve(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x, const double * b, const dim_t sweeps, const bool_t x_is_initial);

void Paso_Preconditioner_LocalSmoother_Sweep(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);
void Paso_Preconditioner_LocalSmoother_Sweep_sequential(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);
void Paso_Preconditioner_LocalSmoother_Sweep_tiled(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);
void Paso_Preconditioner_LocalSmoother_Sweep_colored(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);


/* Local preconditioner */
struct Paso_Preconditioner_AMG {
   dim_t level;
   dim_t n;
   dim_t n_F;
   dim_t n_block;
   Paso_SystemMatrix * A_C;  /* coarse level matrix */
   Paso_SystemMatrix * P;   /* prolongation n x n_C*/ 
   Paso_SystemMatrix * R;   /* restriction  n_C x n */
   
   Paso_Preconditioner_Smoother* Smoother;
   dim_t post_sweeps;
   dim_t pre_sweeps;
   index_t reordering;  /* applied reordering in direct solver */
   dim_t refinements;  /* number of refinements in direct solver (typically =0) */
   double* r;         /* buffer for residual */
   double* x_C;       /* solution of coarse level system */
   double* b_C;       /* right hand side of coarse level system */
   struct Paso_Preconditioner_AMG * AMG_C;
};
typedef struct Paso_Preconditioner_AMG Paso_Preconditioner_AMG;

void Paso_Preconditioner_AMG_free(Paso_Preconditioner_AMG * in);
Paso_Preconditioner_AMG* Paso_Preconditioner_AMG_alloc(Paso_SystemMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Preconditioner_AMG_solve(Paso_SystemMatrix* A, Paso_Preconditioner_AMG * amg, double * x, double * b);

/* Local AMG preconditioner */
struct Paso_Preconditioner_LocalAMG {
   dim_t level;
   dim_t n;
   dim_t n_F;
   dim_t n_block;
   Paso_SparseMatrix * A_C;  /* coarse level matrix */
   Paso_SparseMatrix * P;   /* prolongation n x n_C*/ 
   Paso_SparseMatrix * R;   /* restriction  n_C x n */

   Paso_Preconditioner_LocalSmoother* Smoother;
   dim_t post_sweeps;
   dim_t pre_sweeps;
   index_t reordering;  /* applied reordering in direct solver */
   dim_t refinements;  /* number of refinements in direct solver (typically =0) */
   double* r;         /* buffer for residual */
   double* x_C;       /* solution of coarse level system */
   double* b_C;       /* right hand side of coarse level system */
   struct Paso_Preconditioner_LocalAMG * AMG_C;
};
typedef struct Paso_Preconditioner_LocalAMG Paso_Preconditioner_LocalAMG;

void Paso_Preconditioner_LocalAMG_free(Paso_Preconditioner_LocalAMG * in);
Paso_Preconditioner_LocalAMG* Paso_Preconditioner_LocalAMG_alloc(Paso_SparseMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Preconditioner_LocalAMG_solve(Paso_SparseMatrix* A, Paso_Preconditioner_LocalAMG * amg, double * x, double * b);

void Paso_Preconditioner_LocalAMG_RungeStuebenSearch(const dim_t n, const index_t* offset, const dim_t* degree, const index_t* S, index_t*split_marker, const bool_t usePanel);
void Paso_Preconditioner_LocalAMG_setStrongConnections_Block(Paso_SparseMatrix* A, dim_t *degree, index_t *S, const double theta, const double tau);
void Paso_Preconditioner_LocalAMG_setStrongConnections(Paso_SparseMatrix* A, dim_t *degree, index_t *S, const double theta, const double tau);
Paso_SparseMatrix* Paso_Preconditioner_LocalAMG_getProlongation(Paso_SparseMatrix* A_p, const index_t* offset_S, const dim_t* degree_S, const index_t* S, const dim_t n_C, const index_t* counter_C, const index_t interpolation_method);
void Paso_Preconditioner_LocalAMG_setDirectProlongation_Block(Paso_SparseMatrix* P_p, const Paso_SparseMatrix* A_p, const index_t *counter_C);
void Paso_Preconditioner_LocalAMG_setDirectProlongation(Paso_SparseMatrix* P_p, const Paso_SparseMatrix* A_p, const index_t *counter_C);
void Paso_Preconditioner_LocalAMG_setClassicProlongation(Paso_SparseMatrix* P_p, Paso_SparseMatrix* A_p, const index_t* offset_S, const dim_t* degree_S, const index_t* S, const index_t *counter_C);
void Paso_Preconditioner_LocalAMG_setClassicProlongation_Block(Paso_SparseMatrix* P_p, Paso_SparseMatrix* A_p, const index_t* offset_S, const dim_t* degree_S, const index_t* S, const index_t *counter_C);
index_t Paso_Preconditioner_LocalAMG_getMaxLevel(const Paso_Preconditioner_LocalAMG * in);
double Paso_Preconditioner_LocalAMG_getCoarseLevelSparsity(const Paso_Preconditioner_LocalAMG * in);
dim_t Paso_Preconditioner_LocalAMG_getNumCoarseUnknwons(const Paso_Preconditioner_LocalAMG * in);
void Paso_Preconditioner_LocalAMG_enforceFFConnectivity(const dim_t n, const index_t* offset_S, const dim_t* degree_S, const index_t* S, index_t*split_marker);



/*===============================================*/
/* ILU preconditioner */
struct Paso_Solver_ILU {
  double* factors;
};
typedef struct Paso_Solver_ILU Paso_Solver_ILU;



/* RILU preconditioner */
struct Paso_Solver_RILU {
  dim_t n;
  dim_t n_block;
  dim_t n_F;
  dim_t n_C;
  double* inv_A_FF;
  index_t* A_FF_pivot;
  Paso_SparseMatrix * A_FC;
  Paso_SparseMatrix * A_CF;
  index_t* rows_in_F;
  index_t* rows_in_C;
  index_t* mask_F;
  index_t* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  struct Paso_Solver_RILU * RILU_of_Schur;
};
typedef struct Paso_Solver_RILU Paso_Solver_RILU;



/* general preconditioner interface */

typedef struct Paso_Preconditioner {
  dim_t type;
  dim_t sweeps;
  /* jacobi preconditioner */
  Paso_Preconditioner_Smoother* jacobi;
  /* Gauss-Seidel preconditioner */
  Paso_Preconditioner_Smoother* gs;  
  /* amg preconditioner */
  Paso_Preconditioner_LocalAMG* localamg;
  Paso_Preconditioner_LocalSmoother* localamgsubstitute;
  
  /* ilu preconditioner */
  Paso_Solver_ILU* ilu;
  /* rilu preconditioner */
  Paso_Solver_RILU* rilu;
  
} Paso_Preconditioner;

void Paso_Preconditioner_free(Paso_Preconditioner*);
Paso_Preconditioner* Paso_Preconditioner_alloc(Paso_SystemMatrix* A,Paso_Options* options);
void Paso_Preconditioner_solve(Paso_Preconditioner* prec, Paso_SystemMatrix* A,double*,double*);


/*******************************************/
void Paso_Solver_ILU_free(Paso_Solver_ILU * in);
Paso_Solver_ILU* Paso_Solver_getILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveILU(Paso_SparseMatrix * A, Paso_Solver_ILU * ilu, double * x, const double * b);

void Paso_Solver_RILU_free(Paso_Solver_RILU * in);
Paso_Solver_RILU* Paso_Solver_getRILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveRILU(Paso_Solver_RILU * rilu, double * x, double * b);

void Paso_Solver_updateIncompleteSchurComplement(Paso_SparseMatrix* A_CC, Paso_SparseMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot, Paso_SparseMatrix *A_FC);


#endif /* #ifndef INC_PRECS */
