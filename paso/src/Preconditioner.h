
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

#define MAX_BLOCK_SIZE 3

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

void Paso_Preconditioner_Smoother_solve(Paso_SystemMatrix* A, Paso_Preconditioner_Smoother * gs, double * x, const double * b, const dim_t sweeps);
void Paso_Preconditioner_LocalSmoother_solve(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x, const double * b, const dim_t sweeps);

void Paso_Preconditioner_LocalSmoother_Sweep(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);
void Paso_Preconditioner_LocalSmoother_Sweep_sequential(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);
void Paso_Preconditioner_LocalSmoother_Sweep_tiled(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);
void Paso_Preconditioner_LocalSmoother_Sweep_colored(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * gs, double * x);

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

struct Paso_Solver_Smoother {
  dim_t ID;  
  Paso_Preconditioner_LocalSmoother* Jacobi;
  Paso_Preconditioner_LocalSmoother* GS;
};
typedef struct  Paso_Solver_Smoother  Paso_Solver_Smoother;

/* AMG preconditioner */
struct Paso_Solver_AMG {
  dim_t n;
  dim_t level;
  bool_t coarsest_level;
  dim_t n_block;
  dim_t n_F;
  dim_t n_C;
  
  Paso_SparseMatrix * A_FF;
  Paso_SparseMatrix * A_FC;
  Paso_SparseMatrix * A_CF;
  Paso_SparseMatrix * W_FC;
  
  Paso_SparseMatrix * P;
  Paso_SparseMatrix * R;
  
  index_t* rows_in_F;
  index_t* rows_in_C;
  index_t* mask_F;
  index_t* mask_C;
  double* x_F;
  double* b_F;
  double* x_C;
  double* b_C;
  
  dim_t post_sweeps;
  dim_t pre_sweeps;
  
  Paso_SparseMatrix * A;
  Paso_SparseMatrix * AOffset1;
  Paso_SparseMatrix * AUnrolled;
  void* solver;
  Paso_Solver_Smoother* Smoother;
  struct Paso_Solver_AMG * AMG_of_Coarse;
};
typedef struct Paso_Solver_AMG Paso_Solver_AMG;


/* AMLI preconditioner */
struct Paso_Solver_AMLI {
  dim_t n;
  dim_t level;
  bool_t coarsest_level;
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
  
  dim_t post_sweeps;
  dim_t pre_sweeps;
  
  Paso_SparseMatrix * A;
  Paso_SparseMatrix * AOffset1;
  void* solver;
  Paso_Preconditioner_LocalSmoother* GS;
  struct Paso_Solver_AMLI * AMLI_of_Schur;
};
typedef struct Paso_Solver_AMLI Paso_Solver_AMLI;

/* AMLI preconditioner on blocks*/
struct Paso_Solver_AMLI_System {
    dim_t block_size;
    Paso_SparseMatrix *block[MAX_BLOCK_SIZE];
    Paso_Solver_AMLI *amliblock[MAX_BLOCK_SIZE];
};
typedef struct Paso_Solver_AMLI_System Paso_Solver_AMLI_System;


/* general preconditioner interface */

typedef struct Paso_Solver_Preconditioner {
  dim_t type;
  dim_t sweeps;
  /* jacobi preconditioner */
  Paso_Preconditioner_Smoother* jacobi;
  /* Gauss-Seidel preconditioner */
  Paso_Preconditioner_Smoother* gs;
  
  /* ilu preconditioner */
  Paso_Solver_ILU* ilu;
  /* rilu preconditioner */
  Paso_Solver_RILU* rilu;
  /* amg preconditioner */
  Paso_Solver_AMG* amg;
  /* amg preconditioner */
  Paso_Solver_AMLI* amli;
  /* amg on System */
  Paso_Solver_AMLI_System* amliSystem;
  
} Paso_Solver_Preconditioner;

void Paso_Preconditioner_free(Paso_Solver_Preconditioner*);
Paso_Solver_Preconditioner* Paso_Preconditioner_alloc(Paso_SystemMatrix* A,Paso_Options* options);
void Paso_Preconditioner_solve(Paso_Solver_Preconditioner* prec, Paso_SystemMatrix* A,double*,double*);




/* AMG: */
void Paso_Solver_AMG_free(Paso_Solver_AMG * in);
Paso_Solver_AMG* Paso_Solver_getAMG(Paso_SparseMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Solver_solveAMG(Paso_Solver_AMG * amg, double * x, double * b);

/*******************************************/
void Paso_Solver_ILU_free(Paso_Solver_ILU * in);
Paso_Solver_ILU* Paso_Solver_getILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveILU(Paso_SparseMatrix * A, Paso_Solver_ILU * ilu, double * x, const double * b);

void Paso_Solver_RILU_free(Paso_Solver_RILU * in);
Paso_Solver_RILU* Paso_Solver_getRILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveRILU(Paso_Solver_RILU * rilu, double * x, double * b);



void Paso_Solver_AMLI_System_free(Paso_Solver_AMLI_System * in);
void Paso_Solver_AMLI_free(Paso_Solver_AMLI * in);
Paso_Solver_AMLI* Paso_Solver_getAMLI(Paso_SparseMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Solver_solveAMLI(Paso_Solver_AMLI * amli, double * x, double * b);

void Paso_Solver_updateIncompleteSchurComplement(Paso_SparseMatrix* A_CC, Paso_SparseMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot, Paso_SparseMatrix *A_FC);


#endif /* #ifndef INC_PRECS */
