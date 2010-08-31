
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

/* jacobi  preconditioner */

typedef struct Paso_Preconditioner_Jacobi {
  double* values;
  index_t* pivot;
} Paso_Preconditioner_Jacobi;

/* GS preconditioner */
typedef struct Paso_Preconditioner_LocalGS {
   double* diag;
   double* buffer;
   index_t* pivot;
   dim_t sweeps;
} Paso_Preconditioner_LocalGS;

typedef struct Paso_Preconditioner_GS {
   Paso_Preconditioner_LocalGS* localGS;
   bool_t is_local;
} Paso_Preconditioner_GS;



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
  Paso_Preconditioner_Jacobi* Jacobi;
  Paso_Preconditioner_LocalGS* GS;
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
  Paso_Preconditioner_Jacobi* GS;
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


/* AMG preconditioner on blocks*/
struct Paso_Solver_AMG_System {
    dim_t block_size;
    Paso_SparseMatrix *block[MAX_BLOCK_SIZE];
    Paso_Solver_AMG *amgblock[MAX_BLOCK_SIZE];
};
typedef struct Paso_Solver_AMG_System Paso_Solver_AMG_System;

/* general preconditioner interface */

typedef struct Paso_Solver_Preconditioner {
  dim_t type;

  /* jacobi preconditioner */
  Paso_Preconditioner_Jacobi* jacobi;
  /* Gauss-Seidel preconditioner */
  Paso_Preconditioner_GS* gs;
  
  /* ilu preconditioner */
  Paso_Solver_ILU* ilu;
  /* rilu preconditioner */
  Paso_Solver_RILU* rilu;
  /* amg preconditioner */
  Paso_Solver_AMG* amg;
  /* amg on System */
  Paso_Solver_AMG_System* amgSystem;
  /* amg preconditioner */
  Paso_Solver_AMLI* amli;
  /* amg on System */
  Paso_Solver_AMLI_System* amliSystem;
  
} Paso_Solver_Preconditioner;

void Paso_Preconditioner_free(Paso_Solver_Preconditioner*);
Paso_Solver_Preconditioner* Paso_Preconditioner_alloc(Paso_SystemMatrix* A,Paso_Options* options);
void Paso_Preconditioner_solve(Paso_Solver_Preconditioner* prec, Paso_SystemMatrix* A,double*,double*);

/* JACOBI */
Paso_Preconditioner_Jacobi* Paso_Preconditioner_Jacobi_alloc(Paso_SystemMatrix * A_p);
Paso_Preconditioner_Jacobi* Paso_Preconditioner_LocalJacobi_alloc(Paso_SparseMatrix * A_p);
void Paso_Preconditioner_Jacobi_free(Paso_Preconditioner_Jacobi * in);

void Paso_Preconditioner_Jacobi_solve(Paso_SystemMatrix * A_p, Paso_Preconditioner_Jacobi * prec, double * x, double * b);
void Paso_Preconditioner_LocalJacobi_solve(Paso_SparseMatrix * A_p, Paso_Preconditioner_Jacobi * prec, double * x, double * b);

/* GAUSS SEIDEL */
void Paso_Preconditioner_GS_free(Paso_Preconditioner_GS * in);
void Paso_Preconditioner_LocalGS_free(Paso_Preconditioner_LocalGS * in);
Paso_Preconditioner_GS* Paso_Preconditioner_GS_alloc(Paso_SystemMatrix * A_p, dim_t sweeps, bool_t is_local, bool_t verbose);

Paso_Preconditioner_LocalGS* Paso_Preconditioner_LocalGS_alloc(Paso_SparseMatrix * A_p, dim_t sweeps, bool_t verbose);
void Paso_Preconditioner_GS_solve(Paso_SystemMatrix* A, Paso_Preconditioner_GS * gs, double * x, const double * b);
void Paso_Preconditioner_LocalGS_solve(Paso_SparseMatrix* A, Paso_Preconditioner_LocalGS * gs, double * x, const double * b);

void Paso_Preconditioner_LocalGS_Sweep(Paso_SparseMatrix* A, Paso_Preconditioner_LocalGS * gs, double * x);
void Paso_Preconditioner_LocalGS_Sweep_sequential(Paso_SparseMatrix* A, Paso_Preconditioner_LocalGS * gs, double * x);
void Paso_Preconditioner_LocalGS_Sweep_tiled(Paso_SparseMatrix* A, Paso_Preconditioner_LocalGS * gs, double * x);
void Paso_Preconditioner_LocalGS_Sweep_colored(Paso_SparseMatrix* A, Paso_Preconditioner_LocalGS * gs, double * x);

/*******************************************/
void Paso_Solver_ILU_free(Paso_Solver_ILU * in);
Paso_Solver_ILU* Paso_Solver_getILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveILU(Paso_SparseMatrix * A, Paso_Solver_ILU * ilu, double * x, const double * b);

void Paso_Solver_RILU_free(Paso_Solver_RILU * in);
Paso_Solver_RILU* Paso_Solver_getRILU(Paso_SparseMatrix * A_p,bool_t verbose);
void Paso_Solver_solveRILU(Paso_Solver_RILU * rilu, double * x, double * b);

void Paso_Solver_AMG_System_free(Paso_Solver_AMG_System * in);
void Paso_Solver_AMG_free(Paso_Solver_AMG * in);
Paso_Solver_AMG* Paso_Solver_getAMG(Paso_SparseMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Solver_solveAMG(Paso_Solver_AMG * amg, double * x, double * b);

void Paso_Solver_AMLI_System_free(Paso_Solver_AMLI_System * in);
void Paso_Solver_AMLI_free(Paso_Solver_AMLI * in);
Paso_Solver_AMLI* Paso_Solver_getAMLI(Paso_SparseMatrix * A_p,dim_t level,Paso_Options* options);
void Paso_Solver_solveAMLI(Paso_Solver_AMLI * amli, double * x, double * b);

void Paso_Solver_updateIncompleteSchurComplement(Paso_SparseMatrix* A_CC, Paso_SparseMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot, Paso_SparseMatrix *A_FC);


#endif /* #ifndef INC_PRECS */
