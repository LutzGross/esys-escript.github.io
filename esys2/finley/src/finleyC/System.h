/* $Id$ */

/**************************************************************/

/*   Finley: SystemMatrix and SystemVector */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_SYSTEM
#define INC_FINLEY_SYSTEM


#include "Mesh.h"
#include "escript/Data/DataC.h"

/**************************************************************/

/*  this struct holds a stiffness matrix: */

/* TODO: implement macro for make file */
/* libraries for linear solvers*/
#define NO_LIB 0
#define SGI_SCSL 1
/* TODO: implement macro for make file */

/* TODO: review need for PTR_INDEX/PTR_OFFSET */
/* Typically INDEX_OFFSET and PTR_OFFSET are 1 if the matrix is
   handed over to a library written in FORTRAN. */

#define INDEX_OFFSET 0
#define PTR_OFFSET 0
#define FINLEY_DEFAULT_MATRIX_TYPE CSR

/* matrix type */
#define  CSC 0
#define  CSR 1

typedef int Finley_SystemMatrixType;

typedef struct Finley_SystemMatrix {
  Finley_SystemMatrixType type;
  int reference_counter;
  maybelong total_row_block_size;
  maybelong total_col_block_size;

  maybelong row_block_size;
  maybelong col_block_size;

  maybelong num_rows;
  maybelong num_cols;
  int symmetric;
  maybelong* ptr;
  maybelong* index;
  maybelong lenOfVal;
  double *val;

  void* solve;  /* pointer to data needed by the direct solver */
  void* iterative; /* pointer to data needed by the iterative solver */

} Finley_SystemMatrix;

/* solver options */

#define NO_REORDERING 0
#define MINIMUM_FILL_IN 1
#define NESTED_DISSECTION 2

#define DEFAULT_METHOD 0
#define PCG 1
#define CR 2
#define CGS 3
#define BICGSTAB 4
#define SSOR 5
#define ILU0 6
#define ILUT 7
#define JACOBI 8

typedef struct {
    int verbose;
    int reordering;
    double tolerance;
    double final_residual;
    int iterative_method;
    int preconditioner;
    int iter_max;
    int iter;
    double drop_tolerance;
    double drop_storage;
    int iterative;

} Finley_SolverOptions;

/*  interfaces: */

Finley_SystemMatrix* Finley_SystemMatrix_alloc(Finley_Mesh*,Finley_SystemMatrixType,int,int,int,int,int);
int Finley_comparIndex(const void *,const void *);
void Finley_SystemMatrix_dealloc(Finley_SystemMatrix*);

void Finley_SystemMatrix_setValues(Finley_SystemMatrix*,double);
void Finley_SystemMatrix_copy(Finley_SystemMatrix*,double*);
void Finley_SystemMatrix_add(Finley_SystemMatrix*,int,maybelong*, int,int,maybelong*,int, double*);
void Finley_RawScaledSystemMatrixVector(double alpha, Finley_SystemMatrix* A, double* in, double beta, double* out);
void Finley_SystemMatrixVector(escriptDataC *,Finley_SystemMatrix *, escriptDataC *);
void Finley_ScaledSystemMatrixVector(double, Finley_SystemMatrix *, escriptDataC *,double,escriptDataC *);

void Finley_SystemMatrix_saveMM(Finley_SystemMatrix *, char *);
void Finley_SystemMatrix_solve(Finley_SystemMatrix* A,escriptDataC* out,escriptDataC* in,Finley_SolverOptions* options);
void Finley_SystemMatrix_solve_free(Finley_SystemMatrix* in);
void Finley_SystemMatrix_iterative(Finley_SystemMatrix* A,escriptDataC* out,escriptDataC* in,Finley_SolverOptions* options);
void Finley_SystemMatrix_iterative_free(Finley_SystemMatrix* in);
void Finley_SystemMatrix_nullifyRowsAndCols(Finley_SystemMatrix* A, escriptDataC* mask_row, escriptDataC* mask_col, double main_diagonal_value);
void Finley_SystemMatrixNullify(Finley_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Finley_SystemMatrix_setDefaults(Finley_SolverOptions*);

#endif /* #ifndef INC_FINLEY_SYSTEM */

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.3  2004/08/28 12:58:08  gross
 * SimpleSolve is not running yet: problem with == of functionsspace
 *
 * Revision 1.2  2004/07/20 01:51:59  gross
 * finley CPP is almost compiling now
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
