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


#include "escript/Data/DataC.h"
#include "SystemPattern.h"

/**************************************************************/

/*  this struct holds a stiffness matrix: */

/* TODO: implement macro for make file */
/* libraries for linear solvers*/
#define NO_LIB 0
#define SGI_SCSL 1

/* matrix type */
#define  CSC 0
#define  CSR 1
/* these formats are used in the SCSL context */
#define  CSC_SYM 2
#define  CSR_SYM 3
#define  CSC_BLK1 4
#define  CSR_BLK1 5
#define  CSC_BLK1_SYM 6
#define  CSR_BLK1_SYM 7

typedef int Finley_SystemMatrixType;

typedef struct Finley_SystemMatrix {
  Finley_SystemMatrixType type;
  int reference_counter;

  maybelong logical_row_block_size;
  maybelong logical_col_block_size;
  maybelong logical_block_size;

  maybelong row_block_size;
  maybelong col_block_size;
  maybelong block_size;

  maybelong num_rows;
  maybelong num_cols;

  Finley_SystemMatrixPattern* pattern;

  size_t len;
  double *val;

  void* direct;  /* pointer to data needed by the direct solver */
  void* iterative; /* pointer to data needed by the iterative solver */

} Finley_SystemMatrix;

/* solver options */


#include "escript/Data/UtilC.h"

typedef struct {
    int method;
    int symmetric;
    double tolerance;

    int verbose;
    int reordering;
    double final_residual;
    int preconditioner;
    int iter_max;
    int iter;
    double drop_tolerance;
    double drop_storage;
    int truncation;
    int restart;


} Finley_SolverOptions;

/*  interfaces: */

Finley_SystemMatrix* Finley_SystemMatrix_alloc(Finley_SystemMatrixType,Finley_SystemMatrixPattern*,int,int);
Finley_SystemMatrix* Finley_SystemMatrix_reference(Finley_SystemMatrix*);
void Finley_SystemMatrix_dealloc(Finley_SystemMatrix*);

void Finley_SystemMatrix_setValues(Finley_SystemMatrix*,double);
void Finley_SystemMatrix_copy(Finley_SystemMatrix*,double*);
void Finley_SystemMatrix_add(Finley_SystemMatrix*,int,maybelong*, int,int,maybelong*,int, double*);
void Finley_RawScaledSystemMatrixVector(double alpha, Finley_SystemMatrix* A, double* in, double beta, double* out);
void Finley_SystemMatrixVector(escriptDataC *,Finley_SystemMatrix *, escriptDataC *);
void Finley_ScaledSystemMatrixVector(double, Finley_SystemMatrix *, escriptDataC *,double,escriptDataC *);

void Finley_SystemMatrix_saveMM(Finley_SystemMatrix *, char *);
Finley_SystemMatrix* Finley_SystemMatrix_loadMM_toCSR(char *);
void Finley_SystemMatrix_solve(Finley_SystemMatrix* A,escriptDataC* out,escriptDataC* in,Finley_SolverOptions* options);
void Finley_SystemMatrix_solve_free(Finley_SystemMatrix* in);
void Finley_SystemMatrix_nullifyRowsAndCols(Finley_SystemMatrix* A, escriptDataC* mask_row, escriptDataC* mask_col, double main_diagonal_value);
void Finley_SystemMatrixNullify(Finley_SystemMatrix* A, double* mask_row, double* mask_col, double main_diagonal_value);
void Finley_SystemMatrix_setDefaults(Finley_SolverOptions*);
int Finley_SystemMatrix_getSystemMatrixTypeId(int, int);
Finley_SystemMatrix* Finley_SystemMatrix_getSubmatrix(Finley_SystemMatrix* A,int,maybelong*,maybelong*);

#endif /* #ifndef INC_FINLEY_SYSTEM */
