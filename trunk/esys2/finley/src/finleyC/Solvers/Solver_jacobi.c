/* $Id$ */

/**************************************************************/

/* Finley: jacobi preconditioner: */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Solver.h"

/**************************************************************/

/* frees the Jacobi preconditioner */

void Finley_Solver_Jacobi_free(Finley_Solver_Jacobi * in) {
  if (in!=NULL) {
     MEMFREE(in->values);
     MEMFREE(in->pivot);
     MEMFREE(in);
  }
}

/**************************************************************/

/* Jacobi precondioner set up */

Finley_Solver_Jacobi* Finley_Solver_getJacobi(Finley_SystemMatrix * A_p) {
  Finley_Solver_Jacobi* out=NULL;
  int n = A_p->num_cols;
  int n_block=A_p->row_block_size;
  int block_size=A_p->block_size;
  int i,iPtr;
  double A11,A12,A13,A21,A22,A23,A31,A32,A33,D;
  /* check matrix is square */
  if (A_p->col_block_size !=A_p->row_block_size) {
    Finley_ErrorCode = TYPE_ERROR;
    sprintf(Finley_ErrorMsg, "Jacobi preconditioner square block size.");
    return NULL;
  }
  /* check matrix is square */
  if (n_block>3) {
    Finley_ErrorCode = TYPE_ERROR;
    sprintf(Finley_ErrorMsg, "Right now the Jacobi preconditioner supports block size less than 4 only");
    return NULL;
  }
  /* allocate vector to hold main diagonal entries: */
  out=MEMALLOC(1,Finley_Solver_Jacobi);
  if (! Finley_checkPtr(out)) {
      /* allocate vector to hold main diagonal entries: */
      out->n_block=n_block;
      out->n=n;
      out->values = MEMALLOC( ((size_t) n) * ((size_t) block_size),double);
      out->pivot = NULL; /* later use */
      if (! (Finley_checkPtr(out->values))) {
        if (n_block==1) {
           #pragma omp parallel for private(i, iPtr) schedule(static)
           for (i = 0; i < A_p->pattern->n_ptr; i++) {
              out->values[i]=1.;
              /* find main diagonal */
              for (iPtr = A_p->pattern->ptr[i]; iPtr < A_p->pattern->ptr[i + 1]; iPtr++) {
                  if (A_p->pattern->index[iPtr]==i+INDEX_OFFSET) {
                      if (ABS(A_p->val[iPtr])>0.) out->values[i]=1./A_p->val[iPtr];
                      break;
                  }
              }
           }
        } else if (n_block==2) {
           #pragma omp parallel for private(i, iPtr, A11,A12,A21,A22,D) schedule(static)
           for (i = 0; i < A_p->pattern->n_ptr; i++) {
              out->values[i*4+0]= 1.;
              out->values[i*4+1]= 0.;
              out->values[i*4+2]= 0.;
              out->values[i*4+3]= 1.;
              /* find main diagonal */
              for (iPtr = A_p->pattern->ptr[i]; iPtr < A_p->pattern->ptr[i + 1]; iPtr++) {
                  if (A_p->pattern->index[iPtr]==i+INDEX_OFFSET) {
                     A11=A_p->val[iPtr*4];
                     A12=A_p->val[iPtr*4+2];
                     A21=A_p->val[iPtr*4+1];
                     A22=A_p->val[iPtr*4+3];
                     D = A11*A22-A12*A21;
                     if (ABS(D) > 0 ){
                        D=1./D;
                        out->values[i*4]= A22*D;
                        out->values[i*4+1]=-A21*D;
                        out->values[i*4+2]=-A12*D;
                        out->values[i*4+3]= A11*D;
                     }
                     break;
                  }
              }
           }
        } else if (n_block==3) {
           #pragma omp parallel for private(i, iPtr,A11,A12,A13,A21,A22,A23,A31,A32,A33,D) schedule(static)
           for (i = 0; i < A_p->pattern->n_ptr; i++) {
              out->values[i*9  ]=1.;
              out->values[i*9+1]=0.;
              out->values[i*9+2]=0.;
              out->values[i*9+3]=0.;
              out->values[i*9+4]=1.;
              out->values[i*9+5]=0.;
              out->values[i*9+6]=0.;
              out->values[i*9+7]=0.;
              out->values[i*9+8]=1.;
              /* find main diagonal */
              for (iPtr = A_p->pattern->ptr[i]; iPtr < A_p->pattern->ptr[i + 1]; iPtr++) {
                  if (A_p->pattern->index[iPtr]==i+INDEX_OFFSET) {
                      A11=A_p->val[iPtr*9  ];
                      A21=A_p->val[iPtr*9+1];
                      A31=A_p->val[iPtr*9+2];
                      A12=A_p->val[iPtr*9+3];
                      A22=A_p->val[iPtr*9+4];
                      A32=A_p->val[iPtr*9+5];
                      A13=A_p->val[iPtr*9+6];
                      A23=A_p->val[iPtr*9+7];
                      A33=A_p->val[iPtr*9+8];
                      D  =  A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);
                      if (ABS(D) > 0 ){
                         D=1./D;
                         out->values[i*9  ]=(A22*A33-A23*A32)*D;
                         out->values[i*9+1]=(A31*A23-A21*A33)*D;
                         out->values[i*9+2]=(A21*A32-A31*A22)*D;
                         out->values[i*9+3]=(A13*A32-A12*A33)*D;
                         out->values[i*9+4]=(A11*A33-A31*A13)*D;
                         out->values[i*9+5]=(A12*A31-A11*A32)*D;
                         out->values[i*9+6]=(A12*A23-A13*A22)*D;
                         out->values[i*9+7]=(A13*A21-A11*A23)*D;
                         out->values[i*9+8]=(A11*A22-A12*A21)*D;
                      }
                      break;
                  }
              }
           }
        }
      }
  }
  if (Finley_ErrorCode == NO_ERROR) {
     return out;
  } else {
     Finley_Solver_Jacobi_free(out);
     return NULL;
  }
} 
/**************************************************************/

/* applies Jacobi preconditioner */

/* should be called within a parallel region                                              */
/* barrier synconization should be performed to make sure that the input vector available */

void Finley_Solver_solveJacobi(Finley_Solver_Jacobi * prec, double * x, double * b) {
     Finley_Solver_applyBlockDiagonalMatrix(prec->n_block,prec->n,prec->values,prec->pivot,x,b);
     return;
}
