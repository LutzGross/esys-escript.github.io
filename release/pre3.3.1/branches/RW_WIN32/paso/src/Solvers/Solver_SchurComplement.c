/* $Id$ */

/**************************************************************/

/* Paso: updates A_CC <- ACC-ACF AFF^{-1} AFC                 */

/* no check of consistency of matrices !!!!                   */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: gross@access.edu.au                                */

/**************************************************************/

#include "paso/Paso.h"
#include "paso/SystemMatrix.h"
#include "Solver.h"

/**************************************************************/



void Paso_Solver_updateIncompleteSchurComplement(Paso_SystemMatrix* A_CC,Paso_SystemMatrix *A_CF,double* invA_FF,index_t* A_FF_pivot,Paso_SystemMatrix *A_FC) {
  index_t iPtr_CC,*index_CC,col_CF,col_FC, *where_p,iPtr_CC_2,i,iPtr_CF,iPtr_FC;
  dim_t index_CC_len;
  bool_t set_A;
  dim_t n_rows=A_CC->num_rows;
  dim_t n_block=A_CC->row_block_size;
  register double A_CF_11,A_CF_21,A_CF_31,A_CF_12,A_CF_22,A_CF_32,A_CF_13,A_CF_23,A_CF_33,
         invA_FF_11,invA_FF_21,invA_FF_31,invA_FF_12,invA_FF_22,invA_FF_32,invA_FF_13,invA_FF_23,invA_FF_33,
         A11,A21,A31,A12,A22,A32,A13,A23,A33,A_FC_11,A_FC_21,A_FC_31,A_FC_12,A_FC_22,A_FC_32,A_FC_13,A_FC_23,A_FC_33;
  if (n_block==1) {
     #pragma omp parallel for private(i,iPtr_CC,index_CC,index_CC_len,col_CF,set_A,iPtr_CF,iPtr_FC,col_FC,where_p,A11) schedule(static)
     for (i = 0; i < n_rows;++i) {
        iPtr_CC=A_CC->pattern->ptr[i];
        index_CC=&(A_CC->pattern->index[iPtr_CC]);
        index_CC_len=(size_t)(A_CC->pattern->ptr[i+1]-A_CC->pattern->ptr[i]);
        /* now we run through the columns of A_CF in row  i */
        for (iPtr_CF = A_CF->pattern->ptr[i]; iPtr_CF < A_CF->pattern->ptr[i + 1]; ++iPtr_CF) {
             col_CF=A_CF->pattern->index[iPtr_CF];
             set_A=TRUE;
             for (iPtr_FC = A_FC->pattern->ptr[col_CF]; iPtr_FC < A_FC->pattern->ptr[col_CF + 1]; ++iPtr_FC) {
                col_FC=A_FC->pattern->index[iPtr_FC];
                /* is (i,col_FC) in the shape of A_CC ? */
                where_p=(index_t*)bsearch(&col_FC,index_CC,index_CC_len,sizeof(index_t),Paso_comparIndex);
                if (where_p!=NULL) {
                    if (set_A) {
                       A11=A_CF->val[iPtr_CF]*invA_FF[col_CF];
                       set_A=FALSE;
                   }
                   A_CC->val[iPtr_CC+(index_t)(where_p-index_CC)]-=A11*A_FC->val[iPtr_FC];
                }
            } /* end of iPtr_FC loop */
         } /* end of iPtr_CF loop */
      } /* end of irow loop */
   } else if (n_block==2) {
      #pragma omp parallel for private(i,iPtr_CC,index_CC,index_CC_len,iPtr_CF,col_CF,iPtr_FC,col_FC,where_p,iPtr_CC_2,set_A,A_CF_11,A_CF_21,A_CF_12,A_CF_22,invA_FF_11,invA_FF_21,invA_FF_12,invA_FF_22,A11,A21,A12,A22,A_FC_11,A_FC_21,A_FC_12,A_FC_22) schedule(static)
     for (i = 0; i < n_rows;++i) {
        iPtr_CC=A_CC->pattern->ptr[i];
        index_CC=&(A_CC->pattern->index[iPtr_CC]);
        index_CC_len=(size_t)(A_CC->pattern->ptr[i+1]-A_CC->pattern->ptr[i]);
        /* now we run through the columns of A_CF in row  i */
        for (iPtr_CF = A_CF->pattern->ptr[i]; iPtr_CF < A_CF->pattern->ptr[i + 1]; ++iPtr_CF) {
             col_CF=A_CF->pattern->index[iPtr_CF];
             set_A=TRUE;
             for (iPtr_FC = A_FC->pattern->ptr[col_CF]; iPtr_FC < A_FC->pattern->ptr[col_CF + 1]; ++iPtr_FC) {
                col_FC=A_FC->pattern->index[iPtr_FC];
                /* is (i,col_FC) in the shape of A_CC ? */
                where_p=(index_t*)bsearch(&col_FC,index_CC,index_CC_len,sizeof(index_t),Paso_comparIndex);
                if (where_p!=NULL) {
                    iPtr_CC_2=iPtr_CC+(index_t)(where_p-index_CC);
                    /* this calculutes A_CF*invA_FF(i,col_CF) */
                    if (set_A) {
                       A_CF_11=A_CF->val[4*iPtr_CF  ];
                       A_CF_21=A_CF->val[4*iPtr_CF+1];
                       A_CF_12=A_CF->val[4*iPtr_CF+2];
                       A_CF_22=A_CF->val[4*iPtr_CF+3];

                       invA_FF_11=invA_FF[4*col_CF  ];
                       invA_FF_21=invA_FF[4*col_CF+1];
                       invA_FF_12=invA_FF[4*col_CF+2];
                       invA_FF_22=invA_FF[4*col_CF+3];

                       A11=A_CF_11*invA_FF_11+A_CF_12*invA_FF_21;
                       A21=A_CF_21*invA_FF_11+A_CF_22*invA_FF_21;
                       A12=A_CF_11*invA_FF_12+A_CF_12*invA_FF_22;
                       A22=A_CF_21*invA_FF_12+A_CF_22*invA_FF_22;
 
                       set_A=FALSE;
                   }

                   A_FC_11=A_FC->val[4*iPtr_FC  ];
                   A_FC_21=A_FC->val[4*iPtr_FC+1];
                   A_FC_12=A_FC->val[4*iPtr_FC+2];
                   A_FC_22=A_FC->val[4*iPtr_FC+3];

                   A_CC->val[4*iPtr_CC_2  ]-=A11*A_FC_11+A12*A_FC_21;
                   A_CC->val[4*iPtr_CC_2+1]-=A21*A_FC_11+A22*A_FC_21;
                   A_CC->val[4*iPtr_CC_2+2]-=A11*A_FC_12+A12*A_FC_22;
                   A_CC->val[4*iPtr_CC_2+3]-=A21*A_FC_12+A22*A_FC_22;

                }
            } /* end of iPtr_FC loop */
         } /* end of iPtr_CF loop */
      } /* end of irow loop */
   } else if (n_block==3) {
      #pragma omp parallel for private(i,iPtr_CC,index_CC,index_CC_len,iPtr_CF,col_CF,iPtr_FC,col_FC,where_p,iPtr_CC_2,set_A,A_CF_11,A_CF_21,A_CF_31,A_CF_12,A_CF_22,A_CF_32,A_CF_13,A_CF_23,A_CF_33,invA_FF_11,invA_FF_21,invA_FF_31,invA_FF_12,invA_FF_22,invA_FF_32,invA_FF_13,invA_FF_23,invA_FF_33,A11,A21,A31,A12,A22,A32,A13,A23,A33,A_FC_11,A_FC_21,A_FC_31,A_FC_12,A_FC_22,A_FC_32,A_FC_13,A_FC_23,A_FC_33) schedule(static)
     for (i = 0; i < n_rows;++i) {
        iPtr_CC=A_CC->pattern->ptr[i];
        index_CC=&(A_CC->pattern->index[iPtr_CC]);
        index_CC_len=(size_t)(A_CC->pattern->ptr[i+1]-A_CC->pattern->ptr[i]);
        /* now we run through the columns of A_CF in row  i */
        for (iPtr_CF = A_CF->pattern->ptr[i]; iPtr_CF < A_CF->pattern->ptr[i + 1]; ++iPtr_CF) {
             col_CF=A_CF->pattern->index[iPtr_CF];
             set_A=TRUE;
             for (iPtr_FC = A_FC->pattern->ptr[col_CF]; iPtr_FC < A_FC->pattern->ptr[col_CF + 1]; ++iPtr_FC) {
                col_FC=A_FC->pattern->index[iPtr_FC];
                /* is (i,col_FC) in the shape of A_CC ? */
                where_p=(index_t*)bsearch(&col_FC,index_CC,index_CC_len,sizeof(index_t),Paso_comparIndex);
                if (where_p!=NULL) {
                    iPtr_CC_2=iPtr_CC+(index_t)(where_p-index_CC);
                    /* this calculutes A_CF*invA_FF(i,col_CF) */
                    if (set_A) {
                       A_CF_11=A_CF->val[9*iPtr_CF  ];
                       A_CF_21=A_CF->val[9*iPtr_CF+1];
                       A_CF_31=A_CF->val[9*iPtr_CF+2];
                       A_CF_12=A_CF->val[9*iPtr_CF+3];
                       A_CF_22=A_CF->val[9*iPtr_CF+4];
                       A_CF_32=A_CF->val[9*iPtr_CF+5];
                       A_CF_13=A_CF->val[9*iPtr_CF+6];
                       A_CF_23=A_CF->val[9*iPtr_CF+7];
                       A_CF_33=A_CF->val[9*iPtr_CF+8];

                       invA_FF_11=invA_FF[9*col_CF  ];
                       invA_FF_21=invA_FF[9*col_CF+1];
                       invA_FF_31=invA_FF[9*col_CF+2];
                       invA_FF_12=invA_FF[9*col_CF+3];
                       invA_FF_22=invA_FF[9*col_CF+4];
                       invA_FF_32=invA_FF[9*col_CF+5];
                       invA_FF_13=invA_FF[9*col_CF+6];
                       invA_FF_23=invA_FF[9*col_CF+7];
                       invA_FF_33=invA_FF[9*col_CF+8];

                       A11=A_CF_11*invA_FF_11+A_CF_12*invA_FF_21+A_CF_13*invA_FF_31;
                       A21=A_CF_21*invA_FF_11+A_CF_22*invA_FF_21+A_CF_23*invA_FF_31;
                       A31=A_CF_31*invA_FF_11+A_CF_32*invA_FF_21+A_CF_33*invA_FF_31;
                       A12=A_CF_11*invA_FF_12+A_CF_12*invA_FF_22+A_CF_13*invA_FF_32;
                       A22=A_CF_21*invA_FF_12+A_CF_22*invA_FF_22+A_CF_23*invA_FF_32;
                       A32=A_CF_31*invA_FF_12+A_CF_32*invA_FF_22+A_CF_33*invA_FF_32;
                       A13=A_CF_11*invA_FF_13+A_CF_12*invA_FF_23+A_CF_13*invA_FF_33;
                       A23=A_CF_21*invA_FF_13+A_CF_22*invA_FF_23+A_CF_23*invA_FF_33;
                       A33=A_CF_31*invA_FF_13+A_CF_32*invA_FF_23+A_CF_33*invA_FF_33;
 
                       set_A=FALSE;
                   }

                   A_FC_11=A_FC->val[9*iPtr_FC  ];
                   A_FC_21=A_FC->val[9*iPtr_FC+1];
                   A_FC_31=A_FC->val[9*iPtr_FC+2];
                   A_FC_12=A_FC->val[9*iPtr_FC+3];
                   A_FC_22=A_FC->val[9*iPtr_FC+4];
                   A_FC_32=A_FC->val[9*iPtr_FC+5];
                   A_FC_13=A_FC->val[9*iPtr_FC+6];
                   A_FC_23=A_FC->val[9*iPtr_FC+7];
                   A_FC_33=A_FC->val[9*iPtr_FC+8];

                   A_CC->val[9*iPtr_CC_2  ]-=A11*A_FC_11+A12*A_FC_21+A13*A_FC_31;
                   A_CC->val[9*iPtr_CC_2+1]-=A21*A_FC_11+A22*A_FC_21+A23*A_FC_31;
                   A_CC->val[9*iPtr_CC_2+2]-=A31*A_FC_11+A32*A_FC_21+A33*A_FC_31;
                   A_CC->val[9*iPtr_CC_2+3]-=A11*A_FC_12+A12*A_FC_22+A13*A_FC_32;
                   A_CC->val[9*iPtr_CC_2+4]-=A21*A_FC_12+A22*A_FC_22+A23*A_FC_32;
                   A_CC->val[9*iPtr_CC_2+5]-=A31*A_FC_12+A32*A_FC_22+A33*A_FC_32;
                   A_CC->val[9*iPtr_CC_2+6]-=A11*A_FC_13+A12*A_FC_23+A13*A_FC_33;
                   A_CC->val[9*iPtr_CC_2+7]-=A21*A_FC_13+A22*A_FC_23+A23*A_FC_33;
                   A_CC->val[9*iPtr_CC_2+8]-=A31*A_FC_13+A32*A_FC_23+A33*A_FC_33;

                }
            } /* end of iPtr_FC loop */
         } /* end of iPtr_CF loop */
      } /* end of irow loop */
   }
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:40  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:50  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
