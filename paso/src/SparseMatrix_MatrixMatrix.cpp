
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************

 Paso: Sparse matrix product

*****************************************************************************

   Author: artak@uq.edu.au, l.gross@uq.edu.au

*****************************************************************************/

#include "SparseMatrix.h"
#include "PasoException.h"
#include "PasoUtil.h" // comparIndex

namespace paso {

// forward declarations
void SparseMatrix_MatrixMatrix_DD(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B);
void SparseMatrix_MatrixMatrix_DB(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B);
void SparseMatrix_MatrixMatrix_BD(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B);
void SparseMatrix_MatrixMatrix_BB(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B);

SparseMatrix_ptr<double> SparseMatrix_MatrixMatrix(const_SparseMatrix_ptr<double> A,
                                           const_SparseMatrix_ptr<double> B)
{
    SparseMatrixType C_type;
    SparseMatrix_ptr<double> out;
    if (!( (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) ||
           (A->type & MATRIX_FORMAT_DEFAULT) ||
           (A->type & MATRIX_FORMAT_BLK1) )) {
        throw PasoException("SparseMatrix::MatrixMatrix: Unsupported matrix format of A.");
    }
    if (!( (B->type & MATRIX_FORMAT_DIAGONAL_BLOCK) ||
           (B->type & MATRIX_FORMAT_DEFAULT) ||
           (B->type & MATRIX_FORMAT_BLK1) )) {
        throw PasoException("SparseMatrix::MatrixMatrix: Unsupported matrix format of B.");
    }
    if (A->col_block_size != B->row_block_size) {
        throw PasoException("SparseMatrix::MatrixMatrix: Column block size of A and row block size of B must match.");
    }
    if (A->numCols != B->numRows) {
        throw PasoException("SparseMatrix::MatrixMatrix: number of columns of A and number of rows of B must match.");
    }

    if ( (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) &&
         (B->type & MATRIX_FORMAT_DIAGONAL_BLOCK) ) {
        C_type=MATRIX_FORMAT_DIAGONAL_BLOCK;
    } else {
        C_type=MATRIX_FORMAT_DEFAULT;
    }

    Pattern_ptr outpattern(A->pattern->multiply(MATRIX_FORMAT_DEFAULT, B->pattern));

    out.reset(new SparseMatrix<double>(C_type, outpattern, A->row_block_size, B->col_block_size, false));

    if ( (A->row_block_size == 1) && (B->col_block_size ==1 ) && (A->col_block_size ==1) ) {
        SparseMatrix_MatrixMatrix_DD(out, A, B);
    } else {
        if (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
            if (B->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
                SparseMatrix_MatrixMatrix_DD(out, A, B);
            } else {
                SparseMatrix_MatrixMatrix_DB(out, A, B);
            }
        } else {
            if (B->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
                SparseMatrix_MatrixMatrix_BD(out, A, B);
            } else {
                SparseMatrix_MatrixMatrix_BB(out, A, B);
            }
        }
    }
    return out;
}

// not good for block size 1
void SparseMatrix_MatrixMatrix_BB(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B)
{
    const dim_t n = C->numRows;
    const dim_t row_block_size = C->row_block_size;
    const dim_t col_block_size = C->col_block_size;
    const dim_t A_col_block_size = A->col_block_size;
    const dim_t C_block_size =C->block_size;
    const dim_t B_block_size =B->block_size;
    const dim_t A_block_size =A->block_size;
    double *C_ij, *A_ik, *B_kj;
    double rtmp, C_ij_00, C_ij_10, C_ij_20, C_ij_30, C_ij_01, C_ij_11, C_ij_21, C_ij_31, C_ij_02, C_ij_12, C_ij_22, C_ij_32, C_ij_03, C_ij_13, C_ij_23, C_ij_33;
    dim_t i, ib, irb, icb;
    index_t ij_ptrC, j, ik_ptrA, k, kj_ptrB, *start_p, *where_p;

    if ( (row_block_size == 2) && (col_block_size == 2) && (A_col_block_size == 2) ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_01, C_ij_11)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];

               C_ij_00=0;
               C_ij_10=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij=&C->val[ij_ptrC*4];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                        kj_ptrB += (index_t)(where_p-start_p);
                        A_ik=&(A->val[ik_ptrA*4]);
                        B_kj=&(B->val[kj_ptrB*4]);

                        C_ij_00 +=A_ik[0+2*0]*B_kj[0+2*0]+A_ik[0+2*1]*B_kj[1+2*0];
                        C_ij_10 +=A_ik[1+2*0]*B_kj[0+2*0]+A_ik[1+2*1]*B_kj[1+2*0];

                        C_ij_01 +=A_ik[0+2*0]*B_kj[0+2*1]+A_ik[0+2*1]*B_kj[1+2*1];
                        C_ij_11 +=A_ik[1+2*0]*B_kj[0+2*1]+A_ik[1+2*1]*B_kj[1+2*1];
                  }

               }
               C_ij[0+2*0]=C_ij_00;
               C_ij[1+2*0]=C_ij_10;
               C_ij[0+2*1]=C_ij_01;
               C_ij[1+2*1]=C_ij_11;

            }
         }
      } /* end of parallel region */

   } else if ( (row_block_size == 3) && (col_block_size ==3 ) && (A_col_block_size == 3)  ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_20, C_ij_01, C_ij_11, C_ij_21, C_ij_02, C_ij_12, C_ij_22)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];

               C_ij_00=0;
               C_ij_10=0;
               C_ij_20=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij_21=0;
               C_ij_02=0;
               C_ij_12=0;
               C_ij_22=0;
               C_ij = &C->val[ij_ptrC*9];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                        kj_ptrB += (index_t)(where_p-start_p);
                        A_ik=&(A->val[ik_ptrA*9]);
                        B_kj=&(B->val[kj_ptrB*9]);

                        C_ij_00 +=A_ik[0+3*0]*B_kj[0+3*0]
                                 +A_ik[0+3*1]*B_kj[1+3*0]
                                 +A_ik[0+3*2]*B_kj[2+3*0];
                        C_ij_10 +=A_ik[1+3*0]*B_kj[0+3*0]
                                 +A_ik[1+3*1]*B_kj[1+3*0]
                                 +A_ik[1+3*2]*B_kj[2+3*0];
                        C_ij_20 +=A_ik[2+3*0]*B_kj[0+3*0]
                                 +A_ik[2+3*1]*B_kj[1+3*0]
                                 +A_ik[2+3*2]*B_kj[2+3*0];

                        C_ij_01 +=A_ik[0+3*0]*B_kj[0+3*1]
                                 +A_ik[0+3*1]*B_kj[1+3*1]
                                 +A_ik[0+3*2]*B_kj[2+3*1];
                        C_ij_11 +=A_ik[1+3*0]*B_kj[0+3*1]
                                 +A_ik[1+3*1]*B_kj[1+3*1]
                                 +A_ik[1+3*2]*B_kj[2+3*1];
                        C_ij_21 +=A_ik[2+3*0]*B_kj[0+3*1]
                                 +A_ik[2+3*1]*B_kj[1+3*1]
                                 +A_ik[2+3*2]*B_kj[2+3*1];

                        C_ij_01 +=A_ik[0+3*0]*B_kj[0+3*2]
                                 +A_ik[0+3*1]*B_kj[1+3*2]
                                 +A_ik[0+3*2]*B_kj[2+3*2];
                        C_ij_11 +=A_ik[1+3*0]*B_kj[0+3*2]
                                 +A_ik[1+3*1]*B_kj[1+3*2]
                                 +A_ik[1+3*2]*B_kj[2+3*2];
                        C_ij_21 +=A_ik[2+3*0]*B_kj[0+3*2]
                                 +A_ik[2+3*1]*B_kj[1+3*2]
                                 +A_ik[2+3*2]*B_kj[2+3*2];
                  }

               }
               C_ij[0+3*0]=C_ij_00;
               C_ij[1+3*0]=C_ij_10;
               C_ij[2+3*0]=C_ij_20;
               C_ij[0+3*1]=C_ij_01;
               C_ij[1+3*1]=C_ij_11;
               C_ij[2+3*1]=C_ij_21;
               C_ij[0+3*2]=C_ij_02;
               C_ij[1+3*2]=C_ij_12;
               C_ij[2+3*2]=C_ij_22;
            }
         }
      } /* end of parallel region */
   } else if ( (row_block_size == 4) && (col_block_size ==4 ) && (A_col_block_size == 4)  ){
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_20, C_ij_30, C_ij_01, C_ij_11, C_ij_21, C_ij_31, C_ij_02, C_ij_12, C_ij_22, C_ij_32, C_ij_03, C_ij_13, C_ij_23, C_ij_33)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_20=0;
               C_ij_30=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij_21=0;
               C_ij_31=0;
               C_ij_02=0;
               C_ij_12=0;
               C_ij_22=0;
               C_ij_32=0;
               C_ij_03=0;
               C_ij_13=0;
               C_ij_23=0;
               C_ij_33=0;
               C_ij = &C->val[ij_ptrC*16];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case */
                           kj_ptrB += (index_t)(where_p-start_p);
                           A_ik=&(A->val[ik_ptrA*16]);
                           B_kj=&(B->val[kj_ptrB*16]);

                           C_ij_00 +=A_ik[0+4*0]*B_kj[0+4*0]
                                    +A_ik[0+4*1]*B_kj[1+4*0]
                                    +A_ik[0+4*2]*B_kj[2+4*0]
                                    +A_ik[0+4*3]*B_kj[3+4*0];
                           C_ij_10 +=A_ik[1+4*0]*B_kj[0+4*0]
                                    +A_ik[1+4*1]*B_kj[1+4*0]
                                    +A_ik[1+4*2]*B_kj[2+4*0]
                                    +A_ik[1+4*3]*B_kj[3+4*0];
                           C_ij_20 +=A_ik[2+4*0]*B_kj[0+4*0]
                                    +A_ik[2+4*1]*B_kj[1+4*0]
                                    +A_ik[2+4*2]*B_kj[2+4*0]
                                    +A_ik[2+4*3]*B_kj[3+4*0];
                           C_ij_30 +=A_ik[3+4*0]*B_kj[0+4*0]
                                    +A_ik[3+4*1]*B_kj[1+4*0]
                                    +A_ik[3+4*2]*B_kj[2+4*0]
                                    +A_ik[3+4*3]*B_kj[3+4*0];

                           C_ij_01 +=A_ik[0+4*0]*B_kj[0+4*1]
                                    +A_ik[0+4*1]*B_kj[1+4*1]
                                    +A_ik[0+4*2]*B_kj[2+4*1]
                                    +A_ik[0+4*3]*B_kj[3+4*1];
                           C_ij_11 +=A_ik[1+4*0]*B_kj[0+4*1]
                                    +A_ik[1+4*1]*B_kj[1+4*1]
                                    +A_ik[1+4*2]*B_kj[2+4*1]
                                    +A_ik[1+4*3]*B_kj[3+4*1];
                           C_ij_21 +=A_ik[2+4*0]*B_kj[0+4*1]
                                    +A_ik[2+4*1]*B_kj[1+4*1]
                                    +A_ik[2+4*2]*B_kj[2+4*1]
                                    +A_ik[2+4*3]*B_kj[3+4*1];
                           C_ij_31 +=A_ik[3+4*0]*B_kj[0+4*1]
                                    +A_ik[3+4*1]*B_kj[1+4*1]
                                    +A_ik[3+4*2]*B_kj[2+4*1]
                                    +A_ik[3+4*3]*B_kj[3+4*1];

                           C_ij_02 +=A_ik[0+4*0]*B_kj[0+4*2]
                                    +A_ik[0+4*1]*B_kj[1+4*2]
                                    +A_ik[0+4*2]*B_kj[2+4*2]
                                    +A_ik[0+4*3]*B_kj[3+4*2];
                           C_ij_12 +=A_ik[1+4*0]*B_kj[0+4*2]
                                    +A_ik[1+4*1]*B_kj[1+4*2]
                                    +A_ik[1+4*2]*B_kj[2+4*2]
                                    +A_ik[1+4*3]*B_kj[3+4*2];
                           C_ij_22 +=A_ik[2+4*0]*B_kj[0+4*2]
                                    +A_ik[2+4*1]*B_kj[1+4*2]
                                    +A_ik[2+4*2]*B_kj[2+4*2]
                                    +A_ik[2+4*3]*B_kj[3+4*2];
                           C_ij_32 +=A_ik[3+4*0]*B_kj[0+4*2]
                                    +A_ik[3+4*1]*B_kj[1+4*2]
                                    +A_ik[3+4*2]*B_kj[2+4*2]
                                    +A_ik[3+4*3]*B_kj[3+4*2];

                           C_ij_03 +=A_ik[0+4*0]*B_kj[0+4*3]
                                    +A_ik[0+4*1]*B_kj[1+4*3]
                                    +A_ik[0+4*2]*B_kj[2+4*3]
                                    +A_ik[0+4*3]*B_kj[3+4*3];
                           C_ij_13 +=A_ik[1+4*0]*B_kj[0+4*3]
                                    +A_ik[1+4*1]*B_kj[1+4*3]
                                    +A_ik[1+4*2]*B_kj[2+4*3]
                                    +A_ik[1+4*3]*B_kj[3+4*3];
                           C_ij_23 +=A_ik[2+4*0]*B_kj[0+4*3]
                                    +A_ik[2+4*1]*B_kj[1+4*3]
                                    +A_ik[2+4*2]*B_kj[2+4*3]
                                    +A_ik[2+4*3]*B_kj[3+4*3];
                           C_ij_33 +=A_ik[3+4*0]*B_kj[0+4*3]
                                    +A_ik[3+4*1]*B_kj[1+4*3]
                                    +A_ik[3+4*2]*B_kj[2+4*3]
                                    +A_ik[3+4*3]*B_kj[3+4*3];
                     }
                  }
                  C_ij[0+4*0]=C_ij_00;
                  C_ij[1+4*0]=C_ij_10;
                  C_ij[2+4*0]=C_ij_20;
                  C_ij[3+4*0]=C_ij_30;
                  C_ij[0+4*1]=C_ij_01;
                  C_ij[1+4*1]=C_ij_11;
                  C_ij[2+4*1]=C_ij_21;
                  C_ij[3+4*1]=C_ij_31;
                  C_ij[0+4*2]=C_ij_02;
                  C_ij[1+4*2]=C_ij_12;
                  C_ij[2+4*2]=C_ij_22;
                  C_ij[3+4*2]=C_ij_32;
                  C_ij[0+4*3]=C_ij_03;
                  C_ij[1+4*3]=C_ij_13;
                  C_ij[2+4*3]=C_ij_23;
                  C_ij[3+4*3]=C_ij_33;
            }
         }
      } /* end of parallel region */

   } else { // blocksizes differ or >4
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, rtmp, A_ik,B_kj )
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*C_block_size]);
               for (ib=0; ib<C_block_size; ++ib)  C_ij[ib]=0;
               for(ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*A_block_size]);
                       B_kj=&(B->val[kj_ptrB*B_block_size]);

                       for (irb=0; irb<row_block_size; ++irb) {
                          for (icb=0; icb<col_block_size; ++icb) {
                             rtmp=C_ij[irb+row_block_size*icb];
                             for (ib=0; ib<A_col_block_size; ++ib) {
                                rtmp+=A_ik[irb+row_block_size*ib]*B_kj[ib+A_col_block_size*icb];
                             }
                             C_ij[irb+row_block_size*icb]=rtmp;
                          }
                       }
                  }
               }
            }
         }
      } /* end of parallel region */
   }
}

// not good for block size 1
void SparseMatrix_MatrixMatrix_DB(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B)
{
    const dim_t n = C->numRows;
    const dim_t row_block_size = C->row_block_size;
    const dim_t col_block_size = C->col_block_size;
    const dim_t A_col_block_size = A->col_block_size;
    const dim_t C_block_size =C->block_size;
    const dim_t B_block_size =B->block_size;
    const dim_t A_block_size =A->block_size;
    double *C_ij, *A_ik, *B_kj;
    double rtmp, C_ij_00, C_ij_10, C_ij_20, C_ij_30, C_ij_01, C_ij_11, C_ij_21, C_ij_31, C_ij_02, C_ij_12, C_ij_22, C_ij_32, C_ij_03, C_ij_13, C_ij_23, C_ij_33;
    dim_t i, ib, irb, icb;
    index_t ij_ptrC, j, ik_ptrA, k, kj_ptrB, *start_p, *where_p;

    if ( (row_block_size == 2) && (col_block_size ==2 ) && (A_block_size == 2) ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_01, C_ij_11)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij = &C->val[ij_ptrC*4];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA  < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                        kj_ptrB += (index_t)(where_p-start_p);
                        A_ik=&(A->val[ik_ptrA*2]);
                        B_kj=&(B->val[kj_ptrB*4]);

                        C_ij_00 +=A_ik[0]*B_kj[0+2*0];
                        C_ij_10 +=A_ik[1]*B_kj[1+2*0];

                        C_ij_01 +=A_ik[0]*B_kj[0+2*1];
                        C_ij_11 +=A_ik[1]*B_kj[1+2*1];
                  }

               }
               C_ij[0+2*0]=C_ij_00;
               C_ij[1+2*0]=C_ij_10;
               C_ij[0+2*1]=C_ij_01;
               C_ij[1+2*1]=C_ij_11;

            }
         }
      } /* end of parallel region */

   } else if ( (row_block_size == 3) && (col_block_size ==3 ) && (A_block_size == 3)  ){
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_20, C_ij_01, C_ij_11, C_ij_21, C_ij_02, C_ij_12, C_ij_22)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_20=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij_21=0;
               C_ij_02=0;
               C_ij_12=0;
               C_ij_22=0;
               C_ij = &C->val[ij_ptrC*9];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                        kj_ptrB += (index_t)(where_p-start_p);
                        A_ik=&(A->val[ik_ptrA*3]);
                        B_kj=&(B->val[kj_ptrB*9]);

                        C_ij_00 +=A_ik[0]*B_kj[0+3*0];
                        C_ij_10 +=A_ik[1]*B_kj[1+3*0];
                        C_ij_20 +=A_ik[2]*B_kj[2+3*0];

                        C_ij_01 +=A_ik[0]*B_kj[0+3*1];
                        C_ij_11 +=A_ik[1]*B_kj[1+3*1];
                        C_ij_21 +=A_ik[2]*B_kj[2+3*1];

                        C_ij_02 +=A_ik[0]*B_kj[0+3*2];
                        C_ij_12 +=A_ik[1]*B_kj[1+3*2];
                        C_ij_22 +=A_ik[2]*B_kj[2+3*2];

                  }

               }
               C_ij[0+3*0]=C_ij_00;
               C_ij[1+3*0]=C_ij_10;
               C_ij[2+3*0]=C_ij_20;
               C_ij[0+3*1]=C_ij_01;
               C_ij[1+3*1]=C_ij_11;
               C_ij[2+3*1]=C_ij_21;
               C_ij[0+3*2]=C_ij_02;
               C_ij[1+3*2]=C_ij_12;
               C_ij[2+3*2]=C_ij_22;
            }
         }
      } /* end of parallel region */
   } else if ( (row_block_size == 4) && (col_block_size ==4 ) && (A_block_size == 4)  ){
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_20, C_ij_30, C_ij_01, C_ij_11, C_ij_21, C_ij_31, C_ij_02, C_ij_12, C_ij_22, C_ij_32, C_ij_03, C_ij_13, C_ij_23, C_ij_33)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_20=0;
               C_ij_30=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij_21=0;
               C_ij_31=0;
               C_ij_02=0;
               C_ij_12=0;
               C_ij_22=0;
               C_ij_32=0;
               C_ij_03=0;
               C_ij_13=0;
               C_ij_23=0;
               C_ij_33=0;
               C_ij = &C->val[ij_ptrC*16];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                           kj_ptrB += (index_t)(where_p-start_p);
                           A_ik=&(A->val[ik_ptrA*4]);
                           B_kj=&(B->val[kj_ptrB*16]);

                           C_ij_00 +=A_ik[0]*B_kj[0+4*0];
                           C_ij_10 +=A_ik[1]*B_kj[1+4*0];
                           C_ij_20 +=A_ik[2]*B_kj[2+4*0];
                           C_ij_30 +=A_ik[3]*B_kj[3+4*0];

                           C_ij_01 +=A_ik[0]*B_kj[0+4*1];
                           C_ij_11 +=A_ik[1]*B_kj[1+4*1];
                           C_ij_21 +=A_ik[2]*B_kj[2+4*1];
                           C_ij_31 +=A_ik[3]*B_kj[3+4*1];

                           C_ij_02 +=A_ik[0]*B_kj[0+4*2];
                           C_ij_12 +=A_ik[1]*B_kj[1+4*2];
                           C_ij_22 +=A_ik[2]*B_kj[2+4*2];
                           C_ij_32 +=A_ik[3]*B_kj[3+4*2];

                           C_ij_03 +=A_ik[0]*B_kj[0+4*3];
                           C_ij_13 +=A_ik[1]*B_kj[1+4*3];
                           C_ij_23 +=A_ik[2]*B_kj[2+4*3];
                           C_ij_33 +=A_ik[3]*B_kj[3+4*3];
                  }
               }
               C_ij[0+4*0]=C_ij_00;
               C_ij[1+4*0]=C_ij_10;
               C_ij[2+4*0]=C_ij_20;
               C_ij[3+4*0]=C_ij_30;
               C_ij[0+4*1]=C_ij_01;
               C_ij[1+4*1]=C_ij_11;
               C_ij[2+4*1]=C_ij_21;
               C_ij[3+4*1]=C_ij_31;
               C_ij[0+4*2]=C_ij_02;
               C_ij[1+4*2]=C_ij_12;
               C_ij[2+4*2]=C_ij_22;
               C_ij[3+4*2]=C_ij_32;
               C_ij[0+4*3]=C_ij_03;
               C_ij[1+4*3]=C_ij_13;
               C_ij[2+4*3]=C_ij_23;
               C_ij[3+4*3]=C_ij_33;
            }
         }
      } /* end of parallel region */

   } else { // blocksizes differ or >4
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, rtmp, A_ik,B_kj )
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*C_block_size]);
               for (ib=0; ib<C_block_size; ++ib)  C_ij[ib]=0;
               for(ik_ptrA = A->pattern->ptr[i]; ik_ptrA  < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*A_block_size]);
                       B_kj=&(B->val[kj_ptrB*B_block_size]);

                       for (irb=0; irb<A_block_size; ++irb) {
                          rtmp=A_ik[irb];
                          for (icb=0; icb<col_block_size; ++icb) {
                                C_ij[irb+row_block_size*icb]+=rtmp*B_kj[irb+A_col_block_size*icb];
                          }
                       }
                  }
               }
            }
         }
      } /* end of parallel region */
   }
}

// not good for block size 1
void SparseMatrix_MatrixMatrix_BD(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B)
{
    const dim_t n = C->numRows;
    const dim_t row_block_size = C->row_block_size;
    const dim_t col_block_size = C->col_block_size;
    const dim_t C_block_size =C->block_size;
    const dim_t B_block_size =B->block_size;
    const dim_t A_block_size =A->block_size;
    double *C_ij, *A_ik, *B_kj;
    double rtmp, C_ij_00, C_ij_10, C_ij_20, C_ij_30, C_ij_01, C_ij_11, C_ij_21, C_ij_31, C_ij_02, C_ij_12, C_ij_22, C_ij_32, C_ij_03, C_ij_13, C_ij_23, C_ij_33;
    dim_t i, ib, irb, icb;
    index_t ij_ptrC, j, ik_ptrA, k, kj_ptrB, *start_p, *where_p;

    if ( (row_block_size == 2) && (col_block_size ==2 ) && (B_block_size == 2) ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_01, C_ij_11)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij = &C->val[ij_ptrC*4];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                        kj_ptrB += (index_t)(where_p-start_p);
                        A_ik=&(A->val[ik_ptrA*4]);
                        B_kj=&(B->val[kj_ptrB*2]);

                        C_ij_00 +=A_ik[0+2*0]*B_kj[0];
                        C_ij_10 +=A_ik[1+2*0]*B_kj[0];

                        C_ij_01 +=A_ik[0+2*1]*B_kj[1];
                        C_ij_11 +=A_ik[1+2*1]*B_kj[1];
                  }

               }
               C_ij[0+2*0]=C_ij_00;
               C_ij[1+2*0]=C_ij_10;
               C_ij[0+2*1]=C_ij_01;
               C_ij[1+2*1]=C_ij_11;

            }
         }
      } /* end of parallel region */

   } else if ( (row_block_size == 3) && (col_block_size ==3 ) && (B_block_size == 3)  ){
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_20, C_ij_01, C_ij_11, C_ij_21, C_ij_02, C_ij_12, C_ij_22)
      {
#pragma omp for schedule(static)
         for (i = 0; i < n; i++) {
            for (ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_20=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij_21=0;
               C_ij_02=0;
               C_ij_12=0;
               C_ij_22=0;
               C_ij = &C->val[ij_ptrC*9];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                        kj_ptrB += (index_t)(where_p-start_p);
                        A_ik=&(A->val[ik_ptrA*9]);
                        B_kj=&(B->val[kj_ptrB*3]);

                        C_ij_00 +=A_ik[0+3*0]*B_kj[0];
                        C_ij_10 +=A_ik[1+3*0]*B_kj[0];
                        C_ij_20 +=A_ik[2+3*0]*B_kj[0];

                        C_ij_01 +=A_ik[0+3*1]*B_kj[1];
                        C_ij_11 +=A_ik[1+3*1]*B_kj[1];
                        C_ij_21 +=A_ik[2+3*1]*B_kj[1];

                        C_ij_02 +=A_ik[0+3*2]*B_kj[2];
                        C_ij_12 +=A_ik[1+3*2]*B_kj[2];
                        C_ij_22 +=A_ik[2+3*2]*B_kj[2];
                  }
               }
               C_ij[0+3*0]=C_ij_00;
               C_ij[1+3*0]=C_ij_10;
               C_ij[2+3*0]=C_ij_20;
               C_ij[0+3*1]=C_ij_01;
               C_ij[1+3*1]=C_ij_11;
               C_ij[2+3*1]=C_ij_21;
               C_ij[0+3*2]=C_ij_02;
               C_ij[1+3*2]=C_ij_12;
               C_ij[2+3*2]=C_ij_22;
            }
         }
      } /* end of parallel region */
   } else if ( (row_block_size == 4) && (col_block_size ==4 ) && (B_block_size == 4)  ){
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, A_ik,B_kj,C_ij_00, C_ij_10, C_ij_20, C_ij_30, C_ij_01, C_ij_11, C_ij_21, C_ij_31, C_ij_02, C_ij_12, C_ij_22, C_ij_32, C_ij_03, C_ij_13, C_ij_23, C_ij_33)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_00=0;
               C_ij_10=0;
               C_ij_20=0;
               C_ij_30=0;
               C_ij_01=0;
               C_ij_11=0;
               C_ij_21=0;
               C_ij_31=0;
               C_ij_02=0;
               C_ij_12=0;
               C_ij_22=0;
               C_ij_32=0;
               C_ij_03=0;
               C_ij_13=0;
               C_ij_23=0;
               C_ij_33=0;
               C_ij = &C->val[ij_ptrC*16];

               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                           kj_ptrB += (index_t)(where_p-start_p);
                           A_ik=&(A->val[ik_ptrA*16]);
                           B_kj=&(B->val[kj_ptrB*4]);

                           C_ij_00 +=A_ik[0+4*0]*B_kj[0];
                           C_ij_10 +=A_ik[1+4*0]*B_kj[0];
                           C_ij_20 +=A_ik[2+4*0]*B_kj[0];
                           C_ij_30 +=A_ik[3+4*0]*B_kj[0];

                           C_ij_01 +=A_ik[0+4*1]*B_kj[1];
                           C_ij_11 +=A_ik[1+4*1]*B_kj[1];
                           C_ij_21 +=A_ik[2+4*1]*B_kj[1];
                           C_ij_31 +=A_ik[3+4*1]*B_kj[1];

                           C_ij_02 +=A_ik[0+4*2]*B_kj[2];
                           C_ij_12 +=A_ik[1+4*2]*B_kj[2];
                           C_ij_22 +=A_ik[2+4*2]*B_kj[2];
                           C_ij_32 +=A_ik[3+4*2]*B_kj[2];

                           C_ij_03 +=A_ik[0+4*3]*B_kj[3];
                           C_ij_13 +=A_ik[1+4*3]*B_kj[3];
                           C_ij_23 +=A_ik[2+4*3]*B_kj[3];
                           C_ij_33 +=A_ik[3+4*3]*B_kj[3];
                     }
                  }
                  C_ij[0+4*0]=C_ij_00;
                  C_ij[1+4*0]=C_ij_10;
                  C_ij[2+4*0]=C_ij_20;
                  C_ij[3+4*0]=C_ij_30;
                  C_ij[0+4*1]=C_ij_01;
                  C_ij[1+4*1]=C_ij_11;
                  C_ij[2+4*1]=C_ij_21;
                  C_ij[3+4*1]=C_ij_31;
                  C_ij[0+4*2]=C_ij_02;
                  C_ij[1+4*2]=C_ij_12;
                  C_ij[2+4*2]=C_ij_22;
                  C_ij[3+4*2]=C_ij_32;
                  C_ij[0+4*3]=C_ij_03;
                  C_ij[1+4*3]=C_ij_13;
                  C_ij[2+4*3]=C_ij_23;
                  C_ij[3+4*3]=C_ij_33;
            }
         }
      } /* end of parallel region */

   } else { // blocksizes differ or >4
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, irb, icb, ib, rtmp, A_ik,B_kj )
      {
#pragma omp for schedule(static)
         for (i = 0; i < n; i++) {
            for (ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*C_block_size]);
               for (ib=0; ib<C_block_size; ++ib)  C_ij[ib]=0;
               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*A_block_size]);
                       B_kj=&(B->val[kj_ptrB*B_block_size]);

                       for (icb=0; icb<B_block_size; ++icb) {
                          rtmp=B_kj[icb];
                          for (irb=0; irb<row_block_size; ++irb) {
                             C_ij[irb+row_block_size*icb]+=A_ik[irb+row_block_size*icb]*rtmp;
                          }
                       }
                  }
               }
            }
         }
      } /* end of parallel region */
   }
}

// not good for block size 1
void SparseMatrix_MatrixMatrix_DD(SparseMatrix_ptr<double> C, const_SparseMatrix_ptr<double> A,
                                  const_SparseMatrix_ptr<double> B)
{
    const dim_t n = C->numRows;
    const dim_t C_block_size =C->block_size;
    const dim_t B_block_size =B->block_size;
    const dim_t A_block_size =A->block_size;
    double *C_ij, *A_ik, *B_kj;
    double C_ij_0, C_ij_1, C_ij_2, C_ij_3;
    dim_t i, ib;
    index_t ij_ptrC, j, ik_ptrA, k, kj_ptrB, *start_p, *where_p;

    if ( (A_block_size == 1) && (B_block_size ==1) && (C_block_size == 1) ) {
#pragma omp parallel private(i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, C_ij_0)
      {
#pragma omp for schedule(static)
         for(i = 0; i < n; i++) {
            for(ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij_0=0;
               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       C_ij_0+=A->val[ik_ptrA]*B->val[kj_ptrB];
                  }

               }
               C->val[ij_ptrC]=C_ij_0;
            }
         }
      } /* end of parallel region */
    } else if ( (A_block_size == 2) && (B_block_size ==2) && (C_block_size == 2) ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, A_ik,B_kj, C_ij_0, C_ij_1)
      {
#pragma omp for schedule(static)
         for (i = 0; i < n; i++) {
            for (ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*2]);
               C_ij_0=0;
               C_ij_1=0;
               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*2]);
                       B_kj=&(B->val[kj_ptrB*2]);

                       C_ij_0+=A_ik[0]*B_kj[0];
                       C_ij_1+=A_ik[1]*B_kj[1];
                  }

               }
               C_ij[0]=C_ij_0;
               C_ij[1]=C_ij_1;
            }
         }
      } /* end of parallel region */
   } else if ( (A_block_size == 3) && (B_block_size ==3) && (C_block_size == 3) ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, A_ik,B_kj, C_ij_0, C_ij_1, C_ij_2)
      {
#pragma omp for schedule(static)
         for (i = 0; i < n; i++) {
            for (ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*C_block_size]);
               C_ij_0=0;
               C_ij_1=0;
               C_ij_2=0;
               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*3]);
                       B_kj=&(B->val[kj_ptrB*3]);

                       C_ij_0+=A_ik[0]*B_kj[0];
                       C_ij_1+=A_ik[1]*B_kj[1];
                       C_ij_2+=A_ik[2]*B_kj[2];
                  }

               }
               C_ij[0]=C_ij_0;
               C_ij[1]=C_ij_1;
               C_ij[2]=C_ij_2;
            }
         }
      } /* end of parallel region */
   } else if ( (A_block_size == 4) && (B_block_size ==4 ) && (C_block_size == 4) ) {
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, A_ik,B_kj, C_ij_0, C_ij_1, C_ij_2, C_ij_3 )
      {
#pragma omp for schedule(static)
         for (i = 0; i < n; i++) {
            for (ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*C_block_size]);
               C_ij_0=0;
               C_ij_1=0;
               C_ij_2=0;
               C_ij_3=0;
               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*4]);
                       B_kj=&(B->val[kj_ptrB*4]);

                       C_ij_0+=A_ik[0]*B_kj[0];
                       C_ij_1+=A_ik[1]*B_kj[1];
                       C_ij_2+=A_ik[2]*B_kj[2];
                       C_ij_3+=A_ik[3]*B_kj[3];
                  }

               }
               C_ij[0]=C_ij_0;
               C_ij[1]=C_ij_1;
               C_ij[2]=C_ij_2;
               C_ij[3]=C_ij_3;
            }
         }
      } /* end of parallel region */
   } else { // blocksizes differ or > 4
#pragma omp parallel private(C_ij, i, ij_ptrC, j, ik_ptrA, k, kj_ptrB, start_p, where_p, ib, A_ik,B_kj )
      {
#pragma omp for schedule(static)
         for (i = 0; i < n; i++) {
            for (ij_ptrC = C->pattern->ptr[i]; ij_ptrC < C->pattern->ptr[i+1]; ++ij_ptrC) {
               j = C->pattern->index[ij_ptrC];
               C_ij=&(C->val[ij_ptrC*C_block_size]);
               for (ib=0; ib<C_block_size; ++ib)  C_ij[ib]=0;
               for (ik_ptrA = A->pattern->ptr[i]; ik_ptrA < A->pattern->ptr[i+1]; ++ik_ptrA ) {
                  k=A->pattern->index[ik_ptrA];
                  kj_ptrB=B->pattern->ptr[k];
                  /* find j in B[k,:] */
                  start_p=&(B->pattern->index[kj_ptrB]);
                  where_p=(index_t*)bsearch(&j, start_p,
                                            B->pattern->ptr[k + 1]-kj_ptrB,
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (where_p != NULL) { // this should always be the case
                       kj_ptrB += (index_t)(where_p-start_p);
                       A_ik=&(A->val[ik_ptrA*A_block_size]);
                       B_kj=&(B->val[kj_ptrB*B_block_size]);
                       for (ib=0; ib<std::min(A_block_size, B_block_size); ++ib) C_ij[ib]+=A_ik[ib]*B_kj[ib];
                  }
               }
            }
         }
      } /* end of parallel region */
   }
}

} // namespace paso

