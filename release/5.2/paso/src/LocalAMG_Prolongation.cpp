
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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


/****************************************************************************/

/* Paso: defines AMG prolongation  */

/****************************************************************************/

/* Author: Artak Amirbekyan, artak@uq.edu.au, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Options.h"
#include "PasoUtil.h"
#include "Preconditioner.h"
#include "SparseMatrix.h"

namespace paso {

/****************************************************************************

    Methods necessary for AMG preconditioner

    Construct the n x n_C prolongation matrix P from A_p.

    The columns in A_p to be considered are marked by counter_C[n] where
    an unknown i to be considered in P is marked by 0<= counter_C[i] < n_C
    and counter_C[i]  gives the new column number in P.
    S defines the strong connections.

    The pattern of P is formed as follows:

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise
    If row i is not C, then P[i,j] <> 0 if counter_C[k]==j (k in C) and (i,k) is strong connection.

    Two settings for P are implemented (see below).
*/

SparseMatrix_ptr Preconditioner_LocalAMG_getProlongation(SparseMatrix_ptr A_p,
        const index_t* offset_S, const dim_t* degree_S, const index_t* S,
        dim_t n_C, const index_t* counter_C, index_t interpolation_method)
{
    const dim_t n_block=A_p->row_block_size;
    index_t *ptr=NULL, *index=NULL,j, iptr;
    const dim_t n =A_p->numRows;
    dim_t i,p,z, len_P;

    ptr=new index_t[n+1];

    // count the number of entries per row in the Prolongation matrix

#pragma omp parallel for private(i,z,iptr,j,p) schedule(static)
    for (i=0; i<n; ++i) {
            if (counter_C[i] >= 0) {
                z=1;    /* i is a C unknown */
            } else {
                z=0;
                iptr=offset_S[i];
                for (p=0; p<degree_S[i]; ++p) {
                    j=S[iptr+p];  /* this is a strong connection */
                    if (counter_C[j]>=0) z++; /* and is in C */
                }
            }
            ptr[i]=z;
    }
    len_P=util::cumsum(n,ptr);
    ptr[n]=len_P;

    /* allocate and create index vector for prolongation: */
    index=new index_t[len_P];

#pragma omp parallel for private(i,z,iptr,j,p)  schedule(static)
    for (i=0;i<n;++i) {
        if (counter_C[i]>=0) {
            index[ptr[i]]=counter_C[i];  /* i is a C unknown */
        } else {
            z=0;
            iptr=offset_S[i];
            for (p=0; p<degree_S[i]; ++p) {
                j=S[iptr+p];  /* this is a strong connection */
                if (counter_C[j]>=0) {  /* and is in C */
                    index[ptr[i]+z]=counter_C[j];
                    z++; /* and is in C */
                }
            }
        }
    }
    Pattern_ptr outpattern;
    outpattern.reset(new Pattern(MATRIX_FORMAT_DEFAULT, n, n_C, ptr, index));
    /* now we need to create a matrix and fill it */
    SparseMatrix_ptr out;
    out.reset(new SparseMatrix(MATRIX_FORMAT_DIAGONAL_BLOCK,
                               outpattern, n_block, n_block, false));

    if ( (interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING) || (interpolation_method == PASO_CLASSIC_INTERPOLATION) ) {
        if (n_block == 1) {
            Preconditioner_LocalAMG_setClassicProlongation(
                    out, A_p, offset_S, degree_S, S, counter_C);
        } else {
            Preconditioner_LocalAMG_setClassicProlongation_Block(
                    out, A_p, offset_S, degree_S, S, counter_C);
        }
    } else {
        if (n_block == 1) {
            Preconditioner_LocalAMG_setDirectProlongation(
                    out, A_p, counter_C);
        } else {
            Preconditioner_LocalAMG_setDirectProlongation_Block(
                    out, A_p, counter_C);
        }
    }
    return out;
}

/*
    Direct Prolongation:
    -------------------

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise.
    If row i is not C, then P[i,j] = - a[i] * A[i,k]/A[i,i] with j=counter_C[k]>=0 and k in S

   and    a[i]=
             alpha[i] = sum_s min(A[i,s],0)/(sum_{s in S and C} min(A[i,s],0))   A[i,k]<0
                   or                                                         if
             beta[i] = sum_s max(A[i,s],0)/(sum_{s in S and C} max(A[i,s],0))   A[i,k]>0


*/

void Preconditioner_LocalAMG_setDirectProlongation(SparseMatrix_ptr P_p,
        const_SparseMatrix_ptr A_p, const index_t* counter_C)
{
   dim_t i;
   const dim_t n =A_p->numRows;
   double alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos, A_ij, A_ii, rtmp;
   index_t iPtr, j, offset;
   index_t *where_p, *start_p;

   #pragma omp parallel for private(A_ii, offset, where_p, start_p, i, alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos,iPtr,j, A_ij , rtmp)  schedule(static)
   for (i=0;i<n;++i) {
      if (counter_C[i]>=0) {
            offset = P_p->pattern->ptr[i];
            P_p->val[offset]=1.;  /* i is a C row */
      } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
         /* if i is an F row we first calculate alpha and beta: */
         sum_all_neg=0; /* sum of all negative values in row i of A */
         sum_all_pos=0; /* sum of all positive values in row i of A */
         sum_strong_neg=0; /* sum of all negative values A_ij where j is in C and strongly connected to i*/
         sum_strong_pos=0; /* sum of all positive values A_ij where j is in C and strongly connected to i*/
         A_ii=0;
         for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
            j=A_p->pattern->index[iPtr];
            A_ij=A_p->val[iPtr];
            if(j==i) {
               A_ii=A_ij;
            } else {

               if(A_ij< 0)  {
                  sum_all_neg+=A_ij;
               } else {
                  sum_all_pos+=A_ij;
               }

               if (counter_C[j]>=0) {
                  /* is i strongly connected with j? We search for counter_C[j] in P[i,:] */
                  start_p=&(P_p->pattern->index[P_p->pattern->ptr[i]]);
                  where_p=(index_t*)bsearch(&(counter_C[j]), start_p,
                                            P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i],
                                            sizeof(index_t),
                                            util::comparIndex);
                  if (! (where_p == NULL) ) { /* yes i strongly connected with j */
                        offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
                        P_p->val[offset]=A_ij; /* will be modified later */
                        if (A_ij< 0)  {
                           sum_strong_neg+=A_ij;
                        } else {
                           sum_strong_pos+=A_ij;
                        }
                  }
               }
            }
         }
         if(sum_strong_neg<0) {
            alpha= sum_all_neg/sum_strong_neg;
         } else {
            alpha=0;
         }
         if(sum_strong_pos>0) {
            beta= sum_all_pos/sum_strong_pos;
         } else {
            beta=0;
            A_ii+=sum_all_pos;
         }
         if ( A_ii > 0.) {
            rtmp=(-1.)/A_ii;
            alpha*=rtmp;
            beta*=rtmp;
         }
         for (iPtr=P_p->pattern->ptr[i];iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
            A_ij=P_p->val[iPtr];
            if (A_ij > 0 ) {
               P_p->val[iPtr]*=beta;
            } else {
               P_p->val[iPtr]*=alpha;
            }
         }
      }
   }
}

void Preconditioner_LocalAMG_setDirectProlongation_Block(SparseMatrix_ptr P_p,
        const_SparseMatrix_ptr A_p, const index_t* counter_C)
{
   dim_t i;
   const dim_t n =A_p->numRows;
   const dim_t row_block=A_p->row_block_size;
   const dim_t A_block = A_p->block_size;
   double *alpha, *beta, *sum_all_neg, *sum_all_pos, *sum_strong_neg, *sum_strong_pos, *A_ii;
   double A_ij, rtmp;
   index_t iPtr, j, offset, ib;
   index_t *where_p, *start_p;

   #pragma omp parallel private(ib, rtmp, A_ii, offset, where_p, start_p, i, alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos,iPtr,j, A_ij )
   {
      sum_all_neg=new  double[row_block]; /* sum of all negative values in row i of A */
      sum_all_pos=new double[row_block]; /* sum of all positive values in row i of A */
      sum_strong_neg=new double[row_block]; /* sum of all negative values A_ij where j is in C and strongly connected to i*/
      sum_strong_pos=new double[row_block]; /* sum of all positive values A_ij where j is in C and strongly connected to i*/
      alpha=new double[row_block];
      beta=new double[row_block];
      A_ii=new double[row_block];

      #pragma omp for schedule(static)
      for (i=0;i<n;++i) {
         if (counter_C[i]>=0) {
            offset = P_p->pattern->ptr[i];
            for (ib =0; ib<row_block; ++ib) P_p->val[row_block*offset+ib]=1.;  /* i is a C row */
         } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
            /* if i is an F row we first calculate alpha and beta: */
            for (ib =0; ib<row_block; ++ib) {
               sum_all_neg[ib]=0;
               sum_all_pos[ib]=0;
               sum_strong_neg[ib]=0;
               sum_strong_pos[ib]=0;
               A_ii[ib]=0;
            }
            for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
               j=A_p->pattern->index[iPtr];
               if(j==i) {
                  for (ib =0; ib<row_block; ++ib) A_ii[ib]=A_p->val[A_block*iPtr+ib+row_block*ib];
               } else {
                  for (ib =0; ib<row_block; ++ib) {
                     A_ij=A_p->val[A_block*iPtr+ib+row_block*ib];
                     if(A_ij< 0)  {
                        sum_all_neg[ib]+=A_ij;
                     } else {
                        sum_all_pos[ib]+=A_ij;
                     }
                  }

                  if (counter_C[j]>=0) {
                     /* is i strongly connected with j? We search for counter_C[j] in P[i,:] */
                     start_p=&(P_p->pattern->index[P_p->pattern->ptr[i]]);
                     where_p=(index_t*)bsearch(&(counter_C[j]), start_p,
                                             P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i],
                                             sizeof(index_t),
                                             util::comparIndex);
                     if (! (where_p == NULL) ) { /* yes i strongly connected with j */
                              offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
                              for (ib =0; ib<row_block; ++ib) {
                                 A_ij=A_p->val[A_block*iPtr+ib+row_block*ib];
                                 P_p->val[row_block*offset+ib]=A_ij; /* will be modified later */
                                 if (A_ij< 0)  {
                                     sum_strong_neg[ib]+=A_ij;
                                 } else {
                                    sum_strong_pos[ib]+=A_ij;
                                 }
                              }
                     }
                  }
               }
            }
            for (ib =0; ib<row_block; ++ib) {
               if(sum_strong_neg[ib]<0) {
                  alpha[ib]= sum_all_neg[ib]/sum_strong_neg[ib];
               } else {
                  alpha[ib]=0;
               }
               if(sum_strong_pos[ib]>0) {
                  beta[ib]= sum_all_pos[ib]/sum_strong_pos[ib];
               } else {
                  beta[ib]=0;
                  A_ii[ib]+=sum_all_pos[ib];
               }
               if ( A_ii[ib] > 0.) {
                  rtmp=(-1./A_ii[ib]);
                  alpha[ib]*=rtmp;
                  beta[ib]*=rtmp;
               }
            }

            for (iPtr=P_p->pattern->ptr[i];iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
               for (ib =0; ib<row_block; ++ib) {
                  A_ij=P_p->val[row_block*iPtr+ib];
                  if (A_ij > 0 ) {
                     P_p->val[row_block*iPtr+ib]*=beta[ib];
                  } else {
                     P_p->val[row_block*iPtr+ib]*=alpha[ib];
                  }
               }
            }
         }
      }/* end i loop */
      delete[] sum_all_neg;
      delete[] sum_all_pos;
      delete[] sum_strong_neg;
      delete[] sum_strong_pos;
      delete[] alpha;
      delete[] beta;
      delete[] A_ii;
   } /* end parallel region */
}

/*
    Classic Prolongation:
    -------------------

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise.
    If row i is not C, then P[i,j] = - 1/a[i] * ( A[i,k] + sum_{l} A[i,l]*A+[l,k]/B[i,k])
             where the summation over l is considering columns which are strongly connected
             to i (l in S[i]) and not in C (counter_C[l]<0) and

                B[i,k]=sum_{m in S_i and in C} A+[k,m]
                a[i]=A[i,i]+sum{l not strongly connected to i} A[i,l]

            A+[i,k]=A[i,k] if sign(A[i,k])==sign(A[i,i])  or 0 otherwise.

*/
void Preconditioner_LocalAMG_setClassicProlongation(SparseMatrix_ptr P_p,
        SparseMatrix_ptr A_p, const index_t* offset_S, const dim_t* degree_S,
        const index_t* S, const index_t* counter_C)
{
   dim_t i, q;
   const dim_t n =A_p->numRows;
   double *D_s=NULL;
   index_t *D_s_offset=NULL, iPtr, iPtr_j;
   const dim_t ll = util::iMax(n, degree_S);
   const index_t *ptr_main_A = A_p->borrowMainDiagonalPointer();

   #pragma omp parallel  private(D_s, D_s_offset, iPtr, q, iPtr_j)
   {
        D_s=new double[ll];
        D_s_offset=new index_t[ll];


        #pragma omp for private(i) schedule(static)
        for (i=0;i<n;++i) {
            if (counter_C[i]>=0) {
                P_p->val[P_p->pattern->ptr[i]]=1.;  /* i is a C row */
            } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
               const index_t *start_s = &(S[offset_S[i]]);
               const index_t *start_p = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
               const dim_t degree_P_i   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
              /* this loop sums up the weak connections in A and creates a
               * list of the strong connected columns which are not in C
               * (=no interpolation nodes) */
              const double A_ii = A_p->val[ptr_main_A[i]];
              double a=A_ii;

              for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
                 const index_t j=A_p->pattern->index[iPtr];
                 const double A_ij=A_p->val[iPtr];
                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection? */
                    const index_t *where_s=(index_t*)bsearch(&j, start_s,degree_S[i],sizeof(index_t), util::comparIndex);
                    if (where_s == NULL) { /* weak connections are accumulated */
                        a+=A_ij;
                    } else {   /* yes i strongly connected with j */
                        if  (counter_C[j]>=0)  { /* j is an interpolation point : add A_ij into P */
                               const index_t *where_p=(index_t*)bsearch(&counter_C[j], start_p,degree_P_i, sizeof(index_t), util::comparIndex);
                               if (where_p == NULL)  {
                                   throw PasoException("Preconditioner_LocalAMG_setClassicProlongation: Interpolation point is missing.");
                               } else {
                                    const index_t offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
                                    P_p->val[offset]+=A_ij;
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
                               const index_t *start_p_j = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
                               const dim_t degree_P_j   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
                               double s=0.;
                               dim_t len_D_s=0;
                               for (iPtr_j=A_p->pattern->ptr[j];iPtr_j<A_p->pattern->ptr[j + 1]; ++iPtr_j) {
                                    const double A_jm=A_p->val[iPtr_j];
                                    const index_t m=A_p->pattern->index[iPtr_j];
                                    /* is m an interpolation point? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_j,degree_P_j, sizeof(index_t), util::comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = P_p->pattern->ptr[i]+ (index_t)(where_p_m-start_p_j);
                                         if (!util::samesign(A_ii,A_jm)) {
                                              D_s[len_D_s]=A_jm;
                                         } else {
                                              D_s[len_D_s]=0.;
                                         }
                                         D_s_offset[len_D_s]=offset_m;
                                         len_D_s++;
                                    }
                               }
                               for (q=0;q<len_D_s;++q) s+=D_s[q];
                               if (std::abs(s)>0) {
                                   s=A_ij/s;
                                   for (q=0;q<len_D_s;++q) {
                                        P_p->val[D_s_offset[q]]+=s*D_s[q];
                                   }
                               } else {
                                   a+=A_ij;
                               }
                          }
                     }
                 }
              }  /* i has been processed, now we need to do some rescaling */
              if (std::abs(a)>0.) {
                   a=-1./a;
                   for (iPtr=P_p->pattern->ptr[i]; iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
                        P_p->val[iPtr]*=a;
                   }
              }
          }
        }  /* end of row i loop */
        delete[] D_s;
        delete[] D_s_offset;
     }    /* end of parallel region */
}

void Preconditioner_LocalAMG_setClassicProlongation_Block(SparseMatrix_ptr P_p,
        SparseMatrix_ptr A_p, const index_t* offset_S, const dim_t* degree_S,
        const index_t* S, const index_t* counter_C)
{
   dim_t i, q, ib;
   const dim_t row_block=A_p->row_block_size;
   const dim_t A_block = A_p->block_size;
   const dim_t n =A_p->numRows;
   double *D_s=NULL;
   index_t *D_s_offset=NULL, iPtr, iPtr_j;
   const dim_t ll = util::iMax(n, degree_S);
   const index_t *ptr_main_A = A_p->borrowMainDiagonalPointer();

   #pragma omp parallel  private(D_s, D_s_offset, iPtr, q, iPtr_j,ib)
   {
        double *a=new double[row_block];
        D_s=new double[row_block*ll];
        D_s_offset=new index_t[row_block*ll];

        #pragma omp for private(i) schedule(static)
        for (i=0;i<n;++i) {
            if (counter_C[i]>=0) {
                const index_t offset = P_p->pattern->ptr[i];
                for (ib =0; ib<row_block; ++ib) P_p->val[row_block*offset+ib]=1.;  /* i is a C row */
            } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
               const index_t *start_s = &(S[offset_S[i]]);
               const index_t *start_p = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
               const dim_t degree_P_i   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
              /* this loop sums up the weak connections in A and creates a
               * list of the strong connected columns which are not in C
               * (=no interpolation nodes) */
              const double *A_ii = &(A_p->val[ptr_main_A[i]*A_block]);
              for (ib=0; ib<row_block; ib++) a[ib]=A_ii[(row_block+1)*ib];

              for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
                 const index_t j=A_p->pattern->index[iPtr];
                 const double* A_ij=&(A_p->val[iPtr*A_block]);

                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
                    const index_t *where_s=(index_t*)bsearch(&j, start_s,degree_S[i],sizeof(index_t), util::comparIndex);
                    if (where_s == NULL) { /* weak connections are accumulated */
                        for (ib=0; ib<row_block; ib++) a[ib]+=A_ij[(row_block+1)*ib];
                    } else {   /* yes i strongly connected with j */
                        if  (counter_C[j]>=0)  { /* j is an interpolation point : add A_ij into P */
                               const index_t *where_p=(index_t*)bsearch(&counter_C[j], start_p,degree_P_i, sizeof(index_t), util::comparIndex);
                               if (where_p == NULL)  {
                                   throw PasoException("Preconditioner_LocalAMG_setClassicProlongation_Block: Interpolation point is missing.");
                               } else {
                                    const index_t offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
                                    for (ib=0; ib<row_block; ib++) P_p->val[offset*row_block+ib] +=A_ij[(row_block+1)*ib];
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
                               const index_t *start_p_j = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
                               const dim_t degree_P_j   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
                               dim_t len_D_s=0;
                               for (iPtr_j=A_p->pattern->ptr[j];iPtr_j<A_p->pattern->ptr[j + 1]; ++iPtr_j) {
                                    const double* A_jm=&(A_p->val[iPtr_j*A_block]);
                                    const index_t m=A_p->pattern->index[iPtr_j];
                                    /* is m an interpolation point? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_j,degree_P_j, sizeof(index_t), util::comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = P_p->pattern->ptr[i]+ (index_t)(where_p_m-start_p_j);
                                         for (ib=0; ib<row_block; ib++) {
                                              if (!util::samesign(A_ii[(row_block+1)*ib],A_jm[(row_block+1)*ib]) ) {
                                                   D_s[len_D_s*row_block+ib]=A_jm[(row_block+1)*ib];
                                              } else {
                                                   D_s[len_D_s*row_block+ib]=0.;
                                              }
                                         }
                                         D_s_offset[len_D_s]=offset_m;
                                         len_D_s++;
                                    }
                               }
                               for (ib=0; ib<row_block; ib++) {
                                   double s=0;
                                   for (q=0;q<len_D_s;++q) s+=D_s[q*row_block+ib];

                                   if (std::abs(s)>0) {
                                       s=A_ij[(row_block+1)*ib]/s;
                                       for (q=0;q<len_D_s;++q) {
                                            P_p->val[D_s_offset[q]*row_block+ib]+=s*D_s[q*row_block+ib];
                                       }
                                   } else {
                                       a[ib]+=A_ij[(row_block+1)*ib];
                                   }
                               }
                          }
                     }
                 }
              }  /* i has been processed, now we need to do some rescaling */
              for (ib=0; ib<row_block; ib++) {
                   double a2=a[ib];
                   if (std::abs(a2)>0.) {
                        a2=-1./a2;
                        for (iPtr=P_p->pattern->ptr[i]; iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
                             P_p->val[iPtr*row_block+ib]*=a2;
                        }
                   }
              }
          }
        }  // end of row i loop
        delete[] D_s;
        delete[] D_s_offset;
        delete[] a;
    } // end of parallel region
}

} // namespace paso

