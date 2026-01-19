
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


/****************************************************************************/

/* Paso: RILU preconditioner with reordering                 */

/****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au                      */

/****************************************************************************/

#include "Paso.h"
#include "BlockOps.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

namespace paso {

void Solver_RILU_free(Solver_RILU* in)
{
    if (in!=NULL) {
        Solver_RILU_free(in->RILU_of_Schur);
        delete[] in->inv_A_FF;
        delete[] in->A_FF_pivot;
        delete[] in->rows_in_F;
        delete[] in->rows_in_C;
        delete[] in->mask_F;
        delete[] in->mask_C;
        delete[] in->x_F;
        delete[] in->b_F;
        delete[] in->x_C;
        delete[] in->b_C;
        delete in;
    }
}

/****************************************************************************/

/*   constructs the block-block factorization of

        [ A_FF A_FC ]
   A_p=
        [ A_CF A_FF ]

to

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FF ]
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ]


   where S=A_FF-ACF*invA_FF*A_FC within the shape of S

   then RILU is applied to S again until S becomes empty
*/
Solver_RILU* Solver_getRILU(SparseMatrix_ptr<double> A_p, bool verbose)
{
    Solver_RILU* out=NULL;
    dim_t n=A_p->numRows;
    dim_t n_block=A_p->row_block_size;
    index_t* mis_marker=NULL;
    index_t* counter=NULL;
    index_t iPtr,*index, *where_p;
    dim_t i,k;
    SparseMatrix_ptr<double> schur;
    double A11,A12,A13,A21,A22,A23,A31,A32,A33,D,time0=0,time1=0;/*,time2=0;*/

    mis_marker=new index_t[n];
    counter=new index_t[n];
    out=new Solver_RILU;
    out->RILU_of_Schur=NULL;
    out->inv_A_FF=NULL;
    out->A_FF_pivot=NULL;
    out->rows_in_F=NULL;
    out->rows_in_C=NULL;
    out->mask_F=NULL;
    out->mask_C=NULL;
    out->x_F=NULL;
    out->b_F=NULL;
    out->x_C=NULL;
    out->b_C=NULL;

    /* identify independent set of rows/columns */
    time0=escript::gettime();
    #pragma omp parallel for private(i) schedule(static)
    for (i=0;i<n;++i) mis_marker[i]=-1;
    A_p->pattern->mis(mis_marker);
    /*time2=escript::gettime()-time0;*/
    #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < n; ++i) counter[i]=mis_marker[i];
    out->n=n;
    out->n_block=n_block;
    out->n_F=util::cumsum(n,counter);
    out->mask_F=new index_t[n];
    out->rows_in_F=new index_t[out->n_F];
    out->inv_A_FF=new double[n_block*n_block*out->n_F];
    out->A_FF_pivot=NULL; /* later use for block size>3 */
#pragma omp parallel
    {
      /* creates an index for F from mask */
      #pragma omp for private(i) schedule(static)
      for (i = 0; i < out->n_F; ++i) out->rows_in_F[i]=-1;
      #pragma omp for private(i) schedule(static)
      for (i = 0; i < n; ++i) {
         if  (mis_marker[i]) {
                out->rows_in_F[counter[i]]=i;
                out->mask_F[i]=counter[i];
         } else {
                out->mask_F[i]=-1;
         }
      }
      #pragma omp for private(i, where_p,iPtr,A11,A12,A13,A21,A22,A23,A31,A32,A33,D,index) schedule(static)
      for (i = 0; i < out->n_F; i++) {
        /* find main diagonal */
        iPtr=A_p->pattern->ptr[out->rows_in_F[i]];
        index=&(A_p->pattern->index[iPtr]);
        where_p=(index_t*)bsearch(&out->rows_in_F[i],
                                index,
                                A_p->pattern->ptr[out->rows_in_F[i] + 1]-A_p->pattern->ptr[out->rows_in_F[i]],
                                sizeof(index_t),
                                util::comparIndex);
        if (where_p==NULL) {
            throw PasoException("Solver_getRILU: main diagonal element missing.");
        } else {
            iPtr+=(index_t)(where_p-index);
            /* get inverse of A_FF block: */
            if (n_block==1) {
               if (std::abs(A_p->val[iPtr])>0.) {
                    out->inv_A_FF[i]=1./A_p->val[iPtr];
               } else {
                    throw PasoException("Solver_getRILU: Break-down in RILU decomposition: non-regular main diagonal block.");
               }
            } else if (n_block==2) {
               A11=A_p->val[iPtr*4];
               A21=A_p->val[iPtr*4+1];
               A12=A_p->val[iPtr*4+2];
               A22=A_p->val[iPtr*4+3];
               D = A11*A22-A12*A21;
               if (std::abs(D) > 0 ){
                    D=1./D;
                    out->inv_A_FF[i*4]= A22*D;
                    out->inv_A_FF[i*4+1]=-A21*D;
                    out->inv_A_FF[i*4+2]=-A12*D;
                    out->inv_A_FF[i*4+3]= A11*D;
               } else {
                    throw PasoException("Solver_getRILU: Break-down in RILU decomposition: non-regular main diagonal block.");
               }
            } else if (n_block==3) {
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
               if (std::abs(D) > 0) {
                    D=1./D;
                    out->inv_A_FF[i*9  ]=(A22*A33-A23*A32)*D;
                    out->inv_A_FF[i*9+1]=(A31*A23-A21*A33)*D;
                    out->inv_A_FF[i*9+2]=(A21*A32-A31*A22)*D;
                    out->inv_A_FF[i*9+3]=(A13*A32-A12*A33)*D;
                    out->inv_A_FF[i*9+4]=(A11*A33-A31*A13)*D;
                    out->inv_A_FF[i*9+5]=(A12*A31-A11*A32)*D;
                    out->inv_A_FF[i*9+6]=(A12*A23-A13*A22)*D;
                    out->inv_A_FF[i*9+7]=(A13*A21-A11*A23)*D;
                    out->inv_A_FF[i*9+8]=(A11*A22-A12*A21)*D;
               } else {
                    throw PasoException("Solver_getRILU: Break-down in RILU decomposition: non-regular main diagonal block.");
               }
           }
        }
      }
    } /* end parallel region */

    // if there are no nodes in the coarse level there is no more
    // work to do
    out->n_C=n-out->n_F;
    if (out->n_C > 0) {
        out->rows_in_C = new index_t[out->n_C];
        out->mask_C = new index_t[n];
        /* creates an index for C from mask */
        #pragma omp parallel for private(i) schedule(static)
        for (i = 0; i < n; ++i) counter[i]=! mis_marker[i];
        util::cumsum(n,counter);
        #pragma omp parallel
        {
              #pragma omp for private(i) schedule(static)
              for (i = 0; i < out->n_C; ++i) out->rows_in_C[i]=-1;
              #pragma omp for private(i) schedule(static)
              for (i = 0; i < n; ++i) {
                 if  (! mis_marker[i]) {
                    out->rows_in_C[counter[i]]=i;
                    out->mask_C[i]=counter[i];
                 } else {
                    out->mask_C[i]=-1;
                 }
              }
        } /* end parallel region */
        /* get A_CF block: */
        out->A_CF=A_p->getSubmatrix(out->n_C, out->n_F, out->rows_in_C, out->mask_F);
        /* get A_FC block: */
        out->A_FC=A_p->getSubmatrix(out->n_F, out->n_C, out->rows_in_F, out->mask_C);
        /* get A_FF block: */
        schur = A_p->getSubmatrix(out->n_C, out->n_C, out->rows_in_C, out->mask_C);
        time0=escript::gettime()-time0;
        time1=escript::gettime();
        /* update A_CC block to get Schur complement and then apply RILU to it */
        Solver_updateIncompleteSchurComplement(schur, out->A_CF, out->inv_A_FF, out->A_FF_pivot, out->A_FC);
        time1=escript::gettime()-time1;
        out->RILU_of_Schur = Solver_getRILU(schur, verbose);
        schur.reset();
        /* allocate work arrays for RILU application */
        out->x_F=new double[n_block*out->n_F];
        out->b_F=new double[n_block*out->n_F];
        out->x_C=new double[n_block*out->n_C];
        out->b_C=new double[n_block*out->n_C];
#pragma omp parallel
        {
#pragma omp for private(i,k) schedule(static)
            for (i = 0; i < out->n_F; ++i) {
                for (k=0; k<n_block;++k) {
                    out->x_F[i*n_block+k]=0.;
                    out->b_F[i*n_block+k]=0.;
                }
            }
#pragma omp for private(i,k) schedule(static)
            for (i = 0; i < out->n_C; ++i) {
                for (k=0; k<n_block;++k) {
                    out->x_C[i*n_block+k]=0.;
                    out->b_C[i*n_block+k]=0.;
                }
            }
        } // end parallel region
    }
    delete[] mis_marker;
    delete[] counter;
    //if (verbose) {
    //    printf("RILU: %d unknowns eliminated. %d left.\n",out->n_F,n-out->n_F);
    //    if (out->n_C>0) {
    //        printf("timing: RILU: MIS/reordering/elimination : %e/%e/%e\n",time2,time0,time1);
    //    } else {
    //        printf("timing: RILU: MIS: %e\n",time2);
    //    }
    //}
    return out;
}

/****************************************************************************/

/* apply RILU precondition b-> x

     in fact it solves

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FF ]  [ x_F ]  = [b_F]
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ]  [ x_C ]  = [b_C]

 in the form

   b->[b_F,b_C]
   x_F=invA_FF*b_F
   b_C=b_C-A_CF*x_F
   x_C=RILU(b_C)
   b_F=b_F-A_FC*x_C
   x_F=invA_FF*b_F
   x<-[x_F,x_C]

 Should be called within a parallel region.
 Barrier synchronization should be performed to make sure that the input
 vector is available.
*/
void Solver_solveRILU(Solver_RILU* rilu, double* x, double* b)
{
    dim_t i,k;
    dim_t n_block=rilu->n_block;

    if (rilu->n_C==0) {
        /* x=invA_FF*b  */
        util::copy(n_block*rilu->n_F, x, b);
        BlockOps_solveAll(n_block,rilu->n_F,rilu->inv_A_FF,rilu->A_FF_pivot,x);
    } else {
        /* b->[b_F,b_C]     */
        if (n_block==1) {
            #pragma omp parallel for private(i) schedule(static)
            for (i=0;i<rilu->n_F;++i) rilu->b_F[i]=b[rilu->rows_in_F[i]];
            #pragma omp parallel for private(i) schedule(static)
            for (i=0;i<rilu->n_C;++i) rilu->b_C[i]=b[rilu->rows_in_C[i]];
        } else {
            #pragma omp parallel for private(i,k) schedule(static)
            for (i=0; i<rilu->n_F; i++)
                for (k=0; k<n_block; k++)
                    rilu->b_F[rilu->n_block*i+k]=b[n_block*rilu->rows_in_F[i]+k];
            #pragma omp parallel for private(i,k) schedule(static)
            for (i=0; i<rilu->n_C; i++)
                for (k=0; k<n_block; k++)
                    rilu->b_C[rilu->n_block*i+k]=b[n_block*rilu->rows_in_C[i]+k];
        }
        /* x_F=invA_FF*b_F  */
        util::copy(n_block*rilu->n_F, rilu->x_F,rilu->b_F);
        BlockOps_solveAll(n_block,rilu->n_F,rilu->inv_A_FF,rilu->A_FF_pivot,rilu->x_F);
        /* b_C=b_C-A_CF*x_F */
        SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,rilu->A_CF,rilu->x_F,1.,rilu->b_C);
        /* x_C=RILU(b_C)     */
        Solver_solveRILU(rilu->RILU_of_Schur,rilu->x_C,rilu->b_C);
        /* b_F=b_F-A_FC*x_C */
        SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,rilu->A_FC,rilu->x_C,1.,rilu->b_F);
        /* x_F=invA_FF*b_F  */
        util::copy(n_block*rilu->n_F, rilu->x_F,rilu->b_F);
        BlockOps_solveAll(n_block,rilu->n_F,rilu->inv_A_FF,rilu->A_FF_pivot,rilu->x_F);
        /* x<-[x_F,x_C]     */
        if (n_block==1) {
           #pragma omp parallel for private(i) schedule(static)
           for (i=0; i<rilu->n; i++) {
              if (rilu->mask_C[i] > -1) {
                  x[i]=rilu->x_C[rilu->mask_C[i]];
              } else {
                  x[i]=rilu->x_F[rilu->mask_F[i]];
              }
           }
        } else {
            #pragma omp parallel for private(i,k) schedule(static)
            for (i=0; i<rilu->n; i++) {
                if (rilu->mask_C[i] > -1) {
                     for (k=0; k<n_block; k++)
                         x[n_block*i+k]=rilu->x_C[n_block*rilu->mask_C[i]+k];
                } else {
                     for (k=0; k<n_block; k++)
                         x[n_block*i+k]=rilu->x_F[n_block*rilu->mask_F[i]+k];
                }
            }
        }
    }
}

} // namespace paso

