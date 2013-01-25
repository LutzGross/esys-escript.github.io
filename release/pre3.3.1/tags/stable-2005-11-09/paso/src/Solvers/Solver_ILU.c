/* $Id$ */

/**************************************************************/

/* Paso: ILU preconditioner with reordering                 */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: gross@access.edu.au                                */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "Util.h"

/**************************************************************/

/* free all memory used by ILU                                */

void Paso_Solver_ILU_free(Paso_Solver_ILU * in) {
     if (in!=NULL) {
        Paso_Solver_ILU_free(in->ILU_of_Schur);
        MEMFREE(in->inv_A_FF);
        MEMFREE(in->A_FF_pivot);
        Paso_SystemMatrix_dealloc(in->A_FC);
        Paso_SystemMatrix_dealloc(in->A_CF);
        MEMFREE(in->rows_in_F);
        MEMFREE(in->rows_in_C);
        MEMFREE(in->mask_F);
        MEMFREE(in->mask_C);
        MEMFREE(in->x_F);
        MEMFREE(in->b_F);
        MEMFREE(in->x_C);
        MEMFREE(in->b_C);
        MEMFREE(in);
     }
}

/**************************************************************/

/*   constructs the block-block factorization of 

        [ A_FF A_FC ]
   A_p= 
        [ A_CF A_FF ]

to 

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FF ]  
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ] 


   where S=A_FF-ACF*invA_FF*A_FC within the shape of S

   then ILU is applied to S again until S becomes empty 

*/
Paso_Solver_ILU* Paso_Solver_getILU(Paso_SystemMatrix * A_p,bool_t verbose) {
  Paso_Solver_ILU* out=NULL;
  dim_t n=A_p->num_rows;
  dim_t n_block=A_p->row_block_size;
  index_t* mis_marker=NULL;  
  index_t* counter=NULL;  
  index_t iPtr,*index, *where_p;
  dim_t i,k;
  Paso_SystemMatrix * schur=NULL;
  double A11,A12,A13,A21,A22,A23,A31,A32,A33,D,time0,time1,time2;
   

  /* identify independend set of rows/columns */
  mis_marker=TMPMEMALLOC(n,index_t);
  counter=TMPMEMALLOC(n,index_t);
  out=MEMALLOC(1,Paso_Solver_ILU);
  out->ILU_of_Schur=NULL;
  out->inv_A_FF=NULL;
  out->A_FF_pivot=NULL;
  out->A_FC=NULL;
  out->A_CF=NULL;
  out->rows_in_F=NULL;
  out->rows_in_C=NULL;
  out->mask_F=NULL;
  out->mask_C=NULL;
  out->x_F=NULL;
  out->b_F=NULL;
  out->x_C=NULL;
  out->b_C=NULL;

  if ( !(Paso_checkPtr(mis_marker) || Paso_checkPtr(out) || Paso_checkPtr(counter) ) ) {
     /* identify independend set of rows/columns */
     time0=Paso_timer();
     Paso_SystemMatrixPattern_mis(A_p->pattern,mis_marker);
     time2=Paso_timer()-time0;
     if (Paso_noError()) {
        #pragma omp parallel for private(i) schedule(static)
        for (i = 0; i < n; ++i) counter[i]=mis_marker[i];
        out->n=n;
        out->n_block=n_block;
        out->n_F=Paso_Util_cumsum(n,counter);
        out->mask_F=MEMALLOC(n,index_t);
        out->rows_in_F=MEMALLOC(out->n_F,index_t);
        out->inv_A_FF=MEMALLOC(n_block*n_block*out->n_F,double);
        out->A_FF_pivot=NULL; /* later use for block size>3 */
        if (! (Paso_checkPtr(out->mask_F) || Paso_checkPtr(out->inv_A_FF) || Paso_checkPtr(out->rows_in_F) ) ) {
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
                                        Paso_comparIndex);
                if (where_p==NULL) {
                    Paso_setError(VALUE_ERROR, "__FILE__: main diagonal element missing.");
                } else {
                    iPtr+=(index_t)(where_p-index);
                    /* get inverse of A_FF block: */
                    if (n_block==1) {
                       if (ABS(A_p->val[iPtr])>0.) {
                            out->inv_A_FF[i]=1./A_p->val[iPtr];
                       } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "__FILE__: Break-down in ILU decomposition: non-regular main diagonal block.");
                       }
                    } else if (n_block==2) {
                       A11=A_p->val[iPtr*4];
                       A21=A_p->val[iPtr*4+1];
                       A12=A_p->val[iPtr*4+2];
                       A22=A_p->val[iPtr*4+3];
                       D = A11*A22-A12*A21;
                       if (ABS(D) > 0 ){
                            D=1./D;
                            out->inv_A_FF[i*4]= A22*D;
                            out->inv_A_FF[i*4+1]=-A21*D;
                            out->inv_A_FF[i*4+2]=-A12*D;
                            out->inv_A_FF[i*4+3]= A11*D;
                       } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "__FILE__:Break-down in ILU decomposition: non-regular main diagonal block.");
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
                       if (ABS(D) > 0 ){
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
                            Paso_setError(ZERO_DIVISION_ERROR, "__FILE__:Break-down in ILU decomposition: non-regular main diagonal block.");
                       }
                   }
                }
              }
           } /* end parallel region */

           if( Paso_noError()) {
              /* if there are no nodes in the coarse level there is no more work to do */
              out->n_C=n-out->n_F;
              if (out->n_C>0) {
                   out->rows_in_C=MEMALLOC(out->n_C,index_t);
                   out->mask_C=MEMALLOC(n,index_t);
                   if (! (Paso_checkPtr(out->mask_C) || Paso_checkPtr(out->rows_in_C) ) ) {
                       /* creates an index for C from mask */
                       #pragma omp parallel for private(i) schedule(static)
                       for (i = 0; i < n; ++i) counter[i]=! mis_marker[i];
                       Paso_Util_cumsum(n,counter);
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
                      out->A_CF=Paso_SystemMatrix_getSubmatrix(A_p,out->n_C,out->rows_in_C,out->mask_F);
                      if (Paso_noError()) {
                         /* get A_FC block: */
                         out->A_FC=Paso_SystemMatrix_getSubmatrix(A_p,out->n_F,out->rows_in_F,out->mask_C);
                         /* get A_FF block: */
                         if (Paso_noError()) {
                            schur=Paso_SystemMatrix_getSubmatrix(A_p,out->n_C,out->rows_in_C,out->mask_C);
                            time0=Paso_timer()-time0;
                            if (Paso_noError()) {
                                time1=Paso_timer();
                                /* update A_CC block to get Schur complement and then apply ILU to it */
                                Paso_Solver_updateIncompleteSchurComplement(schur,out->A_CF,out->inv_A_FF,out->A_FF_pivot,out->A_FC);
                                time1=Paso_timer()-time1;
                                out->ILU_of_Schur=Paso_Solver_getILU(schur,verbose);
                                Paso_SystemMatrix_dealloc(schur);
                            }
                            /* allocate work arrays for ILU application */
                            if (Paso_noError()) {
                              out->x_F=MEMALLOC(n_block*out->n_F,double);
                              out->b_F=MEMALLOC(n_block*out->n_F,double);
                              out->x_C=MEMALLOC(n_block*out->n_C,double);
                              out->b_C=MEMALLOC(n_block*out->n_C,double);
                              if (! (Paso_checkPtr(out->x_F) || Paso_checkPtr(out->b_F) || Paso_checkPtr(out->x_C) || Paso_checkPtr(out->b_C) ) ) {
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
                                  } /* end parallel region */
                              }
                            }
                         }
                     }
                 }
              }
           }
        }
     }
  }
  TMPMEMFREE(mis_marker);
  TMPMEMFREE(counter);
  if (Paso_noError()) {
      if (verbose) {
         printf("ILU: %d unknowns eliminated. %d left.\n",out->n_F,n-out->n_F);
         if (out->n_C>0) {
            printf("timing: ILU: MIS/reordering/elemination : %e/%e/%e\n",time2,time0,time1);
         } else {
            printf("timing: ILU: MIS: %e\n",time2); 
         }
     }
     return out;
  } else  {
     Paso_Solver_ILU_free(out);
     return NULL;
  }
}

/**************************************************************/

/* apply ILU precondition b-> x                               

     in fact it solves 

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FF ]  [ x_F ]  = [b_F]
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ]  [ x_C ]  = [b_C]

 in the form 

   b->[b_F,b_C] 
   x_F=invA_FF*b_F
   b_C=b_C-A_CF*x_F
   x_C=ILU(b_C)
   b_F=b_F-A_FC*x_C
   x_F=invA_FF*b_F
   x<-[x_F,x_C]

 should be called within a parallel region                                              
 barrier synconization should be performed to make sure that the input vector available 

*/

void Paso_Solver_solveILU(Paso_Solver_ILU * ilu, double * x, double * b) {
     dim_t i,k;
     dim_t n_block=ilu->n_block;
     
     if (ilu->n_C==0) {
        /* x=invA_FF*b  */
        Paso_Solver_applyBlockDiagonalMatrix(n_block,ilu->n_F,ilu->inv_A_FF,ilu->A_FF_pivot,x,b);
     } else {
        /* b->[b_F,b_C]     */
        if (n_block==1) {
           #pragma omp for private(i) schedule(static)
           for (i=0;i<ilu->n_F;++i) ilu->b_F[i]=b[ilu->rows_in_F[i]];
           #pragma omp for private(i) schedule(static)
           for (i=0;i<ilu->n_C;++i) ilu->b_C[i]=b[ilu->rows_in_C[i]];
        } else {
           #pragma omp for private(i,k) schedule(static)
           for (i=0;i<ilu->n_F;++i) 
                 for (k=0;k<n_block;k++) ilu->b_F[ilu->n_block*i+k]=b[n_block*ilu->rows_in_F[i]+k];
           #pragma omp for private(i,k) schedule(static)
           for (i=0;i<ilu->n_C;++i) 
                 for (k=0;k<n_block;k++) ilu->b_C[ilu->n_block*i+k]=b[n_block*ilu->rows_in_C[i]+k];
        }
        /* x_F=invA_FF*b_F  */
        Paso_Solver_applyBlockDiagonalMatrix(n_block,ilu->n_F,ilu->inv_A_FF,ilu->A_FF_pivot,ilu->x_F,ilu->b_F);
        /* b_C=b_C-A_CF*x_F */
        Paso_SystemMatrix_MatrixVector(-1.,ilu->A_CF,ilu->x_F,1.,ilu->b_C);
        /* x_C=ILU(b_C)     */
        Paso_Solver_solveILU(ilu->ILU_of_Schur,ilu->x_C,ilu->b_C);
        /* b_F=b_F-A_FC*x_C */
        Paso_SystemMatrix_MatrixVector(-1.,ilu->A_FC,ilu->x_C,1.,ilu->b_F);
        /* x_F=invA_FF*b_F  */
        Paso_Solver_applyBlockDiagonalMatrix(n_block,ilu->n_F,ilu->inv_A_FF,ilu->A_FF_pivot,ilu->x_F,ilu->b_F);
        /* x<-[x_F,x_C]     */
        if (n_block==1) {
           #pragma omp for private(i) schedule(static)
           for (i=0;i<ilu->n;++i) {
              if (ilu->mask_C[i]>-1) {
                  x[i]=ilu->x_C[ilu->mask_C[i]];
              } else {
                  x[i]=ilu->x_F[ilu->mask_F[i]];
              }
           }
        } else {
           #pragma omp for private(i,k) schedule(static)
           for (i=0;i<ilu->n;++i) {
                 if (ilu->mask_C[i]>-1) {
                     for (k=0;k<n_block;k++) x[n_block*i+k]=ilu->x_C[n_block*ilu->mask_C[i]+k];
                 } else {
                     for (k=0;k<n_block;k++) x[n_block*i+k]=ilu->x_F[n_block*ilu->mask_F[i]+k];
                 }
           }
        }
        /* all done */
     }
     return;
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
