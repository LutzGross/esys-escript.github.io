/* $Id$ */

/**************************************************************/

/* Finley: ILU preconditioner with reordering                 */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: gross@access.edu.au                                */

/**************************************************************/

#include "Solver.h"
#include "Util.h"

/**************************************************************/

/* free all memory used by ILU                                */

void Finley_Solver_ILU_free(Finley_Solver_ILU * in) {
     if (in!=NULL) {
        Finley_Solver_ILU_free(in->ILU_of_Schur);
        MEMFREE(in->inv_A_FF);
        MEMFREE(in->A_FF_pivot);
        Finley_SystemMatrix_dealloc(in->A_FC);
        Finley_SystemMatrix_dealloc(in->A_CF);
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
Finley_Solver_ILU* Finley_Solver_getILU(Finley_SystemMatrix * A_p,int verbose) {
  Finley_Solver_ILU* out=NULL;
  maybelong n=A_p->num_rows;
  maybelong n_block=A_p->row_block_size;
  maybelong* mis_marker=NULL;  
  maybelong* counter=NULL;  
  maybelong i,iPtr, *index, *where_p,k;
  Finley_SystemMatrix * schur=NULL;
  double A11,A12,A13,A21,A22,A23,A31,A32,A33,D,time0,time1,time2;
   

  /* identify independend set of rows/columns */
  mis_marker=TMPMEMALLOC(n,maybelong);
  counter=TMPMEMALLOC(n,maybelong);
  out=MEMALLOC(1,Finley_Solver_ILU);
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

  if ( !(Finley_checkPtr(mis_marker) || Finley_checkPtr(out) || Finley_checkPtr(counter) ) ) {
     /* identify independend set of rows/columns */
     time0=Finley_timer();
     Finley_SystemMatrixPattern_mis(A_p->pattern,mis_marker);
     time2=Finley_timer()-time0;
     if (Finley_ErrorCode==NO_ERROR) {
        #pragma omp parallel for private(i) schedule(static)
        for (i = 0; i < n; ++i) counter[i]=mis_marker[i];
        out->n=n;
        out->n_block=n_block;
        out->n_F=Finley_Util_cumsum(n,counter);
        if (verbose) {
           printf("ILU: number of vertices to be eliminated = %d\n",out->n_F);
           printf("ILU: number of vertices in coarse level  = %d\n",n-out->n_F);
        }
        out->mask_F=MEMALLOC(n,maybelong);
        out->rows_in_F=MEMALLOC(out->n_F,maybelong);
        out->inv_A_FF=MEMALLOC(n_block*n_block*out->n_F,double);
        out->A_FF_pivot=NULL; /* later use for block size>3 */
        if (! (Finley_checkPtr(out->mask_F) || Finley_checkPtr(out->inv_A_FF) || Finley_checkPtr(out->rows_in_F) ) ) {
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
                where_p=(maybelong*)bsearch(&out->rows_in_F[i],
                                        index,
                                        A_p->pattern->ptr[out->rows_in_F[i] + 1]-A_p->pattern->ptr[out->rows_in_F[i]],
                                        sizeof(maybelong),
                                        Finley_comparIndex);
                if (where_p==NULL) {
                    Finley_ErrorCode = VALUE_ERROR;
                    sprintf(Finley_ErrorMsg, "no main diagonal in row %d",i);
                } else {
                    iPtr+=(maybelong)(where_p-index);
                    /* get inverse of A_FF block: */
                    if (n_block==1) {
                       if (ABS(A_p->val[iPtr])>0.) {
                            out->inv_A_FF[i]=1./A_p->val[iPtr];
                       } else {
                            Finley_ErrorCode = ZERO_DIVISION_ERROR;
                            sprintf(Finley_ErrorMsg, "Break-down in ILU decomposition: non-regular main diagonal block in row %d",i);
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
                            Finley_ErrorCode = ZERO_DIVISION_ERROR;
                            sprintf(Finley_ErrorMsg, "Break-down in ILU decomposition: non-regular main diagonal block in row %d",i);
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
                            Finley_ErrorCode = ZERO_DIVISION_ERROR;
                            sprintf(Finley_ErrorMsg, "Break-down in ILU decomposition: non-regular main diagonal block in row %d",i);
                       }
                   }
                }
              }
           } /* end parallel region */

           if( Finley_ErrorCode == NO_ERROR) {
              /* if there are no nodes in the coarse level there is no more work to do */
              out->n_C=n-out->n_F;
              if (out->n_C>0) {
                   out->rows_in_C=MEMALLOC(out->n_C,maybelong);
                   out->mask_C=MEMALLOC(n,maybelong);
                   if (! (Finley_checkPtr(out->mask_C) || Finley_checkPtr(out->rows_in_C) ) ) {
                       /* creates an index for C from mask */
                       #pragma omp parallel for private(i) schedule(static)
                       for (i = 0; i < n; ++i) counter[i]=! mis_marker[i];
                       Finley_Util_cumsum(n,counter);
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
                      out->A_CF=Finley_SystemMatrix_getSubmatrix(A_p,out->n_C,out->rows_in_C,out->mask_F);
                      if (Finley_ErrorCode==NO_ERROR) {
                         /* get A_FC block: */
                         out->A_FC=Finley_SystemMatrix_getSubmatrix(A_p,out->n_F,out->rows_in_F,out->mask_C);
                         /* get A_FF block: */
                         if (Finley_ErrorCode==NO_ERROR) {
                            schur=Finley_SystemMatrix_getSubmatrix(A_p,out->n_C,out->rows_in_C,out->mask_C);
                            time0=Finley_timer()-time0;
                            if (Finley_ErrorCode==NO_ERROR) {
                                time1=Finley_timer();
                                /* update A_CC block to get Schur complement and then apply ILU to it */
                                Finley_Solver_updateIncompleteSchurComplement(schur,out->A_CF,out->inv_A_FF,out->A_FF_pivot,out->A_FC);
                                time1=Finley_timer()-time1;
                                out->ILU_of_Schur=Finley_Solver_getILU(schur,verbose);
                                Finley_SystemMatrix_dealloc(schur);
                            }
                            /* allocate work arrays for ILU application */
                            if (Finley_ErrorCode==NO_ERROR) {
                              out->x_F=MEMALLOC(n_block*out->n_F,double);
                              out->b_F=MEMALLOC(n_block*out->n_F,double);
                              out->x_C=MEMALLOC(n_block*out->n_C,double);
                              out->b_C=MEMALLOC(n_block*out->n_C,double);
                              if (! (Finley_checkPtr(out->x_F) || Finley_checkPtr(out->b_F) || Finley_checkPtr(out->x_C) || Finley_checkPtr(out->b_C) ) ) {
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
  if (Finley_ErrorCode==NO_ERROR) {
      if (verbose) {
         /* if (out->n_C>0) {
            printf("timing: ILU: reordering : %e\n",time0);
            printf("timing: ILU: elemination : %e\n",time1);
         }
         printf("timing: ILU: mis %e\n",time2); */
     }
     return out;
  } else  {
     Finley_Solver_ILU_free(out);
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

void Finley_Solver_solveILU(Finley_Solver_ILU * ilu, double * x, double * b) {
     maybelong i,k;
     maybelong n_block=ilu->n_block;
     
     if (ilu->n_C==0) {
        /* x=invA_FF*b  */
        Finley_Solver_applyBlockDiagonalMatrix(n_block,ilu->n_F,ilu->inv_A_FF,ilu->A_FF_pivot,x,b);
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
        Finley_Solver_applyBlockDiagonalMatrix(n_block,ilu->n_F,ilu->inv_A_FF,ilu->A_FF_pivot,ilu->x_F,ilu->b_F);
        /* b_C=b_C-A_CF*x_F */
        Finley_RawScaledSystemMatrixVector(-1.,ilu->A_CF,ilu->x_F,1.,ilu->b_C);
        /* x_C=ILU(b_C)     */
        Finley_Solver_solveILU(ilu->ILU_of_Schur,ilu->x_C,ilu->b_C);
        /* b_F=b_F-A_FC*x_C */
        Finley_RawScaledSystemMatrixVector(-1.,ilu->A_FC,ilu->x_C,1.,ilu->b_F);
        /* x_F=invA_FF*b_F  */
        Finley_Solver_applyBlockDiagonalMatrix(n_block,ilu->n_F,ilu->inv_A_FF,ilu->A_FF_pivot,ilu->x_F,ilu->b_F);
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
