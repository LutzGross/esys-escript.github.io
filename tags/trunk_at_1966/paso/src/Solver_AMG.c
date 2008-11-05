
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: AMG preconditioner with reordering                 */

/**************************************************************/

/* Author: artak@access.edu.au                                */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "PasoUtil.h"
#include "UMFPACK.h"
#include "Pattern_coupling.h"

/**************************************************************/

/* free all memory used by AMG                                */

void Paso_Solver_AMG_free(Paso_Solver_AMG * in) {
     if (in!=NULL) {
        Paso_Solver_AMG_free(in->AMG_of_Schur);
        Paso_Solver_GS_free(in->GS);
        MEMFREE(in->inv_A_FF);
        MEMFREE(in->A_FF_pivot);
        Paso_SparseMatrix_free(in->A_FC);
        Paso_SparseMatrix_free(in->A_CF);
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

   then AMG is applied to S again until S becomes empty 

*/
Paso_Solver_AMG* Paso_Solver_getAMG(Paso_SparseMatrix *A_p,bool_t verbose,dim_t level) {
  Paso_Solver_AMG* out=NULL;
  dim_t n=A_p->numRows;
  dim_t n_block=A_p->row_block_size;
  index_t* mis_marker=NULL;  
  index_t* counter=NULL;
  index_t iPtr,*index, *where_p, iPtr_s;
  dim_t i,k,j,j0;
  Paso_SparseMatrix * schur=NULL;
  Paso_SparseMatrix * schur_withFillIn=NULL;
  double time0,time1,time2,S;
  
  /* identify independend set of rows/columns */
  mis_marker=TMPMEMALLOC(n,index_t);
  counter=TMPMEMALLOC(n,index_t);
  out=MEMALLOC(1,Paso_Solver_AMG);
  out->AMG_of_Schur=NULL;
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
  out->A=Paso_SparseMatrix_getReference(A_p);
  out->GS=Paso_Solver_getGS(A_p,verbose);
  out->GS->sweeps=4;
  out->level=level;
  
  if ( !(Paso_checkPtr(mis_marker) || Paso_checkPtr(out) || Paso_checkPtr(counter) ) ) {
     /* identify independend set of rows/columns */
     time0=Paso_timer();
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;++i) mis_marker[i]=-1;
     /*Paso_Pattern_RS(A_p,mis_marker,0.25);*/
     /*Paso_Pattern_Aggregiation(A_p,mis_marker,0.5);*/
     Paso_Pattern_coup(A_p,mis_marker,0.05);
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
              /* Compute row-sum for getting rs(A_FF)*/
              #pragma omp for private(i,iPtr) schedule(static)
              for (i = 0; i < out->n_F; ++i) {
                out->inv_A_FF[i]=0;
                for (iPtr=A_p->pattern->ptr[out->rows_in_F[i]];iPtr<A_p->pattern->ptr[out->rows_in_F[i] + 1]; ++iPtr) {
                 out->inv_A_FF[i]+=A_p->val[iPtr];
                }
              }
              
/*          for (i=0;i<out->n_F;++i) {
          fprintf(stderr,"GET INV_A_FF %g %d of %d (%d)\n",out->inv_A_FF[i],i,out->rows_in_F[i],(out->inv_A_FF[i]==0));
          }
          fprintf(stderr,"END\n"); 
*/
              #pragma omp for private(i, where_p,iPtr,index) schedule(static)
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
                    Paso_setError(VALUE_ERROR, "Paso_Solver_getAMG: main diagonal element missing.");
                } else {
                    iPtr+=(index_t)(where_p-index);
                    /* get inverse of A_FF block: */
                      S=out->inv_A_FF[i];
                      if (ABS(A_p->val[iPtr])>0.) {
                            if(ABS(S)>0.)
                               out->inv_A_FF[i]=1./S;
                            else
                            {
                               fprintf(stderr,"ROWSUM OF ROW %d->%d is 0\n",i,out->rows_in_F[i]);
                               S=0;
                               for (iPtr=A_p->pattern->ptr[out->rows_in_F[i]];iPtr<A_p->pattern->ptr[out->rows_in_F[i] + 1]; ++iPtr) {
                               if(A_p->val[iPtr]!=0) {
                                fprintf(stderr,"A_p[%d,%d]=%f ",out->rows_in_F[i],A_p->pattern->index[iPtr],A_p->val[iPtr]);
                                S+=A_p->val[iPtr];
                               }
                              }
                              fprintf(stderr,"\n SUMMMMMMM %f\n",S);
                            }
                      } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getAMG: Break-down in AMG decomposition: non-regular main diagonal block.");
                      }
                } 
              }
           } /* end parallel region */

           if( Paso_noError()) {
              /* if there are no nodes in the coarse level there is no more work to do */
              out->n_C=n-out->n_F;
              if (level<2) {
              /*if (out->n_C>10) {*/
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
                      out->A_CF=Paso_SparseMatrix_getSubmatrix(A_p,out->n_C,out->n_F,out->rows_in_C,out->mask_F);
                      if (Paso_noError()) {
                         /* get A_FC block: */
                         out->A_FC=Paso_SparseMatrix_getSubmatrix(A_p,out->n_F,out->n_C,out->rows_in_F,out->mask_C);
                         /* get A_CC block: */
                         if (Paso_noError()) {
                            schur=Paso_SparseMatrix_getSubmatrix(A_p,out->n_C,out->n_C,out->rows_in_C,out->mask_C);
                            
                            /*find the pattern of the schur complement with fill in*/
                            schur_withFillIn=Paso_SparseMatrix_alloc(A_p->type,Paso_Pattern_binop(PATTERN_FORMAT_DEFAULT, schur->pattern, Paso_Pattern_multiply(PATTERN_FORMAT_DEFAULT,out->A_CF->pattern,out->A_FC->pattern)),1,1);

                            /* copy values over*/ 
                            #pragma omp for private(i,iPtr,iPtr_s,j,j0) schedule(static)
                            for (i = 0; i < schur_withFillIn->numRows; ++i) {
                              for (iPtr=schur_withFillIn->pattern->ptr[i];iPtr<schur_withFillIn->pattern->ptr[i + 1]; ++iPtr) {
                                j=schur_withFillIn->pattern->index[iPtr];
                                schur_withFillIn->val[iPtr]=0.;
                                for (iPtr_s=schur->pattern->ptr[i];iPtr_s<schur->pattern->ptr[i + 1]; ++iPtr_s){
                                    j0=schur->pattern->index[iPtr_s];
                                    if (j==j0) {
                                      schur_withFillIn->val[iPtr]=schur->val[iPtr_s];
                                      break;
                                    }
                                }
                              }
                            }
                           time0=Paso_timer()-time0;
                           
                            if (Paso_noError()) {
                                time1=Paso_timer();
                                /* update A_CC block to get Schur complement and then apply AMG to it */
                                Paso_Solver_updateIncompleteSchurComplement(schur_withFillIn,out->A_CF,out->inv_A_FF,out->A_FF_pivot,out->A_FC);
                                time1=Paso_timer()-time1;
                                out->AMG_of_Schur=Paso_Solver_getAMG(schur_withFillIn,verbose,level+1);
                                
                                /*Paso_SparseMatrix_free(schur);*/
                                /* Paso_SparseMatrix_free(schur_withFillIn);*/
                            }
                            /* allocate work arrays for AMG application */
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
         printf("AMG: %d unknowns eliminated. %d left.\n",out->n_F,n-out->n_F);
         if (level<2) {
            printf("timing: AMG: MIS/reordering/elemination : %e/%e/%e\n",time2,time0,time1);
         } else {
            printf("timing: AMG: MIS: %e\n",time2); 
         }
     }
     return out;
  } else  {
     Paso_Solver_AMG_free(out);
     return NULL;
  }
}

/**************************************************************/

/* apply AMG precondition b-> x                               

     in fact it solves 

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FF ]  [ x_F ]  = [b_F]
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ]  [ x_C ]  = [b_C]

 in the form 

   b->[b_F,b_C] 
   x_F=invA_FF*b_F
   b_C=b_C-A_CF*x_F
   x_C=AMG(b_C)
   b_F=b_F-A_FC*x_C
   x_F=invA_FF*b_F
   x<-[x_F,x_C]

 should be called within a parallel region                                              
 barrier synconization should be performed to make sure that the input vector available 

*/

void Paso_Solver_solveAMG(Paso_Solver_AMG * amg, double * x, double * b) {
     dim_t i,oldsweeps;
     double *r=MEMALLOC(amg->n,double);
     /*Paso_Solver_GS* GS=NULL;*/
     double *bold=MEMALLOC(amg->n,double);
     double *bnew=MEMALLOC(amg->n,double);
     double *x0=MEMALLOC(amg->n,double);

     if (amg->level==2) {
     /*if (amg->n_C<=10) {*/
        Paso_UMFPACK1(amg->A,x,b,0);
     } else {
     
        /* presmoothing */
         Paso_Solver_solveGS(amg->GS,x,b);
          oldsweeps=amg->GS->sweeps;
         if (amg->GS->sweeps>1) {
           
           #pragma omp parallel for private(i) schedule(static)
           for (i=0;i<amg->GS->n;++i) bold[i]=b[i];
           
           while(amg->GS->sweeps>1) {
               #pragma omp parallel for private(i) schedule(static)
               for (i=0;i<amg->GS->n;++i) bnew[i]=bold[i]+b[i];
               Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1), amg->A, x, DBLE(1), bnew);
               Paso_Solver_solveGS(amg->GS,x,bnew);
               #pragma omp parallel for private(i) schedule(static)
               for (i=0;i<amg->GS->n;++i) bold[i]=bnew[i];
               amg->GS->sweeps=amg->GS->sweeps-1;
           }
           }
           amg->GS->sweeps=oldsweeps;
           
        /* end of presmoothing */
        
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<amg->n;++i) r[i]=b[i];
         
         /*r=b-Ax*/ 
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
       
        /* b->[b_F,b_C]     */
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amg->n_F;++i) amg->b_F[i]=r[amg->rows_in_F[i]];
        
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amg->n_C;++i) amg->b_C[i]=r[amg->rows_in_C[i]];

        /* x_F=invA_FF*b_F  */
        Paso_Solver_applyBlockDiagonalMatrix(1,amg->n_F,amg->inv_A_FF,amg->A_FF_pivot,amg->x_F,amg->b_F);
        
        /* b_C=b_C-A_CF*x_F */
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A_CF,amg->x_F,1.,amg->b_C);
        
        /* x_C=AMG(b_C)     */
        Paso_Solver_solveAMG(amg->AMG_of_Schur,amg->x_C,amg->b_C);
        
        /* b_F=b_F-A_FC*x_C */
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A_FC,amg->x_C,1.,amg->b_F);
        /* x_F=invA_FF*b_F  */
        Paso_Solver_applyBlockDiagonalMatrix(1,amg->n_F,amg->inv_A_FF,amg->A_FF_pivot,amg->x_F,amg->b_F);
        /* x<-[x_F,x_C]     */
        
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amg->n;++i) {
            if (amg->mask_C[i]>-1) {
                 x[i]+=amg->x_C[amg->mask_C[i]];
            } else {
                 x[i]+=amg->x_F[amg->mask_F[i]];
            }
        }
        

     /*postsmoothing*/
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) r[i]=b[i];
     
     /*for (i=0;i<amg->n;++i) fprintf(stderr,"%lf %lf %d %d %d BBB\n",r[i],x[i],amg->n,amg->A->numCols, amg->level);*/
     
     /*
     for (i = 0; i < amg->A->numRows; ++i) {
          for (iPtr=amg->A->pattern->ptr[i];iPtr<amg->A->pattern->ptr[i + 1]; ++iPtr) {
              j=amg->A->pattern->index[iPtr];
              fprintf(stderr,"A[%d,%d]=%lf ",i,j,amg->A->val[iPtr]);
          } 
          fprintf(stderr,"\n");
     }*/
     
     /*r=b-Ax*/ 
     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
     Paso_Solver_solveGS(amg->GS,x0,r);
     
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) {
      /*fprintf(stderr,"%lf %lf %lf\n",x[i],x0[i],x[i]+x0[i]);*/
      x[i]+=x0[i];
     }
     
     /*end of postsmoothing*/
     
     /*Paso_Solver_GS_free(amg->GS);*/
     }
     MEMFREE(r);
     MEMFREE(bold);
     MEMFREE(bnew);
     MEMFREE(x0);
    return;
}

/*
 * $Log$
 *
 */
