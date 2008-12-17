
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
        Paso_Solver_Jacobi_free(in->GS);
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
  dim_t i,k,j;
  Paso_SparseMatrix * schur=NULL;
  Paso_SparseMatrix * schur_withFillIn=NULL;
  double time0=0,time1=0,time2=0,S=0;
  /*Paso_Pattern* test;*/
  
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
  out->GS=NULL;
  out->A=Paso_SparseMatrix_getReference(A_p);
  /*out->GS=Paso_Solver_getGS(A_p,verbose);*/
  out->GS=Paso_Solver_getJacobi(A_p);
  /*out->GS->sweeps=2;*/
  out->level=level;
   
  if ( !(Paso_checkPtr(mis_marker) || Paso_checkPtr(out) || Paso_checkPtr(counter) ) ) {
     /* identify independend set of rows/columns */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;++i) mis_marker[i]=-1;
     /*Paso_Pattern_RS(A_p,mis_marker,0.25);*/
     /*Paso_Pattern_Aggregiation(A_p,mis_marker,0.5);*/
     Paso_Pattern_coup(A_p,mis_marker,0.05);
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
              /* creates an index for F from mask */
              #pragma omp parallel for private(i) schedule(static)
              for (i = 0; i < out->n_F; ++i) out->rows_in_F[i]=-1;
              #pragma omp parallel for private(i) schedule(static)
              for (i = 0; i < n; ++i) {
                 if  (mis_marker[i]) {
                        out->rows_in_F[counter[i]]=i;
                        out->mask_F[i]=counter[i];
                 } else {
                        out->mask_F[i]=-1;
                 }
              }
              
              /* Compute row-sum for getting rs(A_FF)^-1*/
              #pragma omp parallel for private(i,iPtr,j,S) schedule(static)
              for (i = 0; i < out->n_F; ++i) {
                S=0;
                for (iPtr=A_p->pattern->ptr[out->rows_in_F[i]];iPtr<A_p->pattern->ptr[out->rows_in_F[i] + 1]; ++iPtr) {
                 j=A_p->pattern->index[iPtr];
                 if (mis_marker[j])
                     S+=A_p->val[iPtr];
                }
                   out->inv_A_FF[i]=1./S;
              }

           if( Paso_noError()) {
              /* if there are no nodes in the coarse level there is no more work to do */
              out->n_C=n-out->n_F;
              if (level<5) {
               /*if (out->n_F>500) {*/
                   out->rows_in_C=MEMALLOC(out->n_C,index_t);
                   out->mask_C=MEMALLOC(n,index_t);
                   if (! (Paso_checkPtr(out->mask_C) || Paso_checkPtr(out->rows_in_C) ) ) {
                       /* creates an index for C from mask */
                       #pragma omp parallel for private(i) schedule(static)
                       for (i = 0; i < n; ++i) counter[i]=! mis_marker[i];
                       Paso_Util_cumsum(n,counter);
                          #pragma omp parallel for private(i) schedule(static)
                          for (i = 0; i < out->n_C; ++i) out->rows_in_C[i]=-1;
                          #pragma omp parallel for private(i) schedule(static)
                          for (i = 0; i < n; ++i) {
                             if  (! mis_marker[i]) {
                                out->rows_in_C[counter[i]]=i;
                                out->mask_C[i]=counter[i];
                             } else {
                                out->mask_C[i]=-1;
                             }
                          }
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
                            
                            /*fprintf(stderr,"Sparsity of Schure: %dx%d LEN %d Percentage %f\n",schur_withFillIn->pattern->numOutput,schur_withFillIn->pattern->numInput,schur_withFillIn->len,schur_withFillIn->len/(1.*schur_withFillIn->pattern->numOutput*schur_withFillIn->pattern->numInput));*/
                            
                            /* copy values over*/ 
                            #pragma omp parallel for private(i,iPtr,j,iPtr_s,index,where_p) schedule(static)
                            for (i = 0; i < schur_withFillIn->numRows; ++i) {
                              for (iPtr=schur_withFillIn->pattern->ptr[i];iPtr<schur_withFillIn->pattern->ptr[i + 1]; ++iPtr) {
                                j=schur_withFillIn->pattern->index[iPtr];
                                iPtr_s=schur->pattern->ptr[i];
                                schur_withFillIn->val[iPtr]=0.;
                                index=&(schur->pattern->index[iPtr_s]);
                                where_p=(index_t*)bsearch(&j,
                                        index,
                                        schur->pattern->ptr[i + 1]-schur->pattern->ptr[i],
                                        sizeof(index_t),
                                        Paso_comparIndex);
                                if (where_p!=NULL) {
                                    schur_withFillIn->val[iPtr]=schur->val[iPtr_s+(index_t)(where_p-index)];
                                }
                              }
                            }
                                
                            if (Paso_noError()) {
                                Paso_Solver_updateIncompleteSchurComplement(schur_withFillIn,out->A_CF,out->inv_A_FF,out->A_FF_pivot,out->A_FC);
                                out->AMG_of_Schur=Paso_Solver_getAMG(schur_withFillIn,verbose,level+1);
                                Paso_SparseMatrix_free(schur);
                            }
                            /* allocate work arrays for AMG application */
                            if (Paso_noError()) {
                              out->x_F=MEMALLOC(n_block*out->n_F,double);
                              out->b_F=MEMALLOC(n_block*out->n_F,double);
                              out->x_C=MEMALLOC(n_block*out->n_C,double);
                              out->b_C=MEMALLOC(n_block*out->n_C,double);

                              if (! (Paso_checkPtr(out->x_F) || Paso_checkPtr(out->b_F) || Paso_checkPtr(out->x_C) || Paso_checkPtr(out->b_C) ) ) {
                                    #pragma omp parallel for private(i,k) schedule(static)
                                    for (i = 0; i < out->n_F; ++i) {
                                          for (k=0; k<n_block;++k) {
                                             out->x_F[i*n_block+k]=0.;
                                             out->b_F[i*n_block+k]=0.;
                                          }
                                    }
                                    #pragma omp parallel for private(i,k) schedule(static)
                                    for (i = 0; i < out->n_C; ++i) {
                                        for (k=0; k<n_block;++k) {
                                          out->x_C[i*n_block+k]=0.;
                                          out->b_C[i*n_block+k]=0.;
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
     }
  }
  TMPMEMFREE(mis_marker);
  TMPMEMFREE(counter);
  if (Paso_noError()) {
      if (verbose) {
         printf("AMG: %d unknowns eliminated. %d left.\n",out->n_F,n-out->n_F);
         if (level<5) {
         /*if (out->n_F<500) {*/
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
     dim_t i;
     double *r=MEMALLOC(amg->n,double);
     /*Paso_Solver_GS* GS=NULL;*/
     double *x0=MEMALLOC(amg->n,double);
     double time0=0;
     
     if (amg->level==5) {
     /*if (amg->n_F<=500) {*/
      time0=Paso_timer();
        
        Paso_Solver_solveJacobi(amg->GS,x,b);
        
        /* #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amg->n;++i) r[i]=b[i];
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
        Paso_Solver_solveGS(amg->GS,x0,r);
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amg->n;++i) {
         x[i]+=x0[i];
        }
        */
        /*Paso_UMFPACK1(amg->A,x,b,0);*/

       time0=Paso_timer()-time0;
       /*fprintf(stderr,"timing: DIRECT SOLVER: %e/\n",time0);*/
     } else {
     
        /* presmoothing */
         Paso_Solver_solveJacobi(amg->GS,x,b);
         /*
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<amg->n;++i) r[i]=b[i];

          Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
          Paso_Solver_solveGS(amg->GS,x0,r);
          
          #pragma omp parallel for private(i) schedule(static)
          for (i=0;i<amg->n;++i) {
           x[i]+=x0[i];
          }
          */
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
                 x[i]=amg->x_C[amg->mask_C[i]];
            } else {
                 x[i]=amg->x_F[amg->mask_F[i]];
            }
        }
        
     /*postsmoothing*/

     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) r[i]=b[i];
     
     /*r=b-Ax*/ 
     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
     Paso_Solver_solveJacobi(amg->GS,x0,r);
     
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) x[i]+=x0[i];
     
     
     /*#pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) r[i]=b[i];

     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
     Paso_Solver_solveGS(amg->GS,x0,r);
     
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) {
      x[i]+=x0[i];
     }
     */
     /*end of postsmoothing*/
     
     }
     MEMFREE(r);
     MEMFREE(x0);
    return;
}

/*
 * $Log$
 *
 */
