
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: AMG preconditioner                                  */

/**************************************************************/

/* Author: artak@uq.edu.au                                */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "Options.h"
#include "PasoUtil.h"
#include "UMFPACK.h"
#include "MKL.h"
#include "SystemMatrix.h"
#include "Pattern_coupling.h"

/**************************************************************/

/* free all memory used by AMG                                */

void Paso_Solver_AMG_System_free(Paso_Solver_AMG_System * in) {
     dim_t i;
     if (in!=NULL) {
        for (i=0;i<in->block_size;++i) {
          Paso_Solver_AMG_free(in->amgblock[i]);
          Paso_SparseMatrix_free(in->block[i]);
        }
        MEMFREE(in);
     }
}


/* free all memory used by AMG                                */

void Paso_Solver_AMG_free(Paso_Solver_AMG * in) {
     if (in!=NULL) {
        Paso_Solver_Jacobi_free(in->GS);
        Paso_SparseMatrix_free(in->A_FC);
        Paso_SparseMatrix_free(in->A_FF);
        Paso_SparseMatrix_free(in->W_FC);
        Paso_SparseMatrix_free(in->A_CF);
        Paso_SparseMatrix_free(in->P);
        Paso_SparseMatrix_free(in->R);
        
        Paso_SparseMatrix_free(in->A);
        if(in->coarsest_level==TRUE) {
        #ifdef MKL
          Paso_MKL_free1(in->AOffset1);
          Paso_SparseMatrix_free(in->AOffset1);
        #else
          #ifdef UMFPACK
          Paso_UMFPACK1_free((Paso_UMFPACK_Handler*)(in->solver));
          #endif
        #endif
        }
        MEMFREE(in->rows_in_F);
        MEMFREE(in->rows_in_C);
        MEMFREE(in->mask_F);
        MEMFREE(in->mask_C);
        MEMFREE(in->x_F);
        MEMFREE(in->b_F);
        MEMFREE(in->x_C);
        MEMFREE(in->b_C);
        in->solver=NULL;
        Paso_Solver_AMG_free(in->AMG_of_Coarse);
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
Paso_Solver_AMG* Paso_Solver_getAMG(Paso_SparseMatrix *A_p,dim_t level,Paso_Options* options) {
  Paso_Solver_AMG* out=NULL;
  /*
   Paso_Pattern* temp1=NULL;
  Paso_Pattern* temp2=NULL;
  */
  bool_t verbose=options->verbose;
  dim_t n=A_p->numRows;
  dim_t n_block=A_p->row_block_size;
  index_t* mis_marker=NULL;  
  index_t* counter=NULL;
  /*index_t iPtr,*index, *where_p;*/
  dim_t i;
  Paso_SparseMatrix * A_c=NULL;
  
  /*
  double *temp,*temp_1;
  double S;
  index_t iptr;
  */
  
  /*char filename[8];*/
  /*sprintf(filename,"AMGLevel%d",level);
  
  Paso_SparseMatrix_saveMM(A_p,filename);
  */
  
  /*Make sure we have block sizes 1*/
  if (A_p->col_block_size>1) {
     Paso_setError(TYPE_ERROR,"Paso_Solver_getAMG: AMG requires column block size 1.");
     return NULL;
  }
  if (A_p->row_block_size>1) {
     Paso_setError(TYPE_ERROR,"Paso_Solver_getAMG: AMG requires row block size 1.");
     return NULL;
  }
  out=MEMALLOC(1,Paso_Solver_AMG);
  /* identify independend set of rows/columns */
  mis_marker=TMPMEMALLOC(n,index_t);
  counter=TMPMEMALLOC(n,index_t);
  if ( !( Paso_checkPtr(mis_marker) || Paso_checkPtr(counter) || Paso_checkPtr(out)) ) {
     out->AMG_of_Coarse=NULL;
     out->A_FF=NULL;
     out->A_FC=NULL;
     out->A_CF=NULL;
     out->W_FC=NULL;
     out->P=NULL;
     out->R=NULL;
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
     out->GS=NULL;
     out->solver=NULL;
     /*out->GS=Paso_Solver_getGS(A_p,verbose);*/
     out->level=level;
     out->n=n;
     out->n_F=n+1;
     out->n_block=n_block;
     
     if (level==0 || n<=options->min_coarse_matrix_size) {
         out->coarsest_level=TRUE;
         #ifdef MKL
                  out->AOffset1=Paso_SparseMatrix_alloc(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, out->A->pattern,1,1, FALSE);
                  #pragma omp parallel for private(i) schedule(static)
                  for (i=0;i<out->A->len;++i) {
                       out->AOffset1->val[i]=out->A->val[i];
                  }
         #else
            #ifdef UMFPACK 
            #else
                out->GS=Paso_Solver_getJacobi(A_p);
            #endif
         #endif
     } else {
         out->coarsest_level=FALSE;
         out->GS=Paso_Solver_getJacobi(A_p);
 
         /* identify independend set of rows/columns */
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n;++i) mis_marker[i]=-1;

         if (options->coarsening_method == PASO_YAIR_SHAPIRA_COARSENING) {
              Paso_Pattern_YS(A_p,mis_marker,options->coarsening_threshold);
         }
         else if (options->coarsening_method == PASO_RUGE_STUEBEN_COARSENING) {
              Paso_Pattern_RS(A_p,mis_marker,options->coarsening_threshold);
         }
         else if (options->coarsening_method == PASO_AGGREGATION_COARSENING) {
             Paso_Pattern_Aggregiation(A_p,mis_marker,options->coarsening_threshold);
        }
        else {
           /*Default coarseneing*/
            Paso_Pattern_RS_MI(A_p,mis_marker,options->coarsening_threshold);
            /*Paso_Pattern_YS(A_p,mis_marker,options->coarsening_threshold);*/
            /*Paso_Pattern_RS(A_p,mis_marker,options->coarsening_threshold);*/
            /*Paso_Pattern_Aggregiation(A_p,mis_marker,options->coarsening_threshold);*/
            
        }

        #pragma omp parallel for private(i) schedule(static)
        for (i = 0; i < n; ++i) counter[i]=mis_marker[i];
        
        out->n_F=Paso_Util_cumsum(n,counter);
        
        if (out->n_F==0) {
           out->coarsest_level=TRUE;
           level=0;
           if (verbose) { 
               /*printf("AMG coarsening eliminates all unknowns, switching to Jacobi preconditioner.\n");*/
               printf("AMG coarsening does not eliminate any of the unknowns, switching to Jacobi preconditioner.\n");
           }
        }
        else if (out->n_F==n) {
          out->coarsest_level=TRUE;
           level=0;
           if (verbose) { 
               /*printf("AMG coarsening eliminates all unknowns, switching to Jacobi preconditioner.\n");*/
               printf("AMG coarsening eliminates all of the unknowns, switching to Jacobi preconditioner.\n");
          
            }
        } else {
     
              if (Paso_noError()) {
                 
                 /*#pragma omp parallel for private(i) schedule(static)
                 for (i = 0; i < n; ++i) counter[i]=mis_marker[i];
                 out->n_F=Paso_Util_cumsum(n,counter);
                 */
                 
                 out->mask_F=MEMALLOC(n,index_t);
                 out->rows_in_F=MEMALLOC(out->n_F,index_t);
                 if (! (Paso_checkPtr(out->mask_F) || Paso_checkPtr(out->rows_in_F) ) ) {
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
                    
                 }
              }
      
              /*check whether coarsening process actually makes sense to continue.
              if coarse matrix at least smaller by 30% then continue, otherwise we stop.*/
              if ((out->n_F*100/n)<30) {
                    level=1;
                }
             
              if ( Paso_noError()) {
                    /* if there are no nodes in the coarse level there is no more work to do */
                    out->n_C=n-out->n_F;
                   
                   /*if (out->n_F>500) */
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
                    }
              } 
              if ( Paso_noError()) {
                      /* get A_FF block: */
                      /*
                       out->A_FF=Paso_SparseMatrix_getSubmatrix(A_p,out->n_F,out->n_F,out->rows_in_F,out->mask_F);
                      out->A_CF=Paso_SparseMatrix_getSubmatrix(A_p,out->n_C,out->n_F,out->rows_in_C,out->mask_F);
                      out->A_FC=Paso_SparseMatrix_getSubmatrix(A_p,out->n_F,out->n_C,out->rows_in_F,out->mask_C);
                      */
                      
                      /*Compute W_FC*/
                      /*initialy W_FC=A_FC*/
                      out->W_FC=Paso_SparseMatrix_getSubmatrix(A_p,out->n_F,out->n_C,out->rows_in_F,out->mask_C);
                      
                      /*sprintf(filename,"W_FCbefore_%d",level);
                      Paso_SparseMatrix_saveMM(out->W_FC,filename);
                      */
                      Paso_SparseMatrix_updateWeights(A_p,out->W_FC,mis_marker);
                      
                      /*
                       printf("GOT W_FC, but n=%d,n_F=%d,n_C=%d\n",out->n,out->n_F,out->n_C);
                      sprintf(filename,"W_FCafter_%d",level);
                      Paso_SparseMatrix_saveMM(out->W_FC,filename);
                      */
                      /* get Prolongation and Restriction */
                      out->P=Paso_SparseMatrix_getProlongation(out->W_FC,mis_marker);
                      
                      /*
                       printf("GOT Prolongation P->nxc %dx%d\n",out->P->numRows,out->P->numCols);
                      sprintf(filename,"P_%d",level);
                      Paso_SparseMatrix_saveMM(out->P,filename);
                      */
                      
                      out->R=Paso_SparseMatrix_getRestriction(out->P);
                      /*
                      printf("GOT Restriction->cxn %dx%d\n",out->R->numRows,out->R->numCols);
                      sprintf(filename,"R_%d",level);
                      Paso_SparseMatrix_saveMM(out->R,filename);
                      */
                     
              }
              if ( Paso_noError()) {
                    A_c=Paso_Solver_getCoarseMatrix(A_p,out->R,out->P);
                    /*Paso_Solver_getCoarseMatrix(A_c, A_p,out->R,out->P);*/
                    /*
                     sprintf(filename,"A_C_%d",level);
                    Paso_SparseMatrix_saveMM(A_c,filename);
                    */ 
                    out->AMG_of_Coarse=Paso_Solver_getAMG(A_c,level-1,options);
              }

              /* allocate work arrays for AMG application */
              if (Paso_noError()) {
                         out->x_F=MEMALLOC(n_block*out->n_F,double);
                         out->b_F=MEMALLOC(n_block*out->n_F,double);
                         out->x_C=MEMALLOC(n_block*out->n_C,double);
                         out->b_C=MEMALLOC(n_block*out->n_C,double);
      
                         if (! (Paso_checkPtr(out->x_F) || Paso_checkPtr(out->b_F) || Paso_checkPtr(out->x_C) || Paso_checkPtr(out->b_C) ) ) {
                             #pragma omp parallel for private(i) schedule(static)
                             for (i = 0; i < out->n_F; ++i) {
                                         out->x_F[i]=0.;
                                         out->b_F[i]=0.;
                              }
                              #pragma omp parallel for private(i) schedule(static)
                              for (i = 0; i < out->n_C; ++i) {
                                     out->x_C[i]=0.;
                                     out->b_C[i]=0.;
                              }
                         }
              }
            Paso_SparseMatrix_free(A_c);
         }
     }
  }
  TMPMEMFREE(mis_marker);
  TMPMEMFREE(counter);

  if (Paso_noError()) {
      if (verbose && level>0 && !out->coarsest_level) {
         printf("AMG: level: %d: %d unknowns eliminated. %d left.\n",level, out->n_F,out->n_C);
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

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FC ]  [ x_F ]  = [b_F]
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
     double time0=0;
     double *r=NULL, *x0=NULL;
     bool_t verbose=0;
     #ifdef UMFPACK 
          Paso_UMFPACK_Handler * ptr=NULL;
     #endif
     
     r=MEMALLOC(amg->n,double);
     x0=MEMALLOC(amg->n,double);
     
     if (amg->coarsest_level) {
      
      time0=Paso_timer();
      /*If all unknown are eliminated then Jacobi is the best preconditioner*/
      if (amg->n_F==0 || amg->n_F==amg->n) {
        Paso_Solver_solveJacobi(amg->GS,x,b);
      }
       else {
       #ifdef MKL
          Paso_MKL1(amg->AOffset1,x,b,verbose);
       #else
          #ifdef UMFPACK
             ptr=(Paso_UMFPACK_Handler *)(amg->solver);
             Paso_UMFPACK1(&ptr,amg->A,x,b,verbose);
             amg->solver=(void*) ptr;
          #else      
             Paso_Solver_solveJacobi(amg->GS,x,b);
         #endif
       #endif
       }
       time0=Paso_timer()-time0;
       if (verbose) fprintf(stderr,"timing: DIRECT SOLVER: %e\n",time0);
       
     } else {
        /* presmoothing */
         time0=Paso_timer();
         Paso_Solver_solveJacobi(amg->GS,x,b);
         time0=Paso_timer()-time0;
         if (verbose) fprintf(stderr,"timing: Presmooting: %e\n",time0);
        /* end of presmoothing */
        
        
         time0=Paso_timer();
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<amg->n;++i) r[i]=b[i];
         
         /*r=b-Ax*/ 
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
         
         
        /* b_c <- R*r  */
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(1.,amg->R,r,0.,amg->b_C);
        
        time0=Paso_timer()-time0;
        if (verbose) fprintf(stderr,"timing: Before next level: %e\n",time0);
        
        /* x_C=AMG(b_C)     */
        Paso_Solver_solveAMG(amg->AMG_of_Coarse,amg->x_C,amg->b_C);
        
        time0=Paso_timer();
        
        /* x_0 <- P*x_c */
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(1.,amg->P,amg->x_C,0.,x0);
        
        /* x=x+x0 */
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amg->n;++i) x[i]+=x0[i];
        
     /*postsmoothing*/
     time0=Paso_timer();
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i) r[i]=b[i];
     
     /*r=b-Ax */
     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
     Paso_Solver_solveJacobi(amg->GS,x0,r);
     
     
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<amg->n;++i)  {
      x[i]+=x0[i];
      /*printf("x[%d]=%e \n",i,x[i]);*/
     }
     
     time0=Paso_timer()-time0;
     if (verbose) fprintf(stderr,"timing: Postsmoothing: %e\n",time0);

     /*end of postsmoothing*/
     
     }
     MEMFREE(r);
     MEMFREE(x0);

     return;
}
