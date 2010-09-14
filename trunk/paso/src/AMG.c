
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
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
#include "Preconditioner.h"
#include "Options.h"
#include "PasoUtil.h"
#include "UMFPACK.h"
#include "MKL.h"
#include "SystemMatrix.h"
#include "Pattern_coupling.h"

/**************************************************************/

/* free all memory used by AMG                                */

void Paso_Solver_AMG_free(Paso_Solver_AMG * in) {
     if (in!=NULL) {
	Paso_Preconditioner_LocalSmoother_free(in->Smoother);

	
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
          Paso_SparseMatrix_free(in->AUnrolled);
        #else
          #ifdef UMFPACK
          Paso_UMFPACK1_free((Paso_UMFPACK_Handler*)(in->solver));
          Paso_SparseMatrix_free(in->AUnrolled);
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
  bool_t timing=0;
  
  dim_t n=A_p->numRows;
  dim_t n_block=A_p->row_block_size;
  
  
  index_t* mis_marker=NULL;  
  index_t* counter=NULL;
  /*index_t iPtr,*index, *where_p;*/
  dim_t i,j;
  Paso_SparseMatrix * A_c=NULL;
  double time0=0;
  Paso_SparseMatrix * Atemp=NULL;
  double sparsity=0;
  
  /*
  double *temp,*temp_1;
  double S;
  index_t iptr;
  */
  
  /*char filename[8];*/
  
  /*  
  sprintf(filename,"AMGLevel%d",level);
  
  Paso_SparseMatrix_saveMM(A_p,filename);
  */
  
  /*Make sure we have block sizes 1*/
  /*if (A_p->col_block_size>1) {
     Paso_setError(TYPE_ERROR,"Paso_Solver_getAMG: AMG requires column block size 1.");
     return NULL;
  }
  if (A_p->row_block_size>1) {
     Paso_setError(TYPE_ERROR,"Paso_Solver_getAMG: AMG requires row block size 1.");
     return NULL;
  }*/
  out=MEMALLOC(1,Paso_Solver_AMG);

  
  /* identify independend set of rows/columns */
  mis_marker=TMPMEMALLOC(n,index_t);
  counter=TMPMEMALLOC(n,index_t);
  if ( !( Paso_checkPtr(mis_marker) || Paso_checkPtr(counter) || Paso_checkPtr(out)) ) {
     
     
     out->post_sweeps=options->post_sweeps;
     out->pre_sweeps=options->pre_sweeps;
     
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
     out->A=Paso_SparseMatrix_getReference(A_p);
     out->AUnrolled=NULL;
     out->AOffset1=NULL;
     out->solver=NULL;
     out->Smoother=NULL;
     out->level=level;
     out->n=n;
     out->n_F=n+1;
     out->n_block=n_block;

     
     sparsity=(A_p->len*1.)/(1.*A_p->numRows*A_p->numCols);
     
     if (verbose) fprintf(stdout,"Stats: Sparsity of the Coarse Matrix with %d non-zeros (%d,%d) in level %d is %.6f\n",A_p->len,A_p->numRows,A_p->numCols,level,sparsity);
     
    
     if(sparsity>0.05) {
      level=0;
     }
     
         
     if (level==0 || n<=options->min_coarse_matrix_size) {
         out->coarsest_level=TRUE;
         #ifdef MKL
                  out->AUnrolled=Paso_SparseMatrix_unroll(A_p);
                  out->AOffset1=Paso_SparseMatrix_alloc(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, out->AUnrolled->pattern,1,1, FALSE);
                  #pragma omp parallel for private(i) schedule(static)
                  for (i=0;i<out->A->len;++i) {
                       out->AOffset1->val[i]=out->AUnrolled->val[i];
                  }
         #else
            #ifdef UMFPACK
                out->AUnrolled=Paso_SparseMatrix_unroll(A_p);
                /*Paso_SparseMatrix_saveMM(out->AUnrolled,"AUnroll.mat");
                Paso_SparseMatrix_saveMM(A_p,"Aorg.mat");
                */
            #else
              out->Smoother=Paso_Preconditioner_LocalSmoother_alloc(A_p, (options->smoother == PASO_JACOBI), verbose);
            #endif
         #endif
         
     } else {
         out->coarsest_level=FALSE;
	 out->Smoother=Paso_Preconditioner_LocalSmoother_alloc(A_p, (options->smoother == PASO_JACOBI), verbose);
 
         /* identify independend set of rows/columns */
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<n;++i) mis_marker[i]=-1;

          /*mesuring coarsening time */
          time0=Paso_timer();
          
         if (options->coarsening_method == PASO_YAIR_SHAPIRA_COARSENING) {
              Paso_Pattern_YS(A_p,mis_marker,options->coarsening_threshold);
         }
         else if (options->coarsening_method == PASO_RUGE_STUEBEN_COARSENING) {
              Paso_Pattern_RS(A_p,mis_marker,options->coarsening_threshold);
         }
         else if (options->coarsening_method == PASO_AGGREGATION_COARSENING) {
             Paso_Pattern_Aggregiation(A_p,mis_marker,options->coarsening_threshold);
        }
        else if (options->coarsening_method == PASO_STANDARD_COARSENING) {
             Paso_Pattern_Standard_Block(A_p,mis_marker,options->coarsening_threshold);
        }
        else {
           /*Default coarseneing*/
            Paso_Pattern_Standard_Block(A_p,mis_marker,options->coarsening_threshold);
            /*Paso_Pattern_Read("Standard.spl",n,mis_marker);*/
            /*Paso_Pattern_YS(A_p,mis_marker,options->coarsening_threshold);*/
            /*Paso_Pattern_greedy_Agg(A_p,mis_marker,options->coarsening_threshold);*/
            /*Paso_Pattern_greedy(A_p->pattern,mis_marker);*/
            /*Paso_Pattern_Aggregiation(A_p,mis_marker,options->coarsening_threshold);*/
            
        }
        
        if (timing) fprintf(stdout,"timing: Profilining for level %d:\n",level);
        
        time0=Paso_timer()-time0;
        if (timing) fprintf(stdout,"timing: Coarsening: %e\n",time0);

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
              
              /* if(level==1) {
                   printf("##TOTAL: %d, ELIMINATED: %d\n",n,out->n_F);
                   for (i = 0; i < n; ++i) {
                    printf("##%d %d\n",i,!mis_marker[i]);
                   }
                }
              */
              
              /*check whether coarsening process actually makes sense to continue.
              if coarse matrix at least smaller by 30% then continue, otherwise we stop.*/
              if ((out->n_F*100/n)<30) {
                    level=1;
                }
             
              if ( Paso_noError()) {
                    out->n_C=n-out->n_F;
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
                      /* for (i = 0; i < n; ++i) {
                                  printf("##mis_marker[%d]=%d\n",i,mis_marker[i]);
                       }
                      */                   
                      time0=Paso_timer();
                      Paso_SparseMatrix_updateWeights(A_p,out->W_FC,mis_marker);
                      time0=Paso_timer()-time0;
                      if (timing) fprintf(stdout,"timing: updateWeights: %e\n",time0);
        
                      /*
                      sprintf(filename,"W_FCafter_%d",level);
                      Paso_SparseMatrix_saveMM(out->W_FC,filename);
                      */
                      
                      /* get Prolongation and Restriction */
                      time0=Paso_timer();
                      out->P=Paso_SparseMatrix_getProlongation(out->W_FC,mis_marker);
                      time0=Paso_timer()-time0;
                      if (timing) fprintf(stdout,"timing: getProlongation: %e\n",time0);
                      /*out->P=Paso_SparseMatrix_loadMM_toCSR("P1.mtx");*/
                      
                      /*
                      sprintf(filename,"P_%d",level);
                      Paso_SparseMatrix_saveMM(out->P,filename);
                      */
                      
                      time0=Paso_timer();
                      out->R=Paso_SparseMatrix_getRestriction(out->P);
                      time0=Paso_timer()-time0;
                      if (timing) fprintf(stdout,"timing: getRestriction: %e\n",time0);
                      /*out->R=Paso_SparseMatrix_loadMM_toCSR("R1.mtx");*/
                      
                      /*
                      sprintf(filename,"R_%d",level);
                      Paso_SparseMatrix_saveMM(out->R,filename);
                      */
                     
              }
              if ( Paso_noError()) {
                    
                    time0=Paso_timer();
                    
                    Atemp=Paso_SparseMatrix_MatrixMatrix(A_p,out->P);
                    
                    A_c=Paso_SparseMatrix_MatrixMatrix(out->R,Atemp);
                    
                    /*A_c=Paso_SparseMatrix_loadMM_toCSR("A_C1.mtx");*/
                    
                    Paso_SparseMatrix_free(Atemp);
                    
                    /*A_c=Paso_Solver_getCoarseMatrix(A_p,out->R,out->P);*/
                    time0=Paso_timer()-time0;
                    if (timing) fprintf(stdout,"timing: getCoarseMatrix: %e\n",time0);
                    
                                        
                    /*Paso_Solver_getCoarseMatrix(A_c, A_p,out->R,out->P);*/
                    /*                    
                    sprintf(filename,"A_C_%d",level);
                    Paso_SparseMatrix_saveMM(A_c,filename);
                    */
                                         
                    out->AMG_of_Coarse=Paso_Solver_getAMG(A_c,level-1,options);
              }

              /* allocate work arrays for AMG application */
              if (Paso_noError()) {
                         /*
                          out->x_F=MEMALLOC(n_block*out->n_F,double);
                         out->b_F=MEMALLOC(n_block*out->n_F,double);
                         */
                         out->x_C=MEMALLOC(n_block*out->n_C,double);
                         out->b_C=MEMALLOC(n_block*out->n_C,double);
      
                         /*if (! (Paso_checkPtr(out->x_F) || Paso_checkPtr(out->b_F) || Paso_checkPtr(out->x_C) || Paso_checkPtr(out->b_C) ) ) {*/
                         if ( ! ( Paso_checkPtr(out->x_C) || Paso_checkPtr(out->b_C) ) ) {
                             
                             /*
                              #pragma omp parallel for private(i) schedule(static)
                             for (i = 0; i < out->n_F; ++i) {
                                         out->x_F[i]=0.;
                                         out->b_F[i]=0.;
                              }
                             */
                             
                              #pragma omp parallel for private(i,j) schedule(static)
                              for (i = 0; i < out->n_C; ++i) {
                                     for(j=0;j<n_block;++j) {
                                        out->x_C[i*n_block+j]=0.;
                                        out->b_C[i*n_block+j]=0.;
                                     }
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
     dim_t i,j;
     double time0=0;
     double *r=NULL, *x0=NULL;
     bool_t timing=0;
     
     dim_t post_sweeps=amg->post_sweeps;
     dim_t pre_sweeps=amg->pre_sweeps;
     
     #ifdef UMFPACK 
          Paso_UMFPACK_Handler * ptr=NULL;
     #endif
          
     r=MEMALLOC(amg->n*amg->n_block,double);
     x0=MEMALLOC(amg->n*amg->n_block,double);
     
     if (amg->coarsest_level) {
      
      time0=Paso_timer();
      /*If all unknown are eliminated then Jacobi is the best preconditioner*/

      if (amg->n_F==0 || amg->n_F==amg->n) {
	   Paso_Preconditioner_LocalSmoother_solve(amg->A, amg->Smoother,x,b,1, FALSE);
      }
       else {
       #ifdef MKL
          Paso_MKL1(amg->AOffset1,x,b,timing);
       #else
          #ifdef UMFPACK
             ptr=(Paso_UMFPACK_Handler *)(amg->solver);
             Paso_UMFPACK1(&ptr,amg->AUnrolled,x,b,timing);
             amg->solver=(void*) ptr;
          #else      
          Paso_Preconditioner_LocalSmoother_solve(amg->A, amg->Smoother,x,b,1, FALSE);
         #endif
       #endif
       }
      
       time0=Paso_timer()-time0;
       if (timing) fprintf(stdout,"timing: DIRECT SOLVER: %e\n",time0);
       
     } else {
        /* presmoothing */
         time0=Paso_timer();
	 Paso_Preconditioner_LocalSmoother_solve(amg->A, amg->Smoother, x, b, pre_sweeps, FALSE); 
         time0=Paso_timer()-time0;
         if (timing) fprintf(stdout,"timing: Presmooting: %e\n",time0);
	 
         /* end of presmoothing */
        
         time0=Paso_timer();
          #pragma omp parallel for private(i,j) schedule(static)
           for (i=0;i<amg->n;++i) {
            for (j=0;j<amg->n_block;++j) {
              r[i*amg->n_block+j]=b[i*amg->n_block+j];  
            }
           }
         
         /*r=b-Ax*/ 
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amg->A,x,1.,r);
         
        /* b_c = R*r  */
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(1.,amg->R,r,0.,amg->b_C);
        
        time0=Paso_timer()-time0;
        if (timing) fprintf(stdout,"timing: Before next level: %e\n",time0);
        
        /* x_C=AMG(b_C)     */
        Paso_Solver_solveAMG(amg->AMG_of_Coarse,amg->x_C,amg->b_C);
        
        time0=Paso_timer();
        
        /* x_0 = P*x_c */
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(1.,amg->P,amg->x_C,0.,x0);
        
        /* x=x+x0 */
         #pragma omp parallel for private(i,j) schedule(static)
           for (i=0;i<amg->n;++i) {
            for (j=0;j<amg->n_block;++j) {
              x[i*amg->n_block+j]+=x0[i*amg->n_block+j];  
            }
           }
        
      /*postsmoothing*/
      
      /*solve Ax=b with initial guess x */
      time0=Paso_timer();
      Paso_Preconditioner_LocalSmoother_solve(amg->A, amg->Smoother, x, b, post_sweeps, TRUE); 
      time0=Paso_timer()-time0;
      if (timing) fprintf(stdout,"timing: Postsmoothing: %e\n",time0);
 
      /*end of postsmoothing*/
     
     }
     MEMFREE(r);
     MEMFREE(x0);

     return;
}
