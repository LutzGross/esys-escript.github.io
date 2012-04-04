
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

/* Paso: AMG preconditioner  (local version)                  */

/**************************************************************/

/* Author: artak@uq.edu.au, l.gross@uq.edu.au                 */

/**************************************************************/

#define SHOW_TIMING FALSE
#define MY_DEBUG 0
#define MY_DEBUG1 1

#include "Paso.h"
#include "Preconditioner.h"
#include "Options.h"
#include "PasoUtil.h"
#include "UMFPACK.h"
#include "MKL.h"
#include<stdio.h>


/**************************************************************/

/* free all memory used by AMG                                */

void Paso_Preconditioner_AMG_free(Paso_Preconditioner_AMG * in) {
     if (in!=NULL) {
	Paso_Preconditioner_Smoother_free(in->Smoother);
	Paso_SystemMatrix_free(in->P);
	Paso_SystemMatrix_free(in->R);
	Paso_SystemMatrix_free(in->A_C);
	Paso_Preconditioner_AMG_free(in->AMG_C);
	MEMFREE(in->r);
	MEMFREE(in->x_C);
	MEMFREE(in->b_C);

	MEMFREE(in);
     }
}

index_t Paso_Preconditioner_AMG_getMaxLevel(const Paso_Preconditioner_AMG * in) {
   if (in->AMG_C == NULL) {
      return in->level; 
   } else {
      return Paso_Preconditioner_AMG_getMaxLevel(in->AMG_C);
   }
}
double Paso_Preconditioner_AMG_getCoarseLevelSparsity(const Paso_Preconditioner_AMG * in) {
      if (in->AMG_C == NULL) {
	 if (in->A_C == NULL) {
	    return 1.;
	 } else {
	    return Paso_SystemMatrix_getSparsity(in->A_C);
	 }
      } else {
	    return Paso_Preconditioner_AMG_getCoarseLevelSparsity(in->AMG_C);
      }
}
dim_t Paso_Preconditioner_AMG_getNumCoarseUnknwons(const Paso_Preconditioner_AMG * in) {
   if (in->AMG_C == NULL) {
      if (in->A_C == NULL) {
	 return 0;
      } else {
	 return Paso_SystemMatrix_getTotalNumRows(in->A_C);
      }
   } else {
 	 return Paso_Preconditioner_AMG_getNumCoarseUnknwons(in->AMG_C);
   }
}
/*****************************************************************

   constructs AMG
   
******************************************************************/
Paso_Preconditioner_AMG* Paso_Preconditioner_AMG_alloc(Paso_SystemMatrix *A_p,dim_t level,Paso_Options* options) {

  Paso_Preconditioner_AMG* out=NULL;
  Paso_SystemMatrix *A_C=NULL;
  bool_t verbose=options->verbose;

  const dim_t my_n=A_p->mainBlock->numRows;
  const dim_t overlap_n=A_p->row_coupleBlock->numRows;
  
  const dim_t n = my_n + overlap_n;

  const dim_t n_block=A_p->row_block_size;
  index_t* F_marker=NULL, *counter=NULL, *mask_C=NULL, *rows_in_F;
  dim_t i, n_F, n_C, F_flag, *F_set=NULL, total_n_C=0, total_n_F=0;
  double time0=0;
  const double theta = options->coarsening_threshold;
  const double tau = options->diagonal_dominance_threshold;
  const double sparsity=Paso_SystemMatrix_getSparsity(A_p);
  const dim_t total_n=Paso_SystemMatrix_getGlobalTotalNumRows(A_p);


  /*
      is the input matrix A suitable for coarsening?
      
  */
  if ( (sparsity >= options->min_coarse_sparsity) || 
       (total_n <= options->min_coarse_matrix_size) || 
       (level > options->level_max) ) {

        if (verbose) { 
	      /* 
	          print stopping condition:
                      - 'SPAR' = min_coarse_matrix_sparsity exceeded
                      - 'SIZE' = min_coarse_matrix_size exceeded
                      - 'LEVEL' = level_max exceeded
	      */
	      printf("Paso_Preconditioner: AMG: termination of coarsening by "); 

	      if (sparsity >= options->min_coarse_sparsity)
	          printf("SPAR");

	      if (total_n <= options->min_coarse_matrix_size)
	          printf("SIZE");

	      if (level > options->level_max)
	          printf("LEVEL");

	      printf("\n");

        printf("Paso_Preconditioner: AMG level %d (limit = %d) stopped. sparsity = %e (limit = %e), unknowns = %d (limit = %d)\n", 
	       level,  options->level_max, sparsity, options->min_coarse_sparsity, total_n, options->min_coarse_matrix_size);  

       } 

       return NULL;
  }  else {
     /* Start Coarsening : */
     
     /* this is the table for strong connections combining mainBlock, col_coupleBlock and row_coupleBlock */
     const dim_t len_S=A_p->mainBlock->pattern->len + A_p->col_coupleBlock->pattern->len + A_p->row_coupleBlock->pattern->len  + A_p->row_coupleBlock->numRows * A_p->col_coupleBlock->numCols;

     dim_t* degree_S=TMPMEMALLOC(n, dim_t);
     index_t *offset_S=TMPMEMALLOC(n, index_t);
     index_t *S=TMPMEMALLOC(len_S, index_t);
     dim_t* degree_ST=TMPMEMALLOC(n, dim_t);
     index_t *offset_ST=TMPMEMALLOC(n, index_t);
     index_t *ST=TMPMEMALLOC(len_S, index_t);
     
     
     F_marker=TMPMEMALLOC(n,index_t);
     counter=TMPMEMALLOC(n,index_t);

     if ( !( Esys_checkPtr(F_marker) || Esys_checkPtr(counter) || Esys_checkPtr(degree_S) || Esys_checkPtr(offset_S) || Esys_checkPtr(S) 
	|| Esys_checkPtr(degree_ST) || Esys_checkPtr(offset_ST) || Esys_checkPtr(ST) ) ) {
	/*
	       make sure that corresponding values in the row_coupleBlock and col_coupleBlock are identical 
	*/
        Paso_SystemMatrix_copyColCoupleBlock(A_p);
        Paso_SystemMatrix_copyRemoteCoupleBlock(A_p, FALSE); 

	/* 
	      set splitting of unknows:
	    
    	 */
	 time0=Esys_timer();
	 if (n_block>1) {
	       Paso_Preconditioner_AMG_setStrongConnections_Block(A_p, degree_S, offset_S, S, theta,tau);
	 } else {
	       Paso_Preconditioner_AMG_setStrongConnections(A_p, degree_S, offset_S, S, theta,tau);
	 }
	 Paso_Preconditioner_AMG_transposeStrongConnections(n, degree_S, offset_S, S, n, degree_ST, offset_ST, ST);
/*	 Paso_SystemMatrix_extendedRowsForST(A_p, degree_ST, offset_ST, ST);
 */

	 Paso_Preconditioner_AMG_CIJPCoarsening(n,my_n,F_marker,
					        degree_S, offset_S, S, degree_ST, offset_ST, ST,
						A_p->col_coupler->connector,A_p->col_distribution);
      

         /* in BoomerAMG if interpolation is used FF connectivity is required */
/*MPI:
         if (options->interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING) 
                             Paso_Preconditioner_AMG_enforceFFConnectivity(n, A_p->pattern->ptr, degree_S, S, F_marker);  
*/

	 options->coarsening_selection_time=Esys_timer()-time0 + MAX(0, options->coarsening_selection_time);
	 if (Esys_noError() ) {
	    #pragma omp parallel for private(i) schedule(static)
	    for (i = 0; i < n; ++i) F_marker[i]=(F_marker[i] ==  PASO_AMG_IN_F);
	 
	    /*
	       count number of unkowns to be eliminated:
	    */
	    n_F=Paso_Util_cumsum_maskedTrue(n,counter, F_marker);
            /* collect n_F values on all processes, a direct solver should 
                be used if any n_F value is 0 */
            F_set = TMPMEMALLOC(A_p->mpi_info->size, dim_t);
	    #ifdef ESYS_MPI
            MPI_Allgather(&n_F, 1, MPI_INT, F_set, 1, MPI_INT, A_p->mpi_info->comm);
	    #endif
            total_n_F=0;
            F_flag = 1;
            for (i=0; i<A_p->mpi_info->size; i++) {
                total_n_F+=F_set[i];
                if (F_set[i] == 0) {
                  F_flag = 0;
                  break;
                }
            }
            TMPMEMFREE(F_set);

	    n_C=n-n_F;
            total_n_C=total_n-total_n_F;
	    if (verbose) printf("Paso_Preconditioner: AMG (non-local) level %d: %d unknowns are flagged for elimination. %d left.\n",level,total_n_F,total_n_C);

         
/*          if ( n_F == 0 ) {  is a nasty case. a direct solver should be used, return NULL */
            if (F_flag == 0) {
	       out = NULL;
	    } else {
	       out=MEMALLOC(1,Paso_Preconditioner_AMG);
	       if (! Esys_checkPtr(out)) {
		  out->level = level;
		  out->n = n;
		  out->n_F = n_F;
		  out->n_block = n_block;
		  out->A_C = NULL; 
		  out->P = NULL;  
		  out->R = NULL; 	       
		  out->post_sweeps = options->post_sweeps;
		  out->pre_sweeps  = options->pre_sweeps;
		  out->r = NULL;
		  out->x_C = NULL;
		  out->b_C = NULL;
		  out->AMG_C = NULL;
		  out->Smoother=NULL;
	       }
	       mask_C=TMPMEMALLOC(n,index_t);
	       rows_in_F=TMPMEMALLOC(n_F,index_t);
	       Esys_checkPtr(mask_C);
	       Esys_checkPtr(rows_in_F);
	       if ( Esys_noError() ) {

		  out->Smoother = Paso_Preconditioner_Smoother_alloc(A_p, (options->smoother == PASO_JACOBI), 0, verbose);
	  
		  if (total_n_C != 0) {
			   /* if nothing has been removed we have a diagonal dominant matrix and we just run a few steps of the smoother */ 
   
			/* allocate helpers :*/
			out->x_C=MEMALLOC(n_block*n_C,double);
			out->b_C=MEMALLOC(n_block*n_C,double);
			out->r=MEMALLOC(n_block*n,double);
		     
			Esys_checkPtr(out->r);
			Esys_checkPtr(out->x_C);
			Esys_checkPtr(out->b_C);
		     
			if ( Esys_noError() ) {
			   /* creates index for F:*/
			   #pragma omp parallel private(i)
			   {
			      #pragma omp for schedule(static)
			      for (i = 0; i < n; ++i) {
				 if  (F_marker[i]) rows_in_F[counter[i]]=i;
			      }
			   }
			   /*  create mask of C nodes with value >-1, gives new id */
			   i=Paso_Util_cumsum_maskedFalse(n, mask_C, F_marker);
			   /*
			      get Prolongation :	 
			   */					
 
			   time0=Esys_timer();

			   out->P=Paso_Preconditioner_AMG_getProlongation(A_p,offset_S, degree_S,S,n_C,mask_C, options->interpolation_method);

			}

			/*      
			   construct Restriction operator as transposed of Prolongation operator: 
			*/

			if ( Esys_noError()) {
			   time0=Esys_timer();

			   out->R=Paso_Preconditioner_AMG_getRestriction(out->P);

			   if (SHOW_TIMING) printf("timing: level %d: Paso_SystemMatrix_getTranspose: %e\n",level,Esys_timer()-time0);
			}		
			/*
                        construct coarse level matrix:
                        */
                        if ( Esys_noError()) {
                           time0=Esys_timer();

                           A_C = Paso_Preconditioner_AMG_buildInterpolationOperator(A_p, out->P, out->R);

                           if (SHOW_TIMING) printf("timing: level %d : construct coarse matrix: %e\n",level,Esys_timer()-time0);
                        }

			/*
			   constructe courser level:
			   
			*/
			if ( Esys_noError()) {
			   out->AMG_C=Paso_Preconditioner_AMG_alloc(A_C,level+1,options);
			}

			if ( Esys_noError()) {
			  if ( out->AMG_C == NULL ) { 
			      /* merge the system matrix into 1 rank when 
				 it's not suitable coarsening due to the 
				 total number of unknowns are too small */
			      out->A_C=A_C;
			      out->reordering = options->reordering;
                              out->refinements = options->coarse_matrix_refinements;
			      out->verbose = verbose;
			      out->options_smoother = options->smoother;
			  } else {
			      /* finally we set some helpers for the solver step */
			      out->A_C=A_C;
			  }
			}		  
		  }
	       }
	       TMPMEMFREE(mask_C);
	       TMPMEMFREE(rows_in_F);
	    }
	 }

  }
  TMPMEMFREE(counter);
  TMPMEMFREE(F_marker);
  TMPMEMFREE(degree_S);
  TMPMEMFREE(offset_S);
  TMPMEMFREE(S);
  TMPMEMFREE(degree_ST);
  TMPMEMFREE(offset_ST);
  TMPMEMFREE(ST);
  
  }

  if (Esys_noError()) {
     return out;
  } else  {
     Paso_Preconditioner_AMG_free(out);
     return NULL;
  }
}


void Paso_Preconditioner_AMG_solve(Paso_SystemMatrix* A, Paso_Preconditioner_AMG * amg, double * x, double * b) {
     const dim_t n = A->mainBlock->numRows * A->mainBlock->row_block_size;
     double time0=0;
     const dim_t post_sweeps=amg->post_sweeps;
     const dim_t pre_sweeps=amg->pre_sweeps;

     /* presmoothing */
     time0=Esys_timer();
     Paso_Preconditioner_Smoother_solve(A, amg->Smoother, x, b, pre_sweeps, FALSE); 

     time0=Esys_timer()-time0;
     if (SHOW_TIMING) printf("timing: level %d: Presmoothing: %e\n",amg->level, time0); 
     /* end of presmoothing */
	
     if (amg->n_F < amg->n) { /* is there work on the coarse level? */
         time0=Esys_timer();

	 Paso_Copy(n, amg->r, b);                            /*  r <- b */
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(-1.,A,x,1.,amg->r); /*r=r-Ax*/
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(1.,amg->R,amg->r,0.,amg->b_C);  /* b_c = R*r  */

         time0=Esys_timer()-time0;
	 /* coarse level solve */
	 if ( amg->AMG_C == NULL) {
	    time0=Esys_timer();
	    /*  A_C is the coarsest level */
	    Paso_Preconditioner_AMG_mergeSolve(amg); 

	    if (SHOW_TIMING) printf("timing: level %d: DIRECT SOLVER: %e\n",amg->level,Esys_timer()-time0);
	 } else {
	    Paso_Preconditioner_AMG_solve(amg->A_C, amg->AMG_C,amg->x_C,amg->b_C); /* x_C=AMG(b_C)     */
	 }

  	 time0=time0+Esys_timer();
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(1.,amg->P,amg->x_C,1.,x); /* x = x + P*x_c */    

         /*postsmoothing*/
      
        /*solve Ax=b with initial guess x */
        time0=Esys_timer();
        Paso_Preconditioner_Smoother_solve(A, amg->Smoother, x, b, post_sweeps, TRUE); 
        time0=Esys_timer()-time0;
        if (SHOW_TIMING) printf("timing: level %d: Postsmoothing: %e\n",amg->level,time0);
        /*end of postsmoothing*/
     }

     return;
}

/* theta = threshold for strong connections */
/* tau = threshold for diagonal dominance */

/*S_i={j \in N_i; i strongly coupled to j}

in the sense that |A_{ij}| >= theta * max_k |A_{ik}| 
*/

void Paso_Preconditioner_AMG_setStrongConnections(Paso_SystemMatrix* A, 
	    				  dim_t *degree_S, index_t *offset_S, index_t *S,
					  const double theta, const double tau)
{

   const dim_t my_n=A->mainBlock->numRows;
   const dim_t overlap_n=A->row_coupleBlock->numRows;
   
   index_t iptr, i;
   double *threshold_p=NULL; 

   threshold_p = TMPMEMALLOC(2*my_n, double);
   
   #pragma omp parallel for private(i,iptr) schedule(static)
   for (i=0;i<my_n;++i) {        
	 
	 register double max_offdiagonal = 0.;
	 register double sum_row=0;
	 register double main_row=0;
	 register dim_t kdeg=0;
         register const index_t koffset=A->mainBlock->pattern->ptr[i]+A->col_coupleBlock->pattern->ptr[i];

         
	 /* collect information for row i: */
	 #pragma ivdep
	 for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
	    register index_t j=A->mainBlock->pattern->index[iptr];
	    register double fnorm=ABS(A->mainBlock->val[iptr]);
	    if( j != i) {
	       max_offdiagonal = MAX(max_offdiagonal,fnorm);
	       sum_row+=fnorm;
	    } else {
	       main_row=fnorm;
	    }

	 }

	 #pragma ivdep
	 for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	    register double fnorm=ABS(A->col_coupleBlock->val[iptr]);

	    max_offdiagonal = MAX(max_offdiagonal,fnorm);
	    sum_row+=fnorm;
	 }

         /* inspect row i: */
         {
	    const double threshold = theta*max_offdiagonal;
            threshold_p[2*i+1]=threshold;
	    if (tau*main_row < sum_row) { /* no diagonal dominance */
               threshold_p[2*i]=1;
	       #pragma ivdep
	       for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
	          register index_t j=A->mainBlock->pattern->index[iptr];
	          if(ABS(A->mainBlock->val[iptr])>threshold && i!=j) {
		     S[koffset+kdeg] = j;
		     kdeg++;
	          }
	       }
	       #pragma ivdep
	       for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	          register index_t j=A->col_coupleBlock->pattern->index[iptr];
	          if(ABS(A->col_coupleBlock->val[iptr])>threshold) {
		     S[koffset+kdeg] = j + my_n;
		     kdeg++;
	          }
	       }
            } else {
               threshold_p[2*i]=-1;
            }
         }
         offset_S[i]=koffset;
	 degree_S[i]=kdeg;
      }

      /* now we need to distribute the threshold and the diagonal dominance indicator */
      if (A->mpi_info->size > 1) {

          const index_t koffset_0=A->mainBlock->pattern->ptr[my_n]+A->col_coupleBlock->pattern->ptr[my_n]
				 -A->mainBlock->pattern->ptr[0]-A->col_coupleBlock->pattern->ptr[0];
	  
          double *remote_threshold=NULL;
          
	  Paso_Coupler* threshold_coupler=Paso_Coupler_alloc(A->row_coupler->connector  ,2);
          Paso_Coupler_startCollect(threshold_coupler,threshold_p);
          Paso_Coupler_finishCollect(threshold_coupler);
          remote_threshold=threshold_coupler->recv_buffer;

          #pragma omp parallel for private(i,iptr) schedule(static)
          for (i=0; i<overlap_n; i++) {
	      const double threshold = remote_threshold[2*i+1];
	      register dim_t kdeg=0;
              register const index_t koffset=koffset_0+A->row_coupleBlock->pattern->ptr[i]+A->remote_coupleBlock->pattern->ptr[i];
              if (remote_threshold[2*i]>0) {
		#pragma ivdep
		for (iptr=A->row_coupleBlock->pattern->ptr[i];iptr<A->row_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	          register index_t j=A->row_coupleBlock->pattern->index[iptr];
	          if(ABS(A->row_coupleBlock->val[iptr])>threshold) {
		     S[koffset+kdeg] = j ;
		     kdeg++;
	          }
		}

		#pragma ivdep
		for (iptr=A->remote_coupleBlock->pattern->ptr[i];iptr<A->remote_coupleBlock->pattern->ptr[i+1]; iptr++) {
		  register index_t j=A->remote_coupleBlock->pattern->index[iptr];
		  if(ABS(A->remote_coupleBlock->val[iptr])>threshold && i!=j) {
		      S[koffset+kdeg] = j + my_n;
		      kdeg++;
		  }
		}
              }
              offset_S[i+my_n]=koffset;
	      degree_S[i+my_n]=kdeg;
          }

          Paso_Coupler_free(threshold_coupler);
     }
     TMPMEMFREE(threshold_p);
}

/* theta = threshold for strong connections */
/* tau = threshold for diagonal dominance */ 
/*S_i={j \in N_i; i strongly coupled to j}

in the sense that |A_{ij}|_F >= theta * max_k |A_{ik}|_F 
*/
void Paso_Preconditioner_AMG_setStrongConnections_Block(Paso_SystemMatrix* A, 
							dim_t *degree_S, index_t *offset_S, index_t *S,
							const double theta, const double tau)

{
   const dim_t block_size=A->block_size;
   const dim_t my_n=A->mainBlock->numRows;
   const dim_t overlap_n=A->row_coupleBlock->numRows;
   
   index_t iptr, i, bi;
   double *threshold_p=NULL; 
   
   
   threshold_p = TMPMEMALLOC(2*my_n, double);

   #pragma omp parallel private(i,iptr,bi)
   {
   
      dim_t max_deg=0;
      double *rtmp=NULL;

      #pragma omp for schedule(static)
      for (i=0;i<my_n;++i) max_deg=MAX(max_deg, A->mainBlock->pattern->ptr[i+1]-A->mainBlock->pattern->ptr[i]
				     +A->col_coupleBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i]);
      
      rtmp=TMPMEMALLOC(max_deg, double);
      
      #pragma omp for schedule(static) 
      for (i=0;i<my_n;++i) {        
	 register double max_offdiagonal = 0.;
	 register double sum_row=0;
	 register double main_row=0;
	 register index_t rtmp_offset=-A->mainBlock->pattern->ptr[i];
	 register dim_t kdeg=0;
	 register const index_t koffset=A->mainBlock->pattern->ptr[i]+A->col_coupleBlock->pattern->ptr[i];
	 
	 /* collect information for row i: */
	 for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
	    register index_t j=A->mainBlock->pattern->index[iptr];
	    register double fnorm=0;
	    #pragma ivdep
	    for(bi=0;bi<block_size;++bi) {
   	        register double  rtmp2= A->mainBlock->val[iptr*block_size+bi];
	       fnorm+=rtmp2*rtmp2;
	    }
	    fnorm=sqrt(fnorm);
	    rtmp[iptr+rtmp_offset]=fnorm;
	    
	    if( j != i) {
	       max_offdiagonal = MAX(max_offdiagonal,fnorm);
	       sum_row+=fnorm;
	    } else {
	       main_row=fnorm;
	    }
	 
	 }
      
         rtmp_offset+=A->mainBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i];
	 for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	    register double fnorm=0;
	    #pragma ivdep
	    for(bi=0;bi<block_size;++bi) {
	       register double rtmp2 = A->col_coupleBlock->val[iptr*block_size+bi];
	       fnorm+=rtmp2*rtmp2;
	    }
	    fnorm=sqrt(fnorm);
	    
	    rtmp[iptr+rtmp_offset]=fnorm;
	    max_offdiagonal = MAX(max_offdiagonal,fnorm);
	    sum_row+=fnorm;
	 }
	 
      
	 /* inspect row i: */
	 {
	    const double threshold = theta*max_offdiagonal;
	    rtmp_offset=-A->mainBlock->pattern->ptr[i];
	    
	    threshold_p[2*i+1]=threshold;
	    if (tau*main_row < sum_row) { /* no diagonal dominance */
	       threshold_p[2*i]=1;
	       #pragma ivdep
	       for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
		  register index_t j=A->mainBlock->pattern->index[iptr];
		  if(rtmp[iptr+rtmp_offset] > threshold && i!=j) {
		     S[koffset+kdeg] = j;
		     kdeg++;
		  }
	       }
	       rtmp_offset+=A->mainBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i];
	       #pragma ivdep
	       for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
		  register index_t j=A->col_coupleBlock->pattern->index[iptr];
		  if( rtmp[iptr+rtmp_offset] >threshold) {
		     S[koffset+kdeg] = j + my_n;
		     kdeg++;
		  }
	       }
	    } else {
	       threshold_p[2*i]=-1;
	    }
	 }
	 offset_S[i]=koffset;
	 degree_S[i]=kdeg;
      }
      TMPMEMFREE(rtmp);
   }
   /* now we need to distribute the threshold and the diagonal dominance indicator */
   if (A->mpi_info->size > 1) {
      
      const index_t koffset_0=A->mainBlock->pattern->ptr[my_n]+A->col_coupleBlock->pattern->ptr[my_n]
                             -A->mainBlock->pattern->ptr[0]-A->col_coupleBlock->pattern->ptr[0];
      
      double *remote_threshold=NULL;
      
      Paso_Coupler* threshold_coupler=Paso_Coupler_alloc(A->row_coupler->connector  ,2);
      Paso_Coupler_startCollect(threshold_coupler,threshold_p);
      Paso_Coupler_finishCollect(threshold_coupler);
      remote_threshold=threshold_coupler->recv_buffer;
      
      #pragma omp parallel for private(i,iptr) schedule(static)
      for (i=0; i<overlap_n; i++) {
	 
	 const double threshold2 = remote_threshold[2*i+1]*remote_threshold[2*i+1];
	 register dim_t kdeg=0;
	 register const index_t koffset=koffset_0+A->row_coupleBlock->pattern->ptr[i]+A->remote_coupleBlock->pattern->ptr[i];
	 if (remote_threshold[2*i]>0) {
	    #pragma ivdep
	    for (iptr=A->row_coupleBlock->pattern->ptr[i];iptr<A->row_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	       register index_t j=A->row_coupleBlock->pattern->index[iptr];
	       register double fnorm2=0;
	       #pragma ivdepremote_threshold[2*i]
	       for(bi=0;bi<block_size;++bi) {
		  register double rtmp2 = A->row_coupleBlock->val[iptr*block_size+bi];
		  fnorm2+=rtmp2*rtmp2;
	       }
	       
	       if(fnorm2 > threshold2 ) {
		  S[koffset+kdeg] = j ;
		  kdeg++;
	       }
	    }

	    #pragma ivdep
            for (iptr=A->remote_coupleBlock->pattern->ptr[i];iptr<A->remote_coupleBlock->pattern->ptr[i+1]; ++iptr) {
               register index_t j=A->remote_coupleBlock->pattern->index[iptr];
               register double fnorm2=0;
               #pragma ivdepremote_threshold[2*i]
               for(bi=0;bi<block_size;++bi) {
                  register double rtmp2 = A->remote_coupleBlock->val[iptr*block_size+bi];
                  fnorm2+=rtmp2*rtmp2;
               }
               if(fnorm2 > threshold2 && i != j) {
                  S[koffset+kdeg] = j + my_n;
                  kdeg++;
               }
            }
	    
	 }
	 offset_S[i+my_n]=koffset;
	 degree_S[i+my_n]=kdeg;
      }
      Paso_Coupler_free(threshold_coupler);
   }
   TMPMEMFREE(threshold_p);
}

void Paso_Preconditioner_AMG_transposeStrongConnections(const dim_t n, const dim_t* degree_S, const index_t* offset_S, const index_t* S,
							const dim_t nT, dim_t* degree_ST, index_t* offset_ST,index_t* ST)
{
   index_t i, j;
   dim_t p;
   dim_t len=0;
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<nT ;++i) {
      degree_ST[i]=0;
   }
   for (i=0; i<n ;++i) {
      for (p=0; p<degree_S[i]; ++p) degree_ST[ S[offset_S[i]+p] ]++;
   }
   for (i=0; i<nT ;++i) {
      offset_ST[i]=len;
      len+=degree_ST[i];
      degree_ST[i]=0;
   }
   for (i=0; i<n ;++i) {
      for (p=0; p<degree_S[i]; ++p) {
	   j=S[offset_S[i]+p];
	   ST[offset_ST[j]+degree_ST[j]]=i;
	   degree_ST[j]++;
      }
   }
}

int compareindex(const void *a, const void *b)
{
  return (*(int *)a - *(int *)b);
}

void Paso_Preconditioner_AMG_CIJPCoarsening(const dim_t n, const dim_t my_n, index_t*split_marker,
					    const dim_t* degree_S, const index_t* offset_S, const index_t* S,
					    const dim_t* degree_ST, const index_t* offset_ST, const index_t* ST,
					    Paso_Connector* col_connector, Paso_Distribution* col_dist) 
{
   dim_t i, numUndefined,   iter=0;
  index_t iptr, jptr, kptr;
  double *random=NULL, *w=NULL, *Status=NULL;
  index_t * ST_flag=NULL;

  Paso_Coupler* w_coupler=Paso_Coupler_alloc(col_connector  ,1);
   
  w=TMPMEMALLOC(n, double);
  Status=TMPMEMALLOC(n, double);
  random = Paso_Distribution_createRandomVector(col_dist,1);
  ST_flag=TMPMEMALLOC(offset_ST[n-1]+ degree_ST[n-1], index_t);

  #pragma omp parallel for private(i)
  for (i=0; i< my_n; ++i) {
      w[i]=degree_ST[i]+random[i];
      if (degree_ST[i] < 1) {
	   Status[i]=-100; /* F point  */
      } else {
	   Status[i]=1; /* status undefined */
      }
  }

  #pragma omp parallel for private(i, iptr)
  for (i=0; i< n; ++i) {
      for( iptr =0 ; iptr < degree_ST[i]; ++iptr)  {
	 ST_flag[offset_ST[i]+iptr]=1;
      }
  }   

  
  numUndefined = Paso_Distribution_numPositives(Status, col_dist, 1 );
  /* printf(" coarsening loop start: num of undefined rows = %d \n",numUndefined);  */
  iter=0; 
  while (numUndefined > 0) {
     Paso_Coupler_fillOverlap(n, w, w_coupler);

      /* calculate the maximum value of neighbours following active strong connections:
	    w2[i]=MAX(w[k]) with k in ST[i] or S[i] and (i,k) connection is still active  */       
      #pragma omp parallel for private(i, iptr)
      for (i=0; i<my_n; ++i) {
	 if (Status[i]>0) { /* status is still undefined */

	    register bool_t inD=TRUE;
	    const double wi=w[i];

	    for( iptr =0 ; iptr < degree_S[i]; ++iptr) {
	       const index_t k=S[offset_S[i]+iptr];
	       const index_t* start_p = &ST[offset_ST[k]];
	       const index_t* where_p=(index_t*)bsearch(&i, start_p, degree_ST[k], sizeof(index_t), Paso_comparIndex);

	       if (ST_flag[offset_ST[k] + (index_t)(where_p-start_p)]>0) {
		  if (wi <= w[k] ) {
		     inD=FALSE;
		     break;
		  }
	       }
	    
	    }
	    
	    if (inD) {
		  for( iptr =0 ; iptr < degree_ST[i]; ++iptr) {
		     const index_t k=ST[offset_ST[i]+iptr];
		     if ( ST_flag[offset_ST[i]+iptr] > 0 ) {
                       if (wi <= w[k] ) {
			   inD=FALSE;
			   break;
			}
		     }
		  }
	    }    
	    if (inD) { 
	       Status[i]=0.; /* is in D */
	    }
	 }
	 
      }

      Paso_Coupler_fillOverlap(n, Status, w_coupler);


	 /*   remove connection to D points : 
	 
	       for each i in D:
		  for each j in S_i:
		     w[j]--
		     ST_tag[j,i]=-1
		  for each j in ST[i]:
		     ST_tag[i,j]=-1
		     for each k in ST[j]:
			if k in ST[i]:
			   w[j]--;
			ST_tag[j,k]=-1
			
	 */
	 /* w is updated  for local rows only */
	 {
	    #pragma omp parallel for private(i, jptr)
	    for (i=0; i< my_n; ++i) {

	       for (jptr=0; jptr<degree_ST[i]; ++jptr) {
		  const index_t j=ST[offset_ST[i]+jptr];
		  if ( (Status[j] == 0.) && (ST_flag[offset_ST[i]+jptr]>0) ) {
		     w[i]--;
		     ST_flag[offset_ST[i]+jptr]=-1;
		  }
	       }
	       
	    } 
	    #pragma omp parallel for private(i, jptr)
	    for (i=my_n; i< n; ++i) {
	       for (jptr=0; jptr<degree_ST[i]; ++jptr) {
		  const index_t j = ST[offset_ST[i]+jptr];
		  if ( Status[j] == 0. ) ST_flag[offset_ST[i]+jptr]=-1;
	       }
	    }

	    
	    for (i=0; i< n; ++i) {
	       if ( Status[i] == 0. ) {

                     const index_t* start_p = &ST[offset_ST[i]];

		     for (jptr=0; jptr<degree_ST[i]; ++jptr) {
			const index_t j=ST[offset_ST[i]+jptr];
			ST_flag[offset_ST[i]+jptr]=-1;
			for (kptr=0; kptr<degree_ST[j]; ++kptr) {
			   const index_t k=ST[offset_ST[j]+kptr]; 
			   if (NULL != bsearch(&k, start_p, degree_ST[i], sizeof(index_t), Paso_comparIndex) ) { /* k in ST[i] ? */
			      if (ST_flag[offset_ST[j]+kptr] >0) {
				 if (j< my_n ) {
				    w[j]--;
				 }
				 ST_flag[offset_ST[j]+kptr]=-1;
			      }
			   }
			}
		     }
		  }
	    }
	 }

	 /* adjust status */
	 #pragma omp parallel for private(i)
	 for (i=0; i< my_n; ++i) {
	    if ( Status[i] == 0. ) {
	       Status[i] = -10;   /* this is now a C point */
	    } else if (Status[i] == 1. && w[i]<1.) {
	       Status[i] = -100;   /* this is now a F point */  
	    }
	 }
	 
	 i = numUndefined;
	 numUndefined = Paso_Distribution_numPositives(Status, col_dist, 1 );
	 if (numUndefined == i) {
	   Esys_setError(SYSTEM_ERROR, "Can NOT reduce numUndefined."); 
	   return;
	 }

	 iter++;
	 /* printf(" coarsening loop %d: num of undefined rows = %d \n",iter, numUndefined); */

  } /* end of while loop */

  /* map to output :*/
  Paso_Coupler_fillOverlap(n, Status, w_coupler);
  #pragma omp parallel for private(i)
  for (i=0; i< n; ++i) {
	 if (Status[i] > -50.) {
	    split_marker[i]=PASO_AMG_IN_C;
	 } else {
	    split_marker[i]=PASO_AMG_IN_F;
	 }
  }
  /* clean up : */
  Paso_Coupler_free(w_coupler);
  TMPMEMFREE(random);
  TMPMEMFREE(w);
  TMPMEMFREE(Status);
  TMPMEMFREE(ST_flag);
  
  return;
}

/* Merge the system matrix which is distributed on ranks into a complete 
   matrix on rank 0, then solve this matrix on rank 0 only */
Paso_SparseMatrix* Paso_Preconditioner_AMG_mergeSystemMatrix(Paso_SystemMatrix* A) {
  index_t i, iptr, j, n, remote_n, total_n, len, offset, tag;
  index_t row_block_size, col_block_size, block_size;
  index_t size=A->mpi_info->size;
  index_t rank=A->mpi_info->rank;
  index_t *ptr=NULL, *idx=NULL, *ptr_global=NULL, *idx_global=NULL;
  index_t *temp_n=NULL, *temp_len=NULL;
  double  *val=NULL;
  Paso_Pattern *pattern=NULL;
  Paso_SparseMatrix *out=NULL;
  #ifdef ESYS_MPI
    MPI_Request* mpi_requests=NULL;
    MPI_Status* mpi_stati=NULL;
  #else
    int *mpi_requests=NULL, *mpi_stati=NULL;
  #endif

  if (size == 1) {
    n = A->mainBlock->numRows;
    ptr = TMPMEMALLOC(n, index_t); 
    #pragma omp parallel for private(i)
    for (i=0; i<n; i++) ptr[i] = i;
    out = Paso_SparseMatrix_getSubmatrix(A->mainBlock, n, n, ptr, ptr);
    TMPMEMFREE(ptr);
    return out;
  }

  n = A->mainBlock->numRows;
  block_size = A->block_size;

  /* Merge MainBlock and CoupleBlock to get a complete column entries
     for each row allocated to current rank. Output (ptr, idx, val) 
     contains all info needed from current rank to merge a system 
     matrix  */
  Paso_SystemMatrix_mergeMainAndCouple(A, &ptr, &idx, &val);

  #ifdef ESYS_MPI
    mpi_requests=TMPMEMALLOC(size*2,MPI_Request);
    mpi_stati=TMPMEMALLOC(size*2,MPI_Status);
  #else
    mpi_requests=TMPMEMALLOC(size*2,int);
    mpi_stati=TMPMEMALLOC(size*2,int);
  #endif

  /* Now, pass all info to rank 0 and merge them into one sparse 
     matrix */
  if (rank == 0) {
    /* First, copy local ptr values into ptr_global */
    total_n=Paso_SystemMatrix_getGlobalNumRows(A);
    ptr_global = MEMALLOC(total_n+1, index_t);
    memcpy(ptr_global, ptr, (n+1) * sizeof(index_t));
    iptr = n+1;
    MEMFREE(ptr);
    temp_n = TMPMEMALLOC(size, index_t);
    temp_len = TMPMEMALLOC(size, index_t);
    temp_n[0] = iptr;
    
    /* Second, receive ptr values from other ranks */
    for (i=1; i<size; i++) {
      remote_n = A->row_distribution->first_component[i+1] -
		 A->row_distribution->first_component[i];
      #ifdef ESYS_MPI
      MPI_Irecv(&(ptr_global[iptr]), remote_n, MPI_INT, i, 
			A->mpi_info->msg_tag_counter+i,
			A->mpi_info->comm,
			&mpi_requests[i]);
      #endif
      temp_n[i] = remote_n;
      iptr += remote_n;
    }
    #ifdef ESYS_MPI
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter += size;

    /* Then, prepare to receive idx and val from other ranks */
    len = 0;
    offset = -1;
    for (i=0; i<size; i++) {
      if (temp_n[i] > 0) {
	offset += temp_n[i];
	len += ptr_global[offset];
	temp_len[i] = ptr_global[offset];
      }else 
	temp_len[i] = 0;
    }

    idx_global = MEMALLOC(len, index_t);
    iptr = temp_len[0];
    offset = n+1;
    for (i=1; i<size; i++) {
      len = temp_len[i];
      #ifdef ESYS_MPI
      MPI_Irecv(&(idx_global[iptr]), len, MPI_INT, i,
			A->mpi_info->msg_tag_counter+i,
			A->mpi_info->comm,
			&mpi_requests[i]);
      #endif
      remote_n = temp_n[i];
      for (j=0; j<remote_n; j++) {
	ptr_global[j+offset] = ptr_global[j+offset] + iptr;
      }
      offset += remote_n;
      iptr += len;
    }
    memcpy(idx_global, idx, temp_len[0] * sizeof(index_t));
    MEMFREE(idx);
    row_block_size = A->mainBlock->row_block_size;
    col_block_size = A->mainBlock->col_block_size;
    #ifdef ESYS_MPI
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter += size;
    TMPMEMFREE(temp_n);

    /* Then generate the sparse matrix */
    pattern = Paso_Pattern_alloc(A->mainBlock->pattern->type, total_n,
			total_n, ptr_global, idx_global);
    out = Paso_SparseMatrix_alloc(A->mainBlock->type, pattern, 
			row_block_size, col_block_size, FALSE);
    Paso_Pattern_free(pattern);

    /* Finally, receive and copy the value */
    iptr = temp_len[0] * block_size;
    for (i=1; i<size; i++) {
      len = temp_len[i];
      #ifdef ESYS_MPI
      MPI_Irecv(&(out->val[iptr]), len * block_size, MPI_DOUBLE, i,
                        A->mpi_info->msg_tag_counter+i,
                        A->mpi_info->comm,
                        &mpi_requests[i]);
      #endif
      iptr += (len * block_size);
    }
    memcpy(out->val, val, temp_len[0] * sizeof(double) * block_size);
    MEMFREE(val);
    #ifdef ESYS_MPI
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter += size;
    TMPMEMFREE(temp_len);
  } else { /* it's not rank 0 */

    /* First, send out the local ptr */
    tag = A->mpi_info->msg_tag_counter+rank;
    #ifdef ESYS_MPI
    MPI_Issend(&(ptr[1]), n, MPI_INT, 0, tag, A->mpi_info->comm, 
			&mpi_requests[0]);
    #endif

    /* Next, send out the local idx */
    len = ptr[n];
    tag += size;
    #ifdef ESYS_MPI
    MPI_Issend(idx, len, MPI_INT, 0, tag, A->mpi_info->comm, 
			&mpi_requests[1]);
    #endif

    /* At last, send out the local val */
    len *= block_size;
    tag += size;
    #ifdef ESYS_MPI
    MPI_Issend(val, len, MPI_DOUBLE, 0, tag, A->mpi_info->comm, 
			&mpi_requests[2]);

    MPI_Waitall(3, mpi_requests, mpi_stati);
    #endif
    A->mpi_info->msg_tag_counter = tag + size - rank;
    MEMFREE(ptr);
    MEMFREE(idx);
    MEMFREE(val);

    out = NULL;
  }

  TMPMEMFREE(mpi_requests);
  TMPMEMFREE(mpi_stati);
  return out;
}


void Paso_Preconditioner_AMG_mergeSolve(Paso_Preconditioner_AMG * amg) {
  Paso_SystemMatrix *A = amg->A_C;
  Paso_SparseMatrix *A_D, *A_temp;
  double* x=NULL;
  double* b=NULL;
  index_t rank = A->mpi_info->rank;
  index_t size = A->mpi_info->size;
  index_t i, n, p, n_block;
  index_t *counts, *offset, *dist;
  #ifdef ESYS_MPI
  index_t count;
  #endif
  n_block = amg->n_block;
  A_D = Paso_Preconditioner_AMG_mergeSystemMatrix(A); 

  /* First, gather x and b into rank 0 */
  dist = A->pattern->input_distribution->first_component;
  n = Paso_SystemMatrix_getGlobalNumRows(A);
  b = TMPMEMALLOC(n*n_block, double);
  x = TMPMEMALLOC(n*n_block, double);
  counts = TMPMEMALLOC(size, index_t);
  offset = TMPMEMALLOC(size, index_t);

  #pragma omp parallel for private(i,p)
  for (i=0; i<size; i++) {
    p = dist[i];
    counts[i] = (dist[i+1] - p)*n_block;
    offset[i] = p*n_block;
  }
  #ifdef ESYS_MPI
  count = counts[rank];
  MPI_Gatherv(amg->b_C, count, MPI_DOUBLE, b, counts, offset, MPI_DOUBLE, 0, A->mpi_info->comm);
  MPI_Gatherv(amg->x_C, count, MPI_DOUBLE, x, counts, offset, MPI_DOUBLE, 0, A->mpi_info->comm);
  #endif

  if (rank == 0) {
    /* solve locally */
    #ifdef MKL
      A_temp = Paso_SparseMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, A_D);
      A_temp->solver_package = PASO_MKL;
      Paso_SparseMatrix_free(A_D);
      Paso_MKL(A_temp, x, b, amg->reordering, amg->refinements, SHOW_TIMING);
      Paso_SparseMatrix_free(A_temp);
    #else 
      #ifdef UMFPACK
	A_temp = Paso_SparseMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_CSC, A_D);
	A_temp->solver_package = PASO_UMFPACK;
	Paso_SparseMatrix_free(A_D);
	Paso_UMFPACK(A_temp, x, b, amg->refinements, SHOW_TIMING);
	Paso_SparseMatrix_free(A_temp);
      #else
	A_D->solver_p = Paso_Preconditioner_LocalSmoother_alloc(A_D, (amg->options_smoother == PASO_JACOBI), amg->verbose);
	A_D->solver_package = PASO_SMOOTHER;
	Paso_Preconditioner_LocalSmoother_solve(A_D, A_D->solver_p, x, b, amg->pre_sweeps+amg->post_sweeps, FALSE);
	Paso_SparseMatrix_free(A_D);
      #endif
    #endif
  }

  #ifdef ESYS_MPI
  /* now we need to distribute the solution to all ranks */
  MPI_Scatterv(x, counts, offset, MPI_DOUBLE, amg->x_C, count, MPI_DOUBLE, 0, A->mpi_info->comm);
  #endif

  TMPMEMFREE(x);
  TMPMEMFREE(b);
  TMPMEMFREE(counts);
  TMPMEMFREE(offset);
}
