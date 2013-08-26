
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/* Paso: AMG preconditioner  (local version)                  */

/************************************************************************************/

/* Author: artak@uq.edu.au, l.gross@uq.edu.au                 */

/************************************************************************************/

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


/************************************************************************************/

/* free all memory used by AMG                                */

void Paso_Preconditioner_AMG_free(Paso_Preconditioner_AMG * in) {
     if (in!=NULL) {
	Paso_Preconditioner_Smoother_free(in->Smoother);
	Paso_SystemMatrix_free(in->P);
	Paso_SystemMatrix_free(in->R);
	Paso_SystemMatrix_free(in->A_C);
	Paso_Preconditioner_AMG_free(in->AMG_C);
	delete[] in->r;
	delete[] in->x_C;
	delete[] in->b_C;
        Paso_MergedSolver_free(in->merged_solver);

	delete[] in;
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
/***************************************************************************************

   constructs AMG
   
****************************************************************************************/
Paso_Preconditioner_AMG* Paso_Preconditioner_AMG_alloc(Paso_SystemMatrix *A_p,dim_t level,Paso_Options* options) {

  Paso_Preconditioner_AMG* out=NULL;
  Paso_SystemMatrix *A_C=NULL;
  bool verbose=options->verbose;

  const dim_t my_n=A_p->mainBlock->numRows;
  const dim_t overlap_n=A_p->row_coupleBlock->numRows;
  
  const dim_t n = my_n + overlap_n;

  const dim_t n_block=A_p->row_block_size;
  bool* F_marker=NULL;
  index_t *counter=NULL, *mask_C=NULL, *rows_in_F;
  dim_t i, my_n_F, my_n_C, n_C, F_flag, *F_set=NULL, global_n_C=0, global_n_F=0, n_F;
  double time0=0;
  const double theta = options->coarsening_threshold;
  const double tau = options->diagonal_dominance_threshold;
  const double sparsity=Paso_SystemMatrix_getSparsity(A_p);
  const dim_t global_n=Paso_SystemMatrix_getGlobalNumRows(A_p);


  /*
      is the input matrix A suitable for coarsening?
      
  */
  if ( (sparsity >= options->min_coarse_sparsity) || 
       (global_n <= options->min_coarse_matrix_size) || 
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

	      if (global_n <= options->min_coarse_matrix_size)
	          printf("SIZE");

	      if (level > options->level_max)
	          printf("LEVEL");

	      printf("\n");

        printf("Paso_Preconditioner: AMG level %d (limit = %d) stopped. sparsity = %e (limit = %e), unknowns = %d (limit = %d)\n", 
	       level,  options->level_max, sparsity, options->min_coarse_sparsity, global_n, options->min_coarse_matrix_size);  

       } 

       return NULL;
  }  else {
     /* Start Coarsening : */
     
     /* this is the table for strong connections combining mainBlock, col_coupleBlock and row_coupleBlock */
     const dim_t len_S=A_p->mainBlock->pattern->len + A_p->col_coupleBlock->pattern->len + A_p->row_coupleBlock->pattern->len  + A_p->row_coupleBlock->numRows * A_p->col_coupleBlock->numCols;

     dim_t* degree_S=new  dim_t[n];
     index_t *offset_S=new  index_t[n];
     index_t *S=new  index_t[len_S];
     dim_t* degree_ST=new  dim_t[n];
     index_t *offset_ST=new  index_t[n];
     index_t *ST=new  index_t[len_S];
     
     
     F_marker=new bool[n];
     counter=new index_t[n];

     if ( !( Esys_checkPtr(F_marker) || Esys_checkPtr(counter) || Esys_checkPtr(degree_S) || Esys_checkPtr(offset_S) || Esys_checkPtr(S) 
	|| Esys_checkPtr(degree_ST) || Esys_checkPtr(offset_ST) || Esys_checkPtr(ST) ) ) {
	/*
	       make sure that corresponding values in the row_coupleBlock and col_coupleBlock are identical 
	*/
        Paso_SystemMatrix_copyColCoupleBlock(A_p);
        Paso_SystemMatrix_copyRemoteCoupleBlock(A_p, FALSE); 

	/* 
	      set splitting of unknowns:
	    
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

	 Paso_Preconditioner_AMG_CIJPCoarsening(n,my_n, F_marker,
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
	       count number of unknowns to be eliminated:
	    */
	    my_n_F=Paso_Util_cumsum_maskedTrue(my_n,counter, F_marker);
	    n_F=Paso_Util_cumsum_maskedTrue(n,counter, F_marker); 
            /* collect my_n_F values on all processes, a direct solver should 
                be used if any my_n_F value is 0 */
            F_set = new  dim_t[A_p->mpi_info->size];
	    #ifdef ESYS_MPI
            MPI_Allgather(&my_n_F, 1, MPI_INT, F_set, 1, MPI_INT, A_p->mpi_info->comm);
	    #endif
            global_n_F=0;
            F_flag = 1;
            for (i=0; i<A_p->mpi_info->size; i++) {
                global_n_F+=F_set[i];
                if (F_set[i] == 0) F_flag = 0;
            }
            delete[] F_set;

	    my_n_C=my_n-my_n_F;
            global_n_C=global_n-global_n_F;
	    if (verbose) printf("Paso_Preconditioner: AMG (non-local) level %d: %d unknowns are flagged for elimination. %d left.\n",level,global_n_F,global_n_C);

         
/*          if ( n_F == 0 ) {  is a nasty case. a direct solver should be used, return NULL */
            if (F_flag == 0) {
	       out = NULL;
	    } else {
	       out=new Paso_Preconditioner_AMG;
	       if (! Esys_checkPtr(out)) {
		  out->level = level;
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
                  out->merged_solver=NULL;
	       }
	       mask_C=new index_t[n];
	       rows_in_F=new index_t[n_F];
	       Esys_checkPtr(mask_C);
	       Esys_checkPtr(rows_in_F);
	       if ( Esys_noError() ) {

		  out->Smoother = Paso_Preconditioner_Smoother_alloc(A_p, (options->smoother == PASO_JACOBI), 0, verbose);
	  
		  if (global_n_C != 0) {
			/*  create mask of C nodes with value >-1, gives new id */
			n_C=Paso_Util_cumsum_maskedFalse(n, mask_C, F_marker);
			/* if nothing has been removed we have a diagonal dominant matrix and we just run a few steps of the smoother */ 
   
			/* allocate helpers :*/
			out->x_C=new double[n_block*my_n_C];
			out->b_C=new double[n_block*my_n_C];
			out->r=new double[n_block*my_n];
		     
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
			   construct courser level:
			   
			*/
			if ( Esys_noError()) {
			   out->AMG_C=Paso_Preconditioner_AMG_alloc(A_C,level+1,options);
			}

			if ( Esys_noError()) {
			  out->A_C=A_C;
			  if ( out->AMG_C == NULL ) { 
			      /* merge the system matrix into 1 rank when 
				 it's not suitable coarsening due to the 
				 total number of unknowns are too small */
                              out->merged_solver= Paso_MergedSolver_alloc(A_C, options);
			  }
			}		  
		  }
	       }
	       delete[] mask_C;
	       delete[] rows_in_F;
	    }
	 }

  }
  delete[] counter;
  delete[] F_marker;
  delete[] degree_S;
  delete[] offset_S;
  delete[] S;
  delete[] degree_ST;
  delete[] offset_ST;
  delete[] ST;
  
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
	
     time0=Esys_timer();

     Paso_Copy(n, amg->r, b);                            /*  r <- b */
     Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(-1.,A,x,1.,amg->r); /*r=r-Ax*/
     Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(1.,amg->R,amg->r,0.,amg->b_C);  /* b_c = R*r  */

     time0=Esys_timer()-time0;
     /* coarse level solve */
     if ( amg->AMG_C == NULL) {
	    time0=Esys_timer();
	    /*  A_C is the coarsest level */
            Paso_MergedSolver_solve(amg->merged_solver,amg->x_C, amg->b_C);
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

   threshold_p = new  double[2*my_n];
   
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
     delete[] threshold_p;
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
   
   
   threshold_p = new  double[2*my_n];

   #pragma omp parallel private(i,iptr,bi)
   {
   
      dim_t max_deg=0;
      double *rtmp=NULL;

      #pragma omp for schedule(static)
      for (i=0;i<my_n;++i) max_deg=MAX(max_deg, A->mainBlock->pattern->ptr[i+1]-A->mainBlock->pattern->ptr[i]
				     +A->col_coupleBlock->pattern->ptr[i+1]-A->col_coupleBlock->pattern->ptr[i]);
      
      rtmp=new  double[max_deg];
      
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
      delete[] rtmp;
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
   delete[] threshold_p;
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

void Paso_Preconditioner_AMG_CIJPCoarsening(const dim_t n, const dim_t my_n, bool*split_marker,
					    const dim_t* degree_S, const index_t* offset_S, const index_t* S,
					    const dim_t* degree_ST, const index_t* offset_ST, const index_t* ST,
					    Paso_Connector* col_connector, Paso_Distribution* col_dist) 
{
   dim_t i, numUndefined,   iter=0;
  index_t iptr, jptr, kptr;
  double *random=NULL, *w=NULL, *Status=NULL;
  index_t * ST_flag=NULL;

  Paso_Coupler* w_coupler=Paso_Coupler_alloc(col_connector  ,1);
   
  w=new  double[n];
  Status=new  double[n];
  random = Paso_Distribution_createRandomVector(col_dist,1);
  ST_flag=new  index_t[offset_ST[n-1]+ degree_ST[n-1]];

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

	    register bool inD=TRUE;
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
  delete[] random;
  delete[] w;
  delete[] Status;
  delete[] ST_flag;
  
  return;
}

