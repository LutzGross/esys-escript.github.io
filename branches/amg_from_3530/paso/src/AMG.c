
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

/* Author: artak@uq.edu.au, l.gross@uq.edu.au                                */

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
  Paso_SystemMatrix *A_C;
  bool_t verbose=options->verbose;

  const dim_t my_n=A_p->mainBlock->numRows;
  const dim_t overlap_n=A_p->row_coupleBlock->numRows;
  
  const dim_t n = my_n + overlap_n;

  const dim_t n_block=A_p->row_block_size;
//  const dim_t n_block=A_p->block_size;
  index_t* F_marker=NULL, *counter=NULL, *mask_C=NULL, *rows_in_F;
  dim_t i, n_F, n_C, F_flag, *F_set=NULL;
  double time0=0;
  const double theta = options->coarsening_threshold;
  const double tau = options->diagonal_dominance_threshold;
  const double sparsity=Paso_SystemMatrix_getSparsity(A_p);
  const dim_t total_n=Paso_SystemMatrix_getGlobalTotalNumRows(A_p);

  
  /*
      is the input matrix A suitable for coarsening
      
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
	          printf("SPAR ");

	      if (total_n <= options->min_coarse_matrix_size)
	          printf("SIZE ");

	      if (level > options->level_max)
	          printf("LEVEL ");

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
if (MY_DEBUG1) fprintf(stderr, "after copy Row-couple Block, before remote-couple Block for A\n");

        Paso_SystemMatrix_copyRemoteCoupleBlock(A_p, FALSE); 
if (MY_DEBUG1) fprintf(stderr, "after copy Remote-couple Block for A\n");

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
if (MY_DEBUG) fprintf(stderr, "rand %d after set strong connections\n", A_p->mpi_info->rank);
//	 Paso_SystemMatrix_extendedRowsForST(A_p, degree_ST, offset_ST, ST);

if (MY_DEBUG) {
int p;
char *str1, *str2;
str1 = TMPMEMALLOC(10000, char);
str2 = TMPMEMALLOC(15, char);
sprintf(str1, "rank %d n=%d my_n=%d\n", A_p->mpi_info->rank, n, my_n);
for (i=0; i<n; ++i) {
  sprintf(str2, "%d(%d,%d): ", i, offset_S[i], degree_S[i]);
  strcat(str1, str2);
  for (p=0; p<degree_S[i];++p) {
     sprintf(str2, "%d ",S[offset_S[i]+p]);
     strcat(str1, str2);
  }
  sprintf(str1, "%s\n", str1);
}
fprintf(stderr, "%s", str1);
/*for (i=0; i<n; ++i) {
  sprintf(str2, "%d(%d,%d): ", i, offset_ST[i], degree_ST[i]);
  strcat(str1, str2);
  for (p=0; p<degree_ST[i];++p) {
     sprintf(str2, "%d ",ST[offset_ST[i]+p]);
     strcat(str1, str2);
  }
  sprintf(str1, "%s\n", str1);
}
fprintf(stderr, "Now ST: %s", str1);*/

TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

/*if (MY_DEBUG && level ==2) Paso_SystemMatrix_print(A_p); */
	 Paso_Preconditioner_AMG_CIJPCoarsening(n,my_n,F_marker,
					        degree_S, offset_S, S, degree_ST, offset_ST, ST,
						A_p->col_coupler->connector,A_p->col_distribution);
if (MY_DEBUG) fprintf(stderr, "after CIJPCoarsening\n");
      

/*MPI: 
	 Paso_Preconditioner_AMG_RungeStuebenSearch(n, A_p->pattern->ptr, degree_S, S, F_marker, options->usePanel);
*/
	 
         /* in BoomerAMG interpolation is used FF connectiovity is required :*/
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
	    n_C=n-n_F;
	    if (verbose) printf("Paso_Preconditioner: AMG (non-local) level %d: %d unknowns are flagged for elimination. %d left.\n",level,n_F,n-n_F);

            /* collect n_F values on all processes, a direct solver should 
                be used if any n_F value is 0 */
            F_set = TMPMEMALLOC(A_p->mpi_info->size, dim_t);
            MPI_Allgather(&n_F, 1, MPI_INT, F_set, 1, MPI_INT, A_p->mpi_info->comm);
            F_flag = 1;
            for (i=0; i<A_p->mpi_info->size; i++) {
                if (F_set[i] == 0) {
                  F_flag = 0;
                  break;
                }
            }
            TMPMEMFREE(F_set);
if (MY_DEBUG1) fprintf(stderr, "rank %d F_FLAG %d N_F %d\n", A_p->mpi_info->rank, F_flag, n_F);
         
//          if ( n_F == 0 ) {  /*  is a nasty case. a direct solver should be used, return NULL */
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
if (MY_DEBUG1 && Esys_getErrorType())
fprintf(stderr, "rank %d %s\n", A_p->mpi_info->rank, Esys_getErrorMessage());
else 
fprintf(stderr, "rank %d No Zero on Diagonal\n", A_p->mpi_info->rank);
	  
		  if (n_C != 0) {
			   /* if nothing is been removed we have a diagonal dominant matrix and we just run a few steps of the smoother */ 
   
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
			   /*  create mask of C nodes with value >-1 gives new id */
			   i=Paso_Util_cumsum_maskedFalse(n, mask_C, F_marker);
			   /*
			      get Prolongation :	 
			   */					
			   time0=Esys_timer();

if (MY_DEBUG && level == 2) {
int p;
char *str1, *str2;
str1 = TMPMEMALLOC(10000, char);
str2 = TMPMEMALLOC(15, char);
sprintf(str1, "rank %d n=%d my_n=%d i=%d\n", A_p->mpi_info->rank, n, my_n, i);
for (i=0; i<n; ++i) {
  sprintf(str2, "%d(%d,%d,%d): ", i, offset_S[i], degree_S[i], mask_C[i]);
  strcat(str1, str2);
  for (p=0; p<degree_S[i];++p) {
     sprintf(str2, "%d ",S[offset_S[i]+p]);
     strcat(str1, str2);
  }
  sprintf(str1, "%s\n", str1);
}
fprintf(stderr, "%s", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
} 
if (MY_DEBUG1) fprintf(stderr, "Rank %d before get Prolongation level %d\n", A_p->mpi_info->rank, level);
/*if (level == 1){
    int i=0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    fprintf(stderr, "rank %d PID %d on %s ready for attach\n", A_p->mpi_info->rank, getpid(), hostname);
    if (A_p->mpi_info->rank == 0) {
      fflush(stdout);
      while (0 == i)
        sleep(5);
    }
}*/
			   out->P=Paso_Preconditioner_AMG_getProlongation(A_p,offset_S, degree_S,S,n_C,mask_C, options->interpolation_method);
if (MY_DEBUG1) fprintf(stderr, "after get Prolongation %d\n", out->P->mpi_info->rank);
//Esys_setError(SYSTEM_ERROR, "AMG:DONE.");
//return NULL;

			   if (SHOW_TIMING) printf("timing: level %d: getProlongation: %e\n",level, Esys_timer()-time0);
			}
			/*      
			   construct Restriction operator as transposed of Prolongation operator: 
			*/
			if ( Esys_noError()) {
			   time0=Esys_timer();

//			   out->R=Paso_SystemMatrix_getTranspose(out->P);
			   out->R=Paso_Preconditioner_AMG_getRestriction(out->P);
if (MY_DEBUG && level == 3) {
int rank=A_p->mpi_info->rank;
if (rank == 3) fprintf(stderr, "SEND %d TO %d\n", out->R->col_coupler->connector->send->offsetInShared[2] - out->R->col_coupler->connector->send->offsetInShared[1],
out->R->col_coupler->connector->send->neighbor[1]);
if (rank == 4) fprintf(stderr, "RECV %d FM %d\n", out->R->col_coupler->connector->recv->offsetInShared[1] - out->R->col_coupler->connector->recv->offsetInShared[0],
out->R->col_coupler->connector->recv->neighbor[0]);
}


			   if (SHOW_TIMING) printf("timing: level %d: Paso_SystemMatrix_getTranspose: %e\n",level,Esys_timer()-time0);
			}		
			/* 
			construct coarse level matrix:
			*/
			if ( Esys_noError()) {
			   time0=Esys_timer();

if (MY_DEBUG1) fprintf(stderr, "before buildInterpolation\n");
/*if (A_p->mpi_info->rank == 0) {
  Paso_SystemMatrix_print(A_p);
fprintf(stderr, "Now Restriction Matrix ************************\n");
  Paso_SystemMatrix_print(out->R);
}*/
if (MY_DEBUG) {
char *str1, *str2;
int sum, rank, i, iPtr;
sum = A_p->mainBlock->numRows;
rank = A_p->mpi_info->rank;
str1 = TMPMEMALLOC(sum*sum*20+100, char);
str2 = TMPMEMALLOC(20, char);
sprintf(str1, "rank %d A Main Block: %d rows\n-----------\n", rank, sum);
for (i=0; i<sum; i++) {
  sprintf(str2, "Row %d: ",i);
  strcat(str1, str2);
  for (iPtr =A_p->mainBlock->pattern->ptr[i]; iPtr<A_p->mainBlock->pattern->ptr[i+1]; ++iPtr) {
         sprintf(str2, "(%d %f),",A_p->mainBlock->pattern->index[iPtr], A_p->mainBlock->val[iPtr]);
         strcat(str1, str2);
  }
  sprintf(str1, "%s\n", str1);
}
fprintf(stderr, "%s\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

			   A_C = Paso_Preconditioner_AMG_buildInterpolationOperator(A_p, out->P, out->R);

if (MY_DEBUG && A_p->mpi_info->size == 1) {
char *str1, *str2;
int sum, rank, i, iPtr;
Paso_SparseMatrix *Atemp=NULL, *A_C=NULL;
Atemp = Paso_SparseMatrix_MatrixMatrix(A_p->mainBlock,out->P->mainBlock);
A_C = Paso_SparseMatrix_MatrixMatrix(out->R->mainBlock, Atemp);
Paso_SparseMatrix_free(Atemp);
sum = A_C->numRows;
str1 = TMPMEMALLOC(sum*sum*20+100, char);
str2 = TMPMEMALLOC(20, char);
sprintf(str1, "A_C: %d rows\n-----------\n", sum);
for (i=0; i<sum; i++) {
  sprintf(str2, "Row %d: ",i);
  strcat(str1, str2);
  for (iPtr =A_C->pattern->ptr[i]; iPtr<A_C->pattern->ptr[i+1]; ++iPtr) {
         sprintf(str2, "(%d %f),",A_C->pattern->index[iPtr], A_C->val[iPtr]);
         strcat(str1, str2);
  }
  sprintf(str1, "%s\n", str1);
}
fprintf(stderr, "%s\n", str1);
Paso_SparseMatrix_free(A_C);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

//if (A_p->mpi_info->rank == 0) {
if (MY_DEBUG) {
char *str1, *str2;
int sum, rank, i, iPtr;
sum = A_C->mainBlock->numRows;
rank = A_p->mpi_info->rank;
str1 = TMPMEMALLOC(sum*sum*20+100, char);
str2 = TMPMEMALLOC(20, char);
sprintf(str1, "rank %d A_C Main Block: %d rows\n-----------\n", rank, sum);
/*for (i=0; i<sum; i++) {
  sprintf(str2, "Row %d: ",i);
  strcat(str1, str2);
  for (iPtr =A_C->mainBlock->pattern->ptr[i]; iPtr<A_C->mainBlock->pattern->ptr[i+1]; ++iPtr) {
         sprintf(str2, "(%d %f),",A_C->mainBlock->pattern->index[iPtr], A_C->mainBlock->val[iPtr]);
         strcat(str1, str2);
  }
  sprintf(str1, "%s\n", str1);
}*/
fprintf(stderr, "%s\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}
//}
if (MY_DEBUG1) fprintf(stderr, "after buildInterpolation\n");
//Esys_setError(SYSTEM_ERROR, "AMG:DONE.");
//return NULL;
/*			   Atemp=Paso_SystemMatrix_MatrixMatrix(A_p,out->P);
			   A_C=Paso_SystemMatrix_MatrixMatrix(out->R,Atemp);
			   Paso_Preconditioner_AMG_setStrongConnections
			   Paso_SystemMatrix_free(Atemp);
*/

			   if (SHOW_TIMING) printf("timing: level %d : construct coarse matrix: %e\n",level,Esys_timer()-time0);			
			}

			/*
			   constructe courser level:
			   
			*/
			if ( Esys_noError()) {
if (MY_DEBUG1) fprintf(stderr, "before Step into a level deeper\n");
			   out->AMG_C=Paso_Preconditioner_AMG_alloc(A_C,level+1,options);
if (MY_DEBUG1) fprintf(stderr, "after step into a level deeper\n");
			}

			if ( Esys_noError()) {
if (MY_DEBUG) fprintf(stderr, "Now checking whether AMG_C is set:");
			  if ( out->AMG_C == NULL ) { 
if (MY_DEBUG1) fprintf(stderr, " NO, need direct solver!\n");
//			    if (total_n <= options->min_coarse_matrix_size) {
			      /* merge the system matrix into 1 rank when 
				 it's not suitable coarsening due to the 
				 total number of unknowns are too small */
			      out->A_C=A_C;
			      out->reordering = options->reordering;
                              out->refinements = options->coarse_matrix_refinements;
			      out->verbose = verbose;
			      out->options_smoother = options->smoother;
/*			    } else {
fprintf(stderr, " NO, need direct solver!\n");
			      out->reordering = options->reordering;
			      out->refinements = options->coarse_matrix_refinements;
			      // no coarse level matrix has been constructed. use direct solver 
			      #ifdef MKL
				    out->A_C=Paso_SystemMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, A_C);
				    Paso_SystemMatrix_free(A_C);
				    out->A_C->solver_package = PASO_MKL;
				    if (verbose) printf("Paso_Preconditioner: AMG: use MKL direct solver on the coarsest level (number of unknowns = %d).\n",n_C*n_block); 
			      #else
				    #ifdef UMFPACK
				       out->A_C=Paso_SystemMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_CSC, A_C); 
				       Paso_SystemMatrix_free(A_C);
				       out->A_C->solver_package = PASO_UMFPACK;
				       if (verbose) printf("Paso_Preconditioner: AMG: use UMFPACK direct solver on the coarsest level (number of unknowns = %d).\n",n_C*n_block); 
				    #else
				       out->A_C=A_C;
				       out->A_C->solver_p=Paso_Preconditioner_Smoother_alloc(out->A_C, (options->smoother == PASO_JACOBI), verbose);
				       out->A_C->solver_package = PASO_SMOOTHER;
				       if (verbose) printf("Paso_Preconditioner: AMG: use smoother on the coarsest level (number of unknowns = %d).\n",n_C*n_block);
				    #endif
			      #endif
			    }
*/
			  } else {
			      /* finally we set some helpers for the solver step */
if (MY_DEBUG1) fprintf(stderr, "YES\n");
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

if (MY_DEBUG1) fprintf(stderr, "rank %d Now in AMG_solve level %d %d %d\n", A->mpi_info->rank, amg->level, A->is_balanced, Esys_getErrorType());
if (MY_DEBUG) {
int rank=A->mpi_info->rank;
if (rank == 3) fprintf(stderr, "SEND %d TO %d\n", amg->R->col_coupler->connector->send->offsetInShared[2] - amg->R->col_coupler->connector->send->offsetInShared[1],
amg->R->col_coupler->connector->send->neighbor[1]);
if (rank == 4) fprintf(stderr, "RECV %d FM %d\n", amg->R->col_coupler->connector->recv->offsetInShared[1] - amg->R->col_coupler->connector->recv->offsetInShared[0],
amg->R->col_coupler->connector->recv->neighbor[0]);
}

     /* presmoothing */
     time0=Esys_timer();
     Paso_Preconditioner_Smoother_solve(A, amg->Smoother, x, b, pre_sweeps, FALSE); 
//     Paso_Preconditioner_LocalSmoother_solve(A->mainBlock, amg->Smoother->localSmoother, x, b, pre_sweeps, FALSE);
     time0=Esys_timer()-time0;
if (MY_DEBUG1) fprintf(stderr, "rank %d after smoother_solve F%d n%d %d %d\n", A->mpi_info->rank, amg->n_F, amg->n, A->is_balanced, Esys_getErrorType());
     if (SHOW_TIMING) printf("timing: level %d: Presmooting: %e\n",amg->level, time0); 
     /* end of presmoothing */
	
     if (amg->n_F < amg->n) { /* is there work on the coarse level? */
         time0=Esys_timer();
if (MY_DEBUG){
char *str1, *str2;
int sum, rank, i;
str1 = TMPMEMALLOC(2000+100, char);
str2 = TMPMEMALLOC(100, char);
sum = n;
rank = A->mpi_info->rank;
sprintf(str1, "rank %d Level %d (after Smoother) x[%d] = (", rank, amg->level, sum);
for (i=0; i<sum; i++) {
  sprintf(str2, "%f ", x[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

	 Paso_Copy(n, amg->r, b);                            /*  r <- b */
if (MY_DEBUG1) fprintf(stderr, "rank %d after r=b \n", A->mpi_info->rank);
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(-1.,A,x,1.,amg->r); /*r=r-Ax*/

if (MY_DEBUG1) fprintf(stderr, "rank %d after r=r-Ax %d %d\n", A->mpi_info->rank, A->is_balanced, Esys_getErrorType());
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(1.,amg->R,amg->r,0.,amg->b_C);  /* b_c = R*r  */

if (MY_DEBUG) {
  char * str1, *str2;
  int i, q;
  int sum = Paso_SystemMatrix_getTotalNumRows(A);
  int sum1 = Paso_SystemMatrix_getTotalNumRows(amg->A_C);
  str1 = TMPMEMALLOC(sum*sum*30+100, char);
  str2 = TMPMEMALLOC(30, char);
  sprintf(str1, "rank %d Level %d x[%d]=(", A->mpi_info->rank, amg->level, sum);
  for (i=0; i<sum; i++) {
    sprintf(str2, "%f ", x[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  if (amg->level > 1) {
    sprintf(str1, "rank %d A Main Block:\n-----------\n", A->mpi_info->rank);
    for (q=0; q< sum; ++q){
      sprintf(str2, "Row %d: ",q);
      strcat(str1, str2);
      for (i =A->mainBlock->pattern->ptr[q]; i<A->mainBlock->pattern->ptr[q+1]; ++i) {
         sprintf(str2, "(%d %f),",A->mainBlock->pattern->index[i], A->mainBlock->val[i]);
         strcat(str1, str2);
      }
      sprintf(str1, "%s\n", str1);
    }
    fprintf(stderr, "%s)\n", str1);
  }
  sprintf(str1, "rank %d level %d r[%d]=(", A->mpi_info->rank, amg->level, sum);
  for (i=0; i<sum; i++) {
    sprintf(str2, "%g ", amg->r[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  sprintf(str1, "rank %d level %d b_C[%d]=(", A->mpi_info->rank, amg->level, sum1);
  for (i=0; i<sum1; i++) {
    sprintf(str2, "%g ", amg->b_C[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  TMPMEMFREE(str1);
  TMPMEMFREE(str2);
}

         time0=Esys_timer()-time0;
if (MY_DEBUG1) fprintf(stderr, "rank %d after matrix_vector %d %d\n", A->mpi_info->rank, A->is_balanced, Esys_getErrorType());
	 /* coarse level solve */
	 if ( amg->AMG_C == NULL) {
	    time0=Esys_timer();
	    /*  A_C is the coarsest level */
	    Paso_Preconditioner_AMG_mergeSolve(amg); 
if (MY_DEBUG1) fprintf(stderr, "rank %d after merge_solve %d %d\n", A->mpi_info->rank, A->is_balanced, Esys_getErrorType());
	    if (SHOW_TIMING) printf("timing: level %d: DIRECT SOLVER: %e\n",amg->level,Esys_timer()-time0);
	 } else {
	    Paso_Preconditioner_AMG_solve(amg->A_C, amg->AMG_C,amg->x_C,amg->b_C); /* x_C=AMG(b_C)     */
	 }
if (MY_DEBUG) {
  char * str1, *str2;
  int i, sum = Paso_SystemMatrix_getTotalNumRows(amg->A_C);
  if (sum > 20) sum = 20;
  str1 = TMPMEMALLOC(sum*30+100, char);
  str2 = TMPMEMALLOC(30, char);
  sprintf(str1, "rank %d level %d b[%d]=(", A->mpi_info->rank, amg->level, sum);
  for (i=0; i<sum; i++) {
    sprintf(str2, "%g ", amg->b_C[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  sprintf(str1, "rank %d level %d x[%d]=(", A->mpi_info->rank, amg->level, sum);
  for (i=0; i<sum; i++) {
    sprintf(str2, "%g ", amg->x_C[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  TMPMEMFREE(str1);
  TMPMEMFREE(str2);
}

if (MY_DEBUG) fprintf(stderr, "rank %d after AMG_solve %d %d\n", A->mpi_info->rank, A->is_balanced, Esys_getErrorType());  
  	 time0=time0+Esys_timer();
	 Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(1.,amg->P,amg->x_C,1.,x); /* x = x + P*x_c */    
if (MY_DEBUG) fprintf(stderr, "rank %d after MatrixVector %d %d %s\n", A->mpi_info->rank, A->is_balanced, Esys_getErrorType(), Esys_getErrorMessage());
	
         /*postsmoothing*/
      
        /*solve Ax=b with initial guess x */
        time0=Esys_timer();
        Paso_Preconditioner_Smoother_solve(A, amg->Smoother, x, b, post_sweeps, TRUE); 
//	Paso_Preconditioner_LocalSmoother_solve(A->mainBlock, amg->Smoother->localSmoother, x, b, post_sweeps, TRUE);
        time0=Esys_timer()-time0;
        if (SHOW_TIMING) printf("timing: level %d: Postsmoothing: %e\n",amg->level,time0);
        /*end of postsmoothing*/
     
     }
if (MY_DEBUG1) fprintf(stderr, "rank %d Now out of AMG_solve level %d\n", A->mpi_info->rank, amg->level);

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
   index_t rank=A->mpi_info->rank;

   threshold_p = TMPMEMALLOC(2*my_n, double);
   
//   #pragma omp parallel for private(i,iptr) schedule(static)
   for (i=0;i<my_n;++i) {        
	 
	 register double max_offdiagonal = 0.;
	 register double sum_row=0;
	 register double main_row=0;
	 register dim_t kdeg=0;
         register const index_t koffset=A->mainBlock->pattern->ptr[i]+A->col_coupleBlock->pattern->ptr[i];

         
	 /* collect information for row i: */
//	 #pragma ivdep
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

//	 #pragma ivdep
	 for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	    register double fnorm=ABS(A->col_coupleBlock->val[iptr]);

	    max_offdiagonal = MAX(max_offdiagonal,fnorm);
	    sum_row+=fnorm;
	 }

         /* inspect row i: */
         {
	    const double threshold = theta*max_offdiagonal;
            threshold_p[2*i+1]=threshold;
	    if (tau*main_row < sum_row) { /* no diagonal domainance */
               threshold_p[2*i]=1;
//	       #pragma ivdep
	       for (iptr=A->mainBlock->pattern->ptr[i];iptr<A->mainBlock->pattern->ptr[i+1]; ++iptr) {
	          register index_t j=A->mainBlock->pattern->index[iptr];
if (MY_DEBUG && rank == 0 && i == 16) 
fprintf(stderr, "Row 16--debug: j %d, val %f, threshold %f\n", j, ABS(A->mainBlock->val[iptr]), threshold);
	          if(ABS(A->mainBlock->val[iptr])>threshold && i!=j) {
		     S[koffset+kdeg] = j;
		     kdeg++;
	          }
	       }
//	       #pragma ivdep
	       for (iptr=A->col_coupleBlock->pattern->ptr[i];iptr<A->col_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	          register index_t j=A->col_coupleBlock->pattern->index[iptr];
if (MY_DEBUG && rank == 0 && i == 16)
fprintf(stderr, "Row 16--debug: j %d, Cval %f, threshold %f\n", j, ABS(A->col_coupleBlock->val[iptr]), threshold);
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

//          #pragma omp parallel for private(i,iptr) schedule(static)
          for (i=0; i<overlap_n; i++) {
	      const double threshold = remote_threshold[2*i+1];
	      register dim_t kdeg=0;
              register const index_t koffset=koffset_0+A->row_coupleBlock->pattern->ptr[i]+A->remote_coupleBlock->pattern->ptr[i];
              if (remote_threshold[2*i]>0) {
//	         #pragma ivdep
		for (iptr=A->row_coupleBlock->pattern->ptr[i];iptr<A->row_coupleBlock->pattern->ptr[i+1]; ++iptr) {
	          register index_t j=A->row_coupleBlock->pattern->index[iptr];
if (MY_DEBUG && rank == 0 && i == 1)
fprintf(stderr, "Row 9--debug: j %d, Rowval %f, threshold %f\n", j, ABS(A->row_coupleBlock->val[iptr]), threshold);
	          if(ABS(A->row_coupleBlock->val[iptr])>threshold) {
		     S[koffset+kdeg] = j ;
		     kdeg++;
	          }
		}

//		 #pragma ivdep
		for (iptr=A->remote_coupleBlock->pattern->ptr[i];iptr<A->remote_coupleBlock->pattern->ptr[i+1]; iptr++) {
		  register index_t j=A->remote_coupleBlock->pattern->index[iptr];
if (MY_DEBUG && rank == 0 && i == 1)
fprintf(stderr, "Row 9--debug: j %d, Remote-val %f, threshold %f\n", j, ABS(A->remote_coupleBlock->val[iptr]), threshold);
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
   
   #pragma omp parallel private(i,iptr, bi)
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
	    if (tau*main_row < sum_row) { /* no diagonal domainance */
	       threshold_p[2*i]=1;
	       rtmp_offset=-A->mainBlock->pattern->ptr[i];
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
	 register const index_t koffset=koffset_0+A->row_coupleBlock->pattern->ptr[i];
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

               if(fnorm2 > threshold2 ) {
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
  index_t rank=col_connector->mpi_info->rank;

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
  printf(" coarsening loop start: num of undefined rows = %d \n",numUndefined); 

  
  iter=0; 
  while (numUndefined > 0) {
     Paso_Coupler_fillOverlap(n, w, w_coupler);
/*if (MY_DEBUG) {
	int p;
	for (p=0; p<my_n; ++p) {
	   fprintf(stderr, "rank%d %d : %f %f \n",rank, p, w[p], Status[p]);
	}
	for (p=my_n; p<n; ++p) {
	   fprintf(stderr, "rank%d %d : %f \n",rank, p, w[p]);
	}
}*/
     
      /* calculate the maximum value of naigbours following active strong connections:
	    w2[i]=MAX(w[k]) with k in ST[i] or S[i] and (i,k) conenction is still active  */       
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
//if (MY_DEBUG) fprintf(stderr, "rank%d S: %d (%e) -> %d	(%e)\n",rank, i, wi, k, w[k]);	  
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
//if (MY_DEBUG) fprintf(stderr, "rank%d ST: %d (%e) -> %d	(%e)\n",rank, i, wi, k, w[k]);
			
                       if (wi <= w[k] ) {
			   inD=FALSE;
			   break;
			}
		     }
		  }
	    }    
	    if (inD) { 
	       Status[i]=0.; /* is in D */
//if (MY_DEBUG) fprintf(stderr, "rank%d %d is in D\n",rank, i);
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
//if (MY_DEBUG) fprintf(stderr, "rank%d %d reduced by %d\n",rank, i,j);
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
//if (MY_DEBUG) fprintf(stderr, "rank%d check connection: %d %d\n",rank, i,j);
			ST_flag[offset_ST[i]+jptr]=-1;
			for (kptr=0; kptr<degree_ST[j]; ++kptr) {
			   const index_t k=ST[offset_ST[j]+kptr]; 
//if (MY_DEBUG) fprintf(stderr, "rank%d check connection: %d of %d\n",rank, k,j); 
			   if (NULL != bsearch(&k, start_p, degree_ST[i], sizeof(index_t), Paso_comparIndex) ) { /* k in ST[i] ? */
//if (MY_DEBUG) fprintf(stderr, "rank%d found!\n", rank);
			      if (ST_flag[offset_ST[j]+kptr] >0) {
				 if (j< my_n ) {
				    w[j]--;
//if (MY_DEBUG) fprintf(stderr, "rank%d %d reduced by %d and %d \n",rank, j, i,k); 
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
	 printf(" coarsening loop %d: num of undefined rows = %d \n",iter, numUndefined);
  } /* end of while loop */

  /* map to output :*/
  Paso_Coupler_fillOverlap(n, Status, w_coupler);
//  #pragma omp parallel for private(i)
  for (i=0; i< n; ++i) {
	 if (Status[i] > -50.) {
	    split_marker[i]=PASO_AMG_IN_C;
if (MY_DEBUG) fprintf(stderr, "rank%d CIJP C sets: %d\n",rank, i);
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

if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSystemMatrix CP1 %d\n", rank, Esys_getErrorType());
  n = A->mainBlock->numRows;
  block_size = A->block_size;

  /* Merge MainBlock and CoupleBlock to get a complete column entries
     for each row allocated to current rank. Output (ptr, idx, val) 
     contains all info needed from current rank to merge a system 
     matrix  */
  Paso_SystemMatrix_mergeMainAndCouple(A, &ptr, &idx, &val);
if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSystemMatrix CP2 %d %s\n", rank, Esys_getErrorType(), Esys_getErrorMessage());
if (MY_DEBUG && rank == 0) {
char *str1, *str2;
index_t i, j, ib, sum=ptr[n], block=A->block_size;
str1 = TMPMEMALLOC(sum*block*30+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d val[%d X %d] = (", rank, n, block);
for (i=0; i<n; i++){
 sprintf(str2, "Row %d: ", i);
 strcat(str1, str2);
 for (j=ptr[i]; j<ptr[i+1]; j++){
  sprintf(str2, "(%d, ", idx[j]);
  strcat(str1, str2);
  for (ib=0; ib<block; ib++) {
    sprintf(str2, "%f ", val[j*block+ib]);
    strcat(str1, str2);
  }
  sprintf(str1, "%s) ", str1);
 }
 sprintf(str1, "%s\n", str1);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}


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
if (MY_DEBUG) {
char *str1, *str2;
index_t sum=total_n+1;
str1 = TMPMEMALLOC(sum*100+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d ptr_global[%d of %d] = (", rank, n+1, sum);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", ptr_global[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}
    iptr = n+1;
    MEMFREE(ptr);
    temp_n = TMPMEMALLOC(size, index_t);
    temp_len = TMPMEMALLOC(size, index_t);
    temp_n[0] = iptr;
if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSystemMatrix CP3 %d\n", rank, Esys_getErrorType());
    
    /* Second, receive ptr values from other ranks */
    for (i=1; i<size; i++) {
      remote_n = A->row_distribution->first_component[i+1] -
		 A->row_distribution->first_component[i];
if (MY_DEBUG) fprintf(stderr, "rank %d n %d i %d remote_n %d iptr %d\n", rank, n, i, remote_n, iptr);
      MPI_Irecv(&(ptr_global[iptr]), remote_n, MPI_INT, i, 
			A->mpi_info->msg_tag_counter+i,
			A->mpi_info->comm,
			&mpi_requests[i]);
      temp_n[i] = remote_n;
      iptr += remote_n;
    }
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    A->mpi_info->msg_tag_counter += size;
if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSystemMatrix CP4 %d\n", rank, Esys_getErrorType());
if (MY_DEBUG) {
char *str1, *str2;
index_t sum=iptr;
str1 = TMPMEMALLOC(sum*100+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d ptr_global[%d] = (", rank, sum);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", ptr_global[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}


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
if (MY_DEBUG) {
char *str1, *str2;
index_t sum = size;
str1 = TMPMEMALLOC(sum*100+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d temp[%d, len%d] = (", rank, sum, len);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", temp_len[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

    idx_global = MEMALLOC(len, index_t);
    iptr = temp_len[0];
    offset = n+1;
    for (i=1; i<size; i++) {
      len = temp_len[i];
      MPI_Irecv(&(idx_global[iptr]), len, MPI_INT, i,
			A->mpi_info->msg_tag_counter+i,
			A->mpi_info->comm,
			&mpi_requests[i]);
      remote_n = temp_n[i];
if (MY_DEBUG) fprintf(stderr,"RECV %d from %d offset %d tag %d\n", len, i, iptr, A->mpi_info->msg_tag_counter+i);
      for (j=0; j<remote_n; j++) {
	ptr_global[j+offset] = ptr_global[j+offset] + iptr;
if (MY_DEBUG) fprintf(stderr, "rank %d j %d iptr %d offset %d\n", rank, j, iptr, offset);
      }
      offset += remote_n;
      iptr += len;
    }
if (MY_DEBUG) fprintf(stderr, "rank %d copy len %d\n", rank, temp_len[0]);
    memcpy(idx_global, idx, temp_len[0] * sizeof(index_t));
    MEMFREE(idx);
    row_block_size = A->mainBlock->row_block_size;
    col_block_size = A->mainBlock->col_block_size;
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    A->mpi_info->msg_tag_counter += size;
    TMPMEMFREE(temp_n);
if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSystemMatrix CP5 %d\n", rank, Esys_getErrorType());
if (MY_DEBUG) {
char *str1, *str2;
index_t sum= total_n+1;
str1 = TMPMEMALLOC(sum*100+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d ptr_global[%d] = (", rank, sum);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", ptr_global[i+1]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}
if (MY_DEBUG) {
char *str1, *str2;
index_t sum =  ptr_global[total_n];
str1 = TMPMEMALLOC(sum*100+100, char);
str2 = TMPMEMALLOC(30, char);
sprintf(str1, "rank %d idx_global[%d] = (", rank, sum);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", idx_global[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s) -- total_n %d len %d\n", str1, total_n, ptr_global[total_n]);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

    /* Then generate the sparse matrix */
    pattern = Paso_Pattern_alloc(A->mainBlock->pattern->type, total_n,
			total_n, ptr_global, idx_global);
    out = Paso_SparseMatrix_alloc(A->mainBlock->type, pattern, 
			row_block_size, col_block_size, FALSE);
    Paso_Pattern_free(pattern);

if (MY_DEBUG) {
char *str1, *str2;
index_t sum;
str1 = TMPMEMALLOC(size*1000+100, char);
str2 = TMPMEMALLOC(15, char);
sum = out->pattern->ptr[pattern->numOutput];
sprintf(str1, "rank %d out_idx[%d] = (", rank, sum);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", out->pattern->index[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}

    /* Finally, receive and copy the value */
    iptr = temp_len[0] * block_size;
    for (i=1; i<size; i++) {
      len = temp_len[i];
      MPI_Irecv(&(out->val[iptr]), len * block_size, MPI_DOUBLE, i,
                        A->mpi_info->msg_tag_counter+i,
                        A->mpi_info->comm,
                        &mpi_requests[i]);
      iptr += (len * block_size);
    }
    memcpy(out->val, val, temp_len[0] * sizeof(double) * block_size);
    MEMFREE(val);
    MPI_Waitall(size-1, &(mpi_requests[1]), mpi_stati);
    A->mpi_info->msg_tag_counter += size;
    TMPMEMFREE(temp_len);
if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSystemMatrix CP6 %d\n", rank, Esys_getErrorType());
  } else { /* it's not rank 0 */

if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSystemMatrix CP3 n %d\n", rank, n);
    /* First, send out the local ptr */
    tag = A->mpi_info->msg_tag_counter+rank;
    MPI_Issend(&(ptr[1]), n, MPI_INT, 0, tag, A->mpi_info->comm, 
			&mpi_requests[0]);

if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSystemMatrix CP4\n", rank);
if (MY_DEBUG1 && rank == 1) 
fprintf(stderr, "SIZE n %d distribution %d, ptr[n] %d\n", n, A->row_distribution->first_component[rank+1] - A->row_distribution->first_component[rank], ptr[n]);

    /* Next, send out the local idx */
    len = ptr[n];
    tag += size;
    MPI_Issend(idx, len, MPI_INT, 0, tag, A->mpi_info->comm, 
			&mpi_requests[1]);
if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSystemMatrix CP5\n", rank);
if (MY_DEBUG1) fprintf(stderr, "rank %d SEND len%d tag %d \n", rank, len, tag);

    /* At last, send out the local val */
    len *= block_size;
    tag += size;
if (MY_DEBUG) fprintf(stderr, "rank %d send len%d tag %d \n", rank, len, tag);
    MPI_Issend(val, len, MPI_DOUBLE, 0, tag, A->mpi_info->comm, 
			&mpi_requests[2]);

    MPI_Waitall(3, mpi_requests, mpi_stati);
    A->mpi_info->msg_tag_counter = tag + size - rank;
    MEMFREE(ptr);
    MEMFREE(idx);
    MEMFREE(val);
if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSystemMatrix CP6\n", rank);

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
  index_t i, n, p, count, n_block;
  index_t *counts, *offset, *dist;

if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSolve CP1\n", rank);
  n_block = amg->n_block;

  A_D = Paso_Preconditioner_AMG_mergeSystemMatrix(A); 
if (MY_DEBUG && rank == 0) Paso_SystemMatrix_print(A);
if (MY_DEBUG && rank == 0) {
   int q, iPtr, ib, block_size=A_D->block_size, n=A_D->numRows;
   char *str1, *str2;
   str1 = TMPMEMALLOC(n*n*block_size*30+100, char);
   str2 = TMPMEMALLOC(30, char);

   sprintf(str1, "rank %d merged Matrix A_D (N%d X N%d Block%d):\n %d from Local\n", rank, n, n, block_size, A->mainBlock->numRows);
   for (q=0; q< n; ++q){
      sprintf(str2, "Row %d: ",q);
      strcat(str1, str2);
      for (iPtr =A_D->pattern->ptr[q]; iPtr<A_D->pattern->ptr[q+1]; ++iPtr) {
         sprintf(str2, "(%d ",A_D->pattern->index[iPtr]);
         strcat(str1, str2);
         for (ib=0; ib<block_size; ib++){
           sprintf(str2, "%f ", A_D->val[iPtr*block_size+ib]);
           strcat(str1, str2);
         }
         sprintf(str2, "),");
         strcat(str1, str2);
      }
      sprintf(str1, "%s\n", str1);
   }
   fprintf(stderr, "%s", str1);

}

if (MY_DEBUG) {
if (rank == 0) {
char *str1, *str2;
index_t sum;
str1 = TMPMEMALLOC(size*1000+100, char);
str2 = TMPMEMALLOC(15, char);
sum = A_D->pattern->ptr[A_D->pattern->numOutput];
sprintf(str1, "rank %d idx[%d] = (", rank, sum);
for (i=0; i<sum; i++){
  sprintf(str2, "%d ", A_D->pattern->index[i]);
  strcat(str1, str2);
}
fprintf(stderr, "%s)\n", str1);
TMPMEMFREE(str1);
TMPMEMFREE(str2);
}
}

if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSolve CP2 %d\n", rank, Esys_getErrorType());
  /* First, gather x and b into rank 0 */
  dist = A->pattern->input_distribution->first_component;
  n = Paso_SystemMatrix_getGlobalNumRows(A);
  b = TMPMEMALLOC(n*n_block, double);
  x = TMPMEMALLOC(n*n_block, double);
  counts = TMPMEMALLOC(size, index_t);
  offset = TMPMEMALLOC(size, index_t);

  for (i=0; i<size; i++) {
    p = dist[i];
    counts[i] = (dist[i+1] - p)*n_block;
    offset[i] = p*n_block;
  }
  count = counts[rank];
  MPI_Gatherv(amg->b_C, count, MPI_DOUBLE, b, counts, offset, MPI_DOUBLE, 0, A->mpi_info->comm);
  MPI_Gatherv(amg->x_C, count, MPI_DOUBLE, x, counts, offset, MPI_DOUBLE, 0, A->mpi_info->comm);

if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSolve CP3 %d\n", rank, Esys_getErrorType());

if (MY_DEBUG1) {
if (rank == 0) {
  char * str1, *str2;
  int sum = n*n_block;
  if (sum > 20) sum = 20;
  str1 = TMPMEMALLOC(sum*30+100, char);
  str2 = TMPMEMALLOC(30, char);
  sprintf(str1, "b[%d of %d] from all ranks=(", sum, n*n_block);
  for (i=0; i<sum; i++) {
    sprintf(str2, "%g ", b[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  TMPMEMFREE(str1);
  TMPMEMFREE(str2);
}
}

  if (rank == 0) {
    /* solve locally */
    #ifdef MKL
if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSolve CP3_0 MKL\n", rank);
      A_temp = Paso_SparseMatrix_unroll(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, A_D);
      A_temp->solver_package = PASO_MKL;
      Paso_SparseMatrix_free(A_D);
      Paso_MKL(A_temp, x, b, amg->reordering, amg->refinements, SHOW_TIMING);
if (MY_DEBUG) fprintf(stderr, "rank %d In mergeSolve CP3_4\n", rank);
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

if (MY_DEBUG1) {
if (rank == 0) {
  char * str1, *str2;
  int sum = n*n_block;
  if (sum > 20) sum = 20;
  str1 = TMPMEMALLOC(sum*30+100, char);
  str2 = TMPMEMALLOC(30, char);
  sprintf(str1, "Solve x[%d of %d] = ( ", sum, n*n_block);
  for (i=0; i<sum; i++) {
    sprintf(str2, "%g ", x[i]);
    strcat(str1, str2);
  }
  fprintf(stderr, "%s)\n", str1);
  TMPMEMFREE(str1);
  TMPMEMFREE(str2);
}
}

if (MY_DEBUG1) fprintf(stderr, "rank %d In mergeSolve CP4 %d\n", rank, Esys_getErrorType());
  /* now we need to distribute the solution to all ranks */
  MPI_Scatterv(x, counts, offset, MPI_DOUBLE, amg->x_C, count, MPI_DOUBLE, 0, A->mpi_info->comm);

  TMPMEMFREE(x);
  TMPMEMFREE(b);
  TMPMEMFREE(counts);
  TMPMEMFREE(offset);
}
