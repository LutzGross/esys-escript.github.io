/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

/* Paso: defines AMG prolongation  */

/************************************************************************************/

/* Author: Artak Amirbekyan, artak@uq.edu.au, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

/************************************************************************************

    Methods necessary for AMG preconditioner

    construct n x n_C the prolongation matrix P from A_p.
    
    The columns in A_p to be considered are marked by counter_C[n] where
    an unknown i to be considered in P is marked by 0<= counter_C[i] < n_C 
    and counter_C[i]  gives the new column number in P.
    S defines the strong connections.
    
    The pattern of P is formed as follows:

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise.
    If row i is not C, then P[i,j] <> 0 if counter_C[k]==j (k in C) and (i,k) is a strong connection.  
    
    Two settings for P are implemented (see below) 
   
*/

Paso_SystemMatrix* Paso_Preconditioner_AMG_getProlongation(Paso_SystemMatrix* A_p, 
                                                           const index_t* offset_S, const dim_t* degree_S, const index_t* S,
							   const dim_t n_C, index_t* counter_C, const index_t interpolation_method) 
{
   Esys_MPIInfo *mpi_info=Esys_MPIInfo_getReference(A_p->mpi_info);
   Paso_SystemMatrix *out=NULL;
   Paso_SystemMatrixPattern *pattern=NULL;
   Paso_Distribution *input_dist=NULL, *output_dist=NULL;
   Paso_SharedComponents *send =NULL, *recv=NULL;
   Paso_Connector *col_connector=NULL;
   Paso_Pattern *main_pattern=NULL, *couple_pattern=NULL;
   const dim_t row_block_size=A_p->row_block_size;
   const dim_t col_block_size=A_p->col_block_size;
   const dim_t my_n=A_p->mainBlock->numCols;
   const dim_t overlap_n=A_p->col_coupleBlock->numCols;
   const dim_t num_threads=omp_get_max_threads();
   index_t size=mpi_info->size, *dist=NULL;
   index_t *main_p=NULL, *couple_p=NULL, *main_idx=NULL, *couple_idx=NULL;
   index_t *shared=NULL, *offsetInShared=NULL;
   index_t *recv_shared=NULL, *send_shared=NULL;
   index_t sum, i, j, k, l, p, q, iptr;
   index_t my_n_C, global_label, num_neighbors;
   #ifdef ESYS_MPI
   index_t rank=mpi_info->rank;
   #endif
   Esys_MPI_rank *neighbor=NULL;
   #ifdef ESYS_MPI
     MPI_Request* mpi_requests=NULL;
     MPI_Status* mpi_stati=NULL;
   #else
     int *mpi_requests=NULL, *mpi_stati=NULL;
   #endif

   /* number of C points in current distribution */
   my_n_C = 0;
   sum=0;
   if (num_threads>1) {
     #pragma omp parallel private(i,sum)
     {
	sum=0;
	#pragma omp for schedule(static)
	for (i=0;i<my_n;++i) {
	  if (counter_C[i] != -1) {
	    sum++;
	  }
	}
	#pragma omp critical
	{
	    my_n_C += sum;
	}
     }
   } else { /* num_threads=1 */
     for (i=0;i<my_n;++i) {
         if (counter_C[i] != -1) {
            my_n_C++;
         }
      }
   }

   /* create row distribution (output_distribution) and col distribution 
      (input_distribution) */
   /* ??? should I alloc an new Esys_MPIInfo object or reuse the one in
      system matrix A. for now, I'm reuse A->mpi_info ??? */
   dist = A_p->pattern->output_distribution->first_component;
   output_dist=Paso_Distribution_alloc(mpi_info, dist, 1, 0);
   dist = TMPMEMALLOC(size+1, index_t); /* now prepare for col distribution */
   Esys_checkPtr(dist);
   #ifdef ESYS_MPI
   MPI_Allgather(&my_n_C, 1, MPI_INT, dist, 1, MPI_INT, mpi_info->comm);
   #endif
   global_label=0;
   for (i=0; i<size; i++) {
     k = dist[i];
     dist[i] = global_label;
     global_label += k;
   }
   dist[size] = global_label;

   input_dist=Paso_Distribution_alloc(mpi_info, dist, 1, 0);
   TMPMEMFREE(dist);

   /* create pattern for mainBlock and coupleBlock */
   main_p = MEMALLOC(my_n+1, index_t);
   couple_p = MEMALLOC(my_n+1, index_t);
   if (!(Esys_checkPtr(main_p) || Esys_checkPtr(couple_p))) {
     /* count the number of entries per row in the Prolongation matrix :*/
     #pragma omp parallel for private(i,l,k,iptr,j,p) schedule(static)
     for (i=0; i<my_n; i++) {
	l = 0;
	if (counter_C[i]>=0) {
	  k = 1;    /* i is a C unknown */
	} else {
	  k = 0;
	  iptr = offset_S[i];
	  for (p=0; p<degree_S[i]; p++) {
	    j = S[iptr+p];  /* this is a strong connection */
	    if (counter_C[j]>=0) { /* and is in C */
		if (j <my_n) k++;
		else {
		  l++; 
		}
	    }
	  }
	}
	main_p[i] = k;
	couple_p[i] = l;
     }

     /* number of unknowns in the col-coupleBlock of the interplation matrix */
     sum = 0;
     for (i=0;i<overlap_n;++i) {
	if (counter_C[i+my_n] > -1) {
	  counter_C[i+my_n] -= my_n_C;
	  sum++;
	}
     }

     /* allocate and create index vector for prolongation: */
     p = Paso_Util_cumsum(my_n, main_p);
     main_p[my_n] = p;
     main_idx = MEMALLOC(p, index_t);
     p = Paso_Util_cumsum(my_n, couple_p);
     couple_p[my_n] = p;
     couple_idx = MEMALLOC(p, index_t);
     if (!(Esys_checkPtr(main_idx) || Esys_checkPtr(couple_idx))) {
	#pragma omp parallel for private(i,k,l,iptr,j,p)  schedule(static)
	for (i=0; i<my_n; i++) {
	  if (counter_C[i]>=0) {
	    main_idx[main_p[i]]=counter_C[i];  /* i is a C unknown */
	  } else {
	    k = 0;
	    l = 0;
	    iptr = offset_S[i]; 
	    for (p=0; p<degree_S[i]; p++) {
	      j = S[iptr+p]; /* this is a strong connection */
	      if (counter_C[j] >=0) { /* and is in C */
		if (j < my_n) {
		  main_idx[main_p[i]+k] = counter_C[j];
		  k++; 
		} else {
		  couple_idx[couple_p[i]+l] = counter_C[j];
		  l++;
		}
	      }
	    }
	  }
	}
     }
   }

   if (Esys_noError()) {   
     main_pattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, my_n, 
			my_n_C, main_p, main_idx);
     couple_pattern = Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, my_n, 
			sum, couple_p, couple_idx);
   } else {
     MEMFREE(main_p);
     MEMFREE(main_idx);
     MEMFREE(couple_p);
     MEMFREE(couple_idx);
   }

   /* prepare the receiver for the col_connector. 
      Note that the allocation for "shared" assumes the send and receive buffer
      of the interpolation matrix P is no larger than that of matrix A_p. */
   neighbor = TMPMEMALLOC(size, Esys_MPI_rank);
   offsetInShared = TMPMEMALLOC(size+1, index_t);
   recv = A_p->col_coupler->connector->recv;
   send = A_p->col_coupler->connector->send;
   i = recv->numSharedComponents;
   recv_shared = TMPMEMALLOC(i,index_t);
   memset(recv_shared, 0, sizeof(index_t)*i);
   k = send->numSharedComponents;
   send_shared = TMPMEMALLOC(k,index_t);
   if (k > i) i = k;
   shared = TMPMEMALLOC(i, index_t);

   #ifdef ESYS_MPI
     mpi_requests=TMPMEMALLOC(size*2,MPI_Request);
     mpi_stati=TMPMEMALLOC(size*2,MPI_Status);
   #else
     mpi_requests=TMPMEMALLOC(size*2,int);
     mpi_stati=TMPMEMALLOC(size*2,int);
   #endif

   for (p=0; p<send->numNeighbors; p++) {
     i = send->offsetInShared[p];
     #ifdef ESYS_MPI
     MPI_Irecv (&(send_shared[i]), send->offsetInShared[p+1]-i, MPI_INT,
		send->neighbor[p], mpi_info->msg_tag_counter+send->neighbor[p],
		mpi_info->comm, &mpi_requests[p]);
     #endif
   }

   num_neighbors = 0;
   q = 0;
   p = recv->numNeighbors;
   offsetInShared[0]=0;
   for (i=0; i<p; i++) {
     l = 0;
     k = recv->offsetInShared[i+1];
     for (j=recv->offsetInShared[i]; j<k; j++) {
	if (counter_C[recv->shared[j]] > -1) {
	  shared[q] = my_n_C + q;
	  recv_shared[recv->shared[j]-my_n] = 1;
	  q++;
	  l = 1;
	}
     }
     if (l == 1) {
	iptr = recv->neighbor[i];
	neighbor[num_neighbors] = iptr;
	num_neighbors++;
	offsetInShared[num_neighbors] = q;
     }
     #ifdef ESYS_MPI
     MPI_Issend(&(recv_shared[recv->offsetInShared[i]]), 
		k-recv->offsetInShared[i], MPI_INT, recv->neighbor[i],
		mpi_info->msg_tag_counter+rank, mpi_info->comm,
		&mpi_requests[i+send->numNeighbors]);
     #endif
   }
   recv = Paso_SharedComponents_alloc(my_n_C, num_neighbors, neighbor, shared,
                                      offsetInShared, 1, 0, mpi_info);

   /* now we can build the sender */
   #ifdef ESYS_MPI
   MPI_Waitall(recv->numNeighbors+send->numNeighbors, mpi_requests, mpi_stati);
   #endif
   mpi_info->msg_tag_counter += size;
   TMPMEMFREE(mpi_requests);
   TMPMEMFREE(mpi_stati);

   num_neighbors = 0;
   q = 0;
   p = send->numNeighbors;
   offsetInShared[0]=0;
   for (i=0; i<p; i++) {
     l = 0;
     k = send->offsetInShared[i+1];
     for (j=send->offsetInShared[i]; j<k; j++) {
	if (send_shared[j] == 1) {
          shared[q] = counter_C[send->shared[j]];  
          q++;
          l = 1;
        }
     }
     if (l == 1) {
        iptr = send->neighbor[i];
        neighbor[num_neighbors] = iptr;
        num_neighbors++;
        offsetInShared[num_neighbors] = q;
     }
   }

   send = Paso_SharedComponents_alloc(my_n_C, num_neighbors, neighbor, shared,
				      offsetInShared, 1, 0, mpi_info);
   col_connector = Paso_Connector_alloc(send, recv);
   Paso_SharedComponents_free(recv);
   Paso_SharedComponents_free(send);
   TMPMEMFREE(recv_shared);
   TMPMEMFREE(send_shared);
   TMPMEMFREE(neighbor);
   TMPMEMFREE(offsetInShared);
   TMPMEMFREE(shared);

   /* now we need to create the System Matrix 
      TO BE FIXED: at this stage, we only construction col_couple_pattern
      and col_connector for interpolation matrix P. To be completed, 
      row_couple_pattern and row_connector need to be constructed as well */
   if (Esys_noError()) {
     pattern = Paso_SystemMatrixPattern_alloc(MATRIX_FORMAT_DEFAULT, 
		output_dist, input_dist, main_pattern, couple_pattern, 
		couple_pattern, col_connector, col_connector);
     out = Paso_SystemMatrix_alloc(MATRIX_FORMAT_DIAGONAL_BLOCK, pattern,
		row_block_size, col_block_size, FALSE);
   } 

   /* now fill in the matrix */
   if (Esys_noError()) {
     if ((interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING) 
        || ( interpolation_method == PASO_CLASSIC_INTERPOLATION) ) {
        if (row_block_size == 1) {
          Paso_Preconditioner_AMG_setClassicProlongation(out, A_p, offset_S, degree_S, S, counter_C);
        } else {
          Paso_Preconditioner_AMG_setClassicProlongation_Block(out, A_p, offset_S, degree_S, S, counter_C);
        }
     } else {
        if (row_block_size == 1) { 
          Paso_Preconditioner_AMG_setDirectProlongation(out, A_p, offset_S, degree_S, S, counter_C);
        } else {
          Paso_Preconditioner_AMG_setDirectProlongation_Block(out, A_p, offset_S, degree_S, S, counter_C);
        }
     }
   }  

   /* clean up */ 
   Paso_SystemMatrixPattern_free(pattern);
   Paso_Pattern_free(main_pattern);
   Paso_Pattern_free(couple_pattern);
   Paso_Connector_free(col_connector);
   Paso_Distribution_free(output_dist);
   Paso_Distribution_free(input_dist);
   if (Esys_noError()) {
      return out;
   } else {
      Paso_SystemMatrix_free(out);
      return NULL;
   }
   
}

/*
    Direct Prolongation:
    -------------------

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise.
    If row i is not C, then P[i,j] = - a[i] * A[i,k]/A[i,i] with j=counter_C[k]>=0 and k in S
   
   and    a[i]= 
             alpha[i] = sum_s min(A[i,s],0)/(sum_{s in S and C} min(A[i,s],0))   A[i,k]<0
                   or                                                         if
             beta[i] = sum_s max(A[i,s],0)/(sum_{s in S and C} max(A[i,s],0))   A[i,k]>0
              

*/

void Paso_Preconditioner_AMG_setDirectProlongation(Paso_SystemMatrix* P, 
	Paso_SystemMatrix* A, const index_t* offset_S, const dim_t* degree_S,
        const index_t* S, const index_t *counter_C) { 
   Paso_SparseMatrix *main_block=P->mainBlock;
   Paso_SparseMatrix *couple_block=P->col_coupleBlock;
   Paso_Pattern *main_pattern=main_block->pattern;
   Paso_Pattern *couple_pattern=couple_block->pattern;
   const dim_t my_n=A->mainBlock->numRows;
   index_t range;

   dim_t i;
   register double alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos, A_ij, A_ii, rtmp;
   register index_t iPtr, j, offset; 
   index_t *where_p, *start_p;
   
   #pragma omp parallel for private(i,offset,sum_all_neg,sum_all_pos,sum_strong_neg,sum_strong_pos,A_ii,range,iPtr,j,A_ij,start_p,where_p,alpha,beta,rtmp) schedule(static)
   for (i=0; i<my_n; i++) {
      if (counter_C[i]>=0) {
	    offset = main_pattern->ptr[i];
	    main_block->val[offset]=1.;  /* i is a C row */
      } else if ((main_pattern->ptr[i + 1] > main_pattern->ptr[i]) || 
		 (couple_pattern->ptr[i +1] > couple_pattern->ptr[i])) {
	 /* if i is an F row we first calculate alpha and beta: */

	 sum_all_neg=0; /* sum of all negative values in row i of A */
	 sum_all_pos=0; /* sum of all positive values in row i of A */
	 sum_strong_neg=0; /* sum of all negative values A_ij where j is in C and strongly connected to i*/
	 sum_strong_pos=0; /* sum of all positive values A_ij where j is in C and strongly connected to i*/
	 A_ii=0; 

	 /* first check the mainBlock */
	 range = A->mainBlock->pattern->ptr[i + 1];
	 for (iPtr=A->mainBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
	    j=A->mainBlock->pattern->index[iPtr];
	    A_ij=A->mainBlock->val[iPtr];
	    if(j==i) {
	       A_ii=A_ij;
	    } else {
	       
	       if(A_ij< 0)  {
		  sum_all_neg+=A_ij;
	       } else {
		  sum_all_pos+=A_ij;
	       }
	       
	       if (counter_C[j]>=0) {
		  /* is i strongly connected with j? We search for counter_C[j] in P[i,:] */ 
		  start_p=&(main_pattern->index[main_pattern->ptr[i]]);
		  where_p=(index_t*)bsearch(&(counter_C[j]), start_p,
					    main_pattern->ptr[i + 1] - main_pattern->ptr[i],
					    sizeof(index_t),
					    Paso_comparIndex);
		  if (! (where_p == NULL) ) { /* yes i strongly connected with j */
			offset = main_pattern->ptr[i]+ (index_t)(where_p-start_p);
			main_block->val[offset]=A_ij; /* will be modified later */
			if (A_ij< 0)  {
			   sum_strong_neg+=A_ij;
			} else {
			   sum_strong_pos+=A_ij;
			}
		  }
	       } 

	    } 
	 }

         /* now we deal with the col_coupleBlock */
         range = A->col_coupleBlock->pattern->ptr[i + 1];
         for (iPtr=A->col_coupleBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
            j=A->col_coupleBlock->pattern->index[iPtr];
            A_ij=A->col_coupleBlock->val[iPtr];
            if(A_ij< 0)  {
                  sum_all_neg+=A_ij;
            } else {
                  sum_all_pos+=A_ij;
            }

            if (counter_C[j+my_n]>=0) {
                  /* is i stronly connect with j? We serach for counter_C[j] in P[i,:] */
                  start_p=&(couple_pattern->index[couple_pattern->ptr[i]]);
                  where_p=(index_t*)bsearch(&(counter_C[j+my_n]), start_p,
                                            couple_pattern->ptr[i + 1] - couple_pattern->ptr[i],
                                            sizeof(index_t),
                                            Paso_comparIndex);
                  if (! (where_p == NULL) ) { /* yes i stronly connect with j */
                        offset = couple_pattern->ptr[i]+ (index_t)(where_p-start_p);
                        couple_block->val[offset]=A_ij; /* will be modified later */
                        if (A_ij< 0)  {
                           sum_strong_neg+=A_ij;
                        } else {
                           sum_strong_pos+=A_ij;
                        }
                  }
            }
         }

	 if(sum_strong_neg<0) { 
	    alpha= sum_all_neg/sum_strong_neg;
	 } else {
	    alpha=0;
	 }
	 if(sum_strong_pos>0) {
	    beta= sum_all_pos/sum_strong_pos;
	 } else {
	    beta=0;
	    A_ii+=sum_all_pos;
	 }
	 if ( A_ii > 0.) {
	    rtmp=(-1.)/A_ii;
	    alpha*=rtmp;
	    beta*=rtmp;
	 }

	 range = main_pattern->ptr[i + 1];
	 for (iPtr=main_pattern->ptr[i]; iPtr<range; iPtr++) {
	    A_ij=main_block->val[iPtr];
	    if (A_ij > 0 ) {
	       main_block->val[iPtr] = A_ij * beta;
	    } else {
	       main_block->val[iPtr] = A_ij * alpha;
	    }
	 }

         range = couple_pattern->ptr[i + 1];
         for (iPtr=couple_pattern->ptr[i]; iPtr<range; iPtr++) {
            A_ij=couple_block->val[iPtr];
            if (A_ij > 0 ) {
               couple_block->val[iPtr] = A_ij * beta;
            } else {
               couple_block->val[iPtr] = A_ij * alpha;
            }
         }
      }
   } 
}

void Paso_Preconditioner_AMG_setDirectProlongation_Block(Paso_SystemMatrix* P,
        Paso_SystemMatrix* A, const index_t* offset_S, const dim_t* degree_S,
        const index_t* S, const index_t *counter_C) {
   Paso_SparseMatrix *main_block=P->mainBlock;
   Paso_SparseMatrix *couple_block=P->col_coupleBlock;
   Paso_Pattern *main_pattern=main_block->pattern;
   Paso_Pattern *couple_pattern=couple_block->pattern;
   const dim_t row_block_size=A->row_block_size;
   const dim_t A_block = A->block_size;
   const dim_t my_n=A->mainBlock->numRows;
   index_t range;

   dim_t i;
   double *alpha, *beta, *sum_all_neg, *sum_all_pos, *sum_strong_neg, *sum_strong_pos, *A_ii;
   register double A_ij, rtmp;
   register index_t iPtr, j, offset, ib; 
   index_t *where_p, *start_p;
   
   #pragma omp parallel private(i,offset,ib,sum_all_neg,sum_all_pos,sum_strong_neg,sum_strong_pos,A_ii,range,iPtr,j,A_ij,start_p,where_p,alpha,beta,rtmp)
   {
      sum_all_neg=TMPMEMALLOC(row_block_size, double); /* sum of all negative values in row i of A */
      sum_all_pos=TMPMEMALLOC(row_block_size, double); /* sum of all positive values in row i of A */
      sum_strong_neg=TMPMEMALLOC(row_block_size, double); /* sum of all negative values A_ij where j is in C and strongly connected to i*/
      sum_strong_pos=TMPMEMALLOC(row_block_size, double); /* sum of all positive values A_ij where j is in C and strongly connected to i*/
      alpha=TMPMEMALLOC(row_block_size, double);
      beta=TMPMEMALLOC(row_block_size, double);
      A_ii=TMPMEMALLOC(row_block_size, double);
      
      #pragma omp for schedule(static)
      for (i=0;i<my_n;++i) {
	 if (counter_C[i]>=0) { /* i is a C row */
	    offset = main_pattern->ptr[i];
	    for (ib =0; ib<row_block_size; ++ib) 
		main_block->val[row_block_size*offset+ib]=1.; 
         } else if ((main_pattern->ptr[i + 1] > main_pattern->ptr[i]) ||
		    (couple_pattern->ptr[i + 1] > couple_pattern->ptr[i])) {
	    /* if i is an F row we first calculate alpha and beta: */
	    for (ib =0; ib<row_block_size; ++ib) {
	       sum_all_neg[ib]=0;
	       sum_all_pos[ib]=0;
	       sum_strong_neg[ib]=0;
	       sum_strong_pos[ib]=0;
	       A_ii[ib]=0;
	    }

	    /* first check the mainBlock */
	    range = A->mainBlock->pattern->ptr[i + 1];
	    for (iPtr=A->mainBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
	       j=A->mainBlock->pattern->index[iPtr];
	       if(j==i) {
		  for (ib =0; ib<row_block_size; ++ib) 
		    A_ii[ib]=A->mainBlock->val[A_block*iPtr+ib+row_block_size*ib];
	       } else {
		  for (ib =0; ib<row_block_size; ++ib) {
		     A_ij=A->mainBlock->val[A_block*iPtr+ib+row_block_size*ib];
		     if(A_ij< 0)  {
			sum_all_neg[ib]+=A_ij;
		     } else {
			sum_all_pos[ib]+=A_ij;
		     }
		  }
	       
		  if (counter_C[j]>=0) {
		     /* is i strongly connected with j? We search for counter_C[j] in P[i,:] */ 
		     start_p=&(main_pattern->index[main_pattern->ptr[i]]);
		     where_p=(index_t*)bsearch(&(counter_C[j]), start_p,
					     main_pattern->ptr[i + 1]-main_pattern->ptr[i],
					     sizeof(index_t),
					     Paso_comparIndex);
		     if (! (where_p == NULL) ) { /* yes i strongly connected with j */
			      offset = main_pattern->ptr[i]+ (index_t)(where_p-start_p);
			      for (ib =0; ib<row_block_size; ++ib) {
				 A_ij=A->mainBlock->val[A_block*iPtr+ib+row_block_size*ib];
				 main_block->val[row_block_size*offset+ib]=A_ij; /* will be modified later */
			         if (A_ij< 0)  {
				     sum_strong_neg[ib]+=A_ij;
				 } else {
				    sum_strong_pos[ib]+=A_ij;
				 }
			      }
		     }
		  } 
	    
	       } 
	    }

	    /* now we deal with the col_coupleBlock */
	    range = A->col_coupleBlock->pattern->ptr[i + 1];
            for (iPtr=A->col_coupleBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
                  j=A->col_coupleBlock->pattern->index[iPtr];
                  for (ib =0; ib<row_block_size; ++ib) {
                     A_ij=A->col_coupleBlock->val[A_block*iPtr+ib+row_block_size*ib];
                     if(A_ij< 0)  {
                        sum_all_neg[ib]+=A_ij;
                     } else {
                        sum_all_pos[ib]+=A_ij;
                     }
                  }

                  if (counter_C[j+my_n]>=0) {
                     /* is i stronly connect with j? We serach for counter_C[j] in P[i,:] */
                     start_p=&(couple_pattern->index[couple_pattern->ptr[i]]);
                     where_p=(index_t*)bsearch(&(counter_C[j+my_n]), start_p,
                                             couple_pattern->ptr[i + 1]-couple_pattern->ptr[i],
                                             sizeof(index_t),
                                             Paso_comparIndex);
                     if (! (where_p == NULL) ) { /* yes i stronly connect with j */
                              offset = couple_pattern->ptr[i]+ (index_t)(where_p-start_p);
                              for (ib =0; ib<row_block_size; ++ib) {
                                 A_ij=A->col_coupleBlock->val[A_block*iPtr+ib+row_block_size*ib];
                                 couple_block->val[row_block_size*offset+ib]=A_ij; /* will be modified later */

                                 if (A_ij< 0)  {
                                     sum_strong_neg[ib]+=A_ij;
                                 } else {
                                    sum_strong_pos[ib]+=A_ij;
                                 }
                              }
                     }
                  }
            }

	    for (ib =0; ib<row_block_size; ++ib) {
	       if(sum_strong_neg[ib]<0) { 
		  alpha[ib]= sum_all_neg[ib]/sum_strong_neg[ib];
	       } else {
		  alpha[ib]=0;
	       }
	       if(sum_strong_pos[ib]>0) {
		  beta[ib]= sum_all_pos[ib]/sum_strong_pos[ib];
	       } else {
		  beta[ib]=0;
		  A_ii[ib]+=sum_all_pos[ib];
	       }
	       if ( A_ii[ib] > 0.) {
		  rtmp=(-1./A_ii[ib]);
		  alpha[ib]*=rtmp;
		  beta[ib]*=rtmp;
	       }
	    }

	    range = main_pattern->ptr[i + 1];      
	    for (iPtr=main_pattern->ptr[i]; iPtr<range; iPtr++) {
	       for (ib =0; ib<row_block_size; ++ib) {
		  A_ij=main_block->val[row_block_size*iPtr+ib];
		  if (A_ij > 0 ) {
		     main_block->val[row_block_size*iPtr+ib] = A_ij * beta[ib];
		  } else {
		     main_block->val[row_block_size*iPtr+ib] = A_ij * alpha[ib];
		  }
	       }
	    }

	    range = couple_pattern->ptr[i + 1];
	    for (iPtr=couple_pattern->ptr[i]; iPtr<range; iPtr++) {
	      for (ib =0; ib<row_block_size; ++ib) {
                  A_ij=couple_block->val[row_block_size*iPtr+ib];
                  if (A_ij > 0 ) {
                     couple_block->val[row_block_size*iPtr+ib] = A_ij * beta[ib];
                  } else {
                     couple_block->val[row_block_size*iPtr+ib] = A_ij * alpha[ib];
                  }
               }
            }


	 }
      }/* end i loop */
      TMPMEMFREE(sum_all_neg); 
      TMPMEMFREE(sum_all_pos); 
      TMPMEMFREE(sum_strong_neg); 
      TMPMEMFREE(sum_strong_pos); 
      TMPMEMFREE(alpha);
      TMPMEMFREE(beta);
      TMPMEMFREE(A_ii);
   } /* end parallel region */
}

/*
    Classic Prolongation:
    -------------------

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise.
    If row i is not C, then P[i,j] = - 1/a[i] * ( A[i,k] + sum_{l} A[i,l]*A+[l,k]/B[i,k]) 
             where the summation over l is considering columns which are strongly connected 
             to i (l in S[i]) and not in C (counter_C[l]<0) and 

                B[i,k]=sum_{m in S_i and in C} A+[k,m]
                a[i]=A[i,i]+sum{l not strongly connected to i} A[i,l]

            A+[i,k]=A[i,k] if sign(A[i,k])==sign(A[i,i])  or 0 otherwise
              

*/
void Paso_Preconditioner_AMG_setClassicProlongation(Paso_SystemMatrix* P, 
	Paso_SystemMatrix* A, const index_t* offset_S, const dim_t* degree_S,
	const index_t* S, const index_t *counter_C) { 
   Paso_SparseMatrix *main_block=P->mainBlock;
   Paso_SparseMatrix *couple_block=P->col_coupleBlock;
   Paso_Pattern *main_pattern=main_block->pattern;
   Paso_Pattern *couple_pattern=couple_block->pattern;
   const dim_t my_n=A->mainBlock->numRows;
   index_t range, range_j;

   dim_t i, q;
   double *D_s=NULL;
   index_t *D_s_offset=NULL, iPtr, iPtr_j;
   const index_t *ptr_main_A = Paso_SparseMatrix_borrowMainDiagonalPointer(A->mainBlock);
   index_t *start_p_main_i, *start_p_couple_i;
   dim_t degree_p_main_i, degree_p_couple_i;
   dim_t len, ll, main_len, len_D_s;

   len = A->col_coupleBlock->numCols;
   len = MAX(A->remote_coupleBlock->numCols, len);
   ll = len + my_n;
   main_len = main_pattern->len;

   #pragma omp parallel private(D_s,D_s_offset,i,start_p_main_i,start_p_couple_i,degree_p_main_i,degree_p_couple_i,range,iPtr,q,range_j,iPtr_j,len_D_s)
   {
        D_s=TMPMEMALLOC(ll,double);
        D_s_offset=TMPMEMALLOC(ll,index_t);

	#pragma omp for schedule(static)
        for (i=0;i<my_n;++i) {
            if (counter_C[i]>=0) {
	        main_block->val[main_pattern->ptr[i]]=1.;  /* i is a C row */
            } else if ((main_pattern->ptr[i + 1] > main_pattern->ptr[i]) ||
		       (couple_pattern->ptr[i + 1] > couple_pattern->ptr[i])) {
	      /* this loop sums up the weak connections in a and creates
		 a list of the strong connected columns which are not in
		 C (=no interpolation nodes) */
	      const index_t *start_s = &(S[offset_S[i]]);
              const double A_ii = A->mainBlock->val[ptr_main_A[i]];
              double a=A_ii;

	      start_p_main_i = &(main_pattern->index[main_pattern->ptr[i]]);
              start_p_couple_i = &(couple_pattern->index[couple_pattern->ptr[i]]);
              degree_p_main_i = main_pattern->ptr[i+1] - main_pattern->ptr[i];
              degree_p_couple_i = couple_pattern->ptr[i+1] - couple_pattern->ptr[i];

	      /* first, check the mainBlock */
	      range = A->mainBlock->pattern->ptr[i + 1];
	      for (iPtr=A->mainBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
	         const index_t j=A->mainBlock->pattern->index[iPtr];
	         const double A_ij=A->mainBlock->val[iPtr];
                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
	            const index_t *where_s=(index_t*)bsearch(&j, start_s,degree_S[i],sizeof(index_t), Paso_comparIndex);
	            if (where_s == NULL) { /* weak connections are accumulated */
                        a+=A_ij;  
                    } else {   /* yes i strongly connected with j */
                        if  (counter_C[j]>=0)  { /* j is an interpolation point : add A_ij into P */
	                       const index_t *where_p=(index_t*)bsearch(&counter_C[j], start_p_main_i,degree_p_main_i, sizeof(index_t), Paso_comparIndex);
                               if (where_p == NULL)  {
                                       Esys_setError(SYSTEM_ERROR, "Paso_Preconditioner_setClassicProlongation: interpolation point is missing.");
                               } else {
  		                    const index_t offset = main_pattern->ptr[i]+ (index_t)(where_p-start_p_main_i);
  	                            main_block->val[offset]+=A_ij; 
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
                               double s=0.;
                               len_D_s=0;

			       /* first, the mainBlock part */
			       range_j = A->mainBlock->pattern->ptr[j + 1];
	                       for (iPtr_j=A->mainBlock->pattern->ptr[j]; iPtr_j<range_j; iPtr_j++) {
	                            const double A_jm=A->mainBlock->val[iPtr_j];
	                            const index_t m=A->mainBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
	                            const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_main_i,degree_p_main_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
  		                         const index_t offset_m = main_pattern->ptr[i]+ (index_t)(where_p_m-start_p_main_i);
                                         if (! SAMESIGN(A_ii,A_jm)) {
                                              D_s[len_D_s]=A_jm;
                                         } else {
                                              D_s[len_D_s]=0.;
                                         }
                                         D_s_offset[len_D_s]=offset_m;
                                         len_D_s++;
                                    }
                               }

			       /* then the coupleBlock part */
			       if (degree_p_couple_i) {
				 range_j = A->col_coupleBlock->pattern->ptr[j + 1];
				 for (iPtr_j=A->col_coupleBlock->pattern->ptr[j]; iPtr_j<range_j; iPtr_j++) {
                                    const double A_jm=A->col_coupleBlock->val[iPtr_j];
                                    const index_t m=A->col_coupleBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m+my_n], start_p_couple_i,degree_p_couple_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = couple_pattern->ptr[i]+ (index_t)(where_p_m-start_p_couple_i);
                                         if (! SAMESIGN(A_ii,A_jm)) {
                                              D_s[len_D_s]=A_jm;
                                         } else {
                                              D_s[len_D_s]=0.;
                                         }
                                         D_s_offset[len_D_s]=offset_m + main_len;
                                         len_D_s++;
                                    } 
				 }
                               }

                               for (q=0;q<len_D_s;++q) s+=D_s[q];
                               if (ABS(s)>0) {
                                   s=A_ij/s;
                                   for (q=0;q<len_D_s;++q) {
					if (D_s_offset[q] < main_len) 
                                          main_block->val[D_s_offset[q]]+=s*D_s[q];
					else 
					  couple_block->val[D_s_offset[q]-main_len]+=s*D_s[q];
                                   }
                               } else {
                                   a+=A_ij;
                               }
                          }
                     }
                 }
              }  

	      if (A->mpi_info->size > 1) {
              /* now, deal with the coupleBlock */
              range = A->col_coupleBlock->pattern->ptr[i + 1];
              for (iPtr=A->col_coupleBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
                 const index_t j=A->col_coupleBlock->pattern->index[iPtr];
                 const double A_ij=A->col_coupleBlock->val[iPtr];
                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
		    index_t t=j+my_n;
                    const index_t *where_s=(index_t*)bsearch(&t, start_s,degree_S[i],sizeof(index_t), Paso_comparIndex);
                    if (where_s == NULL) { /* weak connections are accummulated */
                        a+=A_ij;
                    } else {   /* yes i stronly connect with j */
                        if  (counter_C[t]>=0)  { /* j is an interpolation point : add A_ij into P */
                               const index_t *where_p=(index_t*)bsearch(&counter_C[t], start_p_couple_i,degree_p_couple_i, sizeof(index_t), Paso_comparIndex);
                               if (where_p == NULL)  {
                                       Esys_setError(SYSTEM_ERROR, "Paso_Preconditioner_AMG_setClassicProlongation: interpolation point is missing.");
                               } else {
                                    const index_t offset = couple_pattern->ptr[i]+ (index_t)(where_p-start_p_couple_i);
                                    couple_block->val[offset]+=A_ij;
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
			       double s=0.;
                               len_D_s=0;

                               /* first, the row_coupleBlock part */
                               range_j = A->row_coupleBlock->pattern->ptr[j + 1];
                               for (iPtr_j=A->row_coupleBlock->pattern->ptr[j]; iPtr_j<range_j; iPtr_j++) {
                                    const double A_jm=A->row_coupleBlock->val[iPtr_j];
                                    const index_t m=A->row_coupleBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_main_i,degree_p_main_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = main_pattern->ptr[i]+ (index_t)(where_p_m-start_p_main_i);
                                         if (! SAMESIGN(A_ii,A_jm)) {
                                              D_s[len_D_s]=A_jm;
                                         } else {
                                              D_s[len_D_s]=0.;
                                         }
                                         D_s_offset[len_D_s]=offset_m;
                                         len_D_s++;
                                    }
                               }

                               /* then the remote_coupleBlock part */
			       range_j = A->remote_coupleBlock->pattern->ptr[j + 1];
			       for (iPtr_j=A->remote_coupleBlock->pattern->ptr[j]; iPtr_j<range_j; iPtr_j++) {
                                    const double A_jm=A->remote_coupleBlock->val[iPtr_j];
                                    const index_t m=A->remote_coupleBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m+my_n], start_p_couple_i, degree_p_couple_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = couple_pattern->ptr[i]+ (index_t)(where_p_m-start_p_couple_i);
                                         if (! SAMESIGN(A_ii,A_jm)) {
                                              D_s[len_D_s]=A_jm;
                                         } else {
                                              D_s[len_D_s]=0.;
                                         }
                                         D_s_offset[len_D_s]=offset_m + main_len;
                                         len_D_s++;
                                    }
                               }

                               for (q=0;q<len_D_s;++q) s+=D_s[q];
                               if (ABS(s)>0) {
                                   s=A_ij/s;
                                   for (q=0;q<len_D_s;++q) {
					if (D_s_offset[q] < main_len)
                                          main_block->val[D_s_offset[q]]+=s*D_s[q];
					else
					  couple_block->val[D_s_offset[q]-main_len]+=s*D_s[q];
                                   }
                               } else {
                                   a+=A_ij;
                               }
                          }
                     }
                 }
              }
	      }

	      /* i has been processed, now we need to do some rescaling */
              if (ABS(a)>0.) {
                   a=-1./a;
		   range = main_pattern->ptr[i + 1];
                   for (iPtr=main_pattern->ptr[i]; iPtr<range; iPtr++) {
                        main_block->val[iPtr]*=a;
                   }

		   range = couple_pattern->ptr[i + 1];
                   for (iPtr=couple_pattern->ptr[i]; iPtr<range; iPtr++) {
                        couple_block->val[iPtr]*=a;
                   }
              }
          }
        }  /* end of row i loop */
        TMPMEMFREE(D_s);
        TMPMEMFREE(D_s_offset);
     }    /* end of parallel region */
}

void Paso_Preconditioner_AMG_setClassicProlongation_Block(
	Paso_SystemMatrix* P, Paso_SystemMatrix* A, const index_t* offset_S,
	const dim_t* degree_S, const index_t* S, const index_t *counter_C) {
   Paso_SparseMatrix *main_block=P->mainBlock;
   Paso_SparseMatrix *couple_block=P->col_coupleBlock;
   Paso_Pattern *main_pattern=main_block->pattern;
   Paso_Pattern *couple_pattern=couple_block->pattern;
   const dim_t row_block=A->row_block_size;
   const dim_t my_n=A->mainBlock->numRows;
   const dim_t A_block = A->block_size;
   index_t range, range_j;

   dim_t i, q, ib;
   double *D_s=NULL;
   index_t *D_s_offset=NULL, iPtr, iPtr_j;
   index_t *start_p_main_i, *start_p_couple_i;
   dim_t degree_p_main_i, degree_p_couple_i;
   dim_t len, ll, main_len, len_D_s;
   const index_t *ptr_main_A = Paso_SparseMatrix_borrowMainDiagonalPointer(A->mainBlock);
   
   len = A->col_coupleBlock->numCols;
   len = MAX(A->remote_coupleBlock->numCols, len);
   ll = len + my_n;
   main_len = main_pattern->len;
   #pragma omp parallel private(D_s,D_s_offset,i,ib,start_p_main_i,start_p_couple_i,degree_p_main_i,degree_p_couple_i,range,iPtr,q,range_j,iPtr_j,len_D_s)
   {
        double *a=TMPMEMALLOC(row_block, double);
        D_s=TMPMEMALLOC(row_block*ll,double);
        D_s_offset=TMPMEMALLOC(row_block*ll,index_t);

   	#pragma omp for private(i) schedule(static)
        for (i=0;i<my_n;++i) {
            if (counter_C[i]>=0) {
	        const index_t offset = main_pattern->ptr[i];
	        for (ib =0; ib<row_block; ++ib) main_block->val[row_block*offset+ib]=1.;  /* i is a C row */
            } else if ((main_pattern->ptr[i + 1] > main_pattern->ptr[i]) ||
		       (couple_pattern->ptr[i + 1] > couple_pattern->ptr[i])) {
	      /* this loop sums up the weak connections in a and creates
		 a list of the strong connected columns which are not in
		 C (=no interpolation nodes) */
	      const index_t *start_s = &(S[offset_S[i]]);
	      const double *A_ii = &(A->mainBlock->val[ptr_main_A[i]*A_block]);
	      start_p_main_i = &(main_pattern->index[main_pattern->ptr[i]]);
	      start_p_couple_i = &(couple_pattern->index[couple_pattern->ptr[i]]);
              degree_p_main_i = main_pattern->ptr[i+1] - main_pattern->ptr[i];
	      degree_p_couple_i = couple_pattern->ptr[i+1] - couple_pattern->ptr[i];
              for (ib=0; ib<row_block; ib++) a[ib]=A_ii[(row_block+1)*ib];

	      /* first, check the mainBlock */
	      range = A->mainBlock->pattern->ptr[i + 1];
	      for (iPtr=A->mainBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
	         const index_t j=A->mainBlock->pattern->index[iPtr];
	         const double* A_ij=&(A->mainBlock->val[iPtr*A_block]);

                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
	            const index_t *where_s=(index_t*)bsearch(&j, start_s,degree_S[i],sizeof(index_t), Paso_comparIndex);
	            if (where_s == NULL) { /* weak connections are accumulated */
                        for (ib=0; ib<row_block; ib++) a[ib]+=A_ij[(row_block+1)*ib];
                    } else {   /* yes i strongly connected with j */
                        if  (counter_C[j]>=0)  { /* j is an interpolation point : add A_ij into P */
	                       const index_t *where_p=(index_t*)bsearch(&counter_C[j], start_p_main_i,degree_p_main_i, sizeof(index_t), Paso_comparIndex);
                               if (where_p == NULL)  {
                                       Esys_setError(SYSTEM_ERROR, "Paso_Preconditioner_AMG_setClassicProlongation_Block: interpolation point is missing.");
                               } else {
  		                    const index_t offset = main_pattern->ptr[i]+ (index_t)(where_p-start_p_main_i);
                                    for (ib=0; ib<row_block; ib++) main_block->val[offset*row_block+ib] +=A_ij[(row_block+1)*ib];
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
                               len_D_s=0;

			       /* first, the mainBlock part */
			       range_j = A->mainBlock->pattern->ptr[j + 1];
	                       for (iPtr_j=A->mainBlock->pattern->ptr[j];iPtr_j<range_j; iPtr_j++) {
	                            const double* A_jm=&(A->mainBlock->val[iPtr_j*A_block]);
	                            const index_t m=A->mainBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
	                            const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_main_i,degree_p_main_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
  		                         const index_t offset_m = main_pattern->ptr[i]+ (index_t)(where_p_m-start_p_main_i);
                                         for (ib=0; ib<row_block; ib++) {
                                              if (! SAMESIGN(A_ii[(row_block+1)*ib],A_jm[(row_block+1)*ib]) ) {
                                                   D_s[len_D_s*row_block+ib]=A_jm[(row_block+1)*ib];
                                              } else {
                                                   D_s[len_D_s*row_block+ib]=0.;
                                              }
                                         }
                                         D_s_offset[len_D_s] = offset_m;
                                         len_D_s++;
                                    }
                               }

			       /* then the coupleBlock part */
			       range_j = A->col_coupleBlock->pattern->ptr[j+1];
			       for (iPtr_j=A->col_coupleBlock->pattern->ptr[j];iPtr_j<range_j; iPtr_j++) {
                                    const double* A_jm=&(A->col_coupleBlock->val[iPtr_j*A_block]);
                                    const index_t m=A->col_coupleBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m+my_n], start_p_couple_i,degree_p_couple_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = couple_pattern->ptr[i]+ (index_t)(where_p_m-start_p_couple_i);
                                         for (ib=0; ib<row_block; ib++) {
                                              if (! SAMESIGN(A_ii[(row_block+1)*ib],A_jm[(row_block+1)*ib]) ) {
                                                   D_s[len_D_s*row_block+ib]=A_jm[(row_block+1)*ib];
                                              } else {
                                                   D_s[len_D_s*row_block+ib]=0.;
                                              }
                                         }
                                         D_s_offset[len_D_s]=offset_m + main_len;
                                         len_D_s++;
                                    }
			       }

                               for (ib=0; ib<row_block; ib++) { 
                                   double s=0;
                                   for (q=0;q<len_D_s;++q) s+=D_s[q*row_block+ib];
                        
                                   if (ABS(s)>0) {
                                     s=A_ij[(row_block+1)*ib]/s;
				     for (q=0; q<len_D_s; q++) { 
				       if (D_s_offset[q] < main_len) 
                                            main_block->val[D_s_offset[q]*row_block+ib]+=s*D_s[q*row_block+ib];
				       else{
					    couple_block->val[(D_s_offset[q]-main_len)*row_block+ib]+=s*D_s[q*row_block+ib];
					}
				     }
                                   } else {
                                       a[ib]+=A_ij[(row_block+1)*ib];
                                   }
                               }
                          }
                     }
                 }
              }  

	      if (A->mpi_info->size > 1) {
	      /* now, deal with the coupleBlock */
	      range = A->col_coupleBlock->pattern->ptr[i + 1];
	      for (iPtr=A->col_coupleBlock->pattern->ptr[i]; iPtr<range; iPtr++) {
	         const index_t j=A->col_coupleBlock->pattern->index[iPtr];
	         const double* A_ij=&(A->col_coupleBlock->val[iPtr*A_block]);

                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
		    index_t t=j+my_n;
	            const index_t *where_s=(index_t*)bsearch(&t, start_s,degree_S[i],sizeof(index_t), Paso_comparIndex);
	            if (where_s == NULL) { /* weak connections are accummulated */
                        for (ib=0; ib<row_block; ib++) a[ib]+=A_ij[(row_block+1)*ib];
                    } else {   /* yes i strongly connected with j */
                        if  (counter_C[t]>=0)  { /* j is an interpolation point : add A_ij into P */
	                       const index_t *where_p=(index_t*)bsearch(&counter_C[t], start_p_couple_i,degree_p_couple_i, sizeof(index_t), Paso_comparIndex);
                               if (where_p == NULL)  {
                                       Esys_setError(SYSTEM_ERROR, "Paso_Preconditioner_AMG_setClassicProlongation_Block: interpolation point is missing.");

                               } else {
  		                    const index_t offset = couple_pattern->ptr[i]+ (index_t)(where_p-start_p_couple_i);
                                    for (ib=0; ib<row_block; ib++) couple_block->val[offset*row_block+ib] +=A_ij[(row_block+1)*ib];
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
                               len_D_s=0;

			       /* first, the row_coupleBlock part */
			       range_j = A->row_coupleBlock->pattern->ptr[j + 1];
	                       for (iPtr_j=A->row_coupleBlock->pattern->ptr[j];iPtr_j<range_j; iPtr_j++) {
                                    /* is m an interpolation point ? */
				    index_t m, *where_p_m;
				    double *A_jm;
				    A_jm=&(A->row_coupleBlock->val[iPtr_j*A_block]);
				    m=A->row_coupleBlock->pattern->index[iPtr_j];

				    where_p_m=(index_t*)bsearch(&counter_C[m], start_p_main_i,degree_p_main_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
  		                         const index_t offset_m = main_pattern->ptr[i]+ (index_t)(where_p_m-start_p_main_i);
                                         for (ib=0; ib<row_block; ib++) {
                                              if (! SAMESIGN(A_ii[(row_block+1)*ib],A_jm[(row_block+1)*ib]) ) {
                                                   D_s[len_D_s*row_block+ib]=A_jm[(row_block+1)*ib];
                                              } else {
                                                   D_s[len_D_s*row_block+ib]=0.;
                                              }
                                         }
                                         D_s_offset[len_D_s]=offset_m;
                                         len_D_s++;
                                    }
                               }

			       /* then the remote_coupleBlock part */
			       if (degree_p_couple_i) {
				 range_j = A->remote_coupleBlock->pattern->ptr[j + 1];
				 for (iPtr_j=A->remote_coupleBlock->pattern->ptr[j];iPtr_j<range_j; iPtr_j++) {
                                    const double* A_jm=&(A->remote_coupleBlock->val[iPtr_j*A_block]);
                                    const index_t m=A->remote_coupleBlock->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
                                    const index_t *where_p_m=(index_t*)bsearch(&counter_C[m+my_n], start_p_couple_i,degree_p_couple_i, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
                                         const index_t offset_m = couple_pattern->ptr[i]+ (index_t)(where_p_m-start_p_couple_i);
                                         for (ib=0; ib<row_block; ib++) {
                                              if (! SAMESIGN(A_ii[(row_block+1)*ib],A_jm[(row_block+1)*ib]) ) {
                                                   D_s[len_D_s*row_block+ib]=A_jm[(row_block+1)*ib];
                                              } else {
                                                   D_s[len_D_s*row_block+ib]=0.;
                                              }
                                         }
                                         D_s_offset[len_D_s]=offset_m + main_len;
                                         len_D_s++;
                                    } 
				 }
                               }

                               for (ib=0; ib<row_block; ib++) { 
                                   double s=0;
                                   for (q=0;q<len_D_s;++q) s+=D_s[q*row_block+ib];
                        
                                   if (ABS(s)>0) {
                                       s=A_ij[(row_block+1)*ib]/s;
				       for (q=0;q<len_D_s;++q) {
					 if (D_s_offset[q] < main_len) 
                                            main_block->val[D_s_offset[q]*row_block+ib]+=s*D_s[q*row_block+ib];
					 else
					    couple_block->val[(D_s_offset[q]-main_len)*row_block+ib]+=s*D_s[q*row_block+ib];
				       }
                                   } else {
                                       a[ib]+=A_ij[(row_block+1)*ib];
                                   }
                               }
                          }
                     }
                 }
	      }
	      }

	      /* i has been processed, now we need to do some rescaling */
              for (ib=0; ib<row_block; ib++) { 
                   register double a2=a[ib];
                   if (ABS(a2)>0.) {
                        a2=-1./a2;
			range = main_pattern->ptr[i + 1];
                        for (iPtr=main_pattern->ptr[i]; iPtr<range; iPtr++) {
                             main_block->val[iPtr*row_block+ib]*=a2;
                        }

			range = couple_pattern->ptr[i + 1];
			for (iPtr=couple_pattern->ptr[i]; iPtr<range; iPtr++) {
                             couple_block->val[iPtr*row_block+ib]*=a2;
                        }
                   }
              }
          }
        }  /* end of row i loop */
        TMPMEMFREE(D_s);
        TMPMEMFREE(D_s_offset);
        TMPMEMFREE(a);
     }    /* end of parallel region */
}
