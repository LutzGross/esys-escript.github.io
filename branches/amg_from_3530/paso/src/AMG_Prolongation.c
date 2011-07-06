
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

/* Paso: defines AMG prolongation  */

/**************************************************************/

/* Author: Artak Amirbekyan, artak@uq.edu.au, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

/**************************************************************

    Methods nessecary for AMG preconditioner

    construct n x n_C the prolongation matrix P from A_p.
    
    the columns in A_p to be considered are marked by counter_C[n] where
    an unknown i is to be considered in P is marked by 0<= counter_C[i] < n_C 
    and counter_C[i]  gives the new column number in P. S defines the strong connections.
    
    The pattern of P is formed as follows:

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise
    If row i is not C, then P[i,j] <> 0 if counter_C[k]==j (k in C) and (i,k) is strong connection.  
    
    two settings for P are implemented (see below) 
   
*/

#define MY_DEBUG 1
 
Paso_SystemMatrix* Paso_Preconditioner_AMG_getProlongation(Paso_SystemMatrix* A_p, 
                                                           const index_t* offset_S, const dim_t* degree_S, const index_t* S,
							   const dim_t n_C, const index_t* counter_C, const index_t interpolation_method) 
{
   Esys_MPIInfo *mpi_info=A_p->mpi_info;
   Paso_SparseMatrix *main_block=NULL, *couple_block=NULL;
   Paso_SystemMatrix *out=NULL;
   Paso_SystemMatrixPattern *pattern=NULL;
   Paso_Distribution *input_dist=NULL, *output_dist=NULL;
   Paso_SharedComponents *send =NULL, *recv=NULL;
   Paso_Connector *col_connector=NULL, *row_connector=NULL;
   Paso_Coupler *coupler=NULL;
   Paso_Pattern *main_pattern=NULL, *couple_pattern=NULL;
   const dim_t row_block_size=A_p->row_block_size;
   const dim_t col_block_size=A_p->col_block_size;
   const dim_t my_n=A_p->mainBlock->numCols;
   const dim_t overlap_n=A_p->col_coupleBlock->numCols;
   const dim_t n = my_n + overlap_n;
   const dim_t num_threads=omp_get_max_threads();
   double *couple_marker=NULL;
   index_t size=mpi_info->size, rank=mpi_info->rank, *dist=NULL;
   index_t *main_p=NULL, *couple_p=NULL, *main_idx=NULL, *couple_idx=NULL;
   index_t *shared=NULL, *offsetInShared=NULL;
   index_t sum, j, iptr;;
   dim_t i, my_n_C, k, l, p, q, global_label, num_neighbors;
   dim_t *recv_len=NULL, *send_len=NULL;
   Esys_MPI_rank *neighbor=NULL;

   if (MY_DEBUG) {
     fprintf(stderr, "size=%d rank=%d n=%d my_n=%d overlap_n=%d\n", 
		size, rank, n, my_n, overlap_n);
   }

   /* number of C points in current distribution */
   my_n_C = 0;
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

   if (MY_DEBUG) {
     fprintf(stderr, "my_n_C=%d\n", my_n_C);
   }

   /* create row distribution (output_distribution) and col distribution 
      (input_distribution) */
   /* ??? should I alloc an new Esys_MPIInfo object or reuse the one in
      system matrix A. for now, I'm reuse A->mpi_info ??? */
   dist = A_p->pattern->output_distribution->first_component;
   output_dist=Paso_Distribution_alloc(mpi_info, dist, 1, 0);
   dist = TMPMEMALLOC(size+1, index_t); /* now prepare for col distribution */
   Esys_checkPtr(dist);
   MPI_Allgather(&my_n_C, 1, MPI_INT, dist, 1, MPI_INT, mpi_info->comm);
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
   couple_marker = MEMALLOC(overlap_n, double);
   #pragma omp parallel for private(i) schedule(static)
   for (i=0; i<overlap_n; i++) couple_marker[i] = 0;
   if (!(Esys_checkPtr(main_p) || Esys_checkPtr(couple_p))) {
     /* count the number of entries per row in the Prolongation matrix :*/
     /* #pragma omp parallel for private(i,k,iptr,j,p) schedule(static) */
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
		  couple_marker[j-my_n] = 1;
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
     if (num_threads>1) {
	#pragma omp parallel private(i,p)
	{
          p=0;
	  #pragma omp for schedule(static)
	  for (i=0;i<overlap_n;++i) {
	    if (couple_marker[i] > 0) {
	      p++;
	    }
	  }
	}
	#pragma omp critical
	{
	  sum += p;
	}
     } else { /* num_threads=1 */
	for (i=0;i<overlap_n;++i) {
	  if (couple_marker[i] > 0) {
            sum++;
	  }
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

   /* now we need to create the Sparse Matrix (mainBlock and coupleBlock) */
   if (Esys_noError()) {
     main_block = Paso_SparseMatrix_alloc(MATRIX_FORMAT_DIAGONAL_BLOCK, 
		couple_pattern, row_block_size, col_block_size, FALSE);
     couple_block = Paso_SparseMatrix_alloc(MATRIX_FORMAT_DIAGONAL_BLOCK, 
		couple_pattern, row_block_size, col_block_size, FALSE);
   }

   /* now fill in the matrix */
/*   if (Esys_noError()) {
     if ((interpolation_method == PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING) 
	|| ( interpolation_method == PASO_CLASSIC_INTERPOLATION) ) {
	if (row_block_size == 1) {
	  Paso_Preconditioner_AMG_setClassicProlongation();
	} else {
	  Paso_Preconditioner_AMG_setClassicProlongation_Block();
	}
     } else {
	if (row_block_size == 1) {
	  Paso_Preconditioner_AMG_setDirectProlongation();
	} else {
	  Paso_Preconditioner_AMG_setDirectProlongation_Block();
	}
     }
   }
*/
   /* prepare the receiver for the col_connector. 
      Note that the allocation for "shared" assumes the send and receive buffer
      of the interpolation matrix P is no larger than that of matrix A_p. */
   coupler = Paso_Coupler_alloc(A_p->row_coupler->connector, 1);
   Paso_Coupler_startCollect(coupler, couple_marker);
   recv_len = TMPMEMALLOC(size,dim_t);
   send_len = TMPMEMALLOC(size,dim_t);
   neighbor = TMPMEMALLOC(size, Esys_MPI_rank);
   offsetInShared = TMPMEMALLOC(size+1, index_t);
   recv = A_p->col_coupler->connector->recv;
   i = recv->numSharedComponents;
   k = A_p->col_coupler->connector->send->numSharedComponents;
   if (k > i) i = k;
   shared = TMPMEMALLOC(i, index_t);
   memset(recv_len, 0, sizeof(dim_t)*size);
   num_neighbors = 0;
   q = 0;
   p = recv->numNeighbors;
   offsetInShared[0]=0;
   for (i=0; i<p; i++) {
     l = 0;
     k = recv->offsetInShared[i+1];
     for (j=recv->offsetInShared[i]; j<k; j++) {
	if (couple_marker[j] == 1) {
	  shared[q] = my_n_C + q;
	  q++;
	  l = 1;
	}
     }
     if (l == 1) {
	iptr = recv->neighbor[i];
	neighbor[num_neighbors] = iptr;
        recv_len[iptr] = q - offsetInShared[num_neighbors];
	num_neighbors++;
	offsetInShared[num_neighbors] = q;
     }
   }
   Paso_Coupler_finishCollect(coupler);
   recv = Paso_SharedComponents_alloc(my_n_C, num_neighbors, neighbor, shared,
                                      offsetInShared, 1, 0, mpi_info);

   /* now we can build the sender */
   #ifdef ESYS_MPI
     MPI_Alltoall(recv_len, 1, MPI_INT, send_len, 1, MPI_INT, mpi_info->comm);
   #else
     for (p=0; p<size; p++) snd_len[p] = rcv_len[p];
   #endif
   send = A_p->col_coupler->connector->send;
   num_neighbors = 0;
   q = 0;
   p = send->numNeighbors;
   offsetInShared[0]=0;
   for (i=0; i<p; i++) {
     l = 0;
     k = send->offsetInShared[i+1];
     for (j=send->offsetInShared[i]; j<k; j++) {
        if (coupler->recv_buffer[j] == 1) {
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
   Paso_Coupler_free(coupler);
   send = Paso_SharedComponents_alloc(my_n_C, num_neighbors, neighbor, shared,
				      offsetInShared, 1, 0, mpi_info);
   col_connector = Paso_Connector_alloc(send, recv);
   Paso_SharedComponents_free(recv);
   Paso_SharedComponents_free(send);
   TMPMEMFREE(recv_len);
   TMPMEMFREE(send_len);
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

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise
    If row i is not C, then P[i,j] = - a[i] * A[i,k]/A[i,i] with j=counter_C[k]>=0 and k in S
   
   and    a[i]= 
             alpha[i] = sum_s min(A[i,s],0)/(sum_{s in S and C} min(A[i,s],0))   A[i,k]<0
                   or                                                         if
             beta[i] = sum_s max(A[i,s],0)/(sum_{s in S and C} max(A[i,s],0))   A[i,k]>0
              

*/

void Paso_Preconditioner_AMG_setDirectProlongation(Paso_SparseMatrix* P_p, 
					           const Paso_SparseMatrix* A_p,
						   const index_t *counter_C) { 
   dim_t i;
   const dim_t n =A_p->numRows;
   register double alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos, A_ij, A_ii, rtmp;
   register index_t iPtr, j, offset; 
   index_t *where_p, *start_p;
   
   #pragma omp parallel for private(A_ii, offset, where_p, start_p, i, alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos,iPtr,j, A_ij , rtmp)  schedule(static)
   for (i=0;i<n;++i) {
      if (counter_C[i]>=0) {
	    offset = P_p->pattern->ptr[i];
	    P_p->val[offset]=1.;  /* i is a C row */
      } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
	 /* if i is an F row we first calculate alpha and beta: */
	 sum_all_neg=0; /* sum of all negative values in row i of A */
	 sum_all_pos=0; /* sum of all positive values in row i of A */
	 sum_strong_neg=0; /* sum of all negative values A_ij where j is in C and strongly connected to i*/
	 sum_strong_pos=0; /* sum of all positive values A_ij where j is in C and strongly connected to i*/
	 A_ii=0; 
	 for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
	    j=A_p->pattern->index[iPtr];
	    A_ij=A_p->val[iPtr];
	    if(j==i) {
	       A_ii=A_ij;
	    } else {
	       
	       if(A_ij< 0)  {
		  sum_all_neg+=A_ij;
	       } else {
		  sum_all_pos+=A_ij;
	       }
	       
	       if (counter_C[j]>=0) {
		  /* is i stronly connect with j? We serach for counter_C[j] in P[i,:] */ 
		  start_p=&(P_p->pattern->index[P_p->pattern->ptr[i]]);
		  where_p=(index_t*)bsearch(&(counter_C[j]), start_p,
					    P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i],
					    sizeof(index_t),
					    Paso_comparIndex);
		  if (! (where_p == NULL) ) { /* yes i stronly connect with j */
			offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
			P_p->val[offset]=A_ij; /* will be modified later */
			if (A_ij< 0)  {
			   sum_strong_neg+=A_ij;
			} else {
			   sum_strong_pos+=A_ij;
			}
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
	 for (iPtr=P_p->pattern->ptr[i];iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
	    A_ij=P_p->val[iPtr];
	    if (A_ij > 0 ) {
	       P_p->val[iPtr]*=beta;
	    } else {
	       P_p->val[iPtr]*=alpha;
	    }
	 }
      }
   } 
}

void Paso_Preconditioner_AMG_setDirectProlongation_Block(Paso_SparseMatrix* P_p, 
						         const Paso_SparseMatrix* A_p,
						         const index_t *counter_C) { 
   dim_t i;
   const dim_t n =A_p->numRows;
   const dim_t row_block=A_p->row_block_size;
   const dim_t A_block = A_p->block_size;
   double *alpha, *beta, *sum_all_neg, *sum_all_pos, *sum_strong_neg, *sum_strong_pos, *A_ii;
   register double A_ij, rtmp;
   register index_t iPtr, j, offset, ib; 
   index_t *where_p, *start_p;
   
   #pragma omp parallel private(ib, rtmp, A_ii, offset, where_p, start_p, i, alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos,iPtr,j, A_ij )  
   {
      sum_all_neg=TMPMEMALLOC(row_block, double); /* sum of all negative values in row i of A */
      sum_all_pos=TMPMEMALLOC(row_block, double); /* sum of all positive values in row i of A */
      sum_strong_neg=TMPMEMALLOC(row_block, double); /* sum of all negative values A_ij where j is in C and strongly connected to i*/
      sum_strong_pos=TMPMEMALLOC(row_block, double); /* sum of all positive values A_ij where j is in C and strongly connected to i*/
      alpha=TMPMEMALLOC(row_block, double);
      beta=TMPMEMALLOC(row_block, double);
      A_ii=TMPMEMALLOC(row_block, double);
      
      #pragma omp for schedule(static)
      for (i=0;i<n;++i) {
	 if (counter_C[i]>=0) {
	    offset = P_p->pattern->ptr[i];
	    for (ib =0; ib<row_block; ++ib) P_p->val[row_block*offset+ib]=1.;  /* i is a C row */
         } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
	    /* if i is an F row we first calculate alpha and beta: */
	    for (ib =0; ib<row_block; ++ib) {
	       sum_all_neg[ib]=0;
	       sum_all_pos[ib]=0;
	       sum_strong_neg[ib]=0;
	       sum_strong_pos[ib]=0;
	       A_ii[ib]=0;
	    }
	    for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
	       j=A_p->pattern->index[iPtr];
	       if(j==i) {
		  for (ib =0; ib<row_block; ++ib) A_ii[ib]=A_p->val[A_block*iPtr+ib+row_block*ib];
	       } else {
		  for (ib =0; ib<row_block; ++ib) {
		     A_ij=A_p->val[A_block*iPtr+ib+row_block*ib];
		     if(A_ij< 0)  {
			sum_all_neg[ib]+=A_ij;
		     } else {
			sum_all_pos[ib]+=A_ij;
		     }
		  }
	       
		  if (counter_C[j]>=0) {
		     /* is i stronly connect with j? We serach for counter_C[j] in P[i,:] */ 
		     start_p=&(P_p->pattern->index[P_p->pattern->ptr[i]]);
		     where_p=(index_t*)bsearch(&(counter_C[j]), start_p,
					     P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i],
					     sizeof(index_t),
					     Paso_comparIndex);
		     if (! (where_p == NULL) ) { /* yes i stronly connect with j */
			      offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
			      for (ib =0; ib<row_block; ++ib) {
				 A_ij=A_p->val[A_block*iPtr+ib+row_block*ib];
				 P_p->val[row_block*offset+ib]=A_ij; /* will be modified later */
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
	    for (ib =0; ib<row_block; ++ib) {
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
      
	    for (iPtr=P_p->pattern->ptr[i];iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
	       for (ib =0; ib<row_block; ++ib) {
		  A_ij=P_p->val[row_block*iPtr+ib];
		  if (A_ij > 0 ) {
		     P_p->val[row_block*iPtr+ib]*=beta[ib];
		  } else {
		     P_p->val[row_block*iPtr+ib]*=alpha[ib];
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

    If row i is in C (counter_C[i]>=0), then P[i,j]=1 if j==counter_C[i] or 0 otherwise
    If row i is not C, then P[i,j] = - 1/a[i] * ( A[i,k] + sum_{l} A[i,l]*A+[l,k]/B[i,k]) 
             where the summation over l is considering columns which are strongly connected 
             to i (l in S[i]) and not in C (counter_C[l]<0) and 

                B[i,k]=sum_{m in S_i and in C} A+[k,m]
                a[i]=A[i,i]+sum{l not strongly connected to i} A[i,l]

            A+[i,k]=A[i,k] if sign(A[i,k])==sign(A[i,i])  or 0 otherwise
              

*/
void Paso_Preconditioner_AMG_setClassicProlongation(Paso_SparseMatrix* P_p, 
				 	           Paso_SparseMatrix* A_p,
                                                   const index_t* offset_S, const dim_t* degree_S, const index_t* S,
						   const index_t *counter_C) { 
   dim_t i, q;
   const dim_t n =A_p->numRows;
   double *D_s=NULL;
   index_t *D_s_offset=NULL, iPtr, iPtr_j;
   const dim_t ll = Paso_Util_iMax(n, degree_S);
   const index_t *ptr_main_A = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   

   #pragma omp parallel  private(D_s, D_s_offset, iPtr, q, iPtr_j)
   {
        D_s=TMPMEMALLOC(ll,double);
        D_s_offset=TMPMEMALLOC(ll,index_t);


   	#pragma omp for private(i) schedule(static)
        for (i=0;i<n;++i) {
            if (counter_C[i]>=0) {
	        P_p->val[P_p->pattern->ptr[i]]=1.;  /* i is a C row */
            } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
	       const index_t *start_s = &(S[offset_S[i]]);
	       const index_t *start_p = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
               const dim_t degree_P_i   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
              /* this loop sums up the weak connections in a and creates a list of the strong connected columns 
                                                                      which are not in C (=no interpolation nodes) */
              const double A_ii = A_p->val[ptr_main_A[i]];
              double a=A_ii;

	      for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
	         const index_t j=A_p->pattern->index[iPtr];
	         const double A_ij=A_p->val[iPtr];
                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
	            const index_t *where_s=(index_t*)bsearch(&j, start_s,degree_S[i],sizeof(index_t), Paso_comparIndex);
	            if (where_s == NULL) { /* weak connections are accummulated */
                        a+=A_ij;  
                    } else {   /* yes i stronly connect with j */
                        if  (counter_C[j]>=0)  { /* j is an interpolation point : add A_ij into P */
	                       const index_t *where_p=(index_t*)bsearch(&counter_C[j], start_p,degree_P_i, sizeof(index_t), Paso_comparIndex);
                               if (where_p == NULL)  {
                                       Esys_setError(SYSTEM_ERROR, "Paso_Preconditioner_AMG_setBoomerProlongation: interpolation point is missing.");
                               } else {
  		                    const index_t offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
  	                            P_p->val[offset]+=A_ij; 
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
	                       const index_t *start_p_j = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
                               const dim_t degree_P_j   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
                               double s=0.;
                               dim_t len_D_s=0;
	                       for (iPtr_j=A_p->pattern->ptr[j];iPtr_j<A_p->pattern->ptr[j + 1]; ++iPtr_j) {
	                            const double A_jm=A_p->val[iPtr_j];
	                            const index_t m=A_p->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
	                            const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_j,degree_P_j, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
  		                         const index_t offset_m = P_p->pattern->ptr[i]+ (index_t)(where_p_m-start_p_j);
                                         if (! SAMESIGN(A_ii,A_jm)) {
                                              D_s[len_D_s]=A_jm;
                                         } else {
                                              D_s[len_D_s]=0.;
                                         }
                                         D_s_offset[len_D_s]=offset_m;
                                         len_D_s++;
                                    }
                               }
                               for (q=0;q<len_D_s;++q) s+=D_s[q];
                               if (ABS(s)>0) {
                                   s=A_ij/s;
                                   for (q=0;q<len_D_s;++q) {
                                        P_p->val[D_s_offset[q]]+=s*D_s[q];
                                   }
                               } else {
                                   a+=A_ij;
                               }
                          }
                     }
                 }
              }  /* i has been processed, now we need to do some rescaling */
              if (ABS(a)>0.) {
                   a=-1./a;
                   for (iPtr=P_p->pattern->ptr[i]; iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
                        P_p->val[iPtr]*=a;
                   }
              }
          }
        }  /* endo of row i loop */
        TMPMEMFREE(D_s);
        TMPMEMFREE(D_s_offset);
     }    /* end of parallel region */
}

void Paso_Preconditioner_AMG_setClassicProlongation_Block(Paso_SparseMatrix* P_p, 
				 	           Paso_SparseMatrix* A_p,
                                                   const index_t* offset_S, const dim_t* degree_S, const index_t* S,
						   const index_t *counter_C) { 
   dim_t i, q, ib;
   const dim_t row_block=A_p->row_block_size;
   const dim_t A_block = A_p->block_size;
   const dim_t n =A_p->numRows;
   double *D_s=NULL;
   index_t *D_s_offset=NULL, iPtr, iPtr_j;
   const dim_t ll = Paso_Util_iMax(n, degree_S);
   const index_t *ptr_main_A = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   

   #pragma omp parallel  private(D_s, D_s_offset, iPtr, q, iPtr_j,ib)
   {
        double *a=TMPMEMALLOC(row_block, double);
        D_s=TMPMEMALLOC(row_block*ll,double);
        D_s_offset=TMPMEMALLOC(row_block*ll,index_t);

   	#pragma omp for private(i) schedule(static)
        for (i=0;i<n;++i) {
            if (counter_C[i]>=0) {
	        const index_t offset = P_p->pattern->ptr[i];
	        for (ib =0; ib<row_block; ++ib) P_p->val[row_block*offset+ib]=1.;  /* i is a C row */
            } else if (P_p->pattern->ptr[i + 1] > P_p->pattern->ptr[i]) {
	       const index_t *start_s = &(S[offset_S[i]]);
	       const index_t *start_p = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
               const dim_t degree_P_i   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
              /* this loop sums up the weak connections in a and creates a list of the strong connected columns 
                                                                      which are not in C (=no interpolation nodes) */
              const double *A_ii = &(A_p->val[ptr_main_A[i]*A_block]);
              for (ib=0; ib<row_block; ib++) a[ib]=A_ii[(row_block+1)*ib];
              

	      for (iPtr=A_p->pattern->ptr[i];iPtr<A_p->pattern->ptr[i + 1]; ++iPtr) {
	         const index_t j=A_p->pattern->index[iPtr];
	         const double* A_ij=&(A_p->val[iPtr*A_block]);

                 if ( (i!=j) && (degree_S[j]>0) ) {
                    /* is (i,j) a strong connection ?*/
	            const index_t *where_s=(index_t*)bsearch(&j, start_s,degree_S[i],sizeof(index_t), Paso_comparIndex);
	            if (where_s == NULL) { /* weak connections are accummulated */
                        for (ib=0; ib<row_block; ib++) a[ib]+=A_ij[(row_block+1)*ib];
                    } else {   /* yes i stronly connect with j */
                        if  (counter_C[j]>=0)  { /* j is an interpolation point : add A_ij into P */
	                       const index_t *where_p=(index_t*)bsearch(&counter_C[j], start_p,degree_P_i, sizeof(index_t), Paso_comparIndex);
                               if (where_p == NULL)  {
                                       Esys_setError(SYSTEM_ERROR, "Paso_Preconditioner_AMG_setBoomerProlongation: interpolation point is missing.");
                               } else {
  		                    const index_t offset = P_p->pattern->ptr[i]+ (index_t)(where_p-start_p);
                                    for (ib=0; ib<row_block; ib++) P_p->val[offset*row_block+ib] +=A_ij[(row_block+1)*ib];
                               }
                          } else {  /* j is not an interpolation point */
                               /* find all interpolation points m of k */
	                       const index_t *start_p_j = &(P_p->pattern->index[P_p->pattern->ptr[i]]);
                               const dim_t degree_P_j   = P_p->pattern->ptr[i + 1]-P_p->pattern->ptr[i];
                               dim_t len_D_s=0;
	                       for (iPtr_j=A_p->pattern->ptr[j];iPtr_j<A_p->pattern->ptr[j + 1]; ++iPtr_j) {
	                            const double* A_jm=&(A_p->val[iPtr_j*A_block]);
	                            const index_t m=A_p->pattern->index[iPtr_j];
                                    /* is m an interpolation point ? */
	                            const index_t *where_p_m=(index_t*)bsearch(&counter_C[m], start_p_j,degree_P_j, sizeof(index_t), Paso_comparIndex);
                                    if (! (where_p_m==NULL)) {
  		                         const index_t offset_m = P_p->pattern->ptr[i]+ (index_t)(where_p_m-start_p_j);
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
                               for (ib=0; ib<row_block; ib++) { 
                                   double s=0;
                                   for (q=0;q<len_D_s;++q) s+=D_s[q*row_block+ib];
                        
                                   if (ABS(s)>0) {
                                       s=A_ij[(row_block+1)*ib]/s;
                                       for (q=0;q<len_D_s;++q) {
                                            P_p->val[D_s_offset[q]*row_block+ib]+=s*D_s[q*row_block+ib];
                                       }
                                   } else {
                                       a[ib]+=A_ij[(row_block+1)*ib];
                                   }
                               }
                          }
                     }
                 }
              }  /* i has been processed, now we need to do some rescaling */
              for (ib=0; ib<row_block; ib++) { 
                   register double a2=a[ib];
                   if (ABS(a2)>0.) {
                        a2=-1./a2;
                        for (iPtr=P_p->pattern->ptr[i]; iPtr<P_p->pattern->ptr[i + 1]; ++iPtr) {
                             P_p->val[iPtr*row_block+ib]*=a2;
                        }
                   }
              }
          }
        }  /* endo of row i loop */
        TMPMEMFREE(D_s);
        TMPMEMFREE(D_s_offset);
        TMPMEMFREE(a);
     }    /* end of parallel region */
}
