
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
    
    the columns in A_p to be considered ar marked by counter_C[n] where
    an unknown i to be considered in P is marked by 0<= counter_C[i] < n_C 
    gives the new column number in P. S defines the strong connections.
    
    If row i is in C, then P[i,j]=counter_C[i] if counter_C[i]>=0 and 0 otherwise
    If row i is not C, then P[i,j] = - a[i] * A[i,k]/A[i,i] with j=counter_C[k]>=0 and k in S
   
   and    a[i]= 
             alpha[i] = sum_s min(A[i,s],0)/(sum_{s in S and C} min(A[i,s],0))   A[i,k]<0
                   or                                                         if
             beta[i] = sum_s max(A[i,s],0)/(sum_{s in S and C} max(A[i,s],0))   A[i,k]>0
              

*/
 
Paso_SparseMatrix* Paso_Preconditioner_AMG_getDirectProlongation(const Paso_SparseMatrix* A_p,
								 const dim_t* degree, const index_t* S,
								 const dim_t n_C, const index_t* counter_C) 
{
   Paso_SparseMatrix* out=NULL;
   Paso_Pattern *outpattern=NULL;
   const dim_t n_block=A_p->row_block_size;
   index_t *ptr=NULL, *index=NULL,j, iptr;
   const dim_t n =A_p->numRows; 
   dim_t i,p,z, len_P;
   
   ptr=MEMALLOC(n+1,index_t);
   if (! Esys_checkPtr(ptr)) {

      
      /* count the number of entries per row in the Prolongation matrix :*/
   
      #pragma omp parallel for private(i,z,iptr,j,p) schedule(static)
      for (i=0;i<n;++i) {
	 if (counter_C[i]>=0) {
	    z=1;    /* i is a C unknown */
	 } else {
	    z=0;
	    iptr=A_p->pattern->ptr[i]; 
	    for (p=0; p<degree[i]; ++p) { 
	       j=S[iptr+p];  /* this is a strong connection */
	       if (counter_C[j]>=0) z++; /* and is in C */
	    }
	 }
	 ptr[i]=z;
      }
      len_P=Paso_Util_cumsum(n,ptr);
      ptr[n]=len_P;
      
      /* allocate and create index vector for prolongation: */
      index=MEMALLOC(len_P,index_t);
   
      if (! Esys_checkPtr(index)) {
	 #pragma omp parallel for private(i,z,iptr,j,p)  schedule(static)
	 for (i=0;i<n;++i) {
	    if (counter_C[i]>=0) {
	       index[ptr[i]]=counter_C[i];  /* i is a C unknown */
	    } else {
	       z=0;
	       iptr=A_p->pattern->ptr[i]; 
	       for (p=0; p<degree[i]; ++p) { 
		  j=S[iptr+p];  /* this is a strong connection */
		  if (counter_C[j]>=0) {  /* and is in C */
		     index[ptr[i]+z]=counter_C[j];
		     z++; /* and is in C */
		  }
	       }
	    }
	 } 
      }
   }   
   if (Esys_noError()) {
	 outpattern=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,n,n_C,ptr,index);
   } else {
      MEMFREE(ptr);
      MEMFREE(index);
   }
   /* now we need to create a matrix and fill it */
   if (Esys_noError()) {
      out=Paso_SparseMatrix_alloc(MATRIX_FORMAT_DIAGONAL_BLOCK,outpattern,n_block,n_block,FALSE);
   } 
   
   if (Esys_noError()) {
      if (n_block == 1) {
	 Paso_Preconditioner_AMG_setDirectProlongation(out, A_p, counter_C);
      } else {
	 Paso_Preconditioner_AMG_setDirectProlongation_Block(out, A_p, counter_C);
      }
   }
   Paso_Pattern_free(outpattern);
   if (Esys_noError()) {
      return out;
   } else {
      Paso_SparseMatrix_free(out);
      return NULL;
   }
   
}

void Paso_Preconditioner_AMG_setDirectProlongation(Paso_SparseMatrix* P_p, 
					           const Paso_SparseMatrix* A_p,
						   const index_t *counter_C) { 
   dim_t i;
   const dim_t n =A_p->numRows;
   register double alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos, A_ij, A_ii;
   register index_t iPtr, j, offset; 
   index_t *where_p, *start_p;
   
   #pragma omp parallel for private(A_ii, offset, where_p, start_p, i, alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos,iPtr,j, A_ij )  schedule(static)
   for (i=0;i<n;++i) {
      if (counter_C[i]>=0) {
	    offset = P_p->pattern->ptr[i];
	    P_p->val[offset]=1.;  /* i is a C row */
      } else {
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
	    alpha*=(-1./A_ii);
	    beta*=(-1./A_ii);
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
   register double A_ij;
   register index_t iPtr, j, offset, ib; 
   index_t *where_p, *start_p;
   
   #pragma omp parallel private(ib, A_ii, offset, where_p, start_p, i, alpha, beta, sum_all_neg, sum_all_pos, sum_strong_neg, sum_strong_pos,iPtr,j, A_ij )  
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
	 } else {
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
		  alpha[ib]*=(-1./A_ii[ib]);
		  beta[ib]*=(-1./A_ii[ib]);
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
