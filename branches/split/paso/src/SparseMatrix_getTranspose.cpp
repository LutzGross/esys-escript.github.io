
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/************************************************************************************

  Paso: transposed matrix 

************************************************************************************

   Author: Lutz Gross, l.gross@uq.edu.au 

************************************************************************************/

#include "SparseMatrix.h"

Paso_SparseMatrix* Paso_SparseMatrix_getTranspose(Paso_SparseMatrix* A)
{
   
   Paso_Pattern *ATpattern=NULL;
   Paso_SparseMatrix *AT=NULL;
   
   const dim_t m=A->numCols;
   const dim_t n=A->numRows;
   const dim_t block_size=A->block_size;
   const dim_t col_block_size_A = A->col_block_size;
   const dim_t row_block_size_A = A->row_block_size;
   register index_t iptr_AT,  jptr_A, *start_p, *where_p, iptr2;
   dim_t i;
   register dim_t j, ib, irb, icb;
   
   Paso_IndexListArray* index_list = Paso_IndexListArray_alloc(m);
   for (i=0;i<n;++i) {
      for (iptr2=A->pattern->ptr[i];iptr2<A->pattern->ptr[i+1]; ++iptr2) {
	 j=A->pattern->index[iptr2];
	 Paso_IndexListArray_insertIndex(index_list,j,i);
       }
   }
	 
   ATpattern=Paso_Pattern_fromIndexListArray(0,index_list,0,n,0);
   Paso_IndexListArray_free(index_list);
   AT=Paso_SparseMatrix_alloc(A->type,ATpattern,col_block_size_A,row_block_size_A,FALSE);
   Paso_Pattern_free(ATpattern);
 
   if (  ( (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) && (block_size == 1 ) ) || 
         ( (row_block_size_A == 1 ) && (col_block_size_A == 1)            )  ) {
	 #pragma omp parallel for private(i, iptr_AT, j, jptr_A, start_p, where_p)
	 for (i=0;i<m;++i) {
	    for (iptr_AT=AT->pattern->ptr[i];iptr_AT<AT->pattern->ptr[i+1]; ++iptr_AT) {
	       j=AT->pattern->index[iptr_AT];
	       jptr_A=A->pattern->ptr[j];
	       start_p=&(A->pattern->index[jptr_A]);
	       where_p=(index_t*)bsearch(&i, start_p,
					 A->pattern->ptr[j + 1]-jptr_A,
					 sizeof(index_t),
					 Paso_comparIndex);
		if (! (where_p == NULL) ) { /* this should always be the case */
		    jptr_A += (index_t)(where_p-start_p);
		    AT->val[iptr_AT]=A->val[jptr_A];
		}
	    }
	 }
   } else {
      if (A->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
	 #pragma omp parallel for private(i, iptr_AT, j, jptr_A, start_p, where_p, ib)
	 for (i=0;i<m;++i) {
	    for (iptr_AT=AT->pattern->ptr[i];iptr_AT<AT->pattern->ptr[i+1]; ++iptr_AT) {
	       j=AT->pattern->index[iptr_AT];
	       jptr_A=A->pattern->ptr[j];
	       start_p=&(A->pattern->index[jptr_A]);
	       where_p=(index_t*)bsearch(&i, start_p,
					 A->pattern->ptr[j + 1]-jptr_A,
					 sizeof(index_t),
					 Paso_comparIndex);
		if (! (where_p == NULL) ) { /* this should always be the case */
		    jptr_A += (index_t)(where_p-start_p);
		    for (ib=0; ib < block_size; ++ib )  AT->val[iptr_AT*block_size+ib]=A->val[jptr_A*block_size+ib];
		}
	    }
	 }
      } else {
	 #pragma omp parallel for private(i, iptr_AT, j, jptr_A, start_p, where_p, irb, icb)
	 for (i=0;i<m;++i) {
	    for (iptr_AT=AT->pattern->ptr[i];iptr_AT<AT->pattern->ptr[i+1]; ++iptr_AT) {
	       j=AT->pattern->index[iptr_AT];
	       jptr_A=A->pattern->ptr[j];
	       start_p=&(A->pattern->index[jptr_A]);
	       where_p=(index_t*)bsearch(&i, start_p,
				       A->pattern->ptr[j + 1]-jptr_A,
				       sizeof(index_t),
				       Paso_comparIndex);
	       if (! (where_p == NULL) ) { /* this should always be the case */
		  jptr_A += (index_t)(where_p-start_p);
		  for (irb =0 ; irb < row_block_size_A; ++irb) {
		     for (icb =0 ; icb < col_block_size_A; ++icb) {
			AT->val[iptr_AT*block_size+icb+col_block_size_A*irb]=A->val[jptr_A*block_size+irb+row_block_size_A*icb];
		     }
		  }
	       }
	    }
	 }
      }
   }
   return AT;
}						      
