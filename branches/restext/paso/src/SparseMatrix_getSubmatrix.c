
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

/* Paso: SparseMatrix */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"

/**************************************************************

    returns the submatrix of A where rows are gathered by index row_list 
    and columns are selected by non-negative values of new_col_index.
    if new_col_index[i]>-1 new_col_index[i] gives the column of i in 
    the returned submatrix 

*/


Paso_SparseMatrix* Paso_SparseMatrix_getSubmatrix(Paso_SparseMatrix* A,int n_row_sub,int n_col_sub, index_t* row_list,index_t* new_col_index){
      Paso_Pattern* sub_pattern=NULL;
      Paso_SparseMatrix* out=NULL;
      index_t index_offset=(A->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
      int i,k,tmp,m,subpattern_row;
      int type=A->type;
      Paso_resetError();
      if (A->type & MATRIX_FORMAT_CSC) {
          Paso_setError(TYPE_ERROR,"gathering submatrices supports CSR matrix format only.");
      } else {
         sub_pattern=Paso_Pattern_getSubpattern(A->pattern,n_row_sub,n_col_sub,row_list,new_col_index);
         if (Paso_noError()) {
            /* create the return object */
            out=Paso_SparseMatrix_alloc(type,sub_pattern,A->row_block_size,A->col_block_size,TRUE);
            if (Paso_noError()) {
                 #pragma omp parallel for private(i,k,m,subpattern_row,tmp) schedule(static)
                 for (i=0;i<n_row_sub;++i) {
                     subpattern_row=row_list[i];
                     for (k=A->pattern->ptr[subpattern_row]-index_offset;k<A->pattern->ptr[subpattern_row+1]-index_offset;++k) {
                        tmp=new_col_index[A->pattern->index[k]-index_offset];
                        if (tmp>-1) {
                           #pragma ivdep
                           for (m=out->pattern->ptr[i]-index_offset;m<out->pattern->ptr[i+1]-index_offset;++m) {
                               if (out->pattern->index[m]==tmp+index_offset) {
                                   Paso_copyShortDouble(A->block_size,&(A->val[k*A->block_size]),&(out->val[m*A->block_size]));
                                   break;
                               }
                           }
                        }
                     }
                 }
            }
         }
         Paso_Pattern_free(sub_pattern);
      }
      return out;
}
