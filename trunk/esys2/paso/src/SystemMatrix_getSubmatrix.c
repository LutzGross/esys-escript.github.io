/* $Id$ */

/**************************************************************/

/* Paso: SystemMatrix */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004,2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Util.h"


/**************************************************************

    returns the submatrix of A where rows are gathered by index row_list 
    and columns are selected by non-negative values of new_col_index.
    if new_col_index[i]>-1 new_col_index[i] gives the column of i in 
    the returned submatrix 

*/


Paso_SystemMatrix* Paso_SystemMatrix_getSubmatrix(Paso_SystemMatrix* A,int n_row_sub,index_t* row_list,index_t* new_col_index){
      Paso_SystemMatrixPattern* sub_pattern=NULL;
      Paso_SystemMatrix* out=NULL;
      Paso_resetError();
      int i,k,tmp,m,subpattern_row;
      int type=A->type;
      if (type!=CSR) {
          Paso_setError(TYPE_ERROR,"gathering submatrices supports CSR matrix format only.");
      } else {
         sub_pattern=Paso_SystemMatrixPattern_getSubpattern(A->pattern,n_row_sub,row_list,new_col_index);
         if (Paso_noError()) {
            /* create the return object */
            out=Paso_SystemMatrix_alloc(type,sub_pattern,A->row_block_size,A->col_block_size);
            if (Paso_noError()) {
                 #pragma omp parallel for private(i,k,m,subpattern_row,tmp) schedule(static)
                 for (i=0;i<n_row_sub;++i) {
                     subpattern_row=row_list[i];
                     for (k=A->pattern->ptr[subpattern_row]-PTR_OFFSET;k<A->pattern->ptr[subpattern_row+1]-PTR_OFFSET;++k) {
                        tmp=new_col_index[A->pattern->index[k]-INDEX_OFFSET];
                        if (tmp>-1) {
                           for (m=out->pattern->ptr[i]-PTR_OFFSET;m<out->pattern->ptr[i+1]-PTR_OFFSET;++m) {
                               if (out->pattern->index[m]==tmp+INDEX_OFFSET) {
                                   Paso_copyDouble(A->block_size,&(A->val[k*A->block_size]),&(out->val[m*A->block_size]));
                                   break;
                               }
                           }
                        }
                     }
                 }
            }
         }
         Paso_SystemMatrixPattern_dealloc(sub_pattern);
      }
      return out;
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:48  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
