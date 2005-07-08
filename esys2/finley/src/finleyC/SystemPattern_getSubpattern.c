/* $Id$ */

/**************************************************************/

/* Finley: SystemMatrixPatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Finley.h"
#include "Common.h"
#include "Util.h"
#include "SystemPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */

Finley_SystemMatrixPattern* Finley_SystemMatrixPattern_getSubpattern(Finley_SystemMatrixPattern* pattern, \
                                           int new_n_rows, index_t* row_list,index_t* new_col_index) {
  Finley_SystemMatrixPattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,k,j,subpattern_row,tmp;
  dim_t i;
  Finley_ErrorCode=NO_ERROR;

  ptr=MEMALLOC(new_n_rows+1,index_t);
  if (! Finley_checkPtr(ptr))  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<new_n_rows+1;++i) ptr[i]=0;
        
        /* find the number column entries in each row */
        #pragma omp for private(i,k,j,subpattern_row) schedule(static)
        for (i=0;i<new_n_rows;++i) {
            j=0;
            subpattern_row=row_list[i];
            for (k=pattern->ptr[subpattern_row]-PTR_OFFSET;k<pattern->ptr[subpattern_row+1]-PTR_OFFSET;++k) 
               if (new_col_index[pattern->index[k]-INDEX_OFFSET]>-1) j++;
            ptr[i]=j;
        }
     }
     /* accummulate ptr */
     ptr[new_n_rows]=Finley_Util_cumsum(new_n_rows,ptr);
     index=MEMALLOC(ptr[new_n_rows],index_t);
     if (Finley_checkPtr(index))  {
        MEMFREE(ptr);
     } else {
        /* find the number column entries in each row */
        #pragma omp parallel for private(i,k,j,subpattern_row,tmp) schedule(static)
        for (i=0;i<new_n_rows;++i) {
             j=ptr[i];
             subpattern_row=row_list[i];
             for (k=pattern->ptr[subpattern_row]-PTR_OFFSET;k<pattern->ptr[subpattern_row+1]-PTR_OFFSET;++k) {
                tmp=new_col_index[pattern->index[k]-INDEX_OFFSET];
                if (tmp>-1) {
                    index[j]=tmp;
                    ++j;
                }
             }
        }
        /* create return value */
        out=Finley_SystemMatrixPattern_alloc(new_n_rows,ptr,index);
        if (Finley_ErrorCode!=NO_ERROR) {
          MEMFREE(index);
          MEMFREE(ptr);
        }
     }
  }
  return out;
}
/*
 * $Log$
 * Revision 1.4  2005/07/08 04:07:57  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.2.3  2005/06/29 02:34:56  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.2.2  2005/03/02 23:35:06  gross
 * reimplementation of the ILU in Finley. block size>1 still needs some testing
 *
 * Revision 1.1.2.1  2005/02/18 02:27:31  gross
 * two function that will be used for a reimplementation of the ILU preconditioner
 *
 *
 */
