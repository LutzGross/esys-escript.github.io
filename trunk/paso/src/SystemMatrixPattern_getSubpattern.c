/* $Id$ */

/**************************************************************/

/* Paso: SystemMatrixPatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Util.h"
#include "SystemMatrixPattern.h"

/**************************************************************/

/* creates SystemMatrixPattern  */

Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_getSubpattern(Paso_SystemMatrixPattern* pattern, \
                                           int new_n_rows, index_t* row_list,index_t* new_col_index) {
  index_t index_offset=(pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  Paso_SystemMatrixPattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,k,j,subpattern_row,tmp;
  dim_t i;
  Paso_resetError();

  ptr=MEMALLOC(new_n_rows+1,index_t);
  if (! Paso_checkPtr(ptr))  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<new_n_rows+1;++i) ptr[i]=0;
        
        /* find the number column entries in each row */
        #pragma omp for private(i,k,j,subpattern_row) schedule(static)
        for (i=0;i<new_n_rows;++i) {
            j=0;
            subpattern_row=row_list[i];
            for (k=pattern->ptr[subpattern_row]-index_offset;k<pattern->ptr[subpattern_row+1]-index_offset;++k) 
               if (new_col_index[pattern->index[k]-index_offset]>-1) j++;
            ptr[i]=j;
        }
     }
     /* accummulate ptr */
     ptr[new_n_rows]=Paso_Util_cumsum(new_n_rows,ptr);
     index=MEMALLOC(ptr[new_n_rows],index_t);
     if (Paso_checkPtr(index))  {
        MEMFREE(ptr);
     } else {
        /* find the number column entries in each row */
        #pragma omp parallel for private(i,k,j,subpattern_row,tmp) schedule(static)
        for (i=0;i<new_n_rows;++i) {
             j=ptr[i];
             subpattern_row=row_list[i];
             for (k=pattern->ptr[subpattern_row]-index_offset;k<pattern->ptr[subpattern_row+1]-index_offset;++k) {
                tmp=new_col_index[pattern->index[k]-index_offset];
                if (tmp>-1) {
                    index[j]=tmp;
                    ++j;
                }
             }
        }
        /* create return value */
        out=Paso_SystemMatrixPattern_alloc(pattern->type,new_n_rows,ptr,index);
        if (! Paso_noError()) {
          MEMFREE(index);
          MEMFREE(ptr);
        }
     }
  }
  return out;
}
/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
