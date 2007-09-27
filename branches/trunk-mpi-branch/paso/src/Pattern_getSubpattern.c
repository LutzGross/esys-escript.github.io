
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/* Paso: PatternPattern */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: gross@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "Pattern.h"
#include "PasoUtil.h"

/**************************************************************/

/* creates Pattern  */

Paso_Pattern* Paso_Pattern_getSubpattern(Paso_Pattern* pattern, \
                                           int newNumRows, int newNumCols, index_t* row_list,index_t* new_col_index) {
  index_t index_offset=(pattern->type & PATTERN_FORMAT_OFFSET1 ? 1:0);
  Paso_Pattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,k,j,subpattern_row,tmp;
  dim_t i;
  Paso_resetError();

  ptr=MEMALLOC(newNumRows+1,index_t);
  if (! Paso_checkPtr(ptr))  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<newNumRows+1;++i) ptr[i]=0;
        
        /* find the number column entries in each row */
        #pragma omp for private(i,k,j,subpattern_row) schedule(static)
        for (i=0;i<newNumRows;++i) {
            j=0;
            subpattern_row=row_list[i];
            #pragma ivdep
            for (k=pattern->ptr[subpattern_row]-index_offset;k<pattern->ptr[subpattern_row+1]-index_offset;++k) 
               if (new_col_index[pattern->index[k]-index_offset]>-1) j++;
            ptr[i]=j;
        }
     }
     /* accummulate ptr */
     ptr[newNumRows]=Paso_Util_cumsum(newNumRows,ptr);
     index=MEMALLOC(ptr[newNumRows],index_t);
     if (Paso_checkPtr(index))  {
        MEMFREE(ptr);
     } else {
        /* find the number column entries in each row */
        #pragma omp parallel for private(i,k,j,subpattern_row,tmp) schedule(static)
        for (i=0;i<newNumRows;++i) {
             j=ptr[i];
             subpattern_row=row_list[i];
             #pragma ivdep
             for (k=pattern->ptr[subpattern_row]-index_offset;k<pattern->ptr[subpattern_row+1]-index_offset;++k) {
                tmp=new_col_index[pattern->index[k]-index_offset];
                if (tmp>-1) {
                    index[j]=tmp;
                    ++j;
                }
             }
        }
        /* create return value */
        out=Paso_Pattern_alloc(pattern->type,pattern->input_block_size,pattern->output_block_size,newNumRows,ptr,index);
        if (! Paso_noError()) {
          MEMFREE(index);
          MEMFREE(ptr);
        }
     }
  }
  return out;
}
