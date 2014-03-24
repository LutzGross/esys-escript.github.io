
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


/************************************************************************************/

/* Paso: Pattern */

/************************************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/

#include "Paso.h"
#include "Pattern.h"
#include "PasoUtil.h"

/************************************************************************************/

/* creates Pattern  */

Paso_Pattern* Paso_Pattern_getSubpattern(Paso_Pattern* pattern,
                                         int newNumRows, int newNumCols,
                                         const index_t* row_list,
                                         const index_t* new_col_index)
{
  index_t index_offset=(pattern->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
  Paso_Pattern*out=NULL;
  index_t *ptr=NULL,*index=NULL,k,j,subpattern_row,tmp;
  dim_t i;
  Esys_resetError();

  ptr=new index_t[newNumRows+1];
  if (! Esys_checkPtr(ptr))  {
     #pragma omp parallel
     {
        #pragma omp for private(i) schedule(static)
        for (i=0;i<newNumRows+1;++i) ptr[i]=0;
        
        /* find the number of column entries in each row */
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
     /* accumulate ptr */
     ptr[newNumRows]=Paso_Util_cumsum(newNumRows,ptr);
     index=new index_t[ptr[newNumRows]];
     if (Esys_checkPtr(index))  {
        delete[] ptr;
     } else {
        /* find the number of column entries in each row */
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
        out=Paso_Pattern_alloc(pattern->type,newNumRows,newNumCols,ptr,index);
        if (! Esys_noError()) {
          delete[] index;
          delete[] ptr;
        }
     }
  }
  return out;
}

