
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


/****************************************************************************/

/* Paso: Pattern */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003, 2004, 2005 */
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Pattern.h"
#include "PasoUtil.h"

namespace paso {

// creates a subpattern
Pattern* Pattern_getSubpattern(Pattern* pattern,
                               int newNumRows, int newNumCols,
                               const index_t* row_list,
                               const index_t* new_col_index)
{
    const index_t index_offset=(pattern->type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    Esys_resetError();

    index_t* ptr = new index_t[newNumRows+1];
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (dim_t i=0; i < newNumRows+1; ++i) ptr[i]=0;
    
        /* find the number of column entries in each row */
#pragma omp for schedule(static)
        for (dim_t i=0; i < newNumRows; ++i) {
            index_t j=0;
            const index_t subpattern_row=row_list[i];
            #pragma ivdep
            for (index_t k=pattern->ptr[subpattern_row]-index_offset; k < pattern->ptr[subpattern_row+1]-index_offset; ++k)
                if (new_col_index[pattern->index[k]-index_offset] > -1) j++;
            ptr[i]=j;
        }
    } // parallel section

    // accumulate ptr
    ptr[newNumRows]=Paso_Util_cumsum(newNumRows, ptr);
    index_t* index = new index_t[ptr[newNumRows]];

    // find the number of column entries in each row
#pragma omp parallel for schedule(static)
    for (dim_t i=0; i < newNumRows; ++i) {
        index_t j=ptr[i];
        const index_t subpattern_row=row_list[i];
        #pragma ivdep
        for (index_t k=pattern->ptr[subpattern_row]-index_offset; k < pattern->ptr[subpattern_row+1]-index_offset; ++k) {
            const index_t tmp=new_col_index[pattern->index[k]-index_offset];
            if (tmp > -1) {
                index[j]=tmp;
                ++j;
            }
        }
    }
    // create return value
    Pattern* out=Pattern_alloc(pattern->type, newNumRows, newNumCols, ptr, index);
    if (!Esys_noError()) {
        delete[] index;
        delete[] ptr;
    }
    return out;
}

} // namespace paso

