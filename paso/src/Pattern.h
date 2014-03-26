
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

/*   Paso: CSC/CSR pattern                                                  */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_PATTERN_H__
#define __PASO_PATTERN_H__

#include "Paso.h"
#include "IndexList.h"

namespace paso {

struct Pattern
{
    int type;
    // Number of rows in the ptr array [CSR] / number of cols for CSC
    dim_t numOutput;
    // Number of cols [CSR]
    dim_t numInput;
    // number of non-zeros
    dim_t len;
    // ptr[n] to ptr[n+1] lists indices (in index) of non-zeros in row n
    index_t* ptr;
    // Non-major indices of non-zeros (in CSR this will be col numbers)
    index_t* index;
    // pointer to main diagonal entry
    index_t *main_iptr;
    // number of colors
    dim_t numColors;
    // coloring index: inputs with the same color are not connected
    index_t* coloring;
    dim_t reference_counter;
};

PASO_DLL_API
Pattern* Pattern_alloc(int type, dim_t numOutput, dim_t numInput, index_t* ptr, index_t* index);

PASO_DLL_API
Pattern* Pattern_getReference(Pattern*);

PASO_DLL_API
void Pattern_free(Pattern*);

PASO_DLL_API
int comparIndex(const void *, const void *);

PASO_DLL_API
Pattern* Pattern_unrollBlocks(Pattern*, int, dim_t, dim_t);

PASO_DLL_API
Pattern* Pattern_getSubpattern(Pattern*, dim_t, dim_t, const index_t*, const index_t*);

PASO_DLL_API
bool Pattern_isEmpty(Pattern* in);

PASO_DLL_API
void Pattern_mis(Pattern* pattern_p, index_t* mis_marker);

PASO_DLL_API
void Pattern_reduceBandwidth(Pattern* self, index_t* oldToNew);

PASO_DLL_API
void Pattern_color(Pattern* patter, index_t* num_colors, index_t* colorOf);
Pattern* Pattern_multiply(int type, Pattern* A, Pattern* B);

PASO_DLL_API
Pattern* Pattern_binop(int type, Pattern* A, Pattern* B);

PASO_DLL_API
index_t* Pattern_borrowMainDiagonalPointer(Pattern* A);

PASO_DLL_API
Pattern* Pattern_fromIndexListArray(dim_t n0, Paso_IndexListArray* index_list_array,index_t range_min,index_t range_max, index_t index_offset);

PASO_DLL_API
dim_t Pattern_getNumColors(Pattern* A);

PASO_DLL_API
index_t* Pattern_borrowColoringPointer(Pattern* A);

PASO_DLL_API
dim_t Pattern_maxDeg(Pattern* A);

} // namespace paso

#endif // __PASO_PATTERN_H__

