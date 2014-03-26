
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
 
/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "Pattern.h"

namespace paso {

/* allocates a Pattern  */
Pattern* Pattern_alloc(int type, dim_t numOutput, dim_t numInput, index_t* ptr,
                       index_t* index)
{
    index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    index_t loc_min_index,loc_max_index,min_index=index_offset,max_index=index_offset-1;
    dim_t i;
    Esys_resetError();

    if (ptr!=NULL && index != NULL) {
#pragma omp parallel private(loc_min_index,loc_max_index,i)
      {
        loc_min_index=index_offset;
        loc_max_index=index_offset-1;
        if (type & MATRIX_FORMAT_OFFSET1) {
#pragma omp for schedule(static) 
            for (i=0; i < numOutput; ++i) {
                if (ptr[i]<ptr[i+1]) {
#ifdef USE_QSORTG
                    qsortG(&(index[ptr[i]-1]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t), comparIndex); 
#else
                    qsort(&(index[ptr[i]-1]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t), comparIndex); 
#endif
                    loc_min_index=MIN(loc_min_index,index[ptr[i]-1]);
                    loc_max_index=MAX(loc_max_index,index[ptr[i+1]-2]);
                }
            }
        } else {
#pragma omp for schedule(static) 
            for (i=0; i < numOutput; ++i) {
                if (ptr[i] < ptr[i+1]) {
#ifdef USE_QSORTG
                    qsortG(&(index[ptr[i]]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t), comparIndex); 
#else
                    qsort(&(index[ptr[i]]),(size_t)(ptr[i+1]-ptr[i]),sizeof(index_t), comparIndex); 
#endif
                    loc_min_index=MIN(loc_min_index,index[ptr[i]]);
                    loc_max_index=MAX(loc_max_index,index[ptr[i+1]-1]);
                }
            }
        }
        #pragma omp critical
        {
            min_index=MIN(loc_min_index,min_index);
            max_index=MAX(loc_max_index,max_index);
        }
      } // parallel section

        if ( (min_index<index_offset) || (max_index>=numInput+index_offset) ) {
            Esys_setError(TYPE_ERROR,"Pattern_alloc: Pattern index out of range.");
            return NULL;
        }
    }
    Pattern* out = new Pattern;
    if (!Esys_checkPtr(out)) {
        out->type = type;
        out->reference_counter = 1;
        out->numOutput = numOutput;
        out->numInput = numInput;
        out->ptr = ptr;
        out->index = index;
        out->main_iptr = NULL;
        out->coloring = NULL;
        out->numColors = -1;

        if (out->ptr == NULL) {
            out->len = 0;
        } else {
            out->len=out->ptr[out->numOutput] - index_offset;
        }
    }
    return out;
}

/* returns a reference to in */
Pattern* Pattern_getReference(Pattern* in)
{
    if (in != NULL) {
        ++(in->reference_counter);
    }
    return in;
}
  
/* deallocates a Pattern */
void Pattern_free(Pattern* in)
{
    if (in != NULL) {
        in->reference_counter--;
        if (in->reference_counter <= 0) {
            delete[] in->ptr;
            delete[] in->index;
            delete[] in->main_iptr;
            delete[] in->coloring;
            delete in;
        }
    }
}
/* **************************************************************************/

/*  some routines which help to get the matrix pattern from elements: */

// this routine is used by qsort called in Pattern_alloc
int comparIndex(const void *index1, const void *index2)
{
    index_t Iindex1,Iindex2;
    Iindex1=*(index_t*)index1;
    Iindex2=*(index_t*)index2;
    if (Iindex1 < Iindex2) {
        return -1;
    } else {
        if (Iindex1 > Iindex2) {
            return 1;
        } else {
            return 0;
        }
    }
}

bool Pattern_isEmpty(Pattern* in)
{
    if (in != NULL) {
        if ((in->ptr != NULL) && (in->index != NULL))
            return FALSE;
    }
    return TRUE;
}

/* creates a pattern from a range of indices */
Pattern* Pattern_fromIndexListArray(dim_t n0,
        Paso_IndexListArray* index_list_array,
        index_t range_min, index_t range_max, index_t index_offset)
{
    Pattern* out=NULL;
    index_t* index=NULL;
    const dim_t n=index_list_array->n;
    dim_t* ptr = new index_t[n+1-n0];

    if (!Esys_checkPtr(ptr)) {
        Paso_IndexList* index_list = index_list_array->index_list;

        // get the number of connections per row
#pragma omp parallel for schedule(static)
        for (dim_t i=n0; i < n; ++i) {
            ptr[i-n0]=Paso_IndexList_count(&index_list[i],range_min,range_max);
        }
        // accumulate ptr
        dim_t s=0;
        for (dim_t i=n0; i < n; ++i) {
            const dim_t itmp=ptr[i-n0];
            ptr[i-n0]=s;
            s+=itmp;
        }
        ptr[n-n0]=s;
        // fill index
        index=new index_t[ptr[n-n0]];
        if (!Esys_checkPtr(index)) {
#pragma omp parallel for schedule(static) 
            for (dim_t i=n0; i < n; ++i) {
                Paso_IndexList_toArray(&index_list[i], &index[ptr[i-n0]],
                                       range_min, range_max, index_offset);
            }
            out=Pattern_alloc(MATRIX_FORMAT_DEFAULT, n-n0,
                              range_max+index_offset, ptr, index);
        }
    }
    if (!Esys_noError()) {
        delete[] ptr;
        delete[] index;
        Pattern_free(out);
    }
    return out;
}

index_t* Pattern_borrowMainDiagonalPointer(Pattern* A) 
{
    if (A->main_iptr == NULL) {
        const dim_t n=A->numOutput;
        A->main_iptr=new index_t[n];
        if (!Esys_checkPtr(A->main_iptr)) {
            bool fail = false;
            // identify the main diagonals
#pragma omp parallel for schedule(static)
            for (index_t i = 0; i < n; ++i) {
                index_t *index=&(A->index[A->ptr[i]]);
                index_t *where_p=reinterpret_cast<index_t*>(bsearch(&i,
                            index, (size_t)(A->ptr[i+1] - A->ptr[i]),
                            sizeof(index_t), comparIndex));
                if (where_p == NULL) {
                    fail = true;
                } else {
                    A->main_iptr[i]=A->ptr[i]+(index_t)(where_p-index);
                }
            }
            if (fail) {
                delete[] A->main_iptr;
                A->main_iptr=NULL;
            }
        }
    }
    return A->main_iptr;
}
              
dim_t Pattern_getNumColors(Pattern* A)
{
    // make sure numColors is defined
    Pattern_borrowColoringPointer(A);
    return A->numColors;
}

index_t* Pattern_borrowColoringPointer(Pattern* A)
{
    // is coloring available?
    if (A->coloring == NULL) {
        const dim_t n = A->numInput;
        A->coloring = new index_t[n];
        if (!Esys_checkPtr(A->coloring)) {
            Pattern_color(A, &(A->numColors), A->coloring);
            if (!Esys_noError()) {
                delete[] A->coloring;
                A->coloring = NULL;
            }
        } 
    }
    return A->coloring;
}

dim_t Pattern_maxDeg(Pattern* A)
{
    dim_t deg = 0;
    const dim_t n=A->numInput;

#pragma omp parallel
    {
        dim_t loc_deg=0;
#pragma omp for schedule(static)
        for (dim_t i = 0; i < n; ++i) {
            loc_deg=MAX(loc_deg, A->ptr[i+1]-A->ptr[i]);
        }
#pragma omp critical
        {
            deg = MAX(deg, loc_deg);
        }
    }
    return deg;
}

} // namespace paso

