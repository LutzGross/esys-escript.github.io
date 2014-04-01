
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

#include "Pattern.h"
#include "PasoUtil.h"

using esysUtils::IndexListArray;

namespace paso {

Pattern::Pattern(int ntype, dim_t numOut, dim_t numIn, index_t* inPtr,
                 index_t* idx) :
    type(ntype),
    numOutput(numOut),
    numInput(numIn),
    len(0),
    ptr(inPtr),
    index(idx),
    main_iptr(NULL),
    numColors(-1),
    coloring(NULL)
{
    const index_t index_offset = (ntype & MATRIX_FORMAT_OFFSET1 ? 1:0);
    index_t min_index = index_offset, max_index = index_offset-1;
    Esys_resetError();

    if (inPtr!=NULL && idx != NULL) {
#pragma omp parallel
      {
        index_t loc_min_index=index_offset;
        index_t loc_max_index=index_offset-1;
        if (ntype & MATRIX_FORMAT_OFFSET1) {
#pragma omp for schedule(static)
            for (dim_t i=0; i < numOut; ++i) {
                if (inPtr[i] < inPtr[i+1]) {
                    qsort(&(idx[inPtr[i]-1]),(size_t)(inPtr[i+1]-inPtr[i]),sizeof(index_t), comparIndex);
                    loc_min_index = MIN(loc_min_index, idx[inPtr[i]-1]);
                    loc_max_index = MAX(loc_max_index, idx[inPtr[i+1]-2]);
                }
            }
        } else {
#pragma omp for schedule(static)
            for (dim_t i=0; i < numOut; ++i) {
                if (inPtr[i] < inPtr[i+1]) {
                    qsort(&(idx[inPtr[i]]),(size_t)(inPtr[i+1]-inPtr[i]),sizeof(index_t), comparIndex);
                    loc_min_index = MIN(loc_min_index, idx[inPtr[i]]);
                    loc_max_index = MAX(loc_max_index, idx[inPtr[i+1]-1]);
                }
            }
        }
        #pragma omp critical
        {
            min_index=MIN(loc_min_index, min_index);
            max_index=MAX(loc_max_index, max_index);
        }
      } // parallel section

      if (min_index < index_offset || max_index >= numIn+index_offset) {
          Esys_setError(TYPE_ERROR, "Pattern: Pattern index out of range.");
      }
      len = ptr[numOutput] - index_offset;
    }
}

/* deallocates a Pattern */
Pattern::~Pattern()
{
    delete[] ptr;
    delete[] index;
    delete[] main_iptr;
    delete[] coloring;
}

/* creates a pattern from a range of indices */
Pattern_ptr Pattern::fromIndexListArray(dim_t n0, dim_t n,
                                        const IndexListArray& index_list_array,
                                        index_t range_min, index_t range_max,
                                        index_t index_offset)
{
    dim_t* ptr = new index_t[n+1-n0];

    // get the number of connections per row
#pragma omp parallel for schedule(static)
    for (dim_t i=n0; i < n; ++i) {
        ptr[i-n0]=index_list_array[i].count(range_min, range_max);
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
    index_t* index = new index_t[ptr[n-n0]];
#pragma omp parallel for schedule(static)
    for (dim_t i=n0; i < n; ++i) {
        index_list_array[i].toArray(&index[ptr[i-n0]], range_min, range_max,
                                    index_offset);
    }
    Pattern_ptr out(new Pattern(MATRIX_FORMAT_DEFAULT, n-n0,
                                range_max+index_offset, ptr, index));

    if (!Esys_noError()) {
        delete[] ptr;
        delete[] index;
        out.reset();
    }
    return out;
}

index_t* Pattern::borrowMainDiagonalPointer()
{
    if (main_iptr == NULL) {
        const dim_t n = numOutput;
        main_iptr=new index_t[n];
        bool fail = false;
        // identify the main diagonals
#pragma omp parallel for schedule(static)
        for (index_t i = 0; i < n; ++i) {
            index_t* idx = &index[ptr[i]];
            index_t* where_p=reinterpret_cast<index_t*>(bsearch(&i, idx,
                        (size_t)(ptr[i+1] - ptr[i]),
                        sizeof(index_t), comparIndex));
            if (where_p == NULL) {
                fail = true;
            } else {
                main_iptr[i] = ptr[i]+(index_t)(where_p-idx);
            }
        }
        if (fail) {
            delete[] main_iptr;
            main_iptr = NULL;
        }
    }
    return main_iptr;
}

index_t* Pattern::borrowColoringPointer()
{
    // is coloring available?
    if (coloring == NULL) {
        coloring = new index_t[numInput];
        index_t out = 0;
        const dim_t n = numOutput;
        index_t* mis_marker = new index_t[n];
        // get coloring
#pragma omp parallel for schedule(static)
        for (index_t i = 0; i < n; i++) {
            coloring[i] = -1;
            mis_marker[i] = -1;
        }

        while (Paso_Util_isAny(n, coloring, -1) && Esys_noError()) {
            mis(mis_marker);

#pragma omp parallel for schedule(static)
            for (index_t i = 0; i < n; ++i) {
                if (mis_marker[i])
                    coloring[i] = out;
                mis_marker[i] = coloring[i];
            }
            ++out;
        }
        delete[] mis_marker;
        numColors = out;
    }
    return coloring;
}

// creates a subpattern
Pattern_ptr Pattern::getSubpattern(int newNumRows, int newNumCols,
                                   const index_t* row_list,
                                   const index_t* new_col_index) const
{
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    Esys_resetError();

    index_t* newPtr = new index_t[newNumRows+1];
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (dim_t i=0; i < newNumRows+1; ++i) newPtr[i]=0;
    
        // find the number of column entries in each row
#pragma omp for schedule(static)
        for (dim_t i=0; i < newNumRows; ++i) {
            index_t j=0;
            const index_t subpattern_row=row_list[i];
            #pragma ivdep
            for (index_t k=ptr[subpattern_row]-index_offset;
                         k < ptr[subpattern_row+1]-index_offset; ++k) {
                if (new_col_index[index[k]-index_offset] > -1)
                    j++;
            }
            newPtr[i]=j;
        }
    } // parallel section

    // accumulate ptr
    newPtr[newNumRows]=Paso_Util_cumsum(newNumRows, newPtr);
    index_t* newIndex = new index_t[newPtr[newNumRows]];

    // find the number of column entries in each row
#pragma omp parallel for schedule(static)
    for (dim_t i=0; i < newNumRows; ++i) {
        index_t j=newPtr[i];
        const index_t subpattern_row=row_list[i];
        #pragma ivdep
        for (index_t k=ptr[subpattern_row]-index_offset; k < ptr[subpattern_row+1]-index_offset; ++k) {
            const index_t tmp=new_col_index[index[k]-index_offset];
            if (tmp > -1) {
                newIndex[j]=tmp;
                ++j;
            }
        }
    }
    // create return value
    Pattern_ptr out(new Pattern(type, newNumRows, newNumCols, newPtr, newIndex));
    if (!Esys_noError()) {
        delete[] newIndex;
        delete[] newPtr;
    }
    return out;
}

Pattern_ptr Pattern::unrollBlocks(int newType, dim_t output_block_size,
                                  dim_t input_block_size)
{
    Pattern_ptr out;
    const index_t index_offset_in=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    const index_t index_offset_out=(newType & MATRIX_FORMAT_OFFSET1 ? 1:0);

    Esys_resetError();
    if (output_block_size == 1 && input_block_size == 1 &&
          (type & MATRIX_FORMAT_OFFSET1) == (newType & MATRIX_FORMAT_OFFSET1)) {
        out = shared_from_this();
    } else {
        const dim_t block_size = output_block_size*input_block_size;
        const dim_t new_len = len * block_size;
        const dim_t new_numOutput = numOutput * output_block_size;
        const dim_t new_numInput = numInput * input_block_size;

        index_t* newPtr = new index_t[new_numOutput+1];
        index_t* newIndex = new index_t[new_len];
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (dim_t i=0; i < new_numOutput+1; ++i)
                newPtr[i]=index_offset_out;
   
            #pragma omp single
            newPtr[new_numOutput]=new_len+index_offset_out;
   
            #pragma omp for schedule(static) 
            for (dim_t i=0; i < numOutput; ++i) {
                for (dim_t k=0; k < output_block_size; ++k) {
                    newPtr[i*output_block_size+k]=(ptr[i]-index_offset_in)*block_size+
                                                       (ptr[i+1]-ptr[i])*input_block_size*k+index_offset_out;
                }
            }
             
#pragma omp for schedule(static) 
            for (dim_t i=0; i < new_numOutput; ++i) {
                #pragma ivdep
                for (index_t iPtr=newPtr[i]-index_offset_out; iPtr<newPtr[i+1]-index_offset_out; ++iPtr) {
                    newIndex[iPtr]=index_offset_out;
                }
            }
   
            #pragma omp for schedule(static) 
            for (dim_t i=0; i < numOutput; ++i) {
                for (index_t iPtr=ptr[i]-index_offset_in; iPtr < ptr[i+1]-index_offset_in; ++iPtr) {
                    for (dim_t k=0; k < output_block_size; ++k) {
                        #pragma ivdep
                        for (dim_t j=0; j < input_block_size; ++j) {
                            newIndex[newPtr[i*output_block_size+k]-index_offset_out+(iPtr-(ptr[i]-index_offset_in))*input_block_size+j] = 
                                (index[iPtr]-index_offset_in)*input_block_size+j+index_offset_out;
                        }
                    }
                }
            }
        } // parallel section
        out.reset(new Pattern(newType, new_numOutput, new_numInput, newPtr, newIndex));
        if (!Esys_noError()) {
            delete[] newIndex;
            delete[] newPtr;
        }
    }
    return out;
}

// computes the pattern coming from a matrix-matrix multiplication
Pattern_ptr Pattern::multiply(int type, const_Pattern_ptr B) const
{
    IndexListArray index_list(numOutput);

#pragma omp parallel for schedule(static)
    for (dim_t i = 0; i < numOutput; i++) {
        for (index_t iptrA = ptr[i]; iptrA < ptr[i+1]; ++iptrA) {
            const dim_t j = index[iptrA];
            for (index_t iptrB = B->ptr[j]; iptrB < B->ptr[j+1]; ++iptrB) {
                const dim_t k = B->index[iptrB];
                index_list[i].insertIndex(k);
            }
        }
    }
    return Pattern::fromIndexListArray(0, numOutput, index_list, 0,
                                       B->numInput, 0);
}

/*
 * Computes the pattern  of C = A binary operation B for CSR matrices A,B
 * Note: we do not check whether A_ij(op)B_ij=0
 */
Pattern_ptr Pattern::binop(int type, const_Pattern_ptr B) const
{
    IndexListArray index_list(numOutput);

#pragma omp parallel for schedule(static)
    for (dim_t i = 0; i < B->numOutput; i++) {
        index_t iptrA = ptr[i];
        index_t iptrB = B->ptr[i];

        while (iptrA < ptr[i+1] && iptrB < B->ptr[i+1]) {
            const dim_t j = index[iptrA];
            const dim_t k = B->index[iptrB];
            if (j < k) {
                index_list[i].insertIndex(j);
                iptrA++;
            } else if (j > k) {
                index_list[i].insertIndex(k);
                iptrB++;
            } else { // (j == k)
                index_list[i].insertIndex(j);
                iptrB++;
                iptrA++;
            }
        }
        while(iptrA < ptr[i+1]) {
            const dim_t j = index[iptrA];
            index_list[i].insertIndex(j);
            iptrA++;
        }
        while(iptrB < B->ptr[i+1]) {
            const dim_t k = B->index[iptrB];
            index_list[i].insertIndex(k);
            iptrB++;
        }
    }

    return Pattern::fromIndexListArray(0, numOutput, index_list, 0,
                                       numInput, 0);
}

} // namespace paso

