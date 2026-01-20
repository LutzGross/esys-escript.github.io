
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************/

/*   Paso: CSC/CSR sparse matrix pattern                                    */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_PATTERN_H__
#define __PASO_PATTERN_H__

#include "Paso.h"
#include "PasoException.h"

#include <escript/IndexList.h>

namespace paso {

struct Pattern;
typedef boost::shared_ptr<Pattern> Pattern_ptr;
typedef boost::shared_ptr<const Pattern> const_Pattern_ptr;

struct PASO_DLL_API Pattern : boost::enable_shared_from_this<Pattern>
{
    Pattern(int type, dim_t numOutput, dim_t numInput, index_t* ptr,
            index_t* index);

    ~Pattern();

    Pattern_ptr unrollBlocks(int newType, dim_t outputBlockSize,
                             dim_t inputBlockSize);

    Pattern_ptr getSubpattern(dim_t newNumRows, dim_t newNumCols,
                              const index_t* rowList,
                              const index_t* newColIndex) const;

    /// Searches for a maximal independent set MIS in the matrix pattern
    void mis(index_t* mis_marker) const;

    void reduceBandwidth(index_t* oldToNew);

    Pattern_ptr multiply(int type, const_Pattern_ptr other) const;

    Pattern_ptr binop(int type, const_Pattern_ptr other) const;

    index_t* borrowMainDiagonalPointer();

    static Pattern_ptr fromIndexListArray(dim_t n0, dim_t n,
            const escript::IndexList* index_list_array,
            index_t range_min, index_t range_max, index_t index_offset);

    index_t* borrowColoringPointer();

    dim_t getBandwidth(index_t* label) const;

    inline bool isEmpty() const
    {
        return (!ptr && !index);
    }

    inline dim_t getNumColors()
    {
        // make sure numColors is defined
        borrowColoringPointer();
        return numColors;
    }

    inline dim_t maxDeg() const
    {
        dim_t deg = 0;

        #pragma omp parallel for reduction(max:deg)
        for (dim_t i = 0; i < numInput; ++i) {
                deg=std::max(deg, ptr[i+1]-ptr[i]);
        }
        return deg;
    }

    // convert csr row ptr and col indices to harwell-boeing format
    inline void csrToHB()
    {
        // TODO: add openmp
        if (! (type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
            throw PasoException(
                "Paso: Harwell-Boeing format requires CSR format with index offset 1 and block size 1.");
        }

        if ( !(hb_row == NULL && hb_col == NULL) ) {
            return;
        }

        hb_row = new index_t[len];
        hb_col = new index_t[len];
        for (dim_t i=0, k=0; i<numOutput; i++)
        {
            for (dim_t j=ptr[i]; j<ptr[i+1]; j++, k++)
            {
                hb_row[k] = i+1;
                hb_col[k] = index[j-1];
            }
        }
    }

    // convert csr to harwell-boeing format using MUMPS_INT type
    // This is more efficient than csrToHB() when using MUMPS
    template<typename T>
    inline void csrToHB_typed()
    {
        // TODO: add openmp
        if (! (type & (MATRIX_FORMAT_OFFSET1 + MATRIX_FORMAT_BLK1)) ) {
            throw PasoException(
                "Paso: Harwell-Boeing format requires CSR format with index offset 1 and block size 1.");
        }

        if ( !(hb_row_typed == NULL && hb_col_typed == NULL) ) {
            return;  // Already created
        }

        hb_row_typed = new T[len];
        hb_col_typed = new T[len];
        for (dim_t i=0, k=0; i<numOutput; i++)
        {
            for (dim_t j=ptr[i]; j<ptr[i+1]; j++, k++)
            {
                static_cast<T*>(hb_row_typed)[k] = static_cast<T>(i+1);
                static_cast<T*>(hb_col_typed)[k] = static_cast<T>(index[j-1]);
            }
        }
    }

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
    index_t* main_iptr;
    // number of colors
    dim_t numColors;
    // coloring index: inputs with the same color are not connected
    index_t* coloring;
    // row indices in harwell-boeing format (index_t type)
    index_t* hb_row;
    // col indices in harwell-boeing format (index_t type)
    index_t* hb_col;
    // row indices in harwell-boeing format (typed version for solvers like MUMPS)
    void* hb_row_typed;
    // col indices in harwell-boeing format (typed version for solvers like MUMPS)
    void* hb_col_typed;
};


} // namespace paso

#endif // __PASO_PATTERN_H__

