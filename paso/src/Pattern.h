
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/*   Paso: CSC/CSR sparse matrix pattern                                    */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_PATTERN_H__
#define __PASO_PATTERN_H__

#include "Paso.h"

#include <escript/IndexList.h>

namespace paso {

struct Pattern;
typedef boost::shared_ptr<Pattern> Pattern_ptr;
typedef boost::shared_ptr<const Pattern> const_Pattern_ptr;

struct Pattern : boost::enable_shared_from_this<Pattern>
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
#pragma omp parallel
        {
            dim_t loc_deg=0;
#pragma omp for
            for (dim_t i = 0; i < numInput; ++i) {
                loc_deg=std::max(loc_deg, ptr[i+1]-ptr[i]);
            }
#pragma omp critical
            {
                deg = std::max(deg, loc_deg);
            }
        }
        return deg;
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
};


} // namespace paso

#endif // __PASO_PATTERN_H__

