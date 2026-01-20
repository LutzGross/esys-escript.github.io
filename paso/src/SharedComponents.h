
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************/

/*   Paso: shared components                                                */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_SHAREDCOMPONENTS_H__
#define __PASO_SHAREDCOMPONENTS_H__

#include "Paso.h"

namespace paso {

struct SharedComponents;
typedef boost::shared_ptr<SharedComponents> SharedComponents_ptr;
typedef boost::shared_ptr<const SharedComponents> const_SharedComponents_ptr;

struct PASO_DLL_API SharedComponents
{
    SharedComponents(dim_t localLength, const std::vector<int>& neighbours,
                     const index_t* sharedArray,
                     const std::vector<index_t>& offset,
                     index_t m = 1, index_t b = 0)
        : local_length(localLength*m),
          neighbour(neighbours),
          offsetInShared(offset)
    {
        if (offset.empty()) {
            numSharedComponents = 0;
        } else {
            numSharedComponents = offset[neighbours.size()] * m;
        }
        shared = new index_t[numSharedComponents];
        if (!neighbours.empty() && !offset.empty()) {
            if (m != 1) {
                for (int i = 0; i < offsetInShared.size(); i++) {
                    offsetInShared[i] *= m;
                }
            }
#pragma omp parallel for
            for (dim_t i = 0; i < offset[neighbours.size()]; i++) {
                const index_t itmp = m * sharedArray[i] + b;
                for (dim_t j = 0; j < m; ++j)
                    shared[m*i+j] = itmp+j;
            }
        } else {
            offsetInShared[neighbours.size()] = 0;
        }
    }

    ~SharedComponents()
    {
        delete[] shared;
    }

    /// local array length shared
    dim_t local_length;

    /// list of the processors sharing values with this processor
    std::vector<int> neighbour;

    /// offsetInShared[i] points to the first input value in array shared
    /// for processor i. Has length numNeighbors+1
    std::vector<index_t> offsetInShared;

    /// list of the (local) components which are shared with other processors.
    /// Has length numSharedComponents
    index_t* shared;

    /// = offsetInShared[numNeighbours]
    dim_t numSharedComponents;
};

} // namespace paso

#endif // __PASO_SHAREDCOMPONENTS_H__

