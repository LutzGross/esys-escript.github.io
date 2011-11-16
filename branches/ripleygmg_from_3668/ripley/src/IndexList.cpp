
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/IndexList.h>

using namespace std;

namespace ripley {

Paso_Pattern *createPatternFromIndexMatrix(const IndexMatrix &matrix,
        dim_t n0, dim_t n, index_t rangeMin, index_t rangeMax, index_t offset)
{
    dim_t* ptr = MEMALLOC(n+1-n0, index_t);

    // get the number of connections per row
#pragma omp parallel for schedule(static)
    for (dim_t i = n0; i < n; ++i) {
        ptr[i-n0] = matrix[i].count(rangeMin, rangeMax);
    }
    // accumulate ptr
    index_t s = 0;
    for (dim_t i = n0; i < n; ++i) {
        index_t itmp = ptr[i - n0];
        ptr[i - n0] = s;
        s += itmp;
    }
    ptr[n - n0] = s;

    // fill index
    index_t* index = MEMALLOC(ptr[n-n0], index_t);
#pragma omp parallel for schedule(static)
    for (dim_t i = n0; i < n; ++i) {
        matrix[i].toArray(&index[ptr[i - n0]], rangeMin, rangeMax, offset);
    }

    return Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, n-n0, rangeMax+offset,
           ptr, index);
}

dim_t IndexList::count(index_t rangeMin, index_t rangeMax) const
{
    dim_t res = 0;
    list<index_t>::const_iterator it;
    for (it=m_index.begin(); it!=m_index.end(); it++)
        if (*it >= rangeMin && rangeMax > *it)
            ++res;
    return res;
}

void IndexList::insert(index_t idx)
{
    // only add idx if it doesn't exist yet
    if (find(m_index.begin(), m_index.end(), idx) == m_index.end())
        m_index.push_back(idx);
}

void IndexList::toArray(index_t *array, index_t rangeMin, index_t rangeMax,
                        index_t offset) const
{
    list<index_t>::const_iterator it;
    for (it=m_index.begin(); it!=m_index.end(); it++)
        if (*it >= rangeMin && rangeMax > *it)
            *array++ = (*it) + offset;
}


} // end of namespace ripley

