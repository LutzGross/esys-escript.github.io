
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

#ifndef __RIPLEY_INDEXLIST_H__
#define __RIPLEY_INDEXLIST_H__

#include <ripley/Ripley.h>

extern "C" {
#include <paso/Pattern.h>
}

namespace ripley {

class IndexList;
typedef std::vector<IndexList> IndexMatrix;

/**
   \brief
   Creates a Paso pattern from a range of indices.
*/
Paso_Pattern *createPatternFromIndexMatrix(const IndexMatrix &matrix,
        dim_t n0, dim_t n, index_t rangeMin, index_t rangeMax, index_t offset);

/**
   \brief
   An IndexList stores a row (DOF) of indices to the non-zero columns of a
   matrix.
*/

class IndexList
{
public:
    /**
       \brief
       Constructor.
    */
    IndexList() {}

    /**
       \brief
       Destructor.
    */
    ~IndexList() {}

    /**
       \brief
       Returns the number of row indices that are within the given range.
    */
    dim_t count(index_t rangeMin, index_t rangeMax) const;

    /**
       \brief
       Inserts a row index into the index list if it does not yet exist.
    */
    void insert(index_t idx);

    /**
       \brief
       Copies the given range of row indices into the array after adding an
       offset to the indices.
    */
    void toArray(index_t *array, index_t rangeMin, index_t rangeMax,
                 index_t offset) const;

private:
    std::list<index_t> m_index;
};

} // end of namespace ripley

#endif // __RIPLEY_INDEXLIST_H__

