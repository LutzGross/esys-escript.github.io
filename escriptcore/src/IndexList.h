
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_INDEXLIST_H__
#define __ESCRIPT_INDEXLIST_H__

#include <escript/DataTypes.h>

// pre-reserving saves time under OpenMP. The 85 is a value taken over
// from revision ~101 by jgs.
#define ESYS_INDEXLIST_LENGTH 85

namespace escript {

struct IndexList {
    IndexList() : n(0), extension(NULL) {}
    ~IndexList() { delete extension; }

    DataTypes::index_t m_list[ESYS_INDEXLIST_LENGTH];
    DataTypes::dim_t n;
    IndexList* extension;

    /// inserts row index into the IndexList in if it does not exist
    inline void insertIndex(DataTypes::index_t index)
    {
        for (DataTypes::dim_t i=0; i<n; i++) {
            if (m_list[i] == index)
                return;
        }
        if (n < ESYS_INDEXLIST_LENGTH) {
            m_list[n++] = index;
        } else {
            if (extension == NULL)
                extension = new IndexList();
            extension->insertIndex(index);
        }
    }

    /// counts the number of row indices in the IndexList in
    inline DataTypes::dim_t count(DataTypes::index_t range_min,
                                  DataTypes::index_t range_max) const
    {
        DataTypes::dim_t out=0;
        for (DataTypes::dim_t i=0; i < n; i++) {
            if (m_list[i] >= range_min && range_max > m_list[i])
                ++out;
        }
        if (extension)
            out += extension->count(range_min, range_max);
        return out;
    }

    /// index list to array
    inline void toArray(DataTypes::index_t* array,
                    DataTypes::index_t range_min, DataTypes::index_t range_max,
                    DataTypes::index_t index_offset) const
    {
        DataTypes::index_t idx = 0;
        for (DataTypes::dim_t i=0; i < n; i++) {
            if (m_list[i] >= range_min && range_max > m_list[i]) {
                array[idx] = m_list[i]+index_offset;
                ++idx;
            }
        }
        if (extension)
            extension->toArray(&array[idx], range_min, range_max, index_offset);
    }
};

} // namespace escript

#endif // __ESCRIPT_INDEXLIST_H__

