
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

/*   esysUtils: IndexList                                                   */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __ESYSUTILS_INDEXLIST_H__
#define __ESYSUTILS_INDEXLIST_H__

#include "types.h"

#include <algorithm>
#include <list>
#include <vector>

namespace esysUtils {

struct IndexList;
typedef std::vector<IndexList> IndexListArray;

struct IndexList {
    std::list<index_t> m_list;

    /// inserts row index into the IndexList in if it does not exist
    inline void insertIndex(index_t index)
    {
        if (std::find(m_list.begin(), m_list.end(), index) == m_list.end())
            m_list.push_back(index);
    }

    /// counts the number of row indices in the IndexList in
    inline dim_t count(index_t range_min, index_t range_max) const
    {
        dim_t out=0;
        std::list<index_t>::const_iterator it;
        for (it=m_list.begin(); it != m_list.end(); it++) {
            if (*it >= range_min && range_max > *it)
                ++out;
        }
        return out;
    }

    /// index list to array
    inline void toArray(index_t* array, index_t range_min, index_t range_max,
                        index_t index_offset) const
    {
        index_t idx = 0;
        std::list<index_t>::const_iterator it;
        for (it=m_list.begin(); it != m_list.end(); it++) {
            if (*it >= range_min && range_max > *it) {
                array[idx] = (*it)+index_offset;
                ++idx;
            }
        }
    }
};

} // namespace esysUtils

#endif // __ESYSUTILS_INDEXLIST_H__

