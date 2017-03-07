
/*****************************************************************************
*
* Copyright (c) 2010-2017 by The University of Queensland
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

#ifndef __DUDLEY_ELEMENTTYPE_H__
#define __DUDLEY_ELEMENTTYPE_H__

#include <string>

namespace dudley {

typedef enum {
    Dudley_Point1 = 0,
    Dudley_Line2 = 1,
    Dudley_Tri3 = 2,
    Dudley_Tet4 = 3,
    Dudley_Line2Face = 4,
    Dudley_Tri3Face = 5,
    Dudley_Tet4Face = 6,
    Dudley_NoRef = 7      // marks end of list
} ElementTypeId;

inline ElementTypeId eltTypeFromString(const std::string& s)
{
    if (s == "Point1")
        return Dudley_Point1;
    else if (s == "Line2")
        return Dudley_Line2;
    else if (s == "Tri3")
        return Dudley_Tri3;
    else if (s == "Tet4")
        return Dudley_Tet4;
    else if (s == "Line2Face")
        return Dudley_Line2Face;
    else if (s == "Tri3Face")
        return Dudley_Tri3Face;
    else if (s == "Tet4Face")
        return Dudley_Tet4Face;
    else
        return Dudley_NoRef;
}

}

#endif // __DUDLEY_ELEMENTTYPE_H__

