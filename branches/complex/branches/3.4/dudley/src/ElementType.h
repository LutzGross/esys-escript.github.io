
/*****************************************************************************
*
* Copyright (c) 2010-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef DUDLEY_TYPEID_H
#define DUDLEY_TYPEID_H

typedef enum {
    Dudley_Point1 = 0,
    Dudley_Line2 = 1,
    Dudley_Tri3 = 2,
    Dudley_Tet4 = 3,
    Dudley_Line2Face = 4,
    Dudley_Tri3Face = 5,
    Dudley_Tet4Face = 6,
    Dudley_NoRef = 7			/* marks end of list */
} Dudley_ElementTypeId;

Dudley_ElementTypeId eltTypeFromString(const char *s);
#endif
