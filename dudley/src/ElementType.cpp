
/*****************************************************************************
*
* Copyright (c) 2010-2014 by University of Queensland
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

#include <string.h>
#include "ElementType.h"

Dudley_ElementTypeId eltTypeFromString(const char *s)
{
    if (strcmp(s, "Point1") == 0)
    {
	return Dudley_Point1;
    }
    else if (strcmp(s, "Line2") == 0)
    {
	return Dudley_Line2;
    }
    else if (strcmp(s, "Tri3") == 0)
    {
	return Dudley_Tri3;
    }
    else if (strcmp(s, "Tet4") == 0)
    {
	return Dudley_Tet4;
    }
    else if (strcmp(s, "Line2Face") == 0)
    {
	return Dudley_Line2Face;
    }
    else if (strcmp(s, "Tri3Face") == 0)
    {
	return Dudley_Tri3Face;
    }
    else if (strcmp(s, "Tet4Face") == 0)
    {
	return Dudley_Tet4Face;
    }
    else
	return Dudley_NoRef;
}
