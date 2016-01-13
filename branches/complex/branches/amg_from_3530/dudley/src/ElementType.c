
/*******************************************************
*
* Copyright (c) 2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

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
