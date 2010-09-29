
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

#ifndef DUDLEY_TYPEID_H
#define DUDLEY_TYPEID_H

typedef enum {
  Point1=0,
  Line2=1,
  Tri3=2,
  Tet4=3,
  Line2Face=4,
  Tri3Face=5,
  Tet4Face=6, 
  NoRef=7   /* marks end of list */
} ElementTypeId;


ElementTypeId eltTypeFromString(const char* s);
#endif