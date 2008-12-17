
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "compareCCall.h"

/*
 * Return true if all of the calls match the expected values.
 * Return false if results don't match.
 */
int compareCCall(struct escriptDataC* data, int typeResult) 
{
  int result=0;
  result+=(getFunctionSpaceType(data)==typeResult);
  return result;
}
