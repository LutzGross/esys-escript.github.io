
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


#include "FunctionSpaceException.h"


using namespace escript;


const std::string 
FunctionSpaceException::exceptionNameValue("FunctionSpaceException");


const std::string &
FunctionSpaceException::exceptionName() const
{
  return exceptionNameValue;
}


