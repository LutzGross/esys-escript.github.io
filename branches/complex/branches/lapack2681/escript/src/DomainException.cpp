
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DomainException.h"

using namespace escript;

const std::string 
DomainException::exceptionNameValue("DomainException");


const std::string &
DomainException::exceptionName() const
{
  return exceptionNameValue;
}


