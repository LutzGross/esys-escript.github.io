
/*******************************************************
*
* Copyright (c) 2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "RDomainException.h"


using namespace refine;


const std::string 
RDomainException::exceptionNameValue("RDomainException");


const std::string &
RDomainException::exceptionName() const
{
  return exceptionNameValue;
}

