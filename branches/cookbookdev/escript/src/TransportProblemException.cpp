
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


#include "TransportProblemException.h"


using namespace escript;


const std::string 
TransportProblemException::exceptionNameValue("TransportProblemException");


const std::string &
TransportProblemException::exceptionName() const
{
  return exceptionNameValue;
}


