
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DudleyAdapterException.h"


using namespace dudley;


const std::string 
DudleyAdapterException::exceptionNameValue("DudleyAdapterException");


const std::string &
DudleyAdapterException::exceptionName() const
{
  return exceptionNameValue;
}


