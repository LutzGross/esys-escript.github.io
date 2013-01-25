
/* $Id: DataException.cpp 1312 2007-09-24 06:18:44Z ksteube $ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "DataException.h"

using namespace escript;

const std::string 
DataException::exceptionNameValue("DataException");

const std::string &
DataException::exceptionName() const
{
  return exceptionNameValue;
}

