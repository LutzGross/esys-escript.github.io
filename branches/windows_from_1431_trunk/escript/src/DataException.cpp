
/* $Id: EsysDataException.cpp 1312 2007-09-24 06:18:44Z ksteube $ */

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

#include "EsysDataException.h"
#include <sstream>


using namespace esysUtils;

const std::string 
EsysDataException::exceptionNameValue("EsysDataException");


const std::string &
EsysDataException::exceptionName() const
{
  return exceptionNameValue;
}


