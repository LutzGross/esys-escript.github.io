
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#include "EsysAssertException.h"
#include <sstream>


using namespace esysUtils;

const std::string 
EsysAssertException::exceptionNameValue("EsysAssertException");


const std::string &
EsysAssertException::exceptionName() const
{
  return exceptionNameValue;
}



void 
EsysAssertException::assertFailure (const std::string& assertion,
                           const std::string& date, const std::string& file,
                           int line, const std::string& errDesc)
{
  std::stringstream message;
 
  message << std::endl
          << "EsysAssert(" << assertion << ") failed with message - " 
          << std::endl
          << "\"" << errDesc << "\"" << std::endl
          << "Assertion is located in File : " << file
          << " at Line: " << line << std::endl
          << "File Compilation Date: " << date << std::endl;
 
  throw EsysAssertException(message.str());
}
