
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#include "EsysException.h"

using namespace esysUtils;

const std::string EsysException::exceptionNameValue("GeneralEsysException");

std::ostream &operator<<(std::ostream &output, EsysException &inException){
  output << inException.toString();
  return output;
}

EsysException::EsysException():
Parent(),
m_reason()
{
  updateMessage();   
}

EsysException::EsysException(const std::string &exceptionReason):
Parent(),
m_reason(exceptionReason)
{
  updateMessage();   
}

// Copy Constructor.
// Do not call the parent copy constructor as it has
// undefined effects. In particular, it mat call what() to
// which will resuly on the parent storing a pointer to
// m_exceptionMessage's storage.... esp on winblows.
EsysException::EsysException(const EsysException &other):
Parent(),
m_reason(other.m_reason)
{
  updateMessage();   
}

EsysException &
EsysException::operator=(const EsysException &other) THROW(NO_ARG) 
{
  m_reason = other.m_reason;
  updateMessage();   
  return *this;
}

EsysException::EsysException( const char *cStr ):
Parent(),
m_reason(cStr) 
{
  updateMessage();   
}

EsysException::~EsysException() THROW(NO_ARG)
{}

const std::string & EsysException::exceptionName() const 
{
  return exceptionNameValue;
}

