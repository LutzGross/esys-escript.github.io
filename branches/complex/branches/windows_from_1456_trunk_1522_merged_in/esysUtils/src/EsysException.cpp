
/* $Id$ */

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
EsysException::operator=(const EsysException &other) THROW_ANY 
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

EsysException::~EsysException() THROW_ANY
{}

const std::string & EsysException::exceptionName() const 
{
  return exceptionNameValue;
}

