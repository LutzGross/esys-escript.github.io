
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

ostream &operator<<(ostream &output, EsysException &inException){
  output << inException.toString();
  return output;
}

EsysException::EsysException():
exception(),
m_reason()
{
  updateMessage();   
}

EsysException::EsysException(const std::string &exceptionReason):
exception(),
m_reason(exceptionReason)
{
  updateMessage();   
}

// Copy Constructor.
EsysException::EsysException(const EsysException &inException):
exception(inException)
{
  *this=inException;
}

EsysException::EsysException( const char *cStr ):
exception(),
m_reason(cStr) 
{
  updateMessage();   
}

EsysException::~EsysException()
{}

const std::string & EsysException::exceptionName() const 
{
  return exceptionNameValue;
}

//
// Overloaded assignment operator.
EsysException& EsysException::operator=(const EsysException &inException)
{
  if (this != &inException)
  {
    //
    // call the base class operator=
    // win32 refactor: parent class operator= shares pointer the result is
    // all classes try to free the same pointer, dies on windows badly.
    // exception::operator=(inException);
    //
    // copy the message buffer into this EsysException
    m_reason = inException.m_reason;
    updateMessage();
  }
  return *this;
}
