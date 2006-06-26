/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "EsysException.h"

using namespace std;
using namespace esysUtils;

ostream &operator<<(ostream &output, EsysException &inException){
  output << inException.toString();
  return output;
}

EsysException::EsysException():
  exception() 
{
}

EsysException::EsysException(const string &exceptionReason):
exception()
{
  reason() << exceptionReason;
}

// Copy Constructor.
EsysException::EsysException(const EsysException &inException):
exception(inException)
{
  *this=inException;
}

EsysException::EsysException( const char *cStr ):
  exception() 
{
  reason() << cStr;
}

EsysException::~EsysException() throw()
{}

string EsysException::exceptionName() const 
{
  return "GeneralEsysException";
}
ostringstream& EsysException::reason() 
{
  return m_reason;
}
//
// Overloaded assignment operator.
EsysException& EsysException::operator=(const EsysException &inException) {
  if (this != &inException) {
	  //
	  // call the base class operator=
	  // win32 refactor: parent class operator= shares pointer the result is
	  // all classes try to free the same pointer, dies on windows badly.
	  // this->exception::operator=(dynamic_cast<const exception&>(inException));
	  //
	  // copy the message buffer into this EsysException
	  m_reason << inException.m_reason.str();    // copy the message buffer into this EsysException
  }
  return *this;
}

// return the message as a string
string EsysException::toString() const {
  return exceptionName() + ": " + m_reason.str();
}

const char*  EsysException::what() const throw() {
//
  m_exceptionMessage=toString();
  return m_exceptionMessage.c_str();
}
