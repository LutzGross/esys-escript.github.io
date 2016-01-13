// $Id$
/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
*/
                                                                           
#if !defined escript_EsysAssertException_20040330_H
#define escript_EsysAssertException_20040330_H

#include "EsysException.h"

#include <string>
#include <sstream>

namespace esysUtils {

/**
   \brief
   EsysAssertException exception class.

   Description:
   EsysAssertException exception class.
   The class provides a public function returning the exception name.
*/
class EsysAssertException:public EsysException {

 public:

  /**
     \brief
     Default constructor for the exception.
  */
  EsysAssertException() : EsysException() {}

  /**
     \brief
     Constructor for the exception.
  */
  EsysAssertException(const char *cstr) : EsysException(cstr) {}

  /**
     \brief
     Constructor for the exception.
  */
  EsysAssertException(const std::string &str) : EsysException(str) {}

  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "EsysAssertException";}

  /**
     \brief
     Builds a formatted message and throws an EsysAssertException.
  */
  static void assertFailure (const std::string& assertion,
			     const std::string& date, const std::string& file,
			     int line, const std::string& errDesc)
  {
   std::stringstream message;
 
   message << std::endl
           << "EsysAssert(" << assertion << ") failed with message - " << std::endl
           << "\"" << errDesc << "\"" << std::endl
           << "Assertion is located in File : " << file
           << " at Line: " << line << std::endl
           << "File Compilation Date: " << date << std::endl;
 
   throw EsysAssertException(message.str());
  }

};

} // end of namespace
 
#endif
