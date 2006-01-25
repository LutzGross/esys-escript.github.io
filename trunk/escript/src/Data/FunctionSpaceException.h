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
                                                                           
#if !defined  escript_FunctionSpaceException_20040602_H
#define escript_FunctionSpaceException_20040602_H

#include "EsysException.h"
#include <string>

namespace escript {

/**
   \brief
   FunctionSpaceException exception class.

   Description:
   FunctionSpaceException exception class.
   The class provides a public function returning the exception name
*/
class FunctionSpaceException:public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  FunctionSpaceException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  FunctionSpaceException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  FunctionSpaceException(const std::string &str) : EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "FunctionSpaceException";}
};

} // end of namespace
#endif
