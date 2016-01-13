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
                                                                           
#if !defined  finley_FinleyAdapterException_20040526_H
#define finley_FinleyAdapterException_20040526_H

#include "esysUtils/EsysException.h"
#include <string>

namespace finley {

/**
   \brief
   FinleyAdapterException exception class.

   Description:
   FinleyAdapterException exception class.
   The class provides a public function returning the exception name
*/
class FinleyAdapterException:public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  FinleyAdapterException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  FinleyAdapterException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  FinleyAdapterException(const std::string &str) : EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "FinleyAdapterException";}
};

} // end of namespace
#endif
