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
                                                                           
#if !defined  escript_DomainException_20040608_H
#define escript_DomainException_20040608_H

#include "esysUtils/EsysException.h"
#include <string>

namespace escript {

/**
   \brief
   DomainException exception class.

   Description:
   DomainException exception class.
   The class provides a public function returning the exception name
*/
class DomainException:public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  DomainException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  DomainException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  DomainException(const std::string &str) : EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "DomainException";}
};

} // end of namespace
#endif
