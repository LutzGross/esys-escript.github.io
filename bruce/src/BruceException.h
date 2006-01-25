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
                                                                           
#if !defined finley_BruceException_20050905_H
#define finley_BruceException_20050905_H

#include "EsysException.h"

#include <string>

namespace bruce {

/**
   \brief
   Bruce exception class.

   Description:
   Bruce exception class.
   The class provides a public function returning the exception name.
*/

class BruceException : public esysUtils::EsysException {

 public:

  /**
     \brief
     Default constructor for the exception.
  */
  BruceException() : EsysException() {}

  /**
     \brief
     Constructor for the exception.
  */
  BruceException(const char *cstr) : EsysException(cstr) {}

  /**
     \brief
     Constructor for the exception.
  */
  BruceException(const std::string &str) : EsysException(str) {}

  /**
     \brief
     Returns the name of the exception.
  */
  virtual
  std::string
  exceptionName() const {return "BruceException";}

};

} // end of namespace

#endif
