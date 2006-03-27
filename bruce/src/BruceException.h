/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
                                                                           
#if !defined finley_BruceException_20050905_H
#define finley_BruceException_20050905_H

#include "esysUtils/EsysException.h"

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
