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
#include "system_dep.h"
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
  BRUCE_DLL_API
  BruceException() : EsysException() {}

  /**
     \brief
     Constructor for the exception.
  */
  BRUCE_DLL_API
  BruceException(const char *cstr) : EsysException(cstr) {}

  /**
     \brief
     Constructor for the exception.
  */
  BRUCE_DLL_API
  BruceException(const std::string &str) : EsysException(str) {}

  /**
     \brief
     Returns the name of the exception.
  */
  BRUCE_DLL_API
  virtual
  std::string
  exceptionName() const {return "BruceException";}

};

} // end of namespace

#endif
