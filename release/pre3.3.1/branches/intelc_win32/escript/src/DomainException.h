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

#if !defined  escript_DomainException_20040608_H
#define escript_DomainException_20040608_H
#include "system_dep.h"

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
  ESCRIPT_DLL_API
  DomainException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  ESCRIPT_DLL_API
  DomainException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  ESCRIPT_DLL_API
  DomainException(const std::string &str) : EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  ESCRIPT_DLL_API
  virtual std::string exceptionName() const {return "DomainException";}
};

} // end of namespace
#endif
