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

#if !defined  escript_SystemMatrixException_20040608_H
#define escript_SystemMatrixException_20040608_H

#include "system_dep.h"
#include "esysUtils/EsysException.h"

#include <string>

namespace escript {

/**
   \brief
   SystemMatrixException exception class.

   Description:
   SystemMatrixException exception class.
   The class provides a public function returning the exception name
*/
class SystemMatrixException : public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  ESCRIPT_DLL_API
  SystemMatrixException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  ESCRIPT_DLL_API
  SystemMatrixException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  ESCRIPT_DLL_API
  SystemMatrixException(const std::string &str) : EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  ESCRIPT_DLL_API
  virtual std::string exceptionName() const {return "SystemMatrixException";}
};

} // end of namespace
#endif
