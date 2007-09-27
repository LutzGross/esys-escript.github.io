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
                                                                           
#if !defined  finley_FinleyAdapterException_20040526_H
#define finley_FinleyAdapterException_20040526_H
#include "system_dep.h"

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
  FINLEY_DLL_API
  FinleyAdapterException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  FINLEY_DLL_API
  FinleyAdapterException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  FINLEY_DLL_API
  FinleyAdapterException(const std::string &str) : EsysException(str) {}

  /// Destructor
  FINLEY_DLL_API
  virtual ~FinleyAdapterException() throw() {}
  /**
     \brief
     Returns the name of the exception.
  */
  FINLEY_DLL_API
  virtual std::string exceptionName() const {return "FinleyAdapterException";}
};

} // end of namespace
#endif
