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
class FunctionSpaceException : public esysUtils::EsysException {

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
