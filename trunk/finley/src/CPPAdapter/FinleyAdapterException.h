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

#include "EsysException.h"

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
