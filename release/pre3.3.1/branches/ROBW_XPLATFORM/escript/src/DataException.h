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

#if !defined  escript_DataException_20040324_H
#define escript_DataException_20040324_H

#include "esysUtils/EsysException.h"

#include <string>

namespace escript {

/**
   \brief
   DataException exception class.

   Description:
   DataException exception class.
   The class provides a public function returning the exception name
*/
class DataException:public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  DataException() : esysUtils::EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  DataException(const char *cstr) : esysUtils::EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  DataException(const std::string &str) : esysUtils::EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "DataException";}
};

} // end of namespace
#endif
