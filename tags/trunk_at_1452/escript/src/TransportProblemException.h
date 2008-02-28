
/* $Id:$ */

/*******************************************************
 *
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined  escript_TransportProblemException_20040608_H
#define escript_TransportProblemException_20040608_H

#include "system_dep.h"
#include "esysUtils/EsysException.h"

#include <string>

namespace escript {

/**
   \brief
   TransportProblemException exception class.

   Description:
   TransportProblemException exception class.
   The class provides a public function returning the exception name
*/
class TransportProblemException : public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  ESCRIPT_DLL_API
  TransportProblemException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  ESCRIPT_DLL_API
  TransportProblemException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  ESCRIPT_DLL_API
  TransportProblemException(const std::string &str) : EsysException(str) {}

  /// Destructor
  ESCRIPT_DLL_API
  virtual ~TransportProblemException() throw() {}
  /**
     \brief
     Returns the name of the exception.
  */
  ESCRIPT_DLL_API
  virtual std::string exceptionName() const {return "TransportProblemException";}
};

} // end of namespace
#endif
