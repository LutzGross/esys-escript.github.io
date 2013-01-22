
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#if !defined  finley_FinleyError_20040528_H
#define finley_FinleyError_20040528_H
#include "system_dep.h"

extern "C" {
#include "finley/Finley.h"
}

#include "FinleyAdapterException.h"

#include <string>

namespace finley {
  /**
     \brief
     Provide a C++ interface to the finley C funcion of the same name.
     Needed because of constness problems.
  */
  FINLEY_DLL_API
  void setFinleyError(Finley_ErrorCodeType errorCode, 
		      const std::string& errMess);
 
  /**
     \brief
     Convert a C finley error into a C++ exception.
  */
  FINLEY_DLL_API
  void checkFinleyError();

} // end of namespace
#endif
