
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  dudley_DudleyError_20040528_H
#define dudley_DudleyError_20040528_H
#include "system_dep.h"

extern "C" {
#include "dudley/Dudley.h"
}

#include "DudleyAdapterException.h"

#include <string>

namespace dudley {
  /**
     \brief
     Provide a C++ interface to the dudley C funcion of the same name.
     Needed because of constness problems.
  */
  DUDLEY_DLL_API
  void setDudleyError(Dudley_ErrorCodeType errorCode, 
		      const std::string& errMess);
 
  /**
     \brief
     Convert a C dudley error into a C++ exception.
  */
  DUDLEY_DLL_API
  void checkDudleyError();
  /**
     \brief
     Convert a C paso  error into a C++ exception.
  */
  DUDLEY_DLL_API
  void checkPasoError();
} // end of namespace
#endif
