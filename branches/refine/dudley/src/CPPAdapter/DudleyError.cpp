
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DudleyError.h"
#include <iostream>

namespace dudley {

  void setDudleyError(Dudley_ErrorCodeType errorCode, 
		      const std::string& errMess) 
  {
    Dudley_setError(errorCode,(__const char*)(errMess.c_str()));
  }

  void checkDudleyError() 
  {
    if (Dudley_noError()) {
      return;
    } else {
      //
      // reset the error code to no error otherwise the next call to
      // this function may resurrect a previous error
      Dudley_resetError();
      throw DudleyAdapterException(Dudley_getErrorMessage());
    }
  }
  void checkPasoError() 
  {
    if (Esys_noError()) {
      return;
    } else {
      //
      // reset the error code to no error otherwise the next call to
      // this function may resurrect a previous error
      Esys_resetError();
      throw DudleyAdapterException(Esys_getErrorMessage());
    }
  }

}  // end of namespace
