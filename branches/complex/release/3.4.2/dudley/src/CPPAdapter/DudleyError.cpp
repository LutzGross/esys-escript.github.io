
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


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
