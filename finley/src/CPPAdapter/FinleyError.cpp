
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifdef PASO_MPI
#include <mpi.h>
#endif
#include "FinleyError.h"
#include <iostream>
#include <string>

namespace finley {

  void setFinleyError(Finley_ErrorCodeType errorCode, 
		      const std::string& errMess) 
  {
    Finley_setError(errorCode,(char*)(errMess.c_str()));
  }

  void checkFinleyError() 
  {
    if (Finley_noError()) {
      return;
    } else {
      //
      // reset the error code to no error otherwise the next call to
      // this function may resurrect a previous error
      Finley_resetError();
      throw FinleyAdapterException(Finley_getErrorMessage());
    }
  }
  void checkPasoError() 
  {
    if (Paso_noError()) {
      return;
    } else {
      //
      // reset the error code to no error otherwise the next call to
      // this function may resurrect a previous error
      Paso_resetError();
      throw FinleyAdapterException(Paso_getErrorMessage());
    }
  }

}  // end of namespace
