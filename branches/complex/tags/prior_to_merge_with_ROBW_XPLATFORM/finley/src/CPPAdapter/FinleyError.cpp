/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "FinleyError.h"

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
