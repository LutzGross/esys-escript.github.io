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

#include "finley/CPPAdapter/FinleyError.h"
#include "finley/CPPAdapter/FinleyAdapterException.h"
extern "C" {
#include "finley/finleyC/setFinleyError.h"
}

namespace finley {

  void setFinleyError(enum Finley_ErrorCodeType errorCode, 
		      const std::string& errMess) 
  {
    char message[errMess.size()+1];
    strcpy(message,errMess.c_str());
    setFinleyError(errorCode,message,errMess.size());
  }

  void checkFinleyError() 
  {
    if (Finley_ErrorCode==NO_ERROR) {
      return;
    } else {
      //
      // reset the error code to no error otherwise the next call to
      // this function may resurrect a previous error
      Finley_ErrorCode=NO_ERROR;
      throw FinleyAdapterException(Finley_ErrorMsg);
    }
  }

}  // end of namespace
