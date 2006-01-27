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
                                                                           
#if !defined  finley_FinleyError_20040528_H
#define finley_FinleyError_20040528_H

extern "C" {
#include "Finley.h"
}

#include "FinleyAdapterException.h"

#include <string>

namespace finley {
  /**
     \brief
     Provide a C++ interface to the finley C funcion of the same name.
     Needed because of constness problems.
  */
  void setFinleyError(Finley_ErrorCodeType errorCode, 
		      const std::string& errMess);
 
  /**
     \brief
     Convert a C finley error into a C++ exception.
  */
  void checkFinleyError();
  /**
     \brief
     Convert a C paso  error into a C++ exception.
  */
  void checkPasoError();
} // end of namespace
#endif
