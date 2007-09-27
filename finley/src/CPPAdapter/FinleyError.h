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
                                                                           
#if !defined  finley_FinleyError_20040528_H
#define finley_FinleyError_20040528_H
#include "system_dep.h"

extern "C" {
#include "../Finley.h"
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
  /**
     \brief
     Convert a C paso  error into a C++ exception.
  */
  FINLEY_DLL_API
  void checkPasoError();
} // end of namespace
#endif
