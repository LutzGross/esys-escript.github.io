/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
*/
                                                                           
#if !defined  escript_SystemMatrixException_20040608_H
#define escript_SystemMatrixException_20040608_H

#ifdef MSVC
#ifdef ESCRIPT_EXPORTS
#define ESCRIPT_DLL __declspec(dllexport)
#else
#define ESCRIPT_DLL __declspec(dllimport)
#endif
#else
#define ESCRIPT_DLL
#endif
#include "esysUtils/EsysException.h"
#include <string>

namespace escript {

/**
   \brief
   SystemMatrixException exception class.

   Description:
   SystemMatrixException exception class.
   The class provides a public function returning the exception name
*/
	class ESCRIPT_DLL SystemMatrixException:public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  SystemMatrixException() : EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  SystemMatrixException(const char *cstr) : EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  SystemMatrixException(const std::string &str) : EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "SystemMatrixException";}
};

} // end of namespace
#endif
