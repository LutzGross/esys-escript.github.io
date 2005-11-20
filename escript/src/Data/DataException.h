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
                                                                           
#if !defined  escript_DataException_20040324_H
#define escript_DataException_20040324_H
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
   DataException exception class.

   Description:
   DataException exception class.
   The class provides a public function returning the exception name
*/
	class ESCRIPT_DLL DataException:public esysUtils::EsysException {

 public:
  /**
     \brief
     Default constructor for the exception.
  */
  DataException() : esysUtils::EsysException() {}
  /**
     \brief
     Constructor for the exception.
  */
  DataException(const char *cstr) : esysUtils::EsysException(cstr) {}
  /**
     \brief
     Constructor for the exception.
  */
  DataException(const std::string &str) : esysUtils::EsysException(str) {}
  /**
     \brief
     Returns the name of the exception.
  */
  virtual std::string exceptionName() const {return "DataException";}
};

} // end of namespace
#endif
