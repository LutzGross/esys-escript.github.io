
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __PASO_EXCEPTION_H__
#define __PASO_EXCEPTION_H__

#include <escript/EsysException.h>

namespace paso {

/**
  \brief
  PasoException exception class.

  Description:
  PasoException exception class.
  The class provides a public function returning the exception name
*/
class PasoException : public escript::EsysException
{
public:
    PasoException(const std::string& str) : EsysException(str) {}
    virtual ~PasoException() throw() {}
};

void checkPasoError(); 


} // end of namespace

#endif // __PASO_EXCEPTION_H__

