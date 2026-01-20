
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


#ifndef __FINLEY_EXCEPTION_H__
#define __FINLEY_EXCEPTION_H__

#include <escript/EsysException.h>

namespace finley {

class FinleyException : public escript::EsysException
{
public:
    FinleyException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace

#endif // __FINLEY_EXCEPTION_H__

