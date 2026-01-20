
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined  escript_FunctionSpaceException_20040602_H
#define escript_FunctionSpaceException_20040602_H
#include "system_dep.h"

#include "EsysException.h"

namespace escript
{

class FunctionSpaceException : public EsysException
{
public:
    FunctionSpaceException(const std::string& str) : EsysException(str) {}
    virtual ~FunctionSpaceException() throw() {}
};

} // end of namespace

#endif

