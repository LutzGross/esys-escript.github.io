
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


#if !defined  escript_DomainException_20040608_H
#define escript_DomainException_20040608_H
#include "system_dep.h"

#include "EsysException.h"

namespace escript {

class DomainException : public EsysException
{
public:
    DomainException(const std::string& str) : EsysException(str) {}
    virtual ~DomainException() throw() {}
};

} // end of namespace

#endif

