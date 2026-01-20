
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

#ifndef __ESYS_TRILINOSADAPTEREXCEPTION_H__
#define __ESYS_TRILINOSADAPTEREXCEPTION_H__

#include <escript/EsysException.h>

namespace esys_trilinos {

class TrilinosAdapterException : public escript::EsysException
{
public:
    TrilinosAdapterException(const std::string& str)
        : escript::EsysException(str) {}
};


} // namespace esys_trilinos

#endif // __ESYS_TRILINOSADAPTEREXCEPTION_H__

