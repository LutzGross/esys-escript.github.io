
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

#ifndef __RIPLEY_EXCEPTION_H__
#define __RIPLEY_EXCEPTION_H__

#include <ripley/system_dep.h>
#include <escript/EsysException.h>

namespace ripley {

/**
   \brief
   RipleyException exception class.
*/
class RipleyException : public escript::EsysException
{
public:
    RipleyException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace ripley

#endif // __RIPLEY_EXCEPTION_H__

