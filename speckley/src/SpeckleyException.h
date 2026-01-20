
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

#ifndef __SPECKLEY_EXCEPTION_H__
#define __SPECKLEY_EXCEPTION_H__

#include <speckley/system_dep.h>
#include <escript/EsysException.h>

namespace speckley {

/**
   \brief
   SpeckleyException exception class.
*/
class SpeckleyException : public escript::EsysException
{
public:
    SpeckleyException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace speckley

#endif // __SPECKLEY_EXCEPTION_H__

