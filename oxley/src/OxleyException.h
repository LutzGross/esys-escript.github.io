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
*
*****************************************************************************/

#ifndef __OXLEY_EXCEPTION_H__
#define __OXLEY_EXCEPTION_H__

#include <escript/EsysException.h>

#include <exception>

namespace oxley {

/**
   \brief
   OxleyException exception class.
*/
class OxleyException : public escript::EsysException
{
public:
    OxleyException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace oxley


/**
  \brief
  An exception class that signals an invalid argument value
*/
class ValueError : public escript::EsysException
{
public:
    ValueError(const std::string& str) : EsysException(str) {}
};

#endif
