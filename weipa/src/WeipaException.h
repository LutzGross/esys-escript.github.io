
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

#ifndef __WEIPA_EXCEPTION_H__
#define __WEIPA_EXCEPTION_H__

#include <escript/EsysException.h>

namespace weipa {
/**
   \brief
   WeipaException exception class.
*/
class WeipaException : public escript::EsysException
{
public:
    WeipaException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace weipa

#endif // __WEIPA_EXCEPTION_H__

