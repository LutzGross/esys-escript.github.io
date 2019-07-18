
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
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

