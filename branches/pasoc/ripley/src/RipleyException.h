
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

