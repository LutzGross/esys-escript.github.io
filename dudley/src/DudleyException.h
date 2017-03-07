
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

#ifndef __DUDLEY_EXCEPTION_H__
#define __DUDLEY_EXCEPTION_H__

#include <escript/EsysException.h>

namespace dudley {

class DudleyException : public escript::EsysException
{
public:
    DudleyException(const std::string& str) : escript::EsysException(str) {}
};

} // end of namespace

#endif // __DUDLEY_EXCEPTION_H__
