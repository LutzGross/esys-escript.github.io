/*****************************************************************************
*
* Copyright (c) 2014-2017 by The University of Queensland
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

#include "AbstractReducer.h"
#include "SplitWorldException.h"

namespace escript {

const int AbstractReducer::PARAMTAG=120567;

bool AbstractReducer::hasValue()
{
    return valueadded;
}

double AbstractReducer::getDouble()
{
    throw SplitWorldException("This reducer is not able to provide a single scalar.");  
}

void AbstractReducer::clear()
{
    valueadded=false;
}

void AbstractReducer::newRunJobs()
{
    had_an_export_this_round=false;
}

bool AbstractReducer::canClash()
{
    return false;
}

} // namespace escript

