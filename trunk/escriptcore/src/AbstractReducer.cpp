/*****************************************************************************
*
* Copyright (c) 2014-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include <sstream>
#include <limits>
#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include "AbstractReducer.h"
#include "SplitWorldException.h"

using namespace boost::python;
using namespace escript;


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
}
