/*****************************************************************************
*
* Copyright (c)2015-2017 by The University of Queensland
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

#ifndef __ESCRIPT_PYERR_H__
#define __ESCRIPT_PYERR_H__

#include <boost/python/errors.hpp>

#include <string>

namespace escript {

void getStringFromPyException(boost::python::error_already_set e, std::string& errormsg);

} // namespace escript

#endif // __ESCRIPT_PYERR_H__

