/*****************************************************************************
*
* Copyright (c) 2015-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_PYERR_H__
#define __ESCRIPT_PYERR_H__

#include <boost/python/errors.hpp>

#include <string>

namespace escript {

void getStringFromPyException(boost::python::error_already_set e, std::string& errormsg);

} // namespace escript

#endif // __ESCRIPT_PYERR_H__

