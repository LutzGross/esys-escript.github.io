/*****************************************************************************
*
* Copyright (c) 2015 by University of Queensland
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

// Function factored out of SubWorld code
#ifndef ESPYERR_H
#define ESPYERR_H

#include <string>
#include "system_dep.h"
#include "types.h"
#include "boost/python/errors.hpp"

ESYSUTILS_DLL_API
void getStringFromPyException(boost::python::error_already_set e, std::string& errormsg);

#endif