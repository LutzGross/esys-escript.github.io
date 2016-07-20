
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
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

#ifndef __ESYS_TRILINOSWRAP_UTIL_H__
#define __ESYS_TRILINOSWRAP_UTIL_H__

#include <escript/EsysException.h>

#include <Teuchos_ParameterList.hpp>

#include <boost/python.hpp>
#include <boost/python/dict.hpp>

#include <string>

namespace esys_trilinos {

namespace util {

template<typename T>
void extractParamIfSet(const std::string& name,
                       const boost::python::dict& pyDict,
                       Teuchos::ParameterList& params)
{
    if (pyDict.has_key(name)) {
        boost::python::object bpo = pyDict.get(name);
        if (boost::python::extract<T>(bpo).check()) {
            T val = boost::python::extract<T>(bpo);
            params.set(name, val);
        } else {
            throw escript::ValueError("Wrong type for option " + name);
        }
    }
}

} // namespace util
} // namespace esys_trilinos

#endif // __ESYS_TRILINOSWRAP_UTIL_H__

