
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __FILESAVERS_H__
#define __FILESAVERS_H__

#include <weipa/weipa.h>
#include <escript/AbstractDomain.h>

#include <boost/python/dict.hpp>
#include <string>

namespace weipa {

WEIPA_DLL_API
void saveSilo(const std::string& filename, int cycle, double time,
              escript::Domain_ptr domain,
              const boost::python::dict& datavars);

WEIPA_DLL_API
void saveVTK(const std::string& filename, int cycle, double time,
             escript::Domain_ptr domain,
             const boost::python::dict& datavars, const std::string& metadata,
             const std::string& metadata_schema);


} // namespace weipa

#endif // __FILESAVERS_H__

