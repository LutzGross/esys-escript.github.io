
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

#ifndef __VISITCONTROL_H__
#define __VISITCONTROL_H__

#include <weipa/weipa.h>
#include <escript/AbstractDomain.h>

#include <boost/python/dict.hpp>
#include <string>

namespace weipa {

namespace VisItControl {

    bool initialize(const std::string& simFile, const std::string& comment);

    void publishData(EscriptDataset_ptr dataset);

} // namespace VisItControl

} // namespace weipa

#endif // __VISITCONTROL_H__

