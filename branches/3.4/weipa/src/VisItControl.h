
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef __WEIPA_VISITCONTROL_H__
#define __WEIPA_VISITCONTROL_H__

#include <weipa/weipa.h>
#include <string>

namespace weipa {

namespace VisItControl {

    WEIPA_DLL_API
    bool initialize(const std::string& simFile, const std::string& comment);

    WEIPA_DLL_API
    bool publishData(EscriptDataset_ptr dataset);

} // namespace VisItControl

} // namespace weipa

#endif // __WEIPA_VISITCONTROL_H__

