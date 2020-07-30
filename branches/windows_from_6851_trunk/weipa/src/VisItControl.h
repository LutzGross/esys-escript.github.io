
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

