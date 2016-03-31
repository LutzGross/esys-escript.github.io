/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include <speckley/AbstractAssembler.h>

namespace speckley {

escript::Data unpackData(std::string target,
        std::map<std::string, escript::Data> mapping){
    if (mapping.find(target) == mapping.end())
        return escript::Data();
    return mapping[target];
}

}
