/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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
