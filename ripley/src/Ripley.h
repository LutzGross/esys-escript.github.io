
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

#ifndef __RIPLEY_RIPLEY_H__
#define __RIPLEY_RIPLEY_H__

/************************************************************************************ 
 *  Ripley is a FE domain library with hexagonal/rectangular
 *  elements
 ************************************************************************************/

#include <ripley/system_dep.h>

extern "C" {
#include <esysUtils/Esys_MPI.h>
}

#include <boost/shared_ptr.hpp>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace ripley {

typedef std::pair<index_t,index_t> IndexPair;
typedef std::vector<index_t> IndexVector;
typedef std::vector<Esys_MPI_rank> RankVector;
typedef std::map<std::string,index_t> TagMap;

enum {
    DegreesOfFreedom=1,
    ReducedDegreesOfFreedom=2,
    Nodes=3,
    ReducedNodes=14,
    Elements=4,
    ReducedElements=10,
    FaceElements=5,
    ReducedFaceElements=11,
    Points=6
};

} // namespace ripley

#endif /* __RIPLEY_RIPLEY_H__ */

