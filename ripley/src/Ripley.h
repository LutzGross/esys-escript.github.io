
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __RIPLEY_RIPLEY_H__
#define __RIPLEY_RIPLEY_H__

/*****************************************************************************
 *  Ripley is a FE domain library with regular hexagonal/rectangular
 *  elements
 ****************************************************************************/

#include <ripley/system_dep.h>

#include <escript/EsysMPI.h>

#include <boost/shared_ptr.hpp>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace ripley {

using escript::DataTypes::dim_t;
using escript::DataTypes::index_t;
using escript::DataTypes::cplx_t;
using escript::DataTypes::real_t;

typedef std::pair<index_t,index_t> IndexPair;
typedef std::vector<index_t> IndexVector;
typedef std::vector<real_t> DoubleVector;
typedef std::vector<int> RankVector;
typedef std::map<std::string,int> TagMap;

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

