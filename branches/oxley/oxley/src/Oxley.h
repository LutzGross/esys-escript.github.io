/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
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

#include <escript/EsysMPI.h>

#include <boost/shared_ptr.hpp>
#include <list>
#include <map>
#include <string>
#include <vector>

#define MAXP4ESTNODES 10000 // Maximum allowed nodes in the p4est / p8est
#define MAXREFINEMENTLEVELS 5 // Default levels of refinement
#define MAXTAGS 100 // Maximum allowed number of tags in a domain

#ifdef P4EST_ENABLE_DEBUG
#define LOG_BACKTRACE 1  // Print a backtrace if p4est aborts prematurely
#define LOG_LEVEL 0  // Everything
// #define LOG_LEVEL 4 // Main things that each functi0n does
// #define LOG_LEVEL 8 // Errors only
#else
#define LOG_BACKTRACE 0
#define LOG_LEVEL 9  // Level of logging used by p4est
#endif
// Acceptable values are
/* log priorities */
// #define SC_LP_DEFAULT   (-1)    /**< this selects the SC default threshold */
// #define SC_LP_ALWAYS      0     /**< this will log everything */
// #define SC_LP_TRACE       1     /**< this will prefix file and line number */
// #define SC_LP_DEBUG       2     /**< any information on the internal state */
// #define SC_LP_VERBOSE     3     /**< information on conditions, decisions */
// #define SC_LP_INFO        4     /**< the main things a function is doing */
// #define SC_LP_STATISTICS  5     /**< important for consistency/performance */
// #define SC_LP_PRODUCTION  6     /**< a few lines for a major api function */
// #define SC_LP_ESSENTIAL   7     /**< this logs a few lines max per program */
// #define SC_LP_ERROR       8     /**< this logs errors only */
// #define SC_LP_SILENT      9     /**< this never logs anything */
/*@}*/


#ifndef OXLEYH
#define OXLEYH

namespace oxley {

using escript::DataTypes::dim_t;
using escript::DataTypes::index_t;
using escript::DataTypes::cplx_t;
using escript::DataTypes::real_t;

typedef std::pair<index_t,index_t> IndexPair;
typedef std::vector<index_t> IndexVector;
typedef std::vector<real_t> DoubleVector;
typedef std::vector<int> RankVector;
typedef std::map<std::string,int> TagMap;

// fs types
// enum {
//     DegreesOfFreedom=1,
//     ReducedDegreesOfFreedom=2,
//     Nodes=3,
//     ReducedNodes=14,
//     Elements=4,
//     ReducedElements=10,
//     FaceElements=5,
//     ReducedFaceElements=11,
//     Points=6
// };

} // namespace oxley

#endif
