/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#ifndef __OXLEY_H__
#define __OXLEY_H__

#include <list>
#include <map>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <p4est/p4est_base.h>

#include <escript/EsysMPI.h>

// #include <oxley/Oxley.h>

#define MAXP4ESTNODES 1024*1024 // Maximum allowed nodes in the p4est / p8est
#define MAXTREES 1024*1024 // Maximum allowed trees
#define MAXREFINEMENTLEVELS 16 // Default levels of refinement
#define MAXTAGS 100 // Maximum allowed number of tags in a domain

#define MARE2DEM_TOL 0.25 
#define TOLERANCE 1e-8

#ifdef OXLEY_ENABLE_DEBUG
#define LOG_BACKTRACE 1  // Print a backtrace if p4est aborts prematurely
#define LOG_LEVEL 0  // Everything
// #define LOG_LEVEL 2  // Debug info
// #define LOG_LEVEL 4 // Main things that each function does
// #define LOG_LEVEL 8 // Errors only
#else
// #define LOG_BACKTRACE 0
#define LOG_LEVEL 9  // Nothing
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

typedef std::pair<long,long> LongPair;
typedef std::pair<double,double> DoublePair;
typedef std::tuple<double,double,double> DoubleTuple;
// struct DoubleTuple
// {
//     double x, y, z;

//     DoubleTuple()
//     {
//         x=-1.;y=-1.;z=-1.;
//     }
//     ~DoubleTuple();

//     double get(int n)
//     {
//         if(n == 0)
//             return x;
//         else if(n == 1)
//             return y;
//         else if(n == 2)
//             return z;
//         else
//             return NAN;
//     }
// };

// bool operator==(DoubleTuple& A, DoubleTuple& B)
// {
//     if(    (std::abs(A.get(1)-B.get(1)) < TOLERANCE)
//         && (std::abs(A.get(2)-B.get(2)) < TOLERANCE)
//         && (std::abs(A.get(3)-B.get(3)) < TOLERANCE))
//         return true;
//     else
//         return false;
// }

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

typedef p4est_topidx_t p8est_topidx_t;
typedef p4est_qcoord_t p8est_qcoord_t;
typedef p4est_locidx_t p8est_locidx_t;

// fs types
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

} // namespace oxley

#endif
