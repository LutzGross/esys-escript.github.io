
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/**
\file escript/src/system_dep.h
\ingroup Other
 */
/*
 @(#) system_dep.h
*/


#ifndef escript_system_dep_h
#define escript_system_dep_h


#ifdef NO_FLOAT_H
#   define DBL_EPSILON 2.2204460492503131E-16
#   define DBL_MAX 1.7976931348623157E+308
#   define DBL_MIN 2.2250738585072014E-308
#else /* for the rest of the world */
#   include <float.h>
#endif
#include <limits.h>

#  include <cmath>

#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif

#ifndef SQRT_DBL_EPSILON
#   define SQRT_DBL_EPSILON   1.4901161193847656e-08
#endif

#ifndef M_LN2
#   define M_LN2  0.69314718055994530942  /* log_e 2 */
#endif

#define ESCRIPT_DLL_API

#ifdef _WIN32
#   ifndef ESCRIPT_STATIC_LIB
#      undef ESCRIPT_DLL_API
#      ifdef ESCRIPT_EXPORTS
#         define ESCRIPT_DLL_API __declspec(dllexport)
#      else
#         define ESCRIPT_DLL_API __declspec(dllimport)
#      endif
#   endif
#endif

#ifndef ESCRIPT_MAX_DATA_RANK
#define ESCRIPT_MAX_DATA_RANK 4
#endif

#include <esysUtils/types.h>

#endif

