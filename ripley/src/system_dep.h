
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

#ifndef __RIPLEY_SYSTEM_DEP_H__
#define __RIPLEY_SYSTEM_DEP_H__

#if defined(_WIN32) && defined(__INTEL_COMPILER)
/*
 * The Intel compiler on windows has an "improved" math library compared to
 * the usual Visual C++ one. In particular it has acosh and other similar
 * functions which aren't implemented in Visual C++ math.h.
 * Note you will get a compile time error if any other header (including
 * system ones) includes math.h whilst mathimf.h has been included.
 * As a result system_dep.h must be included FIRST at all times (this
 * prevents math.h from being included).
 */
#   include <mathimf.h>
#else
#   include <math.h>
#endif

#define RIPLEY_DLL_API

#ifdef _WIN32
#   ifndef RIPLEY_STATIC_LIB
#      undef RIPLEY_DLL_API
#      ifdef RIPLEY_EXPORTS
#         define RIPLEY_DLL_API __declspec(dllexport)
#      else
#         define RIPLEY_DLL_API __declspec(dllimport)
#      endif
#   endif
#endif

#endif // __RIPLEY_SYSTEM_DEP_H__

