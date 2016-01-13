
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**
\file esysUtils/src/system_dep.h
\ingroup Other
 */
/*
* @(#) system_dep.h
*/

#ifndef esysutils_system_dep_h
#define esysutils_system_dep_h

#if defined(_WIN32) && defined(__INTEL_COMPILER)
/* The Intel compiler on windows has an "improved" math library compared to the usual Visual C++ one
* In particular it has a acosh and other similar functions which aren't implemented in Visual C++ math.h
* Note you will get a compile time error if any other header (including system ones) include math.h before mathimf.h
* has been included. As a result system_dep.h must be included FIRST at all times (this prevents math.h from being included).
*/
#include <mathimf.h>
#else
#include <math.h>
#endif

#define ESYSUTILS_DLL_API

#ifdef _WIN32
#   ifndef ESYSUTILS_STATIC_LIB
#      undef ESYSUTILS_DLL_API
#      ifdef ESYSUTILS_EXPORTS
#         define ESYSUTILS_DLL_API __declspec(dllexport)
#      else
#         define ESYSUTILS_DLL_API __declspec(dllimport)
#      endif
#   endif

/* This is because of the different declarations of std::exception mentods
*  on windows.
* Also, putting a "throw" in any declaration on windows causes a warning!!!!!!
* If you wish to generate a throw() on other systems, please use 
* THROW(NO_ARG). This is because windows generates warnings if you say
* THROW(), so the NO_ARG trick must be used to avoid the mass of warnings.
*/

#   define THROW(ARG)
#else
#   define THROW(ARG) throw(ARG)
#endif

#define NO_ARG

/* you'll need this one day. */
#ifndef __const
# if (defined __STDC__ && __STDC__) || defined __cplusplus
#  define __const	const
# else
#  define __const
# endif
#endif

#endif
