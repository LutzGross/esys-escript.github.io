
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**
\file system_dep.h
\ingroup Other
 */
//
// @(#) system_dep.h
//

#ifndef esysutils_system_dep_h
#define esysutils_system_dep_h

#if defined(_WIN32) && defined(__INTEL_COMPILER)
// The Intel compiler on windows has an "improved" math library compared to the usual Visual C++ one
// In particular it has a acosh and other similar functions which aren't implemented in Visual C++ math.h
// Note you will get a compile time error if any other header (including system ones) include math.h before mathimf.h
// has been included. As a result system_dep.h must be included FIRST at all times (this prevents math.h from being included).
#include <mathimf.h>
#else
#include <math.h>
#endif

#ifdef _WIN32

// Un-comment this block if you want it dynamic
// and comment out the line immediately after the block.
//#   ifndef INTERFACE_STATIC_LIB
//#      ifdef ESYSUTILS_EXPORTS
//#         define ESYSUTILS_DLL_API __declspec(dllexport)
//#      else
//#         define ESYSUTILS_DLL_API __declspec(dllimport)
//#      endif
//#   endif
#   define ESYSUTILS_DLL_API

// This is because of the different declarations of std::exception mentods
// on windows.
// Also, putting a "throw" in any declaration on windows causes a warning!!!!!!
// If you wish to generate a throw() on other systems, please use 
// THROW(NO_ARG). This is because windows generates warnings if you say
// THROW(), so the NO_ARG trick must be used to avoid the mass of warnings.

#   define THROW(ARG)
#else
#   define ESYSUTILS_DLL_API
#   define THROW(ARG) throw(ARG)
#endif

#define NO_ARG

#endif
