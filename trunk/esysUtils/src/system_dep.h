
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
\file esysUtils/src/system_dep.h
\ingroup Other
 */
/*
* @(#) system_dep.h
*/

#ifndef esysutils_system_dep_h
#define esysutils_system_dep_h

#include <cmath>

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

#endif

#endif

