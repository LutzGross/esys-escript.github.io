
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


/**
\file dudley/src/CPPAdapter/system_dep.h
\ingroup Other
 */
/*
   @(#) system_dep.h
*/

#ifndef dudley_system_dep_h
#define dudley_system_dep_h

#define DUDLEY_DLL_API

#ifdef _WIN32

#   ifndef DUDLEY_STATIC_LIB
#      undef DUDLEY_DLL_API
#      ifdef DUDLEY_EXPORTS
#         define DUDLEY_DLL_API __declspec(dllexport)
#      else
#         define DUDLEY_DLL_API __declspec(dllimport)
#      endif
#   endif
#endif

#endif

