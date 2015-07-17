
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

