
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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
\file escript/src/system_dep.h
\ingroup Other
 */
/*
 @(#) system_dep.h
*/

#ifndef escript_system_dep_h
#define escript_system_dep_h

#define ESCRIPT_DLL_API

#ifdef _WIN32
# undef ESCRIPT_DLL_API
# ifdef ESCRIPT_EXPORTS
#   define ESCRIPT_DLL_API __declspec(dllexport)
# else
#   define ESCRIPT_DLL_API __declspec(dllimport)
# endif
#endif

#endif

