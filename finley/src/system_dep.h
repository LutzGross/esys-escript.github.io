
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014-2017 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
**
*****************************************************************************/


/**
\file finley/src/system_dep.h
\ingroup Other
 */
/*
 @(#) system_dep.h
*/

#ifndef finley_system_dep_h
#define finley_system_dep_h

#define FINLEY_DLL_API

#ifdef _WIN32
# undef FINLEY_DLL_API
# ifdef FINLEY_EXPORTS
#   define FINLEY_DLL_API __declspec(dllexport)
# else
#   define FINLEY_DLL_API __declspec(dllimport)
# endif
#endif

#endif

