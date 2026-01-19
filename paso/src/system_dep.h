
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/**
\file paso/src/system_dep.h
\ingroup Other
 */
/*
 @(#) system_dep.h
*/

#ifndef paso_system_dep_h
#define paso_system_dep_h

#define PASO_DLL_API

#ifdef _WIN32
# undef PASO_DLL_API
# ifdef PASO_EXPORTS
#   define PASO_DLL_API __declspec(dllexport)
# else
#   define PASO_DLL_API __declspec(dllimport)
# endif
#endif

#endif

