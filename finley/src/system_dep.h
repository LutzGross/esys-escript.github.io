
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
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

