
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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
#define ESCRIPT_INLINE_DLL_API

#ifdef _WIN32
# undef ESCRIPT_DLL_API
# undef ESCRIPT_INLINE_DLL_API
# ifdef ESCRIPT_EXPORTS
#   define ESCRIPT_DLL_API __declspec(dllexport)
#   define ESCRIPT_INLINE_DLL_API __declspec(dllexport)
/* TODO: fix boost numpy
#   define BOOST_PYTHON_STATIC_LIB
#   define BOOST_NUMPY_STATIC_LIB */
# else
#   define ESCRIPT_DLL_API __declspec(dllimport)
#   define ESCRIPT_INLINE_DLL_API
# endif

#endif

#endif

