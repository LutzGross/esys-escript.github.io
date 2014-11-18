/**
\file system_dep.h
\ingroup Other
 */
//
// @(#) system_dep.h
//

#ifndef bruce_system_dep_h
#define bruce_system_dep_h

#ifdef _WIN32

# ifdef __INTEL_COMPILER
// The Intel compiler on windows has an "improved" math library compared to the usual Visual C++ one
// In particular it has a acosh and other similar functions which aren't implemented in Visual C++ math.h
// Note you will get a compile time error if any other header (including system ones) includes math.h whilst mathimf.h
// has been included. As a result system_dep.h must be included FIRST at all times (this prevents math.h from being included).
#  include <mathimf.h>
# else
#  include <math.h>
# endif

#  ifndef INTERFACE_STATIC_LIB
#     ifdef BRUCE_EXPORTS
#        define BRUCE_DLL_API __declspec(dllexport)
#     else
#        define BRUCE_DLL_API __declspec(dllimport)
#     endif
#  endif

#else
# define BRUCE_DLL_API
#endif

#endif