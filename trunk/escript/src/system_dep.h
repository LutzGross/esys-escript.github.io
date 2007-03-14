/**
\file system_dep.h
\ingroup Other
 */
//
// @(#) system_dep.h
//

#ifndef escript_system_dep_h
#define escript_system_dep_h


#ifdef NO_FLOAT_H
#   define DBL_EPSILON 2.2204460492503131E-16
#   define DBL_MAX 1.7976931348623157E+308
#   define DBL_MIN 2.2250738585072014E-308
#else /* for the rest of the world */
#   include <float.h>
#endif
#include <limits.h>

# ifdef __INTEL_COMPILER
// The Intel compiler on windows has an "improved" math library compared to the usual Visual C++ one
// In particular it has a acosh and other similar functions which aren't implemented in Visual C++ math.h
// Note you will get a compile time error if any other header (including system ones) includes math.h whilst mathimf.h
// has been included. As a result system_dep.h must be included FIRST at all times (this prevents math.h from being included).
#  include <mathimf.h>
# else
#  include <math.h>
# endif

#ifndef M_PI
#   define M_PI 3.14159265358979323846
#endif

#ifndef SQRT_DBL_EPSILON
#   define SQRT_DBL_EPSILON   1.4901161193847656e-08
#endif

#ifndef M_LN2
#   define M_LN2  0.69314718055994530942  /* log_e 2 */
#endif

#ifdef NO_INV_HYPER
#   include "invhyp.h"
#   define asinh asinh_substitute
#   define acosh acosh_substitute
#   define atanh atanh_substitute
#endif

#ifdef _WIN32

#ifndef INTERFACE_STATIC_LIB
# ifdef ESCRIPT_EXPORTS
#   define ESCRIPT_DLL_API __declspec(dllexport)
#  else
#   define ESCRIPT_DLL_API __declspec(dllimport)
# endif
#endif

#else
# define ESCRIPT_DLL_API
#endif

#endif
