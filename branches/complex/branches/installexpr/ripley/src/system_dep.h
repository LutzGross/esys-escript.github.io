
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef __RIPLEY_SYSTEM_DEP_H__
#define __RIPLEY_SYSTEM_DEP_H__

#if defined(_WIN32) && defined(__INTEL_COMPILER)
/*
 * The Intel compiler on windows has an "improved" math library compared to
 * the usual Visual C++ one. In particular it has acosh and other similar
 * functions which aren't implemented in Visual C++ math.h.
 * Note you will get a compile time error if any other header (including
 * system ones) includes math.h whilst mathimf.h has been included.
 * As a result system_dep.h must be included FIRST at all times (this
 * prevents math.h from being included).
 */
#   include <mathimf.h>
#else
#   include <cmath>
#endif

#define RIPLEY_DLL_API

#ifdef _WIN32
#   ifndef RIPLEY_STATIC_LIB
#      undef RIPLEY_DLL_API
#      ifdef RIPLEY_EXPORTS
#         define RIPLEY_DLL_API __declspec(dllexport)
#      else
#         define RIPLEY_DLL_API __declspec(dllimport)
#      endif
#   endif
#endif


// byte swapping / endianness:

#include <boost/detail/endian.hpp>

namespace ripley {

enum {
    BYTEORDER_NATIVE = BOOST_BYTE_ORDER,
    BYTEORDER_LITTLE_ENDIAN = 1234,
    BYTEORDER_BIG_ENDIAN = 4321
};

enum {
    DATATYPE_INT32 = 1,
    DATATYPE_FLOAT32,
    DATATYPE_FLOAT64
};

} // namespace

#ifdef _WIN32
#include <stdlib.h>
namespace ripley {
inline char* byte_swap32(char* val)
{
    unsigned long* v = reinterpret_cast<unsigned long*>(val);
    *v = _byteswap_ulong(*v);
    return val;
}
} // namespace

#else

#if HAVE_BYTESWAP_H
#   include <byteswap.h>
#elif HAVE_SYS_ENDIAN_H
#   include <sys/endian.h>
#   ifdef bswap32
#     define bswap_32(D) bswap32((D))    
#   endif
#elif HAVE_OSBYTEORDER_H
#   include <libkern/OSByteOrder.h>
#   define bswap_32 OSSwapInt32
#else // uh oh, we can't swap bytes...
#   define bswap_32(D) D
#endif // header selection

namespace ripley {
inline char* byte_swap32(char* val)
{
    unsigned int* v = reinterpret_cast<unsigned int*>(val);
    *v = bswap_32(*v);
    return val;
}
} // namespace ripley

#endif // WIN32


#endif // __RIPLEY_SYSTEM_DEP_H__

