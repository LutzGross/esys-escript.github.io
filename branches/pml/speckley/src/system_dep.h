
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

#ifndef __SPECKLEY_SYSTEM_DEP_H__
#define __SPECKLEY_SYSTEM_DEP_H__

#include <escript/DataTypes.h>

#define Speckley_DLL_API

#ifdef _WIN32
#   ifndef Speckley_STATIC_LIB
#      undef Speckley_DLL_API
#      ifdef Speckley_EXPORTS
#         define Speckley_DLL_API __declspec(dllexport)
#      else
#         define Speckley_DLL_API __declspec(dllimport)
#      endif
#   endif
#endif


// byte swapping / endianness:

#include <boost/detail/endian.hpp>

namespace speckley {

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
namespace speckley {
inline char* byte_swap32(char* val)
{
    unsigned long* v = reinterpret_cast<unsigned long*>(val);
    *v = _byteswap_ulong(*v);
    return val;
}
inline char* byte_swap64(char* val)
{
    unsigned __int64* v = reinterpret_cast<unsigned __int64*>(val);
    *v = _byteswap_uint64(*v);
    return val;
}
} // namespace

#else

#include <stdint.h> // uint64_t

#if HAVE_BYTESWAP_H
#   include <byteswap.h>
#elif HAVE_SYS_ENDIAN_H
#   include <sys/endian.h>
#   ifdef bswap32
#     define bswap_32(D) bswap32((D))
#   endif
#   ifdef bswap64
#     define bswap_64(D) bswap64((D))
#   endif
#elif HAVE_OSBYTEORDER_H
#   include <libkern/OSByteOrder.h>
#   define bswap_32 OSSwapInt32
#   define bswap_64 OSSwapInt64
#else // uh oh, we can't swap bytes...
#   define bswap_32(D) (D)
#   define bswap_64(D) (D)
#endif // header selection

namespace speckley {
inline char* byte_swap32(char* val)
{
    unsigned int* v = reinterpret_cast<unsigned int*>(val);
    *v = bswap_32(*v);
    return val;
}

inline char* byte_swap64(char* val)
{
    uint64_t* v = reinterpret_cast<uint64_t*>(val);
    *v = bswap_64(*v);
    return val;
}
} // namespace speckley

#endif // WIN32


#endif // __SPECKLEY_SYSTEM_DEP_H__

