
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

#ifndef __RIPLEY_SYSTEM_DEP_H__
#define __RIPLEY_SYSTEM_DEP_H__

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

#include <escript/DataTypes.h>

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

namespace ripley {
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
} // namespace ripley

#endif // WIN32


#endif // __RIPLEY_SYSTEM_DEP_H__

