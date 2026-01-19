/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include <escript/DataTypes.h>

// byte swapping / endianness:

#include <boost/detail/endian.hpp>
#ifdef ESYS_DEPRECATED_BOOST_ENDIAN
#include <boost/predef/other/endian.h>
#endif

namespace oxley {

enum {
#ifndef ESYS_DEPRECATED_BOOST_ENDIAN
    BYTEORDER_NATIVE = BOOST_BYTE_ORDER,
#elif defined(ESYS_DEPRECATED_BOOST_ENDIAN) && BOOST_ENDIAN_BIG_BYTE
    BYTEORDER_NATIVE = 4321,
#elif defined(ESYS_DEPRECATED_BOOST_ENDIAN) && BOOST_ENDIAN_LITTLE_BYTE
    BYTEORDER_NATIVE = 1234,
#elif defined(ESYS_DEPRECATED_BOOST_ENDIAN) && BOOST_ENDIAN_LITTLE_WORD
    BYTEORDER_NATIVE = 2134,
#endif
    BYTEORDER_LITTLE_ENDIAN = 1234,
    BYTEORDER_BIG_ENDIAN = 4321
};

enum {
    DATATYPE_INT32 = 1,
    DATATYPE_FLOAT32,
    DATATYPE_FLOAT64
};

} // namespace
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

namespace oxley {
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
} // namespace oxley




