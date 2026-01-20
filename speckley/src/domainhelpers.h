
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
#ifndef _SPECKLEY_DOMAINHELPERS_H_
#define _SPECKLEY_DOMAINHELPERS_H_

#include <speckley/Speckley.h>
#include <escript/Data.h>

namespace speckley {

typedef std::map<std::string, escript::Data> DataMap;

/// returns the data associated with the string key or an empty data object
/// if the map does not contain the given key
inline const escript::Data unpackData(const std::string target,
                                      const DataMap& mapping)
{
    DataMap::const_iterator i = mapping.find(target);
    return (i == mapping.end() ? escript::Data() : i->second);
}

/**
    returns true if the given string is in the mapping AND the resulting
    coefficient is not empty
*/
inline bool isNotEmpty(const std::string target, const DataMap& mapping)
{
    DataMap::const_iterator i = mapping.find(target);
    return i != mapping.end() && !i->second.isEmpty();
}

/// returns true if the target data object is complex
inline bool isComplexCoef(const std::string target,
                                      const DataMap& mapping)
{
    DataMap::const_iterator i = mapping.find(target);
    return i->second.isComplex();
}

/**
    factorises 'product' and inserts the factors into the vector 'factors'
    in order of smallest to largest
*/
void factorise(std::vector<int>& factors, int product);

#ifdef ESYS_HAVE_BOOST_IO
/**
    converts the given gzip compressed char vector into an uncompressed form 
*/
std::vector<char> unzip(const std::vector<char>& compressed);
#endif // ESYS_HAVE_BOOST_IO

} //namespace speckley

#endif // _SPECKLEY_DOMAINHELPERS_H_
