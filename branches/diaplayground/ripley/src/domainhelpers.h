
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/
#ifndef _RIPLEY_DOMAINHELPERS_H_
#define _RIPLEY_DOMAINHELPERS_H_

#include <ripley/Ripley.h>
#include <escript/Data.h>

namespace ripley {

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

/**
    sets va[a] = b and vb[b] = a, used in constructing CSC and CSR matrix
    formats simultaneously
*/
inline void doublyLink(std::vector<ripley::IndexVector>& va,
                       std::vector<ripley::IndexVector>& vb, int a, int b)
{
    va[a].push_back(b);
    vb[b].push_back(a);
}

/**
    factorises 'product' and inserts the factors into the vector 'factors'
    in order of smallest to largest
*/
void factorise(std::vector<int>& factors, int product);


#ifdef USE_BOOSTIO
/**
    converts the given gzip compressed char vector into an uncompressed form 
*/
std::vector<char> unzip(const std::vector<char>& compressed);
#endif // USE_BOOSTIO

} // namespace ripley

#endif // _RIPLEY_DOMAINHELPERS_H_

