
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

#include <ripley/domainhelpers.h>
#include <ripley/RipleyException.h>
#include <cmath>

#ifdef USE_BOOSTIO
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#endif

namespace ripley {

void factorise(std::vector<int>& factors, int product)
{
    int current = product;
    for (int p = 2; p <= sqrt(product); p++) {
        while (current % p == 0) {
            current /= p;
            factors.push_back(p);
        }
    }
    if (current != 1) {
        factors.push_back(current);
    }
}

#ifdef USE_BOOSTIO
std::vector<char> unzip(const std::vector<char>& compressed)
{
    std::vector<char> decompressed = std::vector<char>();

    boost::iostreams::filtering_ostream os;

    os.push(boost::iostreams::gzip_decompressor());
    os.push(boost::iostreams::back_inserter(decompressed));
    try {
        boost::iostreams::write(os, &compressed[0], compressed.size());
    } catch (boost::iostreams::gzip_error e) {
        switch(e.error()) {
            case boost::iostreams::gzip::zlib_error:
                throw RipleyException("Decompressing failed with: zlib error");
            case boost::iostreams::gzip::bad_crc:
                throw RipleyException("Decompressing failed with: CRC error");
            case boost::iostreams::gzip::bad_length:
                throw RipleyException("Decompressing failed with: bad length");
            case boost::iostreams::gzip::bad_header:
                throw RipleyException("Decompressing failed with: bad header");
            case boost::iostreams::gzip::bad_footer:
                throw RipleyException("Decompressing failed with: bad footer");
        }
    }
    return decompressed;
}
#endif

} // namespace ripley

