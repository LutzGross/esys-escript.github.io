
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <speckley/domainhelpers.h>
#include <speckley/SpeckleyException.h>

#include <cmath>

#ifdef ESYS_HAVE_BOOST_IO
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#endif

namespace speckley {

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

#ifdef ESYS_HAVE_BOOST_IO
std::vector<char> unzip(const std::vector<char>& compressed)
{
    std::vector<char> decompressed = std::vector<char>();

    boost::iostreams::filtering_ostream os;

    os.push(boost::iostreams::gzip_decompressor());
    os.push(boost::iostreams::back_inserter(decompressed));
    try {
        boost::iostreams::write(os, &compressed[0], compressed.size());
    } catch (boost::iostreams::gzip_error& e) {
        switch(e.error()) {
            case boost::iostreams::gzip::zlib_error:
                throw speckley::SpeckleyException("Decompressing failed with: zlib error");
            case boost::iostreams::gzip::bad_crc:
                throw speckley::SpeckleyException("Decompressing failed with: CRC error");
            case boost::iostreams::gzip::bad_length:
                throw speckley::SpeckleyException("Decompressing failed with: bad length");
            case boost::iostreams::gzip::bad_header:
                throw speckley::SpeckleyException("Decompressing failed with: bad header");
            case boost::iostreams::gzip::bad_footer:
                throw speckley::SpeckleyException("Decompressing failed with: bad footer");
        }
    }
    return decompressed;
}
#endif

} // namespace speckley

