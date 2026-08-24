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
*
*****************************************************************************/

#include <oxley/MeshIO.h>
#include <oxley/Brick.h>
#include <oxley/OxleyException.h>
#include <oxley/Rectangle.h>

#include <fstream>
#include <sstream>

namespace oxley {

namespace {
const char* HEADER_TAG = "oxley-mesh";
const int HEADER_VERSION = 1;

std::string headerName(const std::string& filename)
{
    return filename + ".oxley";
}
} // anonymous namespace

void writeMeshHeader(const std::string& filename, const MeshHeader& h)
{
    std::ofstream f(headerName(filename).c_str());
    if(!f)
        throw OxleyException("saveMesh: cannot write " + headerName(filename));
    f.precision(17);
    f << HEADER_TAG << " " << HEADER_VERSION << "\n"
      << "dim " << h.dim << "\n"
      << "order " << h.order << "\n"
      << "blocks " << h.n[0] << " " << h.n[1] << " " << h.n[2] << "\n"
      << "origin " << h.origin[0] << " " << h.origin[1] << " " << h.origin[2] << "\n"
      << "extent " << h.extent[0] << " " << h.extent[1] << " " << h.extent[2] << "\n";
    if(!f)
        throw OxleyException("saveMesh: failed writing " + headerName(filename));
}

MeshHeader readMeshHeader(const std::string& filename)
{
    std::ifstream f(headerName(filename).c_str());
    if(!f)
        throw OxleyException("loadMesh: cannot read " + headerName(filename)
                + ". A mesh saved by an older version has no header and can "
                  "only be read back into a domain of matching geometry.");

    std::string tag;
    int version = 0;
    f >> tag >> version;
    if(tag != HEADER_TAG)
        throw OxleyException("loadMesh: " + headerName(filename)
                + " is not an oxley mesh header.");
    if(version != HEADER_VERSION) {
        std::stringstream ss;
        ss << "loadMesh: " << headerName(filename) << " has version " << version
           << ", this build reads version " << HEADER_VERSION << ".";
        throw OxleyException(ss.str());
    }

    MeshHeader h;
    std::string key;
    while(f >> key) {
        if(key == "dim")
            f >> h.dim;
        else if(key == "order")
            f >> h.order;
        else if(key == "blocks")
            f >> h.n[0] >> h.n[1] >> h.n[2];
        else if(key == "origin")
            f >> h.origin[0] >> h.origin[1] >> h.origin[2];
        else if(key == "extent")
            f >> h.extent[0] >> h.extent[1] >> h.extent[2];
        else
            throw OxleyException("loadMesh: unexpected entry '" + key + "' in "
                    + headerName(filename));
    }
    if(h.dim != 2 && h.dim != 3)
        throw OxleyException("loadMesh: the header names neither 2 nor 3 "
                "dimensions.");
    for(int d = 0; d < h.dim; ++d)
        if(h.n[d] <= 0)
            throw OxleyException("loadMesh: the header has a non-positive "
                    "block count.");
    return h;
}

escript::Domain_ptr loadMesh(const std::string& filename,
                             const boost::python::object& comm)
{
    const MeshHeader h = readMeshHeader(filename);
    escript::JMPI jmpi = escript::makeInfoFromPyComm(comm);

    // Build a domain of the recorded geometry, then give it the saved forest.
    // The domain is replaced BEFORE anyone can hold a Data over it, which is
    // what makes this safe where the old in-place loadMesh was not: there the
    // caller already had a domain, and every Data built on it beforehand was
    // left sized for a mesh that no longer existed.
    const std::vector<double> noPoints;
    const std::vector<int> noTags;
    const TagMap noNames;
    const std::vector<int> level(1, 0);

    if(h.dim == 2) {
        Rectangle* dom = new Rectangle(jmpi, h.order, h.n[0], h.n[1],
                h.origin[0], h.origin[1], h.extent[0], h.extent[1],
                noPoints, noTags, noNames, level);
        escript::Domain_ptr result(dom);
        dom->loadMesh(filename);
        return result;
    }

    Brick* dom = new Brick(jmpi, h.order, h.n[0], h.n[1], h.n[2],
            h.origin[0], h.origin[1], h.origin[2],
            h.extent[0], h.extent[1], h.extent[2],
            noPoints, noTags, noNames, level);
    escript::Domain_ptr result(dom);
    dom->loadMesh(filename);
    return result;
}

} // namespace oxley
