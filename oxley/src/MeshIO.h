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

#ifndef __OXLEY_MESHIO_H__
#define __OXLEY_MESHIO_H__

#include <string>

#include <escript/AbstractDomain.h>
#include <escript/EsysMPI.h>

#include <oxley/OxleyDomain.h>

namespace oxley {

/**
   \brief
   The part of a saved mesh that p4est does not store.

   p4est's own files hold the connectivity and the quadrants, but nothing
   about where the domain sits in space or how many blocks it started from -
   that lives in oxley's forestData. Without it a saved mesh can only be read
   back into a domain that already has the right geometry, which is exactly
   what forced the old loadMesh to overwrite an existing domain in place.
   saveMesh writes this alongside, as <filename>.oxley.
*/
struct MeshHeader
{
    int dim = 0;                    ///< 2 or 3
    int order = 1;
    long n[3] = {0, 0, 0};          ///< blocks per axis
    double origin[3] = {0., 0., 0.};
    double extent[3] = {0., 0., 0.};
};

/// writes <filename>.oxley; throws OxleyException if the file cannot be written
void writeMeshHeader(const std::string& filename, const MeshHeader& h);

/// reads <filename>.oxley; throws OxleyException if it is missing or malformed
MeshHeader readMeshHeader(const std::string& filename);

/**
   \brief
   Reads a mesh written by saveMesh and returns it as a NEW domain.

   A free function rather than a method, for the same reason the domain has no
   refine methods: loading into an existing domain replaces its mesh, and every
   Data already built over that domain is then silently the wrong size.

   \param filename the name given to saveMesh, without extension
   \param comm optional MPI communicator, defaulting to the world
*/
escript::Domain_ptr loadMesh(const std::string& filename,
                             const boost::python::object& comm = boost::python::object());

} // namespace oxley

#endif // __OXLEY_MESHIO_H__
