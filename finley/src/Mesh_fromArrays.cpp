
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

#include "FinleyDomain.h"

#include <escript/index.h>

#include <sstream>

namespace finley {

namespace {

/// the contact element type that goes with a given face element type
ElementTypeId contactTypeFor(ElementTypeId faceType)
{
    switch (faceType) {
        case Point1: return Point1_Contact;
        case Line2:  return Line2_Contact;
        case Line3:  return Line3_Contact;
        case Tri3:   return Tri3_Contact;
        case Tri6:   return Tri6_Contact;
        case Rec4:   return Rec4_Contact;
        case Rec8:   return Rec8_Contact;
        case Rec9:   return Rec9_Contact;
        default: break;
    }
    std::stringstream ss;
    ss << "createFromArrays: no contact element type is defined for face "
          "element type " << faceType;
    throw escript::ValueError(ss.str());
}

/// fills an element table from flat arrays. `nodes` holds global node ids and
/// its length fixes the element count. Ids are made globally unique when the
/// caller does not supply them.
void fillElementTable(ElementFile* ef, const std::vector<index_t>& nodes,
                      const std::vector<index_t>& ids,
                      const std::vector<int>& tags,
                      const char* what, escript::JMPI mpiInfo)
{
    const int NN = ef->numNodes;
    if (NN < 1) {
        std::stringstream ss;
        ss << "createFromArrays: " << what << " reference element has no nodes";
        throw escript::ValueError(ss.str());
    }
    if (nodes.size() % NN != 0) {
        std::stringstream ss;
        ss << "createFromArrays: the " << what << " node table has "
           << nodes.size() << " entries, which is not a multiple of the " << NN
           << " nodes per element";
        throw escript::ValueError(ss.str());
    }
    const dim_t numElements = nodes.size() / NN;
    if (!ids.empty() && (dim_t)ids.size() != numElements) {
        std::stringstream ss;
        ss << "createFromArrays: " << what << " has " << numElements
           << " elements but " << ids.size() << " ids";
        throw escript::ValueError(ss.str());
    }
    if (!tags.empty() && (dim_t)tags.size() != numElements) {
        std::stringstream ss;
        ss << "createFromArrays: " << what << " has " << numElements
           << " elements but " << tags.size() << " tags";
        throw escript::ValueError(ss.str());
    }

    // When ids are not supplied, number the elements consecutively across ranks
    // so that they stay unique once the mesh is distributed.
    //
    // The scan runs UNCONDITIONALLY, and only its result is conditional. It is a
    // collective, and `ids.empty()` is a per-rank test: a rank that simply has
    // none of this kind of element - no boundary faces, say, because it owns
    // only interior cells - supplies an empty id list too, and cannot be told
    // apart from a caller that omitted them. Guarding the collective with that
    // test let such a rank enter the scan alone while the others went on, and
    // the run deadlocked here with the ranks in different collectives.
    index_t idOffset = 0;
#ifdef ESYS_MPI
    if (mpiInfo->size > 1) {
        index_t local = numElements;
        index_t scan = 0;
        MPI_Exscan(&local, &scan, 1, MPI_DIM_T, MPI_SUM, mpiInfo->comm);
        if (ids.empty() && mpiInfo->rank != 0)
            idOffset = scan;
    }
#endif

    ef->allocTable(numElements);
    ef->minColor = 0;
    ef->maxColor = numElements > 0 ? numElements - 1 : -1;

#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++) {
        ef->Id[e] = ids.empty() ? (idOffset + e) : ids[e];
        ef->Tag[e] = tags.empty() ? 0 : tags[e];
        ef->Owner[e] = mpiInfo->rank;
        ef->Color[e] = e;
        for (int k = 0; k < NN; k++)
            ef->Nodes[INDEX2(k, e, NN)] = nodes[e * NN + k];
    }
}

} // anonymous namespace

escript::Domain_ptr FinleyDomain::createFromArrays(const MeshArrays& in,
                                                   const std::string& name,
                                                   int order, int reducedOrder,
                                                   bool optimize,
                                                   escript::JMPI mpiInfo)
{
    if (in.numDim != 2 && in.numDim != 3) {
        std::stringstream ss;
        ss << "createFromArrays: numDim is " << in.numDim << ", must be 2 or 3";
        throw escript::ValueError(ss.str());
    }
    if (in.elementType == NoRef || in.faceElementType == NoRef)
        throw escript::ValueError("createFromArrays: element type and face "
                                  "element type must both be set");

    const dim_t numNodes = in.nodeId.size();
    if (in.nodeCoords.size() != (size_t)numNodes * in.numDim) {
        std::stringstream ss;
        ss << "createFromArrays: " << numNodes << " nodes in " << in.numDim
           << " dimensions need " << (size_t)numNodes * in.numDim
           << " coordinates, got " << in.nodeCoords.size();
        throw escript::ValueError(ss.str());
    }
    if (!in.nodeTag.empty() && (dim_t)in.nodeTag.size() != numNodes) {
        std::stringstream ss;
        ss << "createFromArrays: " << numNodes << " nodes but "
           << in.nodeTag.size() << " node tags";
        throw escript::ValueError(ss.str());
    }

    FinleyDomain* out = new FinleyDomain(name, in.numDim, mpiInfo);

    const_ReferenceElementSet_ptr refElements(
            new ReferenceElementSet(in.elementType, order, reducedOrder));
    const_ReferenceElementSet_ptr refFaceElements(
            new ReferenceElementSet(in.faceElementType, order, reducedOrder));
    const_ReferenceElementSet_ptr refContactElements(
            new ReferenceElementSet(contactTypeFor(in.faceElementType), order,
                                    reducedOrder));
    const_ReferenceElementSet_ptr refPoints(
            new ReferenceElementSet(Point1, order, reducedOrder));

    ElementFile* elements = new ElementFile(refElements, mpiInfo);
    out->setElements(elements);
    ElementFile* faces = new ElementFile(refFaceElements, mpiInfo);
    out->setFaceElements(faces);
    out->setContactElements(new ElementFile(refContactElements, mpiInfo));
    out->setPoints(new ElementFile(refPoints, mpiInfo));

    // node table. The global id doubles as the degree of freedom: unlike the
    // structured generators there is no periodicity to fold away here.
    NodeFile* nodes = out->getNodes();
    nodes->allocTable(numNodes);
#pragma omp parallel for
    for (index_t i = 0; i < numNodes; i++) {
        nodes->Id[i] = in.nodeId[i];
        nodes->Tag[i] = in.nodeTag.empty() ? 0 : in.nodeTag[i];
        nodes->globalDegreesOfFreedom[i] = in.nodeId[i];
        for (int d = 0; d < in.numDim; d++)
            nodes->Coordinates[INDEX2(d, i, in.numDim)] =
                    in.nodeCoords[(size_t)i * in.numDim + d];
    }

    fillElementTable(elements, in.elementNodes, in.elementId, in.elementTag,
                     "element", mpiInfo);
    fillElementTable(faces, in.faceNodes, in.faceId, in.faceTag,
                     "face element", mpiInfo);
    out->getContactElements()->allocTable(0);
    out->getPoints()->allocTable(0);

    for (TagMap::const_iterator it = in.tagMap.begin();
         it != in.tagMap.end(); ++it) {
        out->setTagMap(it->first, it->second);
    }

    // resolve the global node references, then distribute and build the
    // overlap, the mappings and the element colouring
    out->resolveNodeIds();
    out->prepare(optimize);
    return out->getPtr();
}

} // namespace finley
