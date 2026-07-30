
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

#include <weipa/OxleyElements.h>
#include <weipa/NodeData.h>
#include <weipa/WeipaException.h>

#include <iostream>

#ifndef VISIT_PLUGIN
#include <oxley/OxleyDomain.h>
#include <oxley/Brick.h>
#include <oxley/Rectangle.h>
#endif

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#endif

#ifndef VISIT_PLUGIN
using escript::DataTypes::dim_t;
#endif

using namespace std;

namespace weipa {
    
//
// Constructor
//
OxleyElements::OxleyElements(const string& elementName, OxleyNodes_ptr nodeData)
    : originalMesh(nodeData), name(elementName), numElements(0),
      numGhostElements(0), nodesPerElement(0),
      type(ZONETYPE_UNKNOWN)
{
    nodeMesh.reset(new OxleyNodes(name));
}

//
// Copy constructor
//
OxleyElements::OxleyElements(const OxleyElements& e)
{
    name = e.name;
    numElements = e.numElements;
    numGhostElements = e.numGhostElements;
    type = e.type;
    nodesPerElement = e.nodesPerElement;
    originalMesh = e.originalMesh;
    if (e.nodeMesh)
        nodeMesh.reset(new OxleyNodes(*e.nodeMesh));
    else
        nodeMesh.reset(new OxleyNodes(name));

    nodes = e.nodes;
    ID = e.ID;
    //tag = e.tag;
    owner = e.owner;
}

//
//
//
bool OxleyElements::initFromOxley(const oxley::OxleyDomain* dom, int fsType)
{
#ifndef VISIT_PLUGIN
    const std::pair<int,dim_t> shape = dom->getDataShape(fsType);
    const dim_t* faces = dom->getNumFacesPerBoundary();

    numElements = shape.second;

    // For FaceElements the mesh-access boundary walk is authoritative: it finds
    // the boundary from the p4est connectivity, whereas getDataShape() reports
    // the count from updateFaceElementCount(), which compares coordinates against
    // the domain extent. Take the count from the same place as the connectivity
    // so the two cannot disagree.
    if (fsType == oxley::FaceElements || fsType == oxley::ReducedFaceElements)
        numElements = (int) dom->getMeshAccess(true).numFaces;

    if (numElements > 0) {
        nodesPerElement = shape.first;
        switch (nodesPerElement) {
            case 2:
                type = ZONETYPE_BEAM;
                break;
            case 4:
                type = ZONETYPE_QUAD;
                break;
            case 8:
                type = ZONETYPE_HEX;
                break;
            default:
                throw WeipaException("OxleyElements: unsupported element shape");
        }
        owner = dom->getOwnerVector(fsType);
        if ((int) owner.size() != numElements)
            owner.assign(numElements, dom->getMPIRank());   // faces are all local

        // serial element ids; TODO(MPI, A6): global element numbering.
        ID.resize(numElements);
        for (int i = 0; i < numElements; i++)
            ID[i] = i;

        nodes.clear();
        if (dom->getDim() == 2) {
            const oxley::Rectangle * rect = static_cast<const oxley::Rectangle *>(dom);

            if (fsType==oxley::Elements) {
                // volume elements: connectivity from the lnodes mesh-access
                // interface, reordering p4est z-order corners into weipa's quad
                // order (z-order 0,1,3,2 == BL,BR,TR,TL).
                const oxley::MeshAccess m = dom->getMeshAccess(true);
                static const int quadMap[4] = {0, 1, 3, 2};
                nodes.reserve((size_t) m.numElements * 4);
                for (long e = 0; e < m.numElements; ++e)
                    for (int k = 0; k < 4; ++k)
                        nodes.push_back((int) m.elementNodes[(size_t) e*4 + quadMap[k]]);
            } else if (fsType==oxley::FaceElements) {
                // Boundary edges straight from the mesh-access interface, which
                // finds them from the p4est connectivity ("this tree face has no
                // neighbour") and winds each so the outward normal is the
                // clockwise rotation of the tangent. The same arrays feed the
                // finley converter, whose normals test (int n.x dS == dim*volume)
                // and per-tag boundary areas check exactly this data.
                //
                // The previous code collected boundary NODES into four buckets
                // and concatenated them, so the Line2 pairs were arbitrary.
                const oxley::MeshAccess m = dom->getMeshAccess(true);
                nodes.reserve((size_t) m.numFaces * 2);
                for (long f = 0; f < m.numFaces; ++f)
                    for (int k = 0; k < 2; ++k)
                        nodes.push_back((int) m.faceNodes[(size_t) f*2 + k]);
                tag.assign(m.faceTags.begin(), m.faceTags.end());
            }
        } else { //3d
            const oxley::Brick * brick = static_cast<const oxley::Brick *>(dom);
            if (fsType==oxley::Elements) {
                // volume elements: connectivity from the lnodes mesh-access
                // interface, reordering p8est z-order corners into weipa's hex
                // order (z-order 0,4,5,1,2,6,7,3).
                const oxley::MeshAccess m = dom->getMeshAccess(true);
                static const int hexMap[8] = {0, 4, 5, 1, 2, 6, 7, 3};
                nodes.reserve((size_t) m.numElements * 8);
                for (long e = 0; e < m.numElements; ++e)
                    for (int k = 0; k < 8; ++k)
                        nodes.push_back((int) m.elementNodes[(size_t) e*8 + hexMap[k]]);
            } else if (fsType==oxley::FaceElements) {
                // Boundary faces straight from the mesh-access interface; see the
                // 2D branch above. Each is a Rec4 wound counter-clockwise as seen
                // from OUTSIDE, so the right-hand rule gives the outward normal.
                const oxley::MeshAccess m = dom->getMeshAccess(true);
                nodes.reserve((size_t) m.numFaces * 4);
                for (long f = 0; f < m.numFaces; ++f)
                    for (int k = 0; k < 4; ++k)
                        nodes.push_back((int) m.faceNodes[(size_t) f*4 + k]);
                tag.assign(m.faceTags.begin(), m.faceTags.end());

            }
        }

        buildMeshes();
    }
    return true;

#else // VISIT_PLUGIN
    return false;
#endif
    return false;
}

StringVec OxleyElements::getMeshNames() const
{
    StringVec res;
    if (nodeMesh)
        res.push_back(nodeMesh->getName());
    return res;
}

StringVec OxleyElements::getVarNames() const
{
    StringVec res;
    res.push_back(name + string("_Id"));
    res.push_back(name + string("_Owner"));
    //res.push_back(name + string("_Tag"));
    return res;
}

const IntVec& OxleyElements::getVarDataByName(const string varName) const
{
    if (varName == name+string("_Id"))
        return ID;
    if (varName == name+string("_Owner"))
        return owner;
    //if (varName == name+string("_Tag"))
    //    return tag;

    throw "Invalid variable name";
}

void OxleyElements::reorderArray(IntVec& v, const IntVec& idx,
                               int elementsPerIndex)
{
    IntVec newArray(v.size());
    IntVec::iterator arrIt = newArray.begin();
    IntVec::const_iterator idxIt;
    if (elementsPerIndex == 1) {
        for (idxIt=idx.begin(); idxIt!=idx.end(); idxIt++) {
            *arrIt++ = v[*idxIt];
        }
    } else {
        for (idxIt=idx.begin(); idxIt!=idx.end(); idxIt++) {
            int i = *idxIt;
            int* start = &v[i*elementsPerIndex];
            copy(start, start+elementsPerIndex, arrIt);
            arrIt += elementsPerIndex;
        }
    }
    v.swap(newArray);
}

IntVec OxleyElements::prepareGhostIndices(int ownIndex)
{
    IntVec indexArray;
    numGhostElements = 0;
    
    // move indices of "ghost zones" to the end to be able to reorder
    // data accordingly
    for (dim_t i=0; i<numElements; i++) {
        if (owner[i] == ownIndex)
            indexArray.push_back(i);
    }

    for (dim_t i=0; i<numElements; i++) {
        if (owner[i] != ownIndex) {
            numGhostElements++;
            indexArray.push_back(i);
        }
    }
    return indexArray;
}

void OxleyElements::reorderGhostZones(int ownIndex)
{
    IntVec indexArray = prepareGhostIndices(ownIndex);

    // move "ghost data" to the end of the arrays
    if (numGhostElements > 0) {
        reorderArray(nodes, indexArray, nodesPerElement);
        reorderArray(owner, indexArray, 1);
        reorderArray(ID, indexArray, 1);
        //reorderArray(tag, indexArray, 1);
    }
}

void OxleyElements::removeGhostZones(int ownIndex)
{
    reorderGhostZones(ownIndex);

    if (numGhostElements > 0) {
        numElements -= numGhostElements;
        nodes.resize(numElements*nodesPerElement);
        owner.resize(numElements);
        ID.resize(numElements);
        //tag.resize(numElements);
        numGhostElements = 0;
    }
}

void OxleyElements::buildMeshes()
{
    // build a new mesh containing only the required nodes
    if (numElements > 0) {
        if (nodeMesh && nodeMesh->getNumNodes() > 0) {
            OxleyNodes_ptr newMesh(new OxleyNodes(nodeMesh, nodes, name));
            nodeMesh.swap(newMesh);
        } else {
            nodeMesh.reset(new OxleyNodes(originalMesh, nodes, name));
        }
#ifdef _DEBUG
        cout << nodeMesh->getName() << " has " << nodeMesh->getNumNodes()
            << " nodes and " << numElements << " elements" << endl;
#endif
    }
}

void OxleyElements::writeConnectivityVTK(ostream& os)
{
    if (numElements > 0) {
        const IntVec& gNI = nodeMesh->getGlobalNodeIndices();
        IntVec::const_iterator it;
        int count = 1;
        for (it=nodes.begin(); it!=nodes.end(); it++, count++) {
            os << gNI[*it];
            if (count % nodesPerElement == 0)
                os << endl;
            else
                os << " ";
        }
    }
}

#ifdef ESYS_HAVE_SILO
inline int toSiloElementType(int type)
{
    switch (type) {
        case ZONETYPE_BEAM: return DB_ZONETYPE_BEAM;
        case ZONETYPE_HEX: return DB_ZONETYPE_HEX;
        case ZONETYPE_POLYGON: return DB_ZONETYPE_POLYGON;
        case ZONETYPE_QUAD: return DB_ZONETYPE_QUAD;
    }
    return 0;
}
#endif

bool OxleyElements::writeToSilo(DBfile* dbfile, const string& siloPath,
                                 const StringVec& labels,
                                 const StringVec& units, bool writeMeshData)
{
#ifdef ESYS_HAVE_SILO
    if (numElements == 0)
        return true;

    int ret;
    if (siloPath != "") {
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    }

    // write out the full mesh in any case
    nodeMesh->setSiloPath(siloPath);
    string siloMeshNameStr = nodeMesh->getFullSiloName();
    const char* siloMeshName = siloMeshNameStr.c_str();
    int arraylen = numElements * nodesPerElement;
    int eltype = toSiloElementType(type);

    string varName = name + string("_zones");
    ret = DBPutZonelist2(dbfile, varName.c_str(), numElements,
            nodeMesh->getNumDims(), &nodes[0], arraylen, 0, 0,
            numGhostElements, &eltype, &nodesPerElement, &numElements, 1, NULL);

    if (ret == 0) {
        CoordArray& coordbase = const_cast<CoordArray&>(nodeMesh->getCoords());
        DBoptlist* optList = NULL;
        int nOpts = labels.size()+units.size();
        if (nOpts>0) {
            optList = DBMakeOptlist(nOpts);
            if (labels.size()>0)
                DBAddOption(optList, DBOPT_XLABEL, (void*)labels[0].c_str());
            if (labels.size()>1)
                DBAddOption(optList, DBOPT_YLABEL, (void*)labels[1].c_str());
            if (labels.size()>2)
                DBAddOption(optList, DBOPT_ZLABEL, (void*)labels[2].c_str());
            if (units.size()>0)
                DBAddOption(optList, DBOPT_XUNITS, (void*)units[0].c_str());
            if (units.size()>1)
                DBAddOption(optList, DBOPT_YUNITS, (void*)units[1].c_str());
            if (units.size()>2)
                DBAddOption(optList, DBOPT_ZUNITS, (void*)units[2].c_str());
        }
        ret = DBPutUcdmesh(dbfile, siloMeshName,
                nodeMesh->getNumDims(), NULL, &coordbase[0],
                nodeMesh->getNumNodes(), numElements, varName.c_str(),
                /*"facelist"*/NULL, DB_FLOAT, optList);

        if (optList)
            DBFreeOptlist(optList);
    }
    
    if (ret != 0)
        return false;

    // write out the element-centered variables if enabled
    if (writeMeshData) {
        varName = name + string("_Id");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
                (float*)&ID[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
                NULL);
        if (ret == 0) {
            varName = name + string("_Owner");
            ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
                (float*)&owner[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
                NULL);
        }
    }

    // "Elements" is a special case
    if (writeMeshData && name == "Elements") {
        nodeMesh->writeToSilo(dbfile);
    }

    return (ret == 0);

#else // !ESYS_HAVE_SILO
    return false;
#endif
}

} // namespace weipa

