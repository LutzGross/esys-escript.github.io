
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

//
// ElementData.cpp
//
#include <escriptreader/ElementData.h>
#include <netcdf.hh>
#if HAVE_SILO
#include <silo.h>
#endif

using namespace std;

enum {
    ZONETYPE_BEAM=1,
    ZONETYPE_HEX,
    ZONETYPE_POLYGON,
    ZONETYPE_QUAD,
    ZONETYPE_TET,
    ZONETYPE_TRIANGLE,
};

// The following arrays contain indices to convert unsupported element
// types into supported ones (lines, triangles, quads, hexahedrons)

static const size_t line3indices[2*2] = {
    0, 2,
    2, 1
};
static const size_t tri6indices[4*3] = {
    0, 3, 5,
    3, 1, 4,
    5, 4, 2,
    3, 4, 5
};
static const size_t rec8indices[6*3] = {
    0, 4, 7,
    4, 1, 5,
    5, 2, 6,
    6, 3, 7,
    4, 6, 7,
    4, 5, 6
};
static const size_t rec9indices[4*4] = {
    0, 4, 8, 7,
    4, 1, 5, 8,
    7, 8, 6, 3,
    8, 5, 2, 6
};
static const size_t tet10indices[8*4] = {
    6, 0, 7, 4,
    2, 6, 9, 5,
    9, 7, 3, 8,
    4, 1, 8, 5,
    6, 7, 9, 8,
    5, 6, 9, 8,
    6, 5, 4, 8,
    7, 6, 4, 8
};
static const size_t hex20indices[36*3] = {
     0,  8, 12,   8,  1, 13,  13,  5, 16,
    16,  4, 12,   8, 13, 16,   8, 16, 12,
     1,  9, 13,   9,  2, 14,  14,  6, 17,
    17,  5, 13,   9, 14, 17,   9, 17, 13,
     2, 10, 14,  10,  3, 15,  15,  7, 18,
    18, 14,  6,  10, 15, 18,  10, 18, 14,
     3, 11, 15,  11,  0, 12,  12,  4, 19,
    19,  7, 15,  11, 12, 19,  11, 19, 15,
     4, 16, 19,  16,  5, 17,  17,  6, 18,
    18,  7, 19,  16, 17, 18,  16, 18, 19,
     3, 10, 11,  10,  2,  9,   9,  1,  8,
     8,  0, 11,  10,  9,  8,  10,  8, 11
};
static const size_t hex27indices[8*8] = {
     0,  8, 20, 11, 12, 21, 26, 24,
     8,  1,  9, 20, 21, 13, 22, 26,
    11, 20, 10,  3, 24, 26, 23, 15,
    20,  9,  2, 10, 26, 22, 14, 23,
    12, 21, 26, 24,  4, 16, 25, 19,
    21, 13, 22, 26, 16,  5, 17, 25,
    24, 26, 23, 15, 19, 25, 18,  7,
    26, 22, 14, 23, 25, 17,  6, 18
};

//
// Constructor
//
ElementData::ElementData(const string& elementName, const Mesh* mainMesh)
    : name(elementName), count(0), reducedCount(0), numGhostElements(0),
    numReducedGhostElements(0), fullMesh(NULL), reducedMesh(NULL),
    originalMesh(mainMesh), fullMeshIsOriginalMesh(false)
{
}

//
// Copy constructor
//
ElementData::ElementData(const ElementData& e)
{
    name = e.name;
    count = e.count;
    reducedCount = e.reducedCount;
    numGhostElements = e.numGhostElements;
    numReducedGhostElements = e.numReducedGhostElements;
    numDims = e.numDims;
    type = e.type;
    reducedType = e.reducedType;
    nodesPerElement = e.nodesPerElement;
    reducedNodesPerElement = e.reducedNodesPerElement;

    originalMesh = e.originalMesh;
    fullMeshIsOriginalMesh = e.fullMeshIsOriginalMesh;
    if (fullMeshIsOriginalMesh)
        fullMesh = const_cast<Mesh*>(originalMesh);
    else if (e.fullMesh)
        fullMesh = new Mesh(*e.fullMesh);
    else
        fullMesh = NULL;

    if (e.reducedMesh)
        reducedMesh = new Mesh(*e.reducedMesh);
    else
        reducedMesh = NULL;

    nodes = e.nodes;
    reducedNodes = e.reducedNodes;
    ID = e.ID;
    color = e.color;
    tag = e.tag;
    owner = e.owner;
    reducedOwner = e.reducedOwner;
    indexArray = e.indexArray;
    reducedIndexArray = e.reducedIndexArray;
    ID2idx = e.ID2idx;
}

//
// Destructor
//
ElementData::~ElementData()
{
    if (reducedMesh)
        delete reducedMesh;
    if (fullMesh && !fullMeshIsOriginalMesh)
        delete fullMesh;
}

StringVec ElementData::getMeshNames() const
{
    StringVec res;
    if (fullMesh && !fullMeshIsOriginalMesh)
        res.push_back(fullMesh->getName());
    if (reducedMesh)
        res.push_back(reducedMesh->getName());
    return res;
}

StringVec ElementData::getVarNames() const
{
    StringVec res;
    if (count > 0) {
        res.push_back(name + string("_Color"));
        res.push_back(name + string("_Id"));
        res.push_back(name + string("_Owner"));
        res.push_back(name + string("_Tag"));
    }
    return res;
}

const IntVec& ElementData::getVarDataByName(const string varName) const
{
    if (varName == name+string("_Color"))
        return color;
    else if (varName == name+string("_Id"))
        return ID;
    else if (varName == name+string("_Owner")) {
        if (reducedCount > 0)
            return reducedOwner;
        else
            return owner;
    } else if (varName == name+string("_Tag"))
        return tag;
    else
        return *(IntVec*)(NULL);
}

//
// Reads element data from given NetCDF file
//
bool ElementData::readFromNc(NcFile* ncfile)
{
    NcAtt* att;
    NcVar* var;

    string num_str("num_");
    num_str += name;

    att = ncfile->get_att(num_str.c_str());
    count = att->as_int(0);

    // Only attempt to read data if there are any elements.
    // Having no elements is not an error.
    if (count == 0)
        return true;

    att = ncfile->get_att((num_str + string("_numNodes")).c_str());
    nodesPerElement = att->as_int(0);

    nodes.insert(nodes.end(), count*nodesPerElement, 0);
    var = ncfile->get_var((name + string("_Nodes")).c_str());
    var->get(&nodes[0], count, nodesPerElement);

    color.insert(color.end(), count, 0);
    var = ncfile->get_var((name + string("_Color")).c_str());
    var->get(&color[0], count);

    ID.insert(ID.end(), count, 0);
    var = ncfile->get_var((name + string("_Id")).c_str());
    var->get(&ID[0], count);

    owner.insert(owner.end(), count, 0);
    var = ncfile->get_var((name + string("_Owner")).c_str());
    var->get(&owner[0], count);

    tag.insert(tag.end(), count, 0);
    var = ncfile->get_var((name + string("_Tag")).c_str());
    var->get(&tag[0], count);

    att = ncfile->get_att((name + string("_TypeId")).c_str());
    FinleyElementInfo f = getFinleyTypeInfo((ElementTypeId)att->as_int(0));
    type = f.elementType;
    reducedType = f.reducedElementType;

    // build elementID->index map
    buildIndexMap();
    
    if (f.elementFactor > 1 || f.reducedElementSize != nodesPerElement)
        buildReducedElements(f);

    buildMeshes();

    return true; 
}

//
//
//
void ElementData::reorderArray(IntVec& v, const IntVec& idx,
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
            copy(&v[i*elementsPerIndex], &v[(i+1)*elementsPerIndex], arrIt);
            arrIt += elementsPerIndex;
        }
    }
    v.swap(newArray);
}

//
//
//
void ElementData::buildReducedElements(const FinleyElementInfo& f)
{
    reducedNodes.clear();
    reducedNodes.insert(reducedNodes.end(),
            f.reducedElementSize*count, 0);
    IntVec::iterator reducedIt = reducedNodes.begin();
    IntVec::const_iterator origIt;
    for (origIt=nodes.begin(); origIt!=nodes.end();
         origIt+=nodesPerElement)
    {
        std::copy(origIt, origIt+f.reducedElementSize, reducedIt);
        reducedIt += f.reducedElementSize;
    }

    // Remove comment to save disk space - we don't really need the full
    // elements except for the main mesh
    if (f.elementFactor > 1 /*&& name == "Elements"*/) {
        // replace each element by multiple smaller ones
        IntVec fullNodes(f.elementSize*f.elementFactor*count);
        IntVec::iterator cellIt = fullNodes.begin();

        // copy owner data
        owner.swap(reducedOwner);
        owner.clear();
        for (int i=0; i < count; i++) {
            owner.insert(owner.end(), f.elementFactor, reducedOwner[i]);
            for (int j=0; j < f.elementFactor*f.elementSize; j++)
                *cellIt++ = nodes[
                    nodesPerElement*i+f.multiCellIndices[j]];
        }
        nodes.swap(fullNodes);
        reducedCount = count;
        reducedNodesPerElement = f.reducedElementSize;
        nodesPerElement = f.elementSize;
        count *= f.elementFactor;
    } else {
        // we only keep the reduced elements but treat them as regular
        // ones, so replace the data accordingly
        nodes.swap(reducedNodes);
        reducedNodes.clear();
        nodesPerElement = f.reducedElementSize;
        type = f.reducedElementType;
    }
}

//
//
//
void ElementData::prepareGhostIndices(int ownIndex)
{
    indexArray.clear();
    reducedIndexArray.clear();
    numGhostElements = 0;
    numReducedGhostElements = 0;
    
    // move indices of "ghost zones" to the end to be able to reorder
    // data accordingly
    for (int i=0; i<count; i++)
        if (owner[i] == ownIndex)
            indexArray.push_back(i);

    for (int i=0; i<count; i++)
        if (owner[i] != ownIndex) {
            numGhostElements++;
            indexArray.push_back(i);
        }

    for (int i=0; i<reducedCount; i++)
        if (reducedOwner[i] == ownIndex)
            reducedIndexArray.push_back(i);

    for (int i=0; i<reducedCount; i++)
        if (reducedOwner[i] != ownIndex) {
            numReducedGhostElements++;
            reducedIndexArray.push_back(i);
        }
}

//
//
//
void ElementData::handleGhostZones(int ownIndex)
{
    prepareGhostIndices(ownIndex);

    // move "ghost data" to the end of the arrays
    if (numGhostElements > 0) {
        reorderArray(nodes, indexArray, nodesPerElement);
        reorderArray(owner, indexArray, 1);
        if (reducedCount == 0) {
            reorderArray(color, indexArray, 1);
            reorderArray(ID, indexArray, 1);
            reorderArray(tag, indexArray, 1);
        }
    }

    if (numReducedGhostElements > 0) {
        reorderArray(reducedNodes, reducedIndexArray, reducedNodesPerElement);
        reorderArray(reducedOwner, reducedIndexArray, 1);
        reorderArray(color, reducedIndexArray, 1);
        reorderArray(ID, reducedIndexArray, 1);
        reorderArray(tag, reducedIndexArray, 1);
    }
}

//
//
//
void ElementData::removeGhostZones()
{
    if (numGhostElements > 0) {
        count -= numGhostElements;
        nodes.resize(count*nodesPerElement);
        owner.resize(count);
        if (reducedCount == 0) {
            color.resize(count);
            ID.resize(count);
            tag.resize(count);
        }
        numGhostElements = 0;
    }

    if (numReducedGhostElements > 0) {
        reducedCount -= numReducedGhostElements;
        reducedNodes.resize(reducedCount*reducedNodesPerElement);
        reducedOwner.resize(reducedCount);
        color.resize(reducedCount);
        ID.resize(reducedCount);
        tag.resize(reducedCount);
        numReducedGhostElements = 0;
    }
}

//
//
//
void ElementData::buildMeshes()
{
    if (count == 0)
        return;

    // use existing original mesh for Elements but build a new mesh
    // containing only the required nodes for other types
    if (name == "Elements") {
        fullMesh = const_cast<Mesh*>(originalMesh);
        fullMeshIsOriginalMesh = true;
    } else {
        // first: build a map of required IDs while translating
        // the original node IDs at the same time
        IntVec::iterator it;
        IndexMap nodeID2idx;
        IntVec nodeIDs;
        size_t newIdx = 0;

        for (it = nodes.begin(); it != nodes.end(); it++) {
            IndexMap::iterator res = nodeID2idx.find(*it);
            if (res == nodeID2idx.end()) {
                nodeIDs.push_back(*it);
                nodeID2idx[*it] = newIdx;
                *it = newIdx++;
            } else {
                *it = res->second;
            }
        }

        // second: use map to fill coordinates for the nodes
        // newIdx contains the number of nodes
        CoordArray coords;
        const CoordArray& origCoords = originalMesh->getCoords();
        for (int dim=0; dim < originalMesh->getNumDims(); dim++) {
            float* c = new float[newIdx];
            coords.push_back(c);
            IndexMap::const_iterator mIt;
            for (mIt = nodeID2idx.begin(); mIt != nodeID2idx.end(); mIt++)
                c[mIt->second] = origCoords[dim][mIt->first];
        }

        if (fullMesh)
            delete fullMesh;

        fullMesh = new Mesh(coords, originalMesh->getNumDims(), newIdx);
        fullMesh->setName(name);
        fullMesh->setIndexMap(nodeID2idx);
        fullMesh->setNodeIDs(nodeIDs);
    }

#ifdef _DEBUG
    cout << fullMesh->getName() << " has " << fullMesh->getNumNodes()
        << " nodes, " << count << " elements" << endl;
#endif

    // build a reduced mesh if necessary
    if (reducedCount > 0) {
        IntVec::iterator it;
        IndexMap nodeID2idx;
        IntVec reducedNodeIDs;
        size_t newIdx = 0;

        for (it = reducedNodes.begin(); it != reducedNodes.end(); it++) {
            int id = originalMesh->getNodeIDs()[*it];
            IndexMap::iterator res = nodeID2idx.find(id);
            if (res == nodeID2idx.end()) {
                reducedNodeIDs.push_back(id);
                nodeID2idx[id] = newIdx;
                *it = newIdx++;
            } else {
                *it = res->second;
            }
        }

        CoordArray coords;
        const CoordArray& origCoords = originalMesh->getCoords();
        for (int dim=0; dim < originalMesh->getNumDims(); dim++) {
            float* c = new float[newIdx];
            coords.push_back(c);
            IndexMap::const_iterator mIt;
            for (mIt = nodeID2idx.begin(); mIt != nodeID2idx.end(); mIt++) {
                IndexMap::const_iterator idx;
                idx = originalMesh->getIndexMap().find(mIt->first);
                c[mIt->second] = origCoords[dim][idx->second];
            }
        }
        if (reducedMesh)
            delete reducedMesh;
        
        reducedMesh = new Mesh(coords, originalMesh->getNumDims(), newIdx);
        reducedMesh->setName(string("Reduced") + name);
        reducedMesh->setIndexMap(nodeID2idx);
        reducedMesh->setNodeIDs(reducedNodeIDs);

#ifdef _DEBUG
        cout << reducedMesh->getName() << " has " << newIdx
            << " nodes, " << reducedCount << " elements" << endl;
#endif
    }
}

#if HAVE_SILO
//
//
//
inline int toSiloElementType(int type)
{
    switch (type) {
        case ZONETYPE_BEAM: return DB_ZONETYPE_BEAM;
        case ZONETYPE_HEX: return DB_ZONETYPE_HEX;
        case ZONETYPE_POLYGON: return DB_ZONETYPE_POLYGON;
        case ZONETYPE_QUAD: return DB_ZONETYPE_QUAD;
        case ZONETYPE_TET: return DB_ZONETYPE_TET;
        case ZONETYPE_TRIANGLE: return DB_ZONETYPE_TRIANGLE;
    }
    return 0;
}
#endif

//
//
//
bool ElementData::writeToSilo(DBfile* dbfile, const string& siloPath)
{
#if HAVE_SILO
    if (count == 0)
        return true;

    const char* meshName;
    int numCells;
    string varName, siloMeshName;
    int ret;

    if (siloPath != "") {
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    }

    // write out the full mesh in any case
    if (siloPath == "/")
        siloMeshName = siloPath + fullMesh->getName();
    else
        siloMeshName = siloPath + string("/") + fullMesh->getName();

    int arraylen = count * nodesPerElement;
    int eltype = toSiloElementType(type);
    varName = name + string("_zones");
    ret = DBPutZonelist2(dbfile, varName.c_str(), count,
            originalMesh->getNumDims(), &nodes[0], arraylen, 0, 0,
            numGhostElements, &eltype, &nodesPerElement, &count, 1, NULL);
    if (ret == 0)
        ret = DBPutUcdmesh(dbfile, siloMeshName.c_str(),
                originalMesh->getNumDims(), NULL, &fullMesh->coords[0],
                fullMesh->getNumNodes(), count, varName.c_str(),
                /*"facelist"*/NULL, DB_FLOAT, NULL);
    
    // Point mesh is useful for debugging
    //DBPutPointmesh(dbfile, "/pointmesh",
    //        originalMesh->getNumDims(), &fullMesh->coords[0],
    //        fullMesh->getNumNodes(), DB_FLOAT, NULL);

    if (ret != 0)
        return false;

    // decide whether to additionally write out the reduced mesh
    if (reducedCount > 0) {
        if (siloPath == "/")
            siloMeshName = siloPath + reducedMesh->getName();
        else
            siloMeshName = siloPath + string("/") + reducedMesh->getName();
        arraylen = reducedCount * reducedNodesPerElement;
        eltype = toSiloElementType(reducedType);
        varName = string("Reduced") + name + string("_zones");
        ret = DBPutZonelist2(dbfile, varName.c_str(), reducedCount,
                originalMesh->getNumDims(), &reducedNodes[0], arraylen, 0, 0,
                numReducedGhostElements, &eltype, &reducedNodesPerElement,
                &reducedCount, 1, NULL);
        if (ret == 0)
            ret = DBPutUcdmesh(dbfile, siloMeshName.c_str(),
                   originalMesh->getNumDims(), NULL, &reducedMesh->coords[0],
                   reducedMesh->getNumNodes(), reducedCount, varName.c_str(),
                   NULL, DB_FLOAT, NULL);
        if (ret != 0)
            return false;
        numCells = reducedCount;
    } else {
        numCells = count;
    }
    meshName = siloMeshName.c_str();

    // finally, write out the element-centered variables on the correct mesh
    varName = name + string("_Color");
    ret = DBPutUcdvar1(dbfile, varName.c_str(), meshName,
            (float*)&color[0], numCells, NULL, 0, DB_INT, DB_ZONECENT, NULL);
    if (ret == 0) {
        varName = name + string("_Id");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), meshName,
            (float*)&ID[0], numCells, NULL, 0, DB_INT, DB_ZONECENT, NULL);
    }
    if (ret == 0) {
        varName = name + string("_Owner");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), meshName,
            (float*)&owner[0], numCells, NULL, 0, DB_INT, DB_ZONECENT, NULL);
    }
    if (ret == 0) {
        varName = name + string("_Tag");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), meshName,
            (float*)&tag[0], numCells, NULL, 0, DB_INT, DB_ZONECENT, NULL);
    }
    
    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !HAVE_SILO
    return false;
#endif
}

//
//
//
FinleyElementInfo ElementData::getFinleyTypeInfo(ElementTypeId typeId)
{
    FinleyElementInfo ret;
    ret.multiCellIndices = NULL;
    ret.elementFactor = 1;

    switch (typeId) {
        case Point1_Contact://untested
        case Line2Face_Contact://untested
        case Line3Face_Contact://untested
        case Line2Face://untested
        case Line3Face://untested
        case Point1://untested
            cerr << "WARNING: Finley type " <<typeId<< " is untested!" << endl;
            ret.elementSize = 1;
            ret.elementType = ZONETYPE_POLYGON;
            break;

        case Tri3Face_Contact://untested
        case Tri3Face://untested
            cerr << "WARNING: Finley type " <<typeId<< " is untested!" << endl;
        case Line2_Contact:
        case Rec4Face_Contact:
        case Rec4Face:
        case Line2:
            ret.elementSize = ret.reducedElementSize = 2;
            ret.elementType = ret.reducedElementType = ZONETYPE_BEAM;
            break;

        case Line3:
            ret.multiCellIndices = line3indices;
            ret.elementFactor = 2;
            // fall through
        case Line3_Contact:
        case Tri6Face_Contact://untested
        case Rec8Face_Contact:
        case Tri6Face://untested
        case Rec8Face:
            //VTK_QUADRATIC_EDGE
            ret.elementSize = ret.reducedElementSize = 2;
            ret.elementType = ret.reducedElementType = ZONETYPE_BEAM;
            break;

        case Tet4Face_Contact://untested
        case Tet4Face://untested
            cerr << "WARNING: Finley type " <<typeId<< " is untested!" << endl;
        case Tri3_Contact:
        case Tri3:
            ret.elementSize = ret.reducedElementSize = 3;
            ret.elementType = ret.reducedElementType = ZONETYPE_TRIANGLE;
            break;

        case Rec4_Contact:
        case Hex8Face_Contact:
        case Hex8Face:
        case Rec4:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_QUAD;
            break;

        case Rec9:
            ret.multiCellIndices = rec9indices;
            ret.elementFactor = 4;
            // fall through
        case Rec9_Contact:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_QUAD;
            break;

        case Tet4:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_TET;
            break;

        case Tri6:
            ret.multiCellIndices = tri6indices;
            ret.elementFactor = 4;
            // fall through
        case Tri6_Contact:
        case Tet10Face_Contact://untested
        case Tet10Face://untested
            //VTK_QUADRATIC_TRIANGLE
            ret.elementSize = ret.reducedElementSize = 3;
            ret.elementType = ret.reducedElementType = ZONETYPE_TRIANGLE;
            break;

        case Rec8:
            ret.multiCellIndices = rec8indices;
            ret.elementFactor = 6;
            // fall through
        case Hex20Face:
        case Rec8_Contact:
        case Hex20Face_Contact:
            //VTK_QUADRATIC_QUAD
            ret.elementSize = 3;
            ret.elementType = ZONETYPE_TRIANGLE; 
            ret.reducedElementSize = 4;
            ret.reducedElementType = ZONETYPE_QUAD;
            break;

        case Tet10:
            //VTK_QUADRATIC_TETRA
            ret.multiCellIndices = tet10indices;
            ret.elementFactor = 8;
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_TET;
            break;

        case Hex20:
            //VTK_QUADRATIC_HEXAHEDRON
            ret.multiCellIndices = hex20indices;
            ret.elementFactor = 36;
            ret.elementSize = 3;
            ret.elementType = ZONETYPE_TRIANGLE;
            ret.reducedElementSize = 8;
            ret.reducedElementType = ZONETYPE_HEX;
            break;

        case Hex27:
            ret.multiCellIndices = hex27indices;
            ret.elementFactor = 8;
            // fall through
        case Hex8:
            ret.elementSize = ret.reducedElementSize = 8;
            ret.elementType = ret.reducedElementType = ZONETYPE_HEX;
            break;

        default:
            cerr << "WARNING: Unknown Finley Type " << typeId << endl;
            break;
    }
    return ret;
}

