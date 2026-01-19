
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

#include <weipa/RipleyElements.h>
#include <weipa/NodeData.h>

#ifndef VISIT_PLUGIN
#include <ripley/RipleyDomain.h>
#endif

#include <iostream>

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
RipleyElements::RipleyElements(const string& elementName, RipleyNodes_ptr nodeData)
    : originalMesh(nodeData), name(elementName), numElements(0),
      numGhostElements(0), nodesPerElement(0),
      type(ZONETYPE_UNKNOWN)
{
    nodeMesh.reset(new RipleyNodes(name));
}

//
// Copy constructor
//
RipleyElements::RipleyElements(const RipleyElements& e)
{
    name = e.name;
    numElements = e.numElements;
    numGhostElements = e.numGhostElements;
    type = e.type;
    nodesPerElement = e.nodesPerElement;
    originalMesh = e.originalMesh;
    if (e.nodeMesh)
        nodeMesh.reset(new RipleyNodes(*e.nodeMesh));
    else
        nodeMesh.reset(new RipleyNodes(name));

    nodes = e.nodes;
    ID = e.ID;
    //tag = e.tag;
    owner = e.owner;
}

//
//
//
bool RipleyElements::initFromRipley(const ripley::RipleyDomain* dom, int fsType)
{
#ifndef VISIT_PLUGIN
    const pair<int,dim_t> shape = dom->getDataShape(fsType);
    const dim_t* faces = dom->getNumFacesPerBoundary();

    numElements = shape.second;

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
        }
        owner = dom->getOwnerVector(fsType);

        const dim_t* iPtr = dom->borrowSampleReferenceIDs(fsType);
        ID.assign(iPtr, iPtr+numElements);

        //iPtr = dom->borrowListOfTags(fsType);
        //tag.assign(iPtr, iPtr+numElements);

        const dim_t* NE = dom->getNumElementsPerDim();
        const dim_t* NN = dom->getNumNodesPerDim();
        nodes.clear();
        if (dom->getDim() == 2) {
            if (fsType==ripley::Elements) {
                int id=0;
                for (dim_t i=0; i<numElements; i++) {
                    nodes.push_back(id);
                    nodes.push_back(id+1);
                    nodes.push_back(id+1+NN[0]);
                    nodes.push_back(id+NN[0]);
                    id++;
                    if ((i+1)%NE[0]==0)
                        id++;
                }
            } else if (fsType==ripley::FaceElements) {
                int id=0;
                for (dim_t i=0; i<faces[0]; i++) {
                    nodes.push_back(id);
                    nodes.push_back(id+NN[0]);
                    id+=NN[0];
                }
                id=NN[0]-1;
                for (dim_t i=0; i<faces[1]; i++) {
                    nodes.push_back(id);
                    nodes.push_back(id+NN[0]);
                    id+=NN[0];
                }
                id=0;
                for (dim_t i=0; i<faces[2]; i++) {
                    nodes.push_back(id);
                    nodes.push_back(id+1);
                    id++;
                }
                id=NN[0]*(NN[1]-1);
                for (dim_t i=0; i<faces[3]; i++) {
                    nodes.push_back(id);
                    nodes.push_back(id+1);
                    id++;
                }
            }
        } else {
            if (fsType==ripley::Elements) {
                int id=0;
                for (dim_t i=0; i<numElements; i++) {
                    nodes.push_back(id);
                    nodes.push_back(id+NN[0]*NN[1]);
                    nodes.push_back(id+NN[0]*NN[1]+1);
                    nodes.push_back(id+1);

                    nodes.push_back(id+NN[0]);
                    nodes.push_back(id+NN[0]*(NN[1]+1));
                    nodes.push_back(id+NN[0]*(NN[1]+1)+1);
                    nodes.push_back(id+1+NN[0]);
                    id++;
                    if ((i+1)%NE[0]==0)
                        id++;
                    if ((i+1)%(NE[0]*NE[1])==0)
                        id+=NN[0];
                }
            } else if (fsType==ripley::FaceElements) {
                const dim_t* NE = dom->getNumElementsPerDim();
                if (faces[0]>0) {
                    for (dim_t k2=0; k2<NE[2]; k2++) {
                        for (dim_t k1=0; k1<NE[1]; k1++) {
                            const dim_t first=k2*NN[0]*NN[1]+k1*NN[0];
                            nodes.push_back(first+NN[0]*NN[1]);
                            nodes.push_back(first);
                            nodes.push_back(first+NN[0]);
                            nodes.push_back(first+NN[0]*(NN[1]+1));
                        }
                    }
                }
                if (faces[1]>0) {
                    for (dim_t k2=0; k2<NE[2]; k2++) {
                        for (dim_t k1=0; k1<NE[1]; k1++) {
                            const dim_t first=k2*NN[0]*NN[1]+(k1+1)*NN[0]-1;
                            nodes.push_back(first+NN[0]*NN[1]);
                            nodes.push_back(first);
                            nodes.push_back(first+NN[0]);
                            nodes.push_back(first+NN[0]*(NN[1]+1));
                        }
                    }
                }
                if (faces[2]>0) {
                    for (dim_t k2=0; k2<NE[2]; k2++) {
                        for (dim_t k0=0; k0<NE[0]; k0++) {
                            const dim_t first=k2*NN[0]*NN[1]+k0;
                            nodes.push_back(first+NN[0]*NN[1]);
                            nodes.push_back(first);
                            nodes.push_back(first+1);
                            nodes.push_back(first+1+NN[0]*NN[1]);
                        }
                    }
                }
                if (faces[3]>0) {
                    for (dim_t k2=0; k2<NE[2]; k2++) {
                        for (dim_t k0=0; k0<NE[0]; k0++) {
                            const dim_t first=(k2+1)*NN[0]*NN[1]-NN[0]+k0;
                            nodes.push_back(first+NN[0]*NN[1]);
                            nodes.push_back(first);
                            nodes.push_back(first+1);
                            nodes.push_back(first+1+NN[0]*NN[1]);
                        }
                    }
                }
                if (faces[4]>0) {
                    for (dim_t k1=0; k1<NE[1]; k1++) {
                        for (dim_t k0=0; k0<NE[0]; k0++) {
                            const dim_t first=k1*NN[0]+k0;
                            nodes.push_back(first);
                            nodes.push_back(first+1);
                            nodes.push_back(first+NN[0]+1);
                            nodes.push_back(first+NN[0]);
                        }
                    }
                }
                if (faces[5]>0) {
                    for (dim_t k1=0; k1<NE[1]; k1++) {
                        for (dim_t k0=0; k0<NE[0]; k0++) {
                            const dim_t first=NN[0]*NN[1]*(NN[2]-1)+k1*NN[0]+k0;
                            nodes.push_back(first);
                            nodes.push_back(first+1);
                            nodes.push_back(first+NN[0]+1);
                            nodes.push_back(first+NN[0]);
                        }
                    }
                }
            }
        }

        buildMeshes();
    }
    return true;

#else // VISIT_PLUGIN
    return false;
#endif
}

StringVec RipleyElements::getMeshNames() const
{
    StringVec res;
    if (nodeMesh)
        res.push_back(nodeMesh->getName());
    return res;
}

StringVec RipleyElements::getVarNames() const
{
    StringVec res;
    res.push_back(name + string("_Id"));
    res.push_back(name + string("_Owner"));
    //res.push_back(name + string("_Tag"));
    return res;
}

const IntVec& RipleyElements::getVarDataByName(const string varName) const
{
    if (varName == name+string("_Id"))
        return ID;
    if (varName == name+string("_Owner"))
        return owner;
    //if (varName == name+string("_Tag"))
    //    return tag;

    throw "Invalid variable name";
}

void RipleyElements::reorderArray(IntVec& v, const IntVec& idx,
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

IntVec RipleyElements::prepareGhostIndices(int ownIndex)
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

void RipleyElements::reorderGhostZones(int ownIndex)
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

void RipleyElements::removeGhostZones(int ownIndex)
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

void RipleyElements::buildMeshes()
{
    // build a new mesh containing only the required nodes
    if (numElements > 0) {
        if (nodeMesh && nodeMesh->getNumNodes() > 0) {
            RipleyNodes_ptr newMesh(new RipleyNodes(nodeMesh, nodes, name));
            nodeMesh.swap(newMesh);
        } else {
            nodeMesh.reset(new RipleyNodes(originalMesh, nodes, name));
        }
#ifdef _DEBUG
        cout << nodeMesh->getName() << " has " << nodeMesh->getNumNodes()
            << " nodes and " << numElements << " elements" << endl;
#endif
    }
}

void RipleyElements::writeConnectivityVTK(ostream& os)
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

bool RipleyElements::writeToSilo(DBfile* dbfile, const string& siloPath,
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

