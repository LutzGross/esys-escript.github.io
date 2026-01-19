
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

#include <weipa/SpeckleyElements.h>
#include <weipa/NodeData.h>

#ifndef VISIT_PLUGIN
#include <speckley/SpeckleyDomain.h>
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
SpeckleyElements::SpeckleyElements(const string& elementName, SpeckleyNodes_ptr nodeData)
    : originalMesh(nodeData), name(elementName), numElements(0),
      numGhostElements(0), nodesPerElement(0),
      type(ZONETYPE_UNKNOWN)
{
    nodeMesh.reset(new SpeckleyNodes(name));
    numGhostElements = 0;
}

//
// Copy constructor
//
SpeckleyElements::SpeckleyElements(const SpeckleyElements& e)
{
    name = e.name;
    numElements = e.numElements;
    numGhostElements = 0;//e.numGhostElements;
    type = e.type;
    nodesPerElement = e.nodesPerElement;
    originalMesh = e.originalMesh;
    if (e.nodeMesh)
        nodeMesh.reset(new SpeckleyNodes(*e.nodeMesh));
    else
        nodeMesh.reset(new SpeckleyNodes(name));

    nodes = e.nodes;
    ID = e.ID;
    //tag = e.tag;
    owner = e.owner;
}

//
//
//
bool SpeckleyElements::initFromSpeckley(const speckley::SpeckleyDomain* dom, int fsType)
{
    if (fsType != speckley::Elements) {
        std::cerr << "Speckley only supports saving via Element functionspaces"
                << std::endl;
        return false;
    }

#ifndef VISIT_PLUGIN
    const pair<int,int> shape = dom->getDataShape(fsType);
    const dim_t* faces = dom->getNumFacesPerBoundary();
    const int* NS = dom->getNumSubdivisionsPerDim();
    const int order = dom->getOrder();

    numElements = shape.second*order*order;
    int nodesPerOrigElement = order*order;
    if (numElements > 0) {
        nodesPerElement = 4;
        if (dom->getDim() == 3) {
            nodesPerElement = 8;
            numElements *= order;
            nodesPerOrigElement *= order;
        }
        owner.assign(numElements, dom->getMPIRank());

        const dim_t* iPtr = dom->borrowSampleReferenceIDs(fsType);
        ID.resize(numElements);
        for (int i = 0; i < shape.second; i++) {
            for (int n = 0; n < nodesPerOrigElement; n++) {
                ID[i*n + n] = iPtr[i];
            }
        }

        const dim_t* NE = dom->getNumElementsPerDim();
        const dim_t* NN = dom->getNumNodesPerDim();
        nodes.clear();
        if (dom->getDim() == 2) {
            type = ZONETYPE_QUAD;
            if (faces[0]==0) {
                owner[0]=(faces[2]==0 ? dom->getMPIRank()-NS[0]-1
                        : dom->getMPIRank()-1);
                for (int i=1; i<NE[1]; i++)
                    owner[i*NE[0]]=dom->getMPIRank()-1;
            }
            if (faces[2]==0) {
                const int first=(faces[0]==0 ? 1 : 0);
                for (int i=first; i<NE[0]; i++)
                    owner[i]=dom->getMPIRank()-NS[0];
            }
            for (int ey = 0; ey < NE[1]; ey++) {
                for (int ex = 0; ex < NE[0]; ex++) {
                    int start = order*(ex + ey*NN[0]);
                    for (int qy=0; qy < order; qy++) {
                        int rowstart = start + qy*NN[0];
                        for (int qx=0; qx < order; qx++) {
                            nodes.push_back(rowstart + qx);
                            nodes.push_back(rowstart + qx + 1);
                            nodes.push_back(rowstart + qx + NN[0] + 1);
                            nodes.push_back(rowstart + qx + NN[0]);
                        }
                    }
                }
            }
        } else {
            type = ZONETYPE_HEX;
            // ownership is not entirely correct but that is not critical.
            // fix when there is time.
            if (faces[1]==0) {
                for (int k2=0; k2<NE[2]; k2++) {
                    for (int k1=0; k1<NE[1]; k1++) {
                        const int e=k2*NE[0]*NE[1]+(k1+1)*NE[0]-1;
                        owner[e]=dom->getMPIRank()+1;
                    }
                }
            }
            if (faces[3]==0) {
                for (int k2=0; k2<NE[2]; k2++) {
                    for (int k0=0; k0<NE[0]; k0++) {
                        const int e=(k2+1)*NE[0]*NE[1]-NE[0]+k0;
                        owner[e]=dom->getMPIRank()+NS[0];
                    }
                }
            }
            if (faces[5]==0) {
                for (int k1=0; k1<NE[1]; k1++) {
                    for (int k0=0; k0<NE[0]; k0++) {
                        const int e=k1*NE[0]+k0+NE[0]*NE[1]*(NE[2]-1);
                        owner[e]=dom->getMPIRank()+NS[0]*NS[1];
                    }
                }
            }
            
            for (int ez = 0; ez < NE[2]; ez++) {
                for (int ey = 0; ey < NE[1]; ey++) {
                    for (int ex = 0; ex < NE[0]; ex++) {
                        int start = order*(ex + ey*NN[0] + ez*NN[0]*NN[1]);
                        for (int qz = 0; qz < order; qz++) {
                            for (int qy=0; qy < order; qy++) {
                                for (int qx=0; qx < order; qx++) {
                                    int xstart = start + qy*NN[0] + qz*NN[0]*NN[1] + qx;
                                    nodes.push_back(xstart);
                                    nodes.push_back(xstart + NN[0]*NN[1]);
                                    nodes.push_back(xstart + NN[0]*NN[1] + 1);
                                    nodes.push_back(xstart + 1);
                                    
                                    nodes.push_back(xstart + NN[0]);
                                    nodes.push_back(xstart + NN[0]*(NN[1]+1));
                                    nodes.push_back(xstart + NN[0]*(NN[1]+1)+1);
                                    nodes.push_back(xstart + NN[0]+1);
                                }
                            }
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

StringVec SpeckleyElements::getMeshNames() const
{
    StringVec res;
    if (nodeMesh)
        res.push_back(nodeMesh->getName());
    return res;
}

StringVec SpeckleyElements::getVarNames() const
{
    StringVec res;
    res.push_back(name + string("_Id"));
    res.push_back(name + string("_Owner"));
    //res.push_back(name + string("_Tag"));
    return res;
}

const IntVec& SpeckleyElements::getVarDataByName(const string varName) const
{
    if (varName == name+string("_Id"))
        return ID;
    if (varName == name+string("_Owner"))
        return owner;
    //if (varName == name+string("_Tag"))
    //    return tag;

    throw "Invalid variable name";
}

void SpeckleyElements::reorderArray(IntVec& v, const IntVec& idx,
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

IntVec SpeckleyElements::prepareGhostIndices(int ownIndex)
{
    IntVec indexArray;
    numGhostElements = 0;
    for (int i=0; i<numElements; i++) {
        indexArray.push_back(i);
    }
    return indexArray;
}

void SpeckleyElements::reorderGhostZones(int ownIndex)
{
    return;
}

void SpeckleyElements::removeGhostZones(int ownIndex)
{
    return;
}

void SpeckleyElements::buildMeshes()
{
    // build a new mesh containing only the required nodes
    if (numElements > 0) {
        if (nodeMesh && nodeMesh->getNumNodes() > 0) {
            SpeckleyNodes_ptr newMesh(new SpeckleyNodes(nodeMesh, nodes, name));
            nodeMesh.swap(newMesh);
        } else {
            nodeMesh.reset(new SpeckleyNodes(originalMesh, nodes, name));
        }
#ifdef _DEBUG
        cout << nodeMesh->getName() << " has " << nodeMesh->getNumNodes()
            << " nodes and " << numElements << " elements" << endl;
#endif
    }
}

void SpeckleyElements::writeConnectivityVTK(ostream& os)
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

bool SpeckleyElements::writeToSilo(DBfile* dbfile, const string& siloPath,
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

