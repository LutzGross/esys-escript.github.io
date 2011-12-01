
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <weipa/RipleyElements.h>
#include <weipa/NodeData.h>

#ifndef VISIT_PLUGIN
#include <ripley/RipleyDomain.h>
#endif

#include <iostream>

#if USE_SILO
#include <silo.h>
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

    NperDim = e.NperDim;
    nodes = e.nodes;
    ID = e.ID;
    tag = e.tag;
    owner = e.owner;
}

//
//
//
bool RipleyElements::initFromRipley(const ripley::RipleyDomain* ripleyDomain,
                                    int fsType)
{
#ifndef VISIT_PLUGIN
    pair<int,int> shape = ripleyDomain->getDataShape(fsType);
    numElements = shape.second;
    if (fsType==ripley::Elements)
        NperDim = ripleyDomain->getNumElementsPerDim();
    else {
        NperDim = ripleyDomain->getNumFacesPerBoundary();
        numElements=0; // ignore faces for now
    }

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
        owner.assign(numElements, ripleyDomain->getMPIRank());

        const int* iPtr = ripleyDomain->borrowSampleReferenceIDs(fsType);
        ID.assign(iPtr, iPtr+numElements);

        //iPtr = ripleyDomain->borrowListOfTagsInUse(fsType);
        tag.assign(iPtr, iPtr+numElements);

        IntVec NN = ripleyDomain->getNumNodesPerDim();
        nodes.clear();
        if (ripleyDomain->getDim() == 2) {
            int id=0;
            for (int i=0; i<numElements; i++) {
                nodes.push_back(id);
                nodes.push_back(id+1);
                nodes.push_back(id+1+NN[0]);
                nodes.push_back(id+NN[0]);
                id++;
                if ((i+1)%NperDim[0]==0)
                    id++;
            }
        } else {
            int id=0;
            for (int i=0; i<numElements; i++) {
                nodes.push_back(id+NN[0]*NN[1]);
                nodes.push_back(id);
                nodes.push_back(id+1);
                nodes.push_back(id+NN[0]*NN[1]+1);

                nodes.push_back(id+NN[0]*(NN[1]+1));
                nodes.push_back(id+NN[0]);
                nodes.push_back(id+1+NN[0]);
                nodes.push_back(id+NN[0]*(NN[1]+1)+1);
                id++;
                if ((i+1)%NperDim[0]==0)
                    id++;
                if ((i+1)%(NperDim[0]*NperDim[1])==0)
                    id+=NN[0];
            }
        }

        nodeMesh=originalMesh;
        //buildMeshes();
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
    res.push_back(name + string("_Tag"));
    return res;
}

const IntVec& RipleyElements::getVarDataByName(const string varName) const
{
    if (varName == name+string("_Id"))
        return ID;
    if (varName == name+string("_Owner"))
        return owner;
    if (varName == name+string("_Tag"))
        return tag;

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
            copy(&v[i*elementsPerIndex], &v[(i+1)*elementsPerIndex], arrIt);
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
    for (int i=0; i<numElements; i++) {
        if (owner[i] == ownIndex)
            indexArray.push_back(i);
    }

    for (int i=0; i<numElements; i++) {
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
        reorderArray(tag, indexArray, 1);
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
        tag.resize(numElements);
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

#if USE_SILO
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
#if USE_SILO
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

        //ret = DBPutQuadmesh(dbfile, siloMeshName, NULL, &coordbase[0],
        //&NperDim[0], nodeMesh->getNumDims(), DB_FLOAT, DB_COLLINEAR, optList);
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
        //ret = DBPutQuadvar1(dbfile, varName.c_str(), siloMeshName,
        //    (float*)&ID[0], &NperDim[0], nodeMesh->getNumDims(), NULL, 0,
        //    DB_INT, DB_ZONECENT, NULL);
        if (ret == 0) {
            varName = name + string("_Owner");
            ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
                (float*)&owner[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
                NULL);
            //ret = DBPutQuadvar1(dbfile, varName.c_str(), siloMeshName,
            //    (float*)&owner[0], &NperDim[0], nodeMesh->getNumDims(), NULL,
            //    0, DB_INT, DB_ZONECENT, NULL);
        }
        //if (ret == 0) {
        //    varName = name + string("_Tag");
        //    ret = DBPutQuadvar1(dbfile, varName.c_str(), siloMeshName,
        //        (float*)&tag[0], dims, nodeMesh->getNumDims(), NULL, 0,
        //        DB_INT, DB_ZONECENT, NULL);
        //}
    }

    // "Elements" is a special case
    if (writeMeshData && name == "Elements") {
        nodeMesh->writeToSilo(dbfile);
    }

    return (ret == 0);

#else // !USE_SILO
    return false;
#endif
}

} // namespace weipa

