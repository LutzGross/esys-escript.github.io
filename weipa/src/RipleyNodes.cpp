
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

#include <weipa/RipleyNodes.h>

#ifndef VISIT_PLUGIN
#include <ripley/RipleyDomain.h>
#endif

#if USE_SILO
#include <silo.h>
#endif

using namespace std;

namespace weipa {

//
// Constructor with name
//
RipleyNodes::RipleyNodes(const string& meshName) :
    numDims(0), numNodes(0), globalNumNodes(0), name(meshName)
{
}

//
//
//
RipleyNodes::RipleyNodes(RipleyNodes_ptr fullNodes, IntVec& requiredNodes,
                   const string& meshName) :
    name(meshName)
{
    numDims = fullNodes->numDims;
    nodeDist = fullNodes->nodeDist;
    globalNumNodes = fullNodes->globalNumNodes;

    // first: find the unique set of required nodes and their IDs while
    // updating the contents of requiredNodes at the same time
    // requiredNodes contains node indices (not IDs!)
    IntVec::iterator it;
    IndexMap indexMap; // maps old index to new index
    size_t newIndex = 0;

    for (it = requiredNodes.begin(); it != requiredNodes.end(); it++) {
        IndexMap::iterator res = indexMap.find(*it);
        if (res == indexMap.end()) {
            nodeID.push_back(fullNodes->nodeID[*it]);
            nodeTag.push_back(fullNodes->nodeTag[*it]);
            indexMap[*it] = newIndex;
            *it = newIndex++;
        } else {
            *it = res->second;
        }
    }

    // second: now that we know how many nodes we need use the map to fill
    // the coordinates
    numNodes = newIndex;
    for (int dim=0; dim<numDims; dim++) {
        const float* origC = fullNodes->coords[dim];
        float* c = new float[numNodes];
        coords.push_back(c);
        IndexMap::const_iterator mIt;
        for (mIt = indexMap.begin(); mIt != indexMap.end(); mIt++) {
            c[mIt->second] = origC[mIt->first];
        }
    }
}

//
// Copy constructor
//
RipleyNodes::RipleyNodes(const RipleyNodes& m)
{
    numDims = m.numDims;
    numNodes = m.numNodes;
    globalNumNodes = m.globalNumNodes;
    nodeID = m.nodeID;
    nodeTag = m.nodeTag;
    nodeDist = m.nodeDist;
    name = m.name;
    for (int i=0; i<numDims; i++) {
        float* c = new float[numNodes];
        copy(m.coords[i], m.coords[i]+numNodes, c);
        coords.push_back(c);
    }
}

//
//
//
RipleyNodes::~RipleyNodes()
{
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
}

//
//
//
bool RipleyNodes::initFromRipley(const ripley::RipleyDomain* dom)
{
#ifndef VISIT_PLUGIN
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();
    nodeID.clear();
    nodeTag.clear();

    numDims = dom->getDim();
    globalNumNodes = dom->getNumDataPointsGlobal();
    pair<int,int> shape = dom->getDataShape(ripley::Nodes);
    numNodes = shape.second;
    nodeDist = dom->getNodeDistribution();

    if (numNodes > 0) {
        for (int d=0; d<numDims; d++) {
            float* c = new float[numNodes];
            coords.push_back(c);
        }
        IntVec NN = dom->getNumNodesPerDim();
        if (numDims==2) {
            pair<double,double> xx=dom->getFirstCoordAndSpacing(0);
            pair<double,double> yy=dom->getFirstCoordAndSpacing(1);
            for (int i1=0; i1<NN[1]; i1++) {
                for (int i0=0; i0<NN[0]; i0++) {
                    coords[0][i0+NN[0]*i1] = (float)(xx.first+i0*xx.second);
                    coords[1][i0+NN[0]*i1] = (float)(yy.first+i1*yy.second);
                }
            }
        } else {
            pair<double,double> xx=dom->getFirstCoordAndSpacing(0);
            pair<double,double> yy=dom->getFirstCoordAndSpacing(1);
            pair<double,double> zz=dom->getFirstCoordAndSpacing(2);
            for (int i2=0; i2<NN[2]; i2++) {
                for (int i1=0; i1<NN[1]; i1++) {
                    for (int i0=0; i0<NN[0]; i0++) {
                        const int index = i0+NN[0]*i1+NN[0]*NN[1]*i2;
                        coords[0][index] = (float)(xx.first+i0*xx.second);
                        coords[1][index] = (float)(yy.first+i1*yy.second);
                        coords[2][index] = (float)(zz.first+i2*zz.second);
                    }
                }
            }
        }

        const int* iPtr = dom->borrowSampleReferenceIDs(ripley::Nodes);
        nodeID.assign(iPtr, iPtr+numNodes);

        //iPtr = dom->borrowListOfTagsInUse(ripley::Nodes);
        nodeTag.assign(iPtr, iPtr+numNodes);
    }

    return true;
#else // VISIT_PLUGIN
    return false;
#endif
}

//
//
//
const IntVec& RipleyNodes::getVarDataByName(const string& name) const
{
    if (name == "Nodes_Id")
        return nodeID;
    if (name == "Nodes_Tag")
        return nodeTag;
    
    throw "Invalid variable name";
}

//
//
//
StringVec RipleyNodes::getVarNames() const
{
    StringVec res;
    res.push_back("Nodes_Id");
    res.push_back("Nodes_Tag");
    return res;
}

//
//
//
void RipleyNodes::writeCoordinatesVTK(ostream& os, int ownIndex)
{
    if (numNodes > 0) {
        int firstId = nodeDist[ownIndex];
        int lastId = nodeDist[ownIndex+1];
        for (size_t i=0; i<numNodes; i++) {
            if (firstId <= nodeID[i] && nodeID[i] < lastId) {
                os << coords[0][i] << " " << coords[1][i] << " ";
                if (numDims == 3)
                    os << coords[2][i];
                else
                    os << 0.;
                os << endl;
            }
        }
    }
}

//
//
//
bool RipleyNodes::writeToSilo(DBfile* dbfile)
{
#if USE_SILO
    if (numNodes == 0)
        return true;

    int ret;

    if (siloPath != "") {
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    }
    string siloMeshName = getFullSiloName();

    // Write node-centered variables
    ret = DBPutUcdvar1(dbfile, "Nodes_Id", siloMeshName.c_str(),
            (float*)&nodeID[0], numNodes, NULL, 0, DB_INT, DB_NODECENT, NULL);

    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_Tag", siloMeshName.c_str(),
                (float*)&nodeTag[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);

    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !USE_SILO
    return false;
#endif
}

} // namespace weipa

