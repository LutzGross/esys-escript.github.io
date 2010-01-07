
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <escriptexport/NodeData.h>

extern "C" {
#include <finley/Mesh.h>
#include <finley/NodeFile.h>
}

#if USE_NETCDF
#include <netcdf.hh>
#endif

#if USE_SILO
#include <silo.h>
#endif

using namespace std;

namespace escriptexport {

//
// Constructor with name
//
NodeData::NodeData(const string& meshName) :
    numDims(0), numNodes(0), name(meshName)
{
}

//
//
//
NodeData::NodeData(NodeData_ptr fullNodes, IntVec& requiredNodes,
                   const string& meshName) :
    name(meshName)
{
    numDims = fullNodes->numDims;
    nodeDist = fullNodes->nodeDist;

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
            nodeGDOF.push_back(fullNodes->nodeGDOF[*it]);
            nodeGNI.push_back(fullNodes->nodeGNI[*it]);
            nodeGRDFI.push_back(fullNodes->nodeGRDFI[*it]);
            nodeGRNI.push_back(fullNodes->nodeGRNI[*it]);
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
NodeData::NodeData(const NodeData& m)
{
    numDims = m.numDims;
    numNodes = m.numNodes;
    nodeID = m.nodeID;
    nodeTag = m.nodeTag;
    nodeGDOF = m.nodeGDOF;
    nodeGNI = m.nodeGNI;
    nodeGRDFI = m.nodeGRDFI;
    nodeGRNI = m.nodeGRNI;
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
NodeData::~NodeData()
{
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
}

//
//
//
bool NodeData::initFromFinley(const Finley_NodeFile* finleyFile)
{
    numDims = finleyFile->numDim;
    numNodes = finleyFile->numNodes;

    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();

    if (numNodes > 0) {
        for (int i=0; i<numDims; i++) {
            double* srcPtr = finleyFile->Coordinates + i;
            float* c = new float[numNodes];
            coords.push_back(c);
            for (int j=0; j<numNodes; j++, srcPtr+=numDims) {
                *c++ = (float) *srcPtr;
            }
        }

        int* iPtr;
 
        iPtr = finleyFile->Id;
        nodeID.clear();
        nodeID.insert(nodeID.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeID.begin());

        iPtr = finleyFile->Tag;
        nodeTag.clear();
        nodeTag.insert(nodeTag.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeTag.begin());

        iPtr = finleyFile->globalDegreesOfFreedom;
        nodeGDOF.clear();
        nodeGDOF.insert(nodeGDOF.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGDOF.begin());

        iPtr = finleyFile->globalNodesIndex;
        nodeGNI.clear();
        nodeGNI.insert(nodeGNI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGNI.begin());

        iPtr = finleyFile->globalReducedDOFIndex;
        nodeGRDFI.clear();
        nodeGRDFI.insert(nodeGRDFI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGRDFI.begin());

        iPtr = finleyFile->globalReducedNodesIndex;
        nodeGRNI.clear();
        nodeGRNI.insert(nodeGRNI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGRNI.begin());

        int mpisize = finleyFile->MPIInfo->size;
        iPtr = finleyFile->nodesDistribution->first_component;
        nodeDist.clear();
        nodeDist.insert(nodeDist.end(), mpisize+1, 0);
        copy(iPtr, iPtr+mpisize+1, nodeDist.begin());
    }
    return true;
}

//
//
//
bool NodeData::readFromNc(NcFile* ncFile)
{
#if USE_NETCDF
    NcAtt* att;
    NcVar* var;
 
    att = ncFile->get_att("numDim");
    numDims = att->as_int(0);

    att = ncFile->get_att("numNodes");
    numNodes = att->as_int(0);

    att = ncFile->get_att("mpi_size");
    int mpisize = att->as_int(0);

    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();
    var = ncFile->get_var("Nodes_Coordinates");
    for (int i=0; i<numDims; i++) {
        float* c = new float[numNodes];
        var->set_cur(0, i);
        var->get(c, numNodes, 1);
        coords.push_back(c);
    }

    nodeID.clear();
    nodeID.insert(nodeID.end(), numNodes, 0);
    var = ncFile->get_var("Nodes_Id");
    var->get(&nodeID[0], numNodes);

    nodeTag.clear();
    nodeTag.insert(nodeTag.end(), numNodes, 0);
    var = ncFile->get_var("Nodes_Tag");
    var->get(&nodeTag[0], numNodes);

    nodeGDOF.clear();
    nodeGDOF.insert(nodeGDOF.end(), numNodes, 0);
    var = ncFile->get_var("Nodes_gDOF");
    var->get(&nodeGDOF[0], numNodes);

    nodeGNI.clear();
    nodeGNI.insert(nodeGNI.end(), numNodes, 0);
    var = ncFile->get_var("Nodes_gNI");
    var->get(&nodeGNI[0], numNodes);

    nodeGRDFI.clear();
    nodeGRDFI.insert(nodeGRDFI.end(), numNodes, 0);
    var = ncFile->get_var("Nodes_grDfI");
    var->get(&nodeGRDFI[0], numNodes);

    nodeGRNI.clear();
    nodeGRNI.insert(nodeGRNI.end(), numNodes, 0);
    var = ncFile->get_var("Nodes_grNI");
    var->get(&nodeGRNI[0], numNodes);

    nodeDist.clear();
    nodeDist.insert(nodeDist.end(), mpisize+1, 0);
    var = ncFile->get_var("Nodes_NodeDistribution");
    var->get(&nodeDist[0], mpisize+1);

    return true;
#else // !USE_NETCDF
    return false;
#endif
}

//
//
//
const IntVec& NodeData::getVarDataByName(const string& name) const
{
    if (name == "Nodes_Id")
        return nodeID;
    else if (name == "Nodes_Tag")
        return nodeTag;
    else if (name == "Nodes_gDOF")
        return nodeGDOF;
    else if (name == "Nodes_gNI")
        return nodeGNI;
    else if (name == "Nodes_grDfI")
        return nodeGRDFI;
    else if (name == "Nodes_grNI")
        return nodeGRNI;
    else
        throw "Invalid variable name";
}

//
//
//
StringVec NodeData::getVarNames() const
{
    StringVec res;
    if (numNodes > 0) {
        res.push_back("Nodes_Id");
        res.push_back("Nodes_Tag");
        res.push_back("Nodes_gDOF");
        res.push_back("Nodes_gNI");
        res.push_back("Nodes_grDfI");
        res.push_back("Nodes_grNI");
    }
    return res;
}

//
//
//
void NodeData::removeGhostNodes(int ownIndex)
{
    if (nodeDist.empty() || ownIndex > nodeDist.size()-1)
        return;

    int firstId = nodeDist[ownIndex];
    int lastId = nodeDist[ownIndex+1];

    // no ghost nodes
    if (lastId-firstId == numNodes)
        return;

    // we have at most lastId-firstId nodes, it could be less however if
    // nodes were culled already
    numNodes = lastId-firstId;

    CoordArray newCoords;
    CoordArray::iterator it;
    for (int i=0; i<numDims; i++) {
        float* c = new float[numNodes];
        newCoords.push_back(c);
    }

    IntVec newNodeID, newNodeTag;
    IntVec newNodeGDOF, newNodeGNI, newNodeGRDFI, newNodeGRNI;

    IndexMap nodeID2idx = getIndexMap();
    int destIdx = 0;
    for (int i=firstId; i<lastId; i++) {
        IndexMap::iterator it = nodeID2idx.find(i);
        if (it == nodeID2idx.end()) {
            continue;
        }
        int idx = it->second;
        for (int dim=0; dim<numDims; dim++) {
            newCoords[dim][destIdx] = coords[dim][idx];
        }
        destIdx++;
        newNodeID.push_back(i);
        newNodeTag.push_back(nodeTag[idx]);
        newNodeGDOF.push_back(nodeGDOF[idx]);
        newNodeGNI.push_back(nodeGNI[idx]);
        newNodeGRDFI.push_back(nodeGRDFI[idx]);
        newNodeGRNI.push_back(nodeGRNI[idx]);
    }

    numNodes = destIdx;

    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;

    coords = newCoords;
    nodeID = newNodeID;
    nodeTag = newNodeTag;
    nodeGDOF = newNodeGDOF;
    nodeGNI = newNodeGNI;
    nodeGRDFI = newNodeGRDFI;
    nodeGRNI = newNodeGRNI;
}

//
//
//
bool NodeData::writeToSilo(DBfile* dbfile)
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
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_gDOF", siloMeshName.c_str(),
                (float*)&nodeGDOF[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_gNI", siloMeshName.c_str(),
                (float*)&nodeGNI[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_grDfI", siloMeshName.c_str(),
                (float*)&nodeGRDFI[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_grNI", siloMeshName.c_str(),
                (float*)&nodeGRNI[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);

    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !USE_SILO
    return false;
#endif
}

} // namespace escriptexport

