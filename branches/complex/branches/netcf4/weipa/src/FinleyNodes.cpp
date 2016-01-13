
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#include <weipa/FinleyNodes.h>

#ifndef VISIT_PLUGIN
#include <dudley/Mesh.h>
#include <dudley/NodeFile.h>
#include <finley/Mesh.h>
#include <finley/NodeFile.h>
#endif

#if USE_NETCDF
#include "esysUtils/netcdf.h"
#endif

#if USE_SILO
#include <silo.h>
#endif

using namespace std;

namespace weipa {

//
// Constructor with name
//
FinleyNodes::FinleyNodes(const string& meshName) :
    numDims(0), numNodes(0), name(meshName)
{
}

//
//
//
FinleyNodes::FinleyNodes(FinleyNodes_ptr fullNodes, IntVec& requiredNodes,
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
FinleyNodes::FinleyNodes(const FinleyNodes& m)
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
FinleyNodes::~FinleyNodes()
{
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
}

//
//
//
bool FinleyNodes::initFromDudley(const Dudley_NodeFile* dudleyFile)
{
#ifndef VISIT_PLUGIN
    numDims = dudleyFile->numDim;
    numNodes = dudleyFile->numNodes;

    int mpisize = dudleyFile->MPIInfo->size;
    int* iPtr = dudleyFile->nodesDistribution->first_component;
    nodeDist.clear();
    nodeDist.insert(nodeDist.end(), mpisize+1, 0);
    copy(iPtr, iPtr+mpisize+1, nodeDist.begin());

    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();
    nodeID.clear();
    nodeTag.clear();
    nodeGDOF.clear();
    nodeGNI.clear();
    nodeGRDFI.clear();
    nodeGRNI.clear();

    if (numNodes > 0) {
        for (int i=0; i<numDims; i++) {
            double* srcPtr = dudleyFile->Coordinates + i;
            float* c = new float[numNodes];
            coords.push_back(c);
            for (int j=0; j<numNodes; j++, srcPtr+=numDims) {
                *c++ = (float) *srcPtr;
            }
        }

        iPtr = dudleyFile->Id;
        nodeID.insert(nodeID.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeID.begin());

        iPtr = dudleyFile->Tag;
        nodeTag.insert(nodeTag.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeTag.begin());

        iPtr = dudleyFile->globalDegreesOfFreedom;
        nodeGDOF.insert(nodeGDOF.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGDOF.begin());

        iPtr = dudleyFile->globalNodesIndex;
        nodeGNI.insert(nodeGNI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGNI.begin());

        iPtr = dudleyFile->globalReducedDOFIndex;
        nodeGRDFI.insert(nodeGRDFI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGRDFI.begin());

        iPtr = dudleyFile->globalReducedNodesIndex;
        nodeGRNI.insert(nodeGRNI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGRNI.begin());

    }
    return true;
#else // VISIT_PLUGIN
    return false;
#endif
}

//
//
//
bool FinleyNodes::initFromFinley(const Finley_NodeFile* finleyFile)
{
#ifndef VISIT_PLUGIN
    numDims = finleyFile->numDim;
    numNodes = finleyFile->numNodes;

    int mpisize = finleyFile->MPIInfo->size;
    int* iPtr = finleyFile->nodesDistribution->first_component;
    nodeDist.clear();
    nodeDist.insert(nodeDist.end(), mpisize+1, 0);
    copy(iPtr, iPtr+mpisize+1, nodeDist.begin());

    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();
    nodeID.clear();
    nodeTag.clear();
    nodeGDOF.clear();
    nodeGNI.clear();
    nodeGRDFI.clear();
    nodeGRNI.clear();

    if (numNodes > 0) {
        for (int i=0; i<numDims; i++) {
            double* srcPtr = finleyFile->Coordinates + i;
            float* c = new float[numNodes];
            coords.push_back(c);
            for (int j=0; j<numNodes; j++, srcPtr+=numDims) {
                *c++ = (float) *srcPtr;
            }
        }

        iPtr = finleyFile->Id;
        nodeID.insert(nodeID.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeID.begin());

        iPtr = finleyFile->Tag;
        nodeTag.insert(nodeTag.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeTag.begin());

        iPtr = finleyFile->globalDegreesOfFreedom;
        nodeGDOF.insert(nodeGDOF.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGDOF.begin());

        iPtr = finleyFile->globalNodesIndex;
        nodeGNI.insert(nodeGNI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGNI.begin());

        iPtr = finleyFile->globalReducedDOFIndex;
        nodeGRDFI.insert(nodeGRDFI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGRDFI.begin());

        iPtr = finleyFile->globalReducedNodesIndex;
        nodeGRNI.insert(nodeGRNI.end(), numNodes, 0);
        copy(iPtr, iPtr+numNodes, nodeGRNI.begin());

    }
    return true;
#else // VISIT_PLUGIN
    return false;
#endif
}

//
//
//
bool FinleyNodes::readFromNc(NcFile* ncFile)
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

    nodeDist.clear();
    nodeDist.insert(nodeDist.end(), mpisize+1, 0);
    var = ncFile->get_var("Nodes_NodeDistribution");
    var->get(&nodeDist[0], mpisize+1);

    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();
    nodeID.clear();
    nodeTag.clear();
    nodeGDOF.clear();
    nodeGNI.clear();
    nodeGRDFI.clear();
    nodeGRNI.clear();

    // Only attempt to read further if there are any nodes.
    // Having no nodes is not an error.
    if (numNodes > 0) {
        var = ncFile->get_var("Nodes_Coordinates");
        for (int i=0; i<numDims; i++) {
            float* c = new float[numNodes];
            var->set_cur(0, i);
            var->get(c, numNodes, 1);
            coords.push_back(c);
        }

        nodeID.insert(nodeID.end(), numNodes, 0);
        var = ncFile->get_var("Nodes_Id");
        var->get(&nodeID[0], numNodes);

        nodeTag.insert(nodeTag.end(), numNodes, 0);
        var = ncFile->get_var("Nodes_Tag");
        var->get(&nodeTag[0], numNodes);

        nodeGDOF.insert(nodeGDOF.end(), numNodes, 0);
        var = ncFile->get_var("Nodes_gDOF");
        var->get(&nodeGDOF[0], numNodes);

        nodeGNI.insert(nodeGNI.end(), numNodes, 0);
        var = ncFile->get_var("Nodes_gNI");
        var->get(&nodeGNI[0], numNodes);

        nodeGRDFI.insert(nodeGRDFI.end(), numNodes, 0);
        var = ncFile->get_var("Nodes_grDfI");
        var->get(&nodeGRDFI[0], numNodes);

        nodeGRNI.insert(nodeGRNI.end(), numNodes, 0);
        var = ncFile->get_var("Nodes_grNI");
        var->get(&nodeGRNI[0], numNodes);
    }

    return true;
#else // !USE_NETCDF
    return false;
#endif
}

//
//
//
const IntVec& FinleyNodes::getVarDataByName(const string& name) const
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
StringVec FinleyNodes::getVarNames() const
{
    StringVec res;
    res.push_back("Nodes_Id");
    res.push_back("Nodes_Tag");
    res.push_back("Nodes_gDOF");
    res.push_back("Nodes_gNI");
    res.push_back("Nodes_grDfI");
    res.push_back("Nodes_grNI");
    return res;
}

//
//
//
int FinleyNodes::getGlobalNumNodes() const
{
    int ret=0;
    if (!nodeDist.empty())
        ret = nodeDist[nodeDist.size()-1];
    return ret;
}

//
//
//
void FinleyNodes::writeCoordinatesVTK(ostream& os, int ownIndex)
{
    if (numNodes > 0) {
        int firstId = nodeDist[ownIndex];
        int lastId = nodeDist[ownIndex+1];
        for (size_t i=0; i<numNodes; i++) {
            if (firstId <= nodeGNI[i] && nodeGNI[i] < lastId) {
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
bool FinleyNodes::writeToSilo(DBfile* dbfile)
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

} // namespace weipa

