
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

#include <weipa/OxleyNodes.h>

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#endif

#ifndef VISIT_PLUGIN
#include <oxley/OxleyDomain.h>
#include <oxley/Brick.h>
#include <oxley/Rectangle.h>
#endif

#ifndef VISIT_PLUGIN
using escript::DataTypes::dim_t;
#endif

#include "p4est/p4est.h"
#include "p4est/p8est.h"

using namespace std;

namespace weipa {

//
// Constructor with name
//
OxleyNodes::OxleyNodes(const string& meshName) :
    numDims(0), numNodes(0), globalNumNodes(0), name(meshName)
{
}

//
//
//
OxleyNodes::OxleyNodes(OxleyNodes_ptr fullNodes, IntVec& requiredNodes,
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
            nodeGNI.push_back(fullNodes->nodeGNI[*it]);
            nodeTag.push_back(fullNodes->nodeTag[*it]);
            indexMap[*it] = newIndex;
            *it = newIndex++;
        } else {
            *it = res->second;
        }
    }

    // carry over the constraints of the nodes that survived, remapped. A master
    // may not be required by any element of this mesh - it is a corner of the
    // coarse neighbour, whose element can be on another rank - and without it the
    // value at the hanging node cannot be computed, so pull it in. It ends up an
    // unreferenced point, which costs a point and no cell.
    const NodeConstraints& fullConstraints = fullNodes->nodeConstraints;
    for (size_t i = 0; i < fullConstraints.size(); i++) {
        const NodeConstraint& src = fullConstraints[i];
        IndexMap::const_iterator res = indexMap.find(src.node);
        if (res == indexMap.end())
            continue;                   // this hanging node is not in this mesh
        NodeConstraint nc;
        nc.node = (int) res->second;
        nc.numMasters = src.numMasters;
        for (int k = 0; k < src.numMasters; k++) {
            IndexMap::const_iterator m = indexMap.find(src.master[k]);
            if (m == indexMap.end()) {
                nodeID.push_back(fullNodes->nodeID[src.master[k]]);
                nodeGNI.push_back(fullNodes->nodeGNI[src.master[k]]);
                nodeTag.push_back(fullNodes->nodeTag[src.master[k]]);
                indexMap[src.master[k]] = newIndex;
                nc.master[k] = (int) newIndex++;
            } else {
                nc.master[k] = (int) m->second;
            }
            nc.weight[k] = src.weight[k];
        }
        nodeConstraints.push_back(nc);
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
OxleyNodes::OxleyNodes(const OxleyNodes& m)
{
    numDims = m.numDims;
    numNodes = m.numNodes;
    globalNumNodes = m.globalNumNodes;
    nodeID = m.nodeID;
    nodeGNI = m.nodeGNI;
    nodeTag = m.nodeTag;
    nodeDist = m.nodeDist;
    nodeConstraints = m.nodeConstraints;
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
OxleyNodes::~OxleyNodes()
{
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
}

//
//
//
bool OxleyNodes::initFromOxley(const oxley::OxleyDomain* dom)
{
#ifndef VISIT_PLUGIN
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
    coords.clear();
    nodeID.clear();
    nodeGNI.clear();
    nodeTag.clear();
    nodeConstraints.clear();

    // Consume the domain's public lnodes-based mesh view; the node numbering
    // scheme stays inside the domain (no p4est / coordinate-hash access here).
    // Everything comes from the clean lnodes mesh-access interface; no p4est
    // or coordinate-hash access here. Node ids are the mesh-access global ids
    // (0..N-1 in serial). weipa maps data to nodes by matching these ids
    // against the data's sample reference ids, so this is consistent whenever
    // the domain's node numbering is the lnodes numbering (true in 2D today;
    // the 3D data path still uses the old numbering until A4/A5).
    // materializeHanging: a hanging position has no lnodes node, and the slot
    // holding a master instead would draw the cell with a corner in the wrong
    // place. The materialised nodes carry no sample, hence nodeConstraints below.
    const oxley::MeshAccess m = dom->getMeshAccess(true);
    numDims = m.numDim;
    numNodes = (int) m.numNodes;

    // The output numbering the domain built: rank r owns exactly
    // [nodeDist[r], nodeDist[r+1]), so a writer emitting one shared point list
    // (VTK) has every node written by exactly one rank, and every rank's
    // connectivity indexes the same list.
    nodeDist.assign(m.outputDistribution.begin(), m.outputDistribution.end());
    globalNumNodes = nodeDist.empty() ? numNodes : (int) nodeDist.back();

    if (numNodes > 0) {
        for (int d = 0; d < numDims; d++) {
            float* c = new float[numNodes];
            for (int i = 0; i < numNodes; i++)
                c[i] = (float) m.nodeCoords[(size_t) i * numDims + d];
            coords.push_back(c);
        }
        // node id labels must match the id space of the escript Data (which
        // uses the domain's node sample ids, a permutation of 0..N-1), so
        // weipa maps data to the correct nodes. Consistent because the data
        // sample order equals the mesh-access lnodes order.
        // ... except for the materialised nodes, which are not samples at all:
        // they take the ids the domain gave them, a block above every lnodes id.
        const dim_t* iPtr = dom->borrowSampleReferenceIDs(oxley::Nodes);
        nodeID.assign(iPtr, iPtr + m.numRealNodes);
        for (long i = m.numRealNodes; i < m.numNodes; i++)
            nodeID.push_back((int) m.nodeGlobalId[i]);
        // node tags are not part of the mesh-access interface yet
        nodeTag.assign(numNodes, 0);
        nodeGNI.assign(m.nodeOutputIndex.begin(), m.nodeOutputIndex.end());

        const int mpc = m.mastersPerConstrainedNode;
        for (size_t i = 0; i < m.constrainedNodes.size(); i++) {
            NodeConstraint nc;
            nc.node = (int) m.constrainedNodes[i];
            nc.numMasters = 0;
            for (int k = 0; k < mpc && k < 4; k++) {
                const long master = m.constraintMasters[i * mpc + k];
                if (master < 0)
                    continue;
                nc.master[nc.numMasters] = (int) master;
                nc.weight[nc.numMasters] = (float) m.constraintWeights[i * mpc + k];
                nc.numMasters++;
            }
            nodeConstraints.push_back(nc);
        }
    }
    return true;
#else // VISIT_PLUGIN
    return false;
#endif
}

//
//
//
const IntVec& OxleyNodes::getVarDataByName(const string& name) const
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
StringVec OxleyNodes::getVarNames() const
{
    StringVec res;
    res.push_back("Nodes_Id");
    res.push_back("Nodes_Tag");
    return res;
}

//
//
//
void OxleyNodes::writeCoordinatesVTK(ostream& os, int ownIndex)
{
    if (numNodes > 0) {
        // by output index, not nodeID: this decides which rank writes which
        // point, so it must use the numbering nodeDist describes
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
bool OxleyNodes::writeToSilo(DBfile* dbfile)
{
#ifdef ESYS_HAVE_SILO
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

#else // !ESYS_HAVE_SILO
    return false;
#endif
}

} // namespace weipa

