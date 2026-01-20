
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
OxleyNodes::OxleyNodes(const OxleyNodes& m)
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
    nodeTag.clear();

    numDims = dom->getDim();
    globalNumNodes = dom->getNumDataPointsGlobal();
    pair<int,dim_t> shape = dom->getDataShape(oxley::Nodes);
    numNodes = dom->getNumNodes();
    // oxley::IndexVector dist = dom->getNodeDistribution();
    // nodeDist.assign(dist.begin(), dist.end());

    if (numNodes > 0) {
        for (int d=0; d<numDims; d++) {
            float* c = new float[numNodes];
            coords.push_back(c);
        }

        if (numDims==2) {
            const oxley::Rectangle * rect = static_cast<const oxley::Rectangle *>(dom);
        #pragma omp parallel for
            for(p4est_topidx_t treeid = rect->p4est->first_local_tree; treeid <= rect->p4est->last_local_tree; ++treeid) {
                p4est_tree_t * tree = p4est_tree_array_index(rect->p4est->trees, treeid);
                sc_array_t * tquadrants = &tree->quadrants;
                p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
                for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
                    p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                    int l = (int) P4EST_QUADRANT_LEN(quad->level);
                    double xy[3];
                    int lxy[4][2]={{0,0},{0,l},{l,0},{l,l}};
                    for(int n = 0; n < 4; n++)
                    {
                        p4est_qcoord_to_vertex(rect->p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy);
                        long nodeid=rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                        coords[0][nodeid]=xy[0];
                        coords[1][nodeid]=xy[1];
                    #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                        std::cout << "coords (" << coords[0][nodeid] << ", " << coords[1][nodeid] << ") " << std::endl;
                    #endif
                    }
                }
            }
        } else {
            const oxley::Brick * brick = static_cast<const oxley::Brick *>(dom);
        #pragma omp parallel for
            for(p4est_topidx_t treeid = brick->p8est->first_local_tree; treeid <= brick->p8est->last_local_tree; ++treeid) {
                p8est_tree_t * tree = p8est_tree_array_index(brick->p8est->trees, treeid);
                sc_array_t * tquadrants = &tree->quadrants;
                oxley::p8est_locidx_t Q = (oxley::p8est_locidx_t) tquadrants->elem_count;
                for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
                    p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
                    // p8est_qcoord_t length = P8EST_QUADRANT_LEN(quad->level);
                    double xy[3];
                    p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, quad->x, quad->y, quad->z, xy);
                    long nodeid=brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
                    coords[0][nodeid]=xy[0];
                    coords[1][nodeid]=xy[1];
                    coords[2][nodeid]=xy[2];
                }
            }
        }
        const dim_t* iPtr = dom->borrowSampleReferenceIDs(oxley::Nodes);
        nodeID.assign(iPtr, iPtr+numNodes);
        // iPtr = dom->borrowListOfTags(oxley::Nodes); //todo
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

