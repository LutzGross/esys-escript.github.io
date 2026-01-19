
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
                throw WeipaException("Unknown shape type");
                break;
        }
        owner = dom->getOwnerVector(fsType);

        const dim_t* iPtr = dom->borrowSampleReferenceIDs(fsType);
        ID.assign(iPtr, iPtr+numElements);

        nodes.clear();
        if (dom->getDim() == 2) {
            const oxley::Rectangle * rect = static_cast<const oxley::Rectangle *>(dom);

            if (fsType==oxley::Elements) {
                for(p4est_topidx_t treeid = rect->p4est->first_local_tree; treeid <= rect->p4est->last_local_tree; ++treeid) {
                    p4est_tree_t * tree = p4est_tree_array_index(rect->p4est->trees, treeid);
                    sc_array_t * tquadrants = &tree->quadrants;
                    p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            // #pragma omp parallel for
                    for(int q = 0; q < Q; ++q) {
                        p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                        long ids[4]={0};
                        rect->getNeighouringNodeIDs(quad->level, quad->x, quad->y, treeid, ids);
                        // Convert the node ordering
                        nodes.push_back(ids[0]);
                        nodes.push_back(ids[2]);
                        nodes.push_back(ids[3]);
                        nodes.push_back(ids[1]);

                        #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                            std::cout << "Element " << 
                                        ids[0] << "-" <<
                                        ids[2] << "-" <<
                                        ids[3] << "-" <<
                                        ids[1] <<  std::endl;
                        #endif
                    }
                }
            } else if (fsType==oxley::FaceElements) {
                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                    long counter = 0;
                #endif
                std::vector<int> lnodes,rnodes,bnodes,tnodes;

                for(p4est_topidx_t treeid = rect->p4est->first_local_tree; treeid <= rect->p4est->last_local_tree; ++treeid) 
                {
                    p4est_tree_t * tree = p4est_tree_array_index(rect->p4est->trees, treeid);
                    sc_array_t * tquadrants = &tree->quadrants;
                    p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
                    for(int q = 0; q < Q; ++q) 
                    {
                        p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                        p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
                        // int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
                        p4est_qcoord_t lxy[4][2] = {{0,0},{l,0},{0,l},{l,l}};
                        double xy[3] = {0};
                        int nodeids[4]={-1};
                        bool do_check_yes_no[4]={false};
                        for(int n = 0; n < 4; n++)
                        {
                            if(rect->isLeftBoundaryNode(quad, n, treeid, l))
                            {
                                p4est_qcoord_to_vertex(rect->p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy);
                                lnodes.push_back(rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes L " << counter++ << ": " << rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second << std::endl;
                                #endif
                            }

                            if(rect->isRightBoundaryNode(quad, n, treeid, l))
                            {
                                p4est_qcoord_to_vertex(rect->p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy);
                                rnodes.push_back(rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes R " << counter++ << ": " << rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second << std::endl;
                                #endif
                            }
                                
                            if(rect->isBottomBoundaryNode(quad, n, treeid, l))
                            {
                                p4est_qcoord_to_vertex(rect->p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy);
                                bnodes.push_back(rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes B " << counter++ << ": " << rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second << std::endl;
                                #endif
                            }
                                
                            if(rect->isTopBoundaryNode(quad, n, treeid, l))
                            {
                                p4est_qcoord_to_vertex(rect->p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy);
                                tnodes.push_back(rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes T " << counter++ << ": " << rect->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second << std::endl;
                                #endif
                            }
                        }
                    }
                }

                for(int i = 0; i < lnodes.size(); i++)
                    nodes.push_back(lnodes[i]);
                for(int i = 0; i < rnodes.size(); i++)
                    nodes.push_back(rnodes[i]);
                for(int i = 0; i < bnodes.size(); i++)
                    nodes.push_back(bnodes[i]);
                for(int i = 0; i < tnodes.size(); i++)
                    nodes.push_back(tnodes[i]);

            }
        } else { //3d
            const oxley::Brick * brick = static_cast<const oxley::Brick *>(dom);
            if (fsType==oxley::Elements) {
                for(oxley::p8est_topidx_t treeid = brick->p8est->first_local_tree; treeid <= brick->p8est->last_local_tree; ++treeid) {
                    p8est_tree_t * tree = p8est_tree_array_index(brick->p8est->trees, treeid);
                    sc_array_t * tquadrants = &tree->quadrants;
                    oxley::p8est_locidx_t Q = (oxley::p8est_locidx_t) tquadrants->elem_count;
            #pragma omp parallel for
                    for(int q = 0; q < Q; ++q) {
                        p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
                        long ids[8]={-1};
                        brick->getNeighouringNodeIDs(quad->level, quad->x, quad->y, quad->z, treeid, ids);

                        nodes.push_back(ids[2]);
                        nodes.push_back(ids[0]);
                        nodes.push_back(ids[1]);
                        nodes.push_back(ids[3]);
                        nodes.push_back(ids[6]);
                        nodes.push_back(ids[4]);
                        nodes.push_back(ids[5]);
                        nodes.push_back(ids[7]);
                    }
                }
            } else if (fsType==oxley::FaceElements) {
                // const dim_t* NE = dom->getNumElements();
                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                    long counter = 0;
                #endif
                std::vector<int> anodes,bnodes,cnodes,dnodes,enodes,fnodes;

                for(oxley::p8est_topidx_t treeid = brick->p8est->first_local_tree; treeid <= brick->p8est->last_local_tree; ++treeid) 
                {
                    p8est_tree_t * tree = p8est_tree_array_index(brick->p8est->trees, treeid);
                    sc_array_t * tquadrants = &tree->quadrants;
                    oxley::p8est_locidx_t Q = (oxley::p8est_locidx_t) tquadrants->elem_count;
                    for(int q = 0; q < Q; ++q) 
                    {
                        p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
                        oxley::p8est_qcoord_t l = P8EST_QUADRANT_LEN(quad->level);
                        oxley::p8est_qcoord_t lxy[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},
                                                    {0,0,l},{l,0,l},{0,l,l},{l,l,l}};
                        double xy[3] = {0};
                        int nodeids[8]={-1};
                        bool do_check_yes_no[8]={false};
                        for(int n = 0; n < 8; n++)
                        {
                            if(brick->isLeftBoundaryNode(quad, n, treeid, l))
                            {
                                p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, 
                                            quad->x+lxy[n][0], quad->y+lxy[n][1], quad->z+lxy[n][2], xy);
                                anodes.push_back(brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes L " << counter++ << ": " 
                                    << brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second << std::endl;
                                #endif
                            }

                            if(brick->isRightBoundaryNode(quad, n, treeid, l))
                            {
                                p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, 
                                            quad->x+lxy[n][0], quad->y+lxy[n][1], quad->z+lxy[n][2], xy);
                                bnodes.push_back(brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes R " << counter++ << ": " 
                                    << brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second << std::endl;
                                #endif
                            }

                            if(brick->isBottomBoundaryNode(quad, n, treeid, l))
                            {
                                p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, 
                                            quad->x+lxy[n][0], quad->y+lxy[n][1], quad->z+lxy[n][2], xy);
                                cnodes.push_back(brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes Bo " << counter++ << ": " 
                                    << brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second << std::endl;
                                #endif
                            }

                            if(brick->isTopBoundaryNode(quad, n, treeid, l))
                            {
                                p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, 
                                            quad->x+lxy[n][0], quad->y+lxy[n][1], quad->z+lxy[n][2], xy);
                                dnodes.push_back(brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes T " << counter++ << ": " 
                                    << brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second << std::endl;
                                #endif
                            }

                            if(brick->isAboveBoundaryNode(quad, n, treeid, l))
                            {
                                p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, 
                                            quad->x+lxy[n][0], quad->y+lxy[n][1], quad->z+lxy[n][2], xy);
                                enodes.push_back(brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes A " << counter++ << ": " 
                                    << brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second << std::endl;
                                #endif
                            }

                            if(brick->isBelowBoundaryNode(quad, n, treeid, l))
                            {
                                p8est_qcoord_to_vertex(brick->p8est->connectivity, treeid, 
                                            quad->x+lxy[n][0], quad->y+lxy[n][1], quad->z+lxy[n][2], xy);
                                fnodes.push_back(brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second);
                                #ifdef OXLEY_ENABLE_DEBUG_WEIPA
                                std::cout << "nodes Be " << counter++ << ": " 
                                    << brick->NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second << std::endl;
                                #endif
                            }                            
                        }
                    }
                }

                for(int i = 0; i < anodes.size(); i++)
                    nodes.push_back(anodes[i]);
                for(int i = 0; i < bnodes.size(); i++)
                    nodes.push_back(bnodes[i]);
                for(int i = 0; i < cnodes.size(); i++)
                    nodes.push_back(cnodes[i]);
                for(int i = 0; i < dnodes.size(); i++)
                    nodes.push_back(dnodes[i]);
                for(int i = 0; i < enodes.size(); i++)
                    nodes.push_back(enodes[i]);
                for(int i = 0; i < fnodes.size(); i++)
                    nodes.push_back(fnodes[i]);
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

