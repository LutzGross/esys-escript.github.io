
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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
            nodeGNI.push_back(fullNodes->nodeGNI[*it]);
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
    nodeGNI.assign(numNodes, 0);

    if (numNodes > 0) {
        for (int d=0; d<numDims; d++) {
            float* c = new float[numNodes];
            coords.push_back(c);
        }
        const IntVec NN = dom->getNumNodesPerDim();
        const IntVec faces = dom->getNumFacesPerBoundary();
        const IntVec NS = dom->getNumSubdivisionsPerDim();

        if (numDims==2) {
            pair<double,double> xx=dom->getFirstCoordAndSpacing(0);
            pair<double,double> yy=dom->getFirstCoordAndSpacing(1);
            const int left=(faces[0]==0 ? 1 : 0);
            const int right=(faces[1]==0 ? NN[0]-2 : NN[0]-1);
            const int bottom=(faces[2]==0 ? 1 : 0);
            const int top=(faces[3]==0 ? NN[1]-2 : NN[1]-1);
            const int ownN0=right-left+1;
            const int ownN1=top-bottom+1;
            const int rx=dom->getMPIRank()%NS[0];
            const int ry=dom->getMPIRank()/NS[0];
            for (int i1=0; i1<NN[1]; i1++) {
                for (int i0=0; i0<NN[0]; i0++) {
                    coords[0][i0+NN[0]*i1] = (float)(xx.first+i0*xx.second);
                    coords[1][i0+NN[0]*i1] = (float)(yy.first+i1*yy.second);
                }
            }
            for (int i1=-1; i1<2; i1++) {
                for (int i0=-1; i0<2; i0++) {
                    // location of neighbour rank
                    const int nrx=rx+i0;
                    const int nry=ry+i1;
                    if (nrx>=0 && nry>=0 && nrx<NS[0] && nry<NS[1]) {
                        // index of first node on neighbour rank
                        const int first=nodeDist[nry*NS[0]+nrx];
                        if (i0==0 && i1==0) {
                            // own nodes
                            for (int y=bottom; y<=top; y++)
                                for (int x=left; x<=right; x++)
                                    nodeGNI[x+y*NN[0]]=first+x-left+(y-bottom)*ownN0;
                        } else if (i0==0) {
                            // top or bottom
                            const int myFirst=(i1==-1 ? 0 : NN[0]*(NN[1]-1));
                            const int nFirst=(i1==-1 ? first+ownN0*(ownN1-1) : first);
                            for (int i=left; i<=right; i++)
                                nodeGNI[myFirst+i]=nFirst+i-left;
                        } else if (i1==0) {
                            // left or right
                            const int myFirst=(i0==-1 ? 0 : NN[0]-1);
                            const int nFirst=(i0==-1 ? first+ownN0-1 : first);
                            for (int i=bottom; i<=top; i++)
                                nodeGNI[myFirst+i*NN[0]]= nFirst+(i-bottom)*ownN0;
                        } else {
                            // corner
                            const int mIdx=(i0+1)/2*(NN[0]-1)+(i1+1)/2*(NN[0]*(NN[1]-1));
                            const int nIdx=(1-i0)/2*(ownN0-1)+(1-i1)/2*(ownN0*(ownN1-1));
                            nodeGNI[mIdx]=first+nIdx;
                        }
                    }
                }
            }

        } else {
            pair<double,double> xx=dom->getFirstCoordAndSpacing(0);
            pair<double,double> yy=dom->getFirstCoordAndSpacing(1);
            pair<double,double> zz=dom->getFirstCoordAndSpacing(2);
            const int left=(faces[0]==0 ? 1 : 0);
            const int right=(faces[1]==0 ? NN[0]-2 : NN[0]-1);
            const int bottom=(faces[2]==0 ? 1 : 0);
            const int top=(faces[3]==0 ? NN[1]-2 : NN[1]-1);
            const int front=(faces[4]==0 ? 1 : 0);
            const int back=(faces[5]==0 ? NN[2]-2 : NN[2]-1);
            const int ownN0=right-left+1;
            const int ownN1=top-bottom+1;
            const int ownN2=back-front+1;
            const int rx=dom->getMPIRank()%NS[0];
            const int ry=dom->getMPIRank()%(NS[0]*NS[1])/NS[0];
            const int rz=dom->getMPIRank()/(NS[0]*NS[1]);
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
            for (int i2=-1; i2<2; i2++) {
                for (int i1=-1; i1<2; i1++) {
                    for (int i0=-1; i0<2; i0++) {
                        // location of neighbour rank
                        const int nrx=rx+i0;
                        const int nry=ry+i1;
                        const int nrz=rz+i2;
                        if (nrx>=0 && nry>=0 && nrz>=0
                                && nrx<NS[0] && nry<NS[1] && nrz<NS[2]) {
                            // index of first node on neighbour rank
                            const int first=nodeDist[nrz*NS[0]*NS[1]+nry*NS[0]+nrx];
                            if (i0==0 && i1==0 && i2==0) {
                                // own nodes
                                for (int z=front; z<=back; z++)
                                    for (int y=bottom; y<=top; y++)
                                        for (int x=left; x<=right; x++)
                                            nodeGNI[x+y*NN[0]+z*NN[0]*NN[1]]=
                                                first+x-left+(y-bottom)*ownN0
                                                +(z-front)*ownN0*ownN1;
                            } else if (i0==0 && i1==0) {
                                // front or back plane
                                for (int y=bottom; y<=top; y++) {
                                    const int myFirst=(i2==-1 ? y*NN[0] : NN[0]*NN[1]*(NN[2]-1)+y*NN[0]);
                                    const int nFirst=(i2==-1 ?
                                            first+ownN0*ownN1*(ownN2-1)+(y-bottom)*ownN0 : first+(y-bottom)*ownN0);
                                    for (int x=left; x<=right; x++)
                                        nodeGNI[myFirst+x]=nFirst+x-left;
                                }
                            } else if (i0==0 && i2==0) {
                                // top or bottom plane
                                for (int z=front; z<=back; z++) {
                                    const int myFirst=(i1==-1 ? z*NN[0]*NN[1] : NN[0]*((z+1)*NN[1]-1));
                                    const int nFirst=(i1==-1 ?
                                            first+ownN0*((z-front+1)*ownN1-1) : first+(z-front)*ownN0*ownN1);
                                    for (int x=left; x<=right; x++)
                                        nodeGNI[myFirst+x]=nFirst+x-left;
                                }
                            } else if (i1==0 && i2==0) {
                                // left or right plane
                                for (int z=front; z<=back; z++) {
                                    const int myFirst=(i0==-1 ? z*NN[0]*NN[1] : NN[0]*(1+z*NN[1])-1);
                                    const int nFirst=(i0==-1 ?
                                            first+ownN0*(1+(z-front)*ownN1)-1 : first+(z-front)*ownN0*ownN1);
                                    for (int y=bottom; y<=top; y++)
                                        nodeGNI[myFirst+y*NN[0]]=nFirst+(y-bottom)*ownN0;
                                }
                            } else if (i0==0) {
                                // edge in x direction
                                const int myFirst=(i1+1)/2*NN[0]*(NN[1]-1)
                                                +(i2+1)/2*NN[0]*NN[1]*(NN[2]-1);
                                const int nFirst=first+(1-i1)/2*ownN0*(ownN1-1)
                                                +(1-i2)/2*ownN0*ownN1*(ownN2-1);
                                for (int i=left; i<=right; i++)
                                    nodeGNI[myFirst+i]=nFirst+i-left;
                            } else if (i1==0) {
                                // edge in y direction
                                const int myFirst=(i0+1)/2*(NN[0]-1)
                                                +(i2+1)/2*NN[0]*NN[1]*(NN[2]-1);
                                const int nFirst=first+(1-i0)/2*(ownN0-1)
                                                +(1-i2)/2*ownN0*ownN1*(ownN2-1);
                                for (int i=bottom; i<=top; i++)
                                    nodeGNI[myFirst+i*NN[0]]=nFirst+(i-bottom)*ownN0;
                            } else if (i2==0) {
                                // edge in z direction
                                const int myFirst=(i0+1)/2*(NN[0]-1)
                                                  +(i1+1)/2*NN[0]*(NN[1]-1);
                                const int nFirst=first+(1-i0)/2*(ownN0-1)
                                                 +(1-i1)/2*ownN0*(ownN1-1);
                                for (int i=front; i<=back; i++)
                                    nodeGNI[myFirst+i*NN[0]*NN[1]]=nFirst+(i-front)*ownN0*ownN1;
                            } else {
                                // corner
                                const int mIdx=(i0+1)/2*(NN[0]-1)
                                               +(i1+1)/2*NN[0]*(NN[1]-1)
                                               +(i2+1)/2*NN[0]*NN[1]*(NN[2]-1);
                                const int nIdx=(1-i0)/2*(ownN0-1)
                                               +(1-i1)/2*ownN0*(ownN1-1)
                                               +(1-i2)/2*ownN0*ownN1*(ownN2-1);
                                nodeGNI[mIdx]=first+nIdx;
                            }
                        }
                    }
                }
            }
        }

        const int* iPtr = dom->borrowSampleReferenceIDs(ripley::Nodes);
        nodeID.assign(iPtr, iPtr+numNodes);

        //iPtr = dom->borrowListOfTags(ripley::Nodes);
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

