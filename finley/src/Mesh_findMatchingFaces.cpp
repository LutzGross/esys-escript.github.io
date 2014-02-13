
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************

  Finley: Mesh

  searches for faces in the mesh which are matching.

*****************************************************************************/

#include "Util.h"
#include "Mesh.h"

namespace finley {

static double lockingGridSize=0.;

/// comparison function for findMatchingFaces
bool FaceCenterCompare(const FaceCenter& e1, const FaceCenter& e2)
{
    for (int i=0; i<e1.x.size(); i++) {
        bool l=(e1.x[i] < e2.x[i]+lockingGridSize);
        bool g=(e2.x[i] < e1.x[i]+lockingGridSize);
        if (! (l && g)) {
            if (l) return true;
            if (g) return false;
        }
    }
    if (e1.refId < e2.refId) {
        return true;
    } else if (e1.refId > e2.refId) {
        return false;
    }
    return false;	// strict order to e1==e2 is false
}

inline double getDist(int e0, int i0, int e1, int i1, int numDim, int NN,
                      const double* X)
{
    double dist=0.;
    for (int i=0; i<numDim; i++) {
        dist=std::max(dist, std::abs(X[INDEX3(i, i0, e0, numDim, NN)]
                    - X[INDEX3(i, i1, e1, numDim, NN)]));
    }
    return dist;
}

void Mesh::findMatchingFaces(double safety_factor, double tolerance,
                             int* numPairs, int* elem0, int* elem1,
                             int* matching_nodes_in_elem1)
{
    const_ReferenceElement_ptr refElement(FaceElements->referenceElementSet->
                                            borrowReferenceElement(false));
    const int numDim=Nodes->numDim;
    const int NN=FaceElements->numNodes;
    const int numNodesOnFace=refElement->Type->numNodesOnFace;
    const int* faceNodes=refElement->Type->faceNodes;
    const int* shiftNodes=refElement->Type->shiftNodes;
    const int* reverseNodes=refElement->Type->reverseNodes;

    if (numNodesOnFace <= 0) {
        char error_msg[LenErrorMsg_MAX];
        sprintf(error_msg, "Mesh::findMatchingFaces: matching faces cannot be applied to face elements of type %s",refElement->Type->Name);
        setError(TYPE_ERROR, error_msg);
        return;
    }
    double* X = new double[NN*numDim*FaceElements->numElements];
    std::vector<FaceCenter> center(FaceElements->numElements);
    int* a1=new int[NN];
    int* a2=new int[NN];
    double h=std::numeric_limits<double>::max();

    // TODO: OMP
    for (int e=0; e<FaceElements->numElements; e++) {
        // get the coordinates of the nodes
        util::gather(NN, &(FaceElements->Nodes[INDEX2(0,e,NN)]), numDim,
                     Nodes->Coordinates, &X[INDEX3(0,0,e,numDim,NN)]);
        // get the element center
        center[e].refId=e;
        center[e].x.assign(numDim, 0);
        for (int i0=0; i0<numNodesOnFace; i0++) {
            for (int i=0; i<numDim; i++)
                center[e].x[i] += X[INDEX3(i,faceNodes[i0],e,numDim,NN)];
        }
        for (int i=0; i<numDim; i++)
            center[e].x[i]/=numNodesOnFace;
        // get the minimum distance between nodes in the element
        for (int i0=0; i0<numNodesOnFace; i0++) {
            for (int i1=i0+1; i1<numNodesOnFace; i1++) {
                double h_local=getDist(e, faceNodes[i0], e, faceNodes[i1], numDim, NN, X);
                h=std::min(h, h_local);
            }
        }
    }
    lockingGridSize=h*std::max(safety_factor, 0.);
#ifdef Finley_TRACE
    printf("locking grid size is %e\n", lockingGridSize);
    printf("absolute tolerance is %e.\n", h * tolerance);
#endif
    // sort the elements by center coordinates (lexicographical)
    std::sort(center.begin(), center.end(), FaceCenterCompare);
    // find elements with matching center
    *numPairs=0;

    // TODO: OMP
    for (int e=0; e<FaceElements->numElements-1 && noError(); e++) {
        double dist=0.;
        for (int i=0; i<numDim; i++)
            dist=std::max(dist, std::abs(center[e].x[i]-center[e+1].x[i]));
        if (dist < h * tolerance) {
            const int e_0=center[e].refId;
            const int e_1=center[e+1].refId;
            elem0[*numPairs]=e_0;
            elem1[*numPairs]=e_1;
            // now the element e_1 is rotated such that the first node in
            // element e_0 and e_1 have the same coordinates
            int* perm=a1;
            int* perm_tmp=a2;
            for (int i=0; i<NN; i++)
                perm[i]=i;
            while (noError()) {
                // if node 0 and perm[0] are the same we are ready
                dist=getDist(e_0, 0, e_1, perm[0], numDim, NN, X);
                if (dist <= h*tolerance)
                    break;
                if (shiftNodes[0]>=0) {
                    // rotate the nodes
                    int* itmp_ptr=perm;
                    perm=perm_tmp;
                    perm_tmp=itmp_ptr;
                    #pragma ivdep
                    for (int i=0; i<NN; i++)
                        perm[i]=perm_tmp[shiftNodes[i]];
                }
                // if the permutation is back at the identity, i.e. perm[0]=0,
                // the faces don't match:
                if (perm[0]==0) {
                    char error_msg[LenErrorMsg_MAX];
                    sprintf(error_msg, "Mesh_findMatchingFaces: couldn't match first node of element %d to touching element %d", e_0, e_1);
                    setError(VALUE_ERROR, error_msg);
                }
            }
            // now we check if the second nodes match
            if (noError()) {
                if (numNodesOnFace > 1) {
                    dist=getDist(e_0, 1, e_1, perm[faceNodes[1]], numDim, NN, X);
                    // if the second node does not match we reverse the
                    // direction of the nodes
                    if (dist > h*tolerance) {
                        // rotate the nodes
                        if (reverseNodes[0] < 0) {
                            char error_msg[LenErrorMsg_MAX];
                            sprintf(error_msg, "Mesh_findMatchingFaces: couldn't match the second node of element %d to touching element %d", e_0, e_1);
                            setError(VALUE_ERROR, error_msg);
                        } else {
                            int* itmp_ptr=perm;
                            perm=perm_tmp;
                            perm_tmp=itmp_ptr;
                            #pragma ivdep
                            for (int i=0; i<NN; i++)
                                perm[i]=perm_tmp[reverseNodes[i]];
                            dist=getDist(e_0, 1, e_1, perm[faceNodes[1]], numDim, NN, X);
                            if (dist > h*tolerance) {
                                char error_msg[LenErrorMsg_MAX];
                                sprintf(error_msg, "Mesh_findMatchingFaces: couldn't match the second node of element %d to touching element %d", e_0, e_1);
                                setError(VALUE_ERROR, error_msg);
                            }
                        }
                    }
                }
            }
            // we check if the rest of the face nodes match
            if (noError()) {
                for (int i=2; i<numNodesOnFace; i++) {
                    const int n=faceNodes[i];
                    dist=getDist(e_0, n, e_1, perm[n], numDim, NN, X);
                    if (dist > h*tolerance) {
                        char error_msg[LenErrorMsg_MAX];
                        sprintf(error_msg, "Mesh_findMatchingFaces: couldn't match the %d-th node of element %d to touching element %d", i, e_0, e_1);
                        setError(VALUE_ERROR, error_msg);
                        break;
                    }
                }
            }
            // copy over the permuted nodes of e_1 into matching_nodes_in_elem1
            if (noError()) {
                for (int i=0; i<NN; i++)
                    matching_nodes_in_elem1[INDEX2(i,*numPairs,NN)]=FaceElements->Nodes[INDEX2(perm[i],e_1,NN)];
            }
            (*numPairs)++;
        }
    }
#ifdef Finley_TRACE
    printf("number of pairs of matching faces %d\n",*numPairs);
#endif

    /* clean up */
    delete[] X;
    delete[] a1;
    delete[] a2;
}

} // namespace finley

