
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


/****************************************************************************

  Finley: Domain

  searches for faces in the mesh which are matching.

*****************************************************************************/

#include "FinleyDomain.h"
#include "Util.h"

#include <escript/index.h>

//#define Finley_TRACE

namespace finley {

static double lockingGridSize = 0.;

// this structure is used for matching surface elements
struct FaceCenter
{
   int refId;
   std::vector<double> x;
};

/// comparison function for findMatchingFaces
bool FaceCenterCompare(const FaceCenter& e1, const FaceCenter& e2)
{
    for (int i = 0; i < e1.x.size(); i++) {
        bool l = (e1.x[i] < e2.x[i]+lockingGridSize);
        bool g = (e2.x[i] < e1.x[i]+lockingGridSize);
        if (! (l && g)) {
            if (l) return true;
            if (g) return false;
        }
    }
    return e1.refId < e2.refId; // strict order to e1==e2 is false
}

inline double getDist(int e0, int i0, int e1, int i1, int numDim, int NN,
                      const double* X)
{
    double dist = 0.;
    for (int i = 0; i < numDim; i++) {
        dist = std::max(dist, std::abs(X[INDEX3(i, i0, e0, numDim, NN)]
                    - X[INDEX3(i, i1, e1, numDim, NN)]));
    }
    return dist;
}

void FinleyDomain::findMatchingFaces(double safety_factor, double tolerance,
                                     int* numPairs, int* elem0, int* elem1,
                                     int* matching_nodes_in_elem1) const
{
    const_ReferenceElement_ptr refElement(m_faceElements->referenceElementSet->
                                            borrowReferenceElement(false));
    const int numDim = m_nodes->numDim;
    const int NN = m_faceElements->numNodes;
    const int numNodesOnFace = refElement->Type->numNodesOnFace;
    const int* faceNodes = refElement->Type->faceNodes;
    const int* shiftNodes = refElement->Type->shiftNodes;
    const int* reverseNodes = refElement->Type->reverseNodes;

    if (numNodesOnFace <= 0) {
        std::stringstream ss;
        ss << "Mesh::findMatchingFaces: matching faces cannot be applied to "
            "face elements of type " << refElement->Type->Name;
        throw escript::ValueError(ss.str());
    }
    double* X = new double[NN * numDim * m_faceElements->numElements];
    std::vector<FaceCenter> center(m_faceElements->numElements);
    int* a1 = new int[NN];
    int* a2 = new int[NN];
    double h = std::numeric_limits<double>::max();

    // TODO: OMP
    for (index_t e = 0; e < m_faceElements->numElements; e++) {
        // get the coordinates of the nodes
        util::gather(NN, &(m_faceElements->Nodes[INDEX2(0,e,NN)]), numDim,
                     m_nodes->Coordinates, &X[INDEX3(0,0,e,numDim,NN)]);
        // get the element center
        center[e].refId = e;
        center[e].x.assign(numDim, 0);
        for (int i0 = 0; i0 < numNodesOnFace; i0++) {
            for (int i = 0; i < numDim; i++)
                center[e].x[i] += X[INDEX3(i,faceNodes[i0],e,numDim,NN)];
        }
        for (int i = 0; i < numDim; i++)
            center[e].x[i] /= numNodesOnFace;
        // get the minimum distance between nodes in the element
        for (int i0 = 0; i0 < numNodesOnFace; i0++) {
            for (int i1 = i0+1; i1 < numNodesOnFace; i1++) {
                double h_local = getDist(e, faceNodes[i0], e, faceNodes[i1], numDim, NN, X);
                h = std::min(h, h_local);
            }
        }
    }
    lockingGridSize = h*std::max(safety_factor, 0.);
#ifdef Finley_TRACE
    std::cout << "locking grid size is " << lockingGridSize << std::endl;
    std::cout << "absolute tolerance is " << h * tolerance << std::endl;
    std::cout << "number of face elements is " << m_faceElements->numElements << std::endl;
#endif
    // sort the elements by center coordinates (lexicographical)
    std::sort(center.begin(), center.end(), FaceCenterCompare);
    // find elements with matching center
    *numPairs = 0;

    // TODO: OMP
    for (index_t e = 0; e < m_faceElements->numElements-1; e++) {
        double dist = 0.;
        for (int i = 0; i < numDim; i++)
            dist = std::max(dist, std::abs(center[e].x[i]-center[e+1].x[i]));
        if (dist < h * tolerance) {
            const int e_0 = center[e].refId;
            const int e_1 = center[e+1].refId;
            elem0[*numPairs] = e_0;
            elem1[*numPairs] = e_1;
            // now the element e_1 is rotated such that the first node in
            // element e_0 and e_1 have the same coordinates
            int* perm = a1;
            int* perm_tmp = a2;
            for (int i = 0; i < NN; i++)
                perm[i] = i;
            while (1) {
                // if node 0 and perm[0] are the same we are ready
                dist = getDist(e_0, 0, e_1, perm[0], numDim, NN, X);
                if (dist <= h*tolerance)
                    break;
                if (shiftNodes[0] >= 0) {
                    // rotate the nodes
                    int* itmp_ptr = perm;
                    perm = perm_tmp;
                    perm_tmp = itmp_ptr;
                    #pragma ivdep
                    for (int i = 0; i < NN; i++)
                        perm[i] = perm_tmp[shiftNodes[i]];
                }
                // if the permutation is back at the identity, i.e. perm[0]=0,
                // the faces don't match:
                if (perm[0] == 0) {
                    std::stringstream ss;
                    ss << "Mesh::findMatchingFaces: couldn't match first node "
                        "of element " << e_0 << " to touching element " << e_1;
                    throw escript::ValueError(ss.str());
                }
            }
            // now we check if the second nodes match
            if (numNodesOnFace > 1) {
                dist = getDist(e_0, 1, e_1, perm[faceNodes[1]], numDim, NN, X);
                // if the second node does not match we reverse the
                // direction of the nodes
                if (dist > h*tolerance) {
                    // rotate the nodes
                    if (reverseNodes[0] < 0) {
                        std::stringstream ss;
                        ss << "Mesh::findMatchingFaces: couldn't match the"
                            " second node of element " << e_0
                            << " to touching element " << e_1;
                        throw escript::ValueError(ss.str());
                    } else {
                        int* itmp_ptr = perm;
                        perm = perm_tmp;
                        perm_tmp = itmp_ptr;
                        #pragma ivdep
                        for (int i = 0; i < NN; i++)
                            perm[i] = perm_tmp[reverseNodes[i]];
                        dist = getDist(e_0, 1, e_1, perm[faceNodes[1]], numDim, NN, X);
                        if (dist > h*tolerance) {
                            std::stringstream ss;
                            ss << "Mesh::findMatchingFaces: couldn't match the"
                                " second node of element " << e_0
                                << " to touching element " << e_1;
                            throw escript::ValueError(ss.str());
                        }
                    }
                }
            }
            // we check if the rest of the face nodes match
            for (int i = 2; i < numNodesOnFace; i++) {
                const int n = faceNodes[i];
                dist = getDist(e_0, n, e_1, perm[n], numDim, NN, X);
                if (dist > h*tolerance) {
                    std::stringstream ss;
                    ss << "Mesh::findMatchingFaces: couldn't match the "
                        << i << "-th node of element " << e_0
                        << " to touching element " << e_1;
                    throw escript::ValueError(ss.str());
                }
            }
            // copy over the permuted nodes of e_1 into matching_nodes_in_elem1
            for (int i = 0; i < NN; i++)
                matching_nodes_in_elem1[INDEX2(i,*numPairs,NN)] =
                    m_faceElements->Nodes[INDEX2(perm[i],e_1,NN)];
            (*numPairs)++;
        }
    }
#ifdef Finley_TRACE
    std::cout << "number of pairs of matching faces " << *numPairs << std::endl;
#endif

    delete[] X;
    delete[] a1;
    delete[] a2;
}

} // namespace finley

