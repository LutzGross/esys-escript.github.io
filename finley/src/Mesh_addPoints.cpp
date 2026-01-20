
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

  Adds points to the Points table of the mesh.

*****************************************************************************/

#include "FinleyDomain.h"

#include <escript/index.h>

using escript::ValueError;

namespace finley {

#ifdef ESYS_MPI
void MPI_minimizeDistance(void* invec_p, void* inoutvec_p, int *len,
                          MPI_Datatype* dtype)
{
    const int numPoints = (*len) / 2;
    double* invec = reinterpret_cast<double*>(invec_p);
    double* inoutvec = reinterpret_cast<double*>(inoutvec_p);
    for (int i = 0; i < numPoints; i++) {
        if (invec[2*i] < inoutvec[2*i]) {
            inoutvec[2*i] = invec[2*i];
            inoutvec[2*i+1] = invec[2*i+1];
        }
    }
}
#endif

void FinleyDomain::addDiracPoints(const std::vector<double>& points,
                                  const std::vector<int>& tags)
{
    const int numDim = getDim();
    if (points.size() % numDim != 0) {
        std::stringstream ss;
        unsigned long size = points.size();
        ss << "Number of coordinates in dirac tags is " << size
           << " - this should be a multiple of the domain's dimensionality "
           << numDim;
        throw ValueError(ss.str());
    }

    // points will be flattened
    int numPoints = points.size() / numDim;
    int numTags = tags.size();

    if (numPoints != numTags)
	   throw ValueError("Number of dirac tags must match number of dirac points.");

    if (numPoints == 0)
        return;

    const escript::DataTypes::real_t LARGE_POSITIVE_FLOAT =
                                         escript::DataTypes::real_t_max();
    ElementFile* oldPoints = m_points;
    const_ReferenceElementSet_ptr refPoints;
    int numOldPoints;
    if (oldPoints == NULL) {
        refPoints.reset(new ReferenceElementSet(Point1, integrationOrder,
                    reducedIntegrationOrder));
        numOldPoints = 0;
    } else {
        refPoints = oldPoints->referenceElementSet;
        numOldPoints = oldPoints->numElements;
    }
    ElementFile* newPoints = new ElementFile(refPoints, m_mpiInfo);

    // first we find the node which is the closest on this processor:
    double *dist_p = new double[numPoints];
    int *node_id_p = new int[numPoints];
    int *point_index_p = new int[numPoints];
    // the code below does actually initialise this before using it

    for (int i = 0; i < numPoints; ++i) {
        dist_p[i] = LARGE_POSITIVE_FLOAT;
        node_id_p[i] = -1;
    }

    const double* coords = m_nodes->Coordinates;
    if (numDim == 3) {
#pragma omp parallel
        {
            for (int i = 0; i < numPoints; ++i) {
                const double X0 = points[INDEX2(0,i,numDim)];
                const double X1 = points[INDEX2(1,i,numDim)];
                const double X2 = points[INDEX2(2,i,numDim)];
                double dist_local = LARGE_POSITIVE_FLOAT;
                int node_id_local = -1;
#pragma omp for
                for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
                    const double D0 = coords[INDEX2(0,n,numDim)] - X0;
                    const double D1 = coords[INDEX2(1,n,numDim)] - X1;
                    const double D2 = coords[INDEX2(2,n,numDim)] - X2;
                    const double d = D0*D0 + D1*D1 + D2*D2;
                    if (d < dist_local) {
                        dist_local = d;
                        node_id_local = n;
                    }
                }
#pragma omp critical
                {
                    if (dist_local < dist_p[i] || (dist_local == dist_p[i] && node_id_p[i]>node_id_local)) {
                        dist_p[i] = dist_local;
                        node_id_p[i] = node_id_local;
                    }
                }
            }
        } // end parallel section
    } else if (numDim == 2) {
#pragma omp parallel
        {
            for (int i = 0; i < numPoints; ++i) {
                const double X0 = points[INDEX2(0,i,numDim)];
                const double X1 = points[INDEX2(1,i,numDim)];
                double dist_local = LARGE_POSITIVE_FLOAT;
                int node_id_local = -1;
#pragma omp for
                for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
                    const double D0 = coords[INDEX2(0,n,numDim)] - X0;
                    const double D1 = coords[INDEX2(1,n,numDim)] - X1;
                    const double d = D0*D0 + D1*D1;
                    if (d < dist_local) {
                        dist_local = d;
                        node_id_local = n;
                    }
                }
#pragma omp critical
                {
                  if (dist_local < dist_p[i] || (dist_local == dist_p[i] && node_id_p[i]>node_id_local)) {
                      dist_p[i] = dist_local;
                      node_id_p[i] = node_id_local;
                  }
                }
            }
        } // end parallel section
    } else { // numDim==1
#pragma omp parallel
        {
            for (int i = 0; i < numPoints; ++i) {
                const double X0 = points[INDEX2(0,i,numDim)];
                double dist_local = LARGE_POSITIVE_FLOAT;
                int node_id_local = -1;
#pragma omp for
                for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
                    const double D0 = coords[INDEX2(0,n,numDim)] - X0;
                    const double d = D0*D0;
                    if (d < dist_local) {
                        dist_local = d;
                        node_id_local = n;
                    }
                }
#pragma omp critical
                {
                    if (dist_local < dist_p[i] || (dist_local == dist_p[i] && node_id_p[i]>node_id_local)) {
                        dist_p[i] = dist_local;
                        node_id_p[i] = node_id_local;
                    }
                }
            }
        } // end parallel section
    } //numDim

#ifdef ESYS_MPI
    // now we need to reduce this across all processors
    const int count = 2*numPoints;
    double *sendbuf = new double[count];
    double *recvbuf = new double[count];

    for (int i = 0; i < numPoints; ++i) {
        sendbuf[2*i  ] = dist_p[i];
        sendbuf[2*i+1] = static_cast<double>(m_nodes->Id[node_id_p[i]]);
    }
    MPI_Op op;
    MPI_Op_create(MPI_minimizeDistance, true, &op);
    MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, m_mpiInfo->comm);
    MPI_Op_free(&op);
    // if the node id has changed we found another node which is closer
    // elsewhere
    for (int i = 0; i < numPoints; ++i) {
        const int best_fit_Id = static_cast<int>(recvbuf[2*i+1]+0.5);
        if (best_fit_Id != m_nodes->Id[node_id_p[i]]) {
            node_id_p[i] = -1;
        }
    }
    delete[] sendbuf;
    delete[] recvbuf;
#endif
    delete[] dist_p;

    // we pick the points to be used on this processor
    int numNewPoints = 0;
    const index_t firstDOF = m_nodes->degreesOfFreedomDistribution->getFirstComponent();
    const index_t lastDOF = m_nodes->degreesOfFreedomDistribution->getLastComponent();

    for (int i = 0; i < numPoints; ++i) {
        if (node_id_p[i] > -1) {
            // this processor uses a node which is identical to point i
            if (m_nodes->globalReducedDOFIndex[node_id_p[i]] > -1) {
                // the point is also used in the reduced mesh
                const index_t global_id = m_nodes->globalDegreesOfFreedom[node_id_p[i]];
                if (firstDOF <= global_id && global_id < lastDOF) {
                    // is this point actually relevant
                    bool notAnOldPoint = true;
                    if (numOldPoints > 0) {
                        // is this point already in the Point table?
                        for (int k = 0; k < numOldPoints; ++k) {
                            if (global_id == oldPoints->Nodes[k]) {
                                notAnOldPoint=false;
                                break;
                            }
                        }
                    }
                    if (notAnOldPoint) {
                        // is this point unique in the new list of points?
                        bool notANewPoint = true;
                        for (int k = 0; k < numNewPoints; ++k) {
                            if (global_id == m_nodes->globalDegreesOfFreedom[node_id_p[point_index_p[k]]]) {
                                notANewPoint = false;
                                break;
                            }
                        }
                        if (notANewPoint) {
                            point_index_p[numNewPoints] = i;
                            numNewPoints++;
                        }
                    }
                }
            }
        }
    }

    // now we are ready to create the new Point table
    newPoints->allocTable(numOldPoints+numNewPoints);
    if (numOldPoints > 0) {
#pragma omp parallel for schedule(static)
        for (int n = 0; n < numOldPoints; n++) {
            newPoints->Owner[n] = oldPoints->Owner[n];
            newPoints->Id[n]    = oldPoints->Id[n];
            newPoints->Tag[n]   = oldPoints->Tag[n];
            newPoints->Nodes[n] = oldPoints->Nodes[n];
            newPoints->Color[n] = 0;
        }
    }
#pragma omp parallel for schedule(static)
    for (int n = 0; n < numNewPoints; n++) {
        const int idx = point_index_p[n];
        newPoints->Owner[numOldPoints+n] = m_mpiInfo->rank;
        newPoints->Id[numOldPoints+n]    = 0;
        newPoints->Tag[numOldPoints+n]   = tags[idx];
        newPoints->Nodes[numOldPoints+n] = node_id_p[idx];
        newPoints->Color[numOldPoints+n] = 0;
    }
    newPoints->minColor = 0;
    newPoints->maxColor = 0;

    // all done, clean up
    delete[] node_id_p;
    delete[] point_index_p;
    delete oldPoints;
    m_points = newPoints;
}

} // namespace finley

