
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


/****************************************************************************

  Adds points to the Points table of the mesh.

*****************************************************************************/

#include "Mesh.h"

#ifdef ESYS_MPI
void Finley_Mesh_MPI_minimizeDistance(void *invec_p, void *inoutvec_p,
                                      int *len, MPI_Datatype *dtype)
{
    const int numPoints = (*len)/2;
    double *invec = reinterpret_cast<double*>(invec_p);
    double *inoutvec = reinterpret_cast<double*>(inoutvec_p);
    for (int i=0; i<numPoints; i++) {
        if (invec[2*i] < inoutvec[2*i]) {
            inoutvec[2*i] = invec[2*i];
            inoutvec[2*i+1] = invec[2*i+1];
        }
    }
}
#endif

void Finley_Mesh_addPoints(Finley_Mesh* mesh, const dim_t numPoints,
                           const double* points_ptr, const index_t* tags_ptr)
{
    Esys_MPIInfo *mpi_info = Esys_MPIInfo_getReference(mesh->MPIInfo);
    ElementFile *oldPoints=mesh->Points;
    Finley_ReferenceElementSet *refPoints=NULL;
    dim_t numOldPoints;
    if (oldPoints == NULL) {
        refPoints=Finley_ReferenceElementSet_alloc(Finley_Point1,
                        mesh->integrationOrder, mesh->reducedIntegrationOrder);
        numOldPoints=0;
    } else {
        refPoints=Finley_ReferenceElementSet_reference(oldPoints->referenceElementSet);
        numOldPoints=mesh->Points->numElements;
    }
    ElementFile *newPoints=new ElementFile(refPoints, mpi_info);

    // first we find the node which is the closest on this processor:
    double *dist_p = new double[numPoints];
    index_t *node_id_p = new index_t[numPoints];
    index_t *point_index_p = new index_t[numPoints];

    for (dim_t i=0; i<numPoints; ++i) {
        dist_p[i]=LARGE_POSITIVE_FLOAT;
        node_id_p[i]=-1;
        node_id_p[i]=-1;
    }

    const double *coords = mesh->Nodes->Coordinates;
    const dim_t numDim=mesh->Nodes->numDim;
    if (numDim == 3) {
#pragma omp parallel
        {
            for (dim_t i=0; i<numPoints; ++i) {
                const double X0=points_ptr[INDEX2(0,i,numDim)];
                const double X1=points_ptr[INDEX2(1,i,numDim)];
                const double X2=points_ptr[INDEX2(2,i,numDim)];
                double dist_local=LARGE_POSITIVE_FLOAT;
                index_t node_id_local=-1;
#pragma omp for
                for (dim_t n=0; n<mesh->Nodes->numNodes; n++) {
                    const double D0=coords[INDEX2(0,n,numDim)] - X0;
                    const double D1=coords[INDEX2(1,n,numDim)] - X1;
                    const double D2=coords[INDEX2(2,n,numDim)] - X2;
                    const double d = D0*D0 + D1*D1 + D2*D2;
                    if (d < dist_local) {
                        dist_local = d;
                        node_id_local = n;
                    }
                }
#pragma omp critical
                {
                    if (dist_local < dist_p[i]) {
                        dist_p[i] = dist_local;
                        node_id_p[i] = node_id_local;
                    }
                }
            }
        } // end parallel section
    } else if (numDim == 2) {
#pragma omp parallel
        {
            for (dim_t i=0; i<numPoints; ++i) {
                const double X0=points_ptr[INDEX2(0,i,numDim)];
                const double X1=points_ptr[INDEX2(1,i,numDim)];
                double dist_local=LARGE_POSITIVE_FLOAT;
                index_t node_id_local=-1;
#pragma omp for
                for (dim_t n=0; n<mesh->Nodes->numNodes; n++) {
                    const double D0=coords[INDEX2(0,n,numDim)] - X0;
                    const double D1=coords[INDEX2(1,n,numDim)] - X1;
                    const double d = D0*D0 + D1*D1;
                    if (d < dist_local) {
                        dist_local = d;
                        node_id_local = n;
                    }
                }
#pragma omp critical
                {
                    if (dist_local < dist_p[i]) {
                        dist_p[i] = dist_local;
                        node_id_p[i] = node_id_local;
                    }
                }
            }
        } // end parallel section
    } else { // numDim==1
#pragma omp parallel
        {
            for (dim_t i=0; i<numPoints; ++i) {
                const double X0=points_ptr[INDEX2(0,i,numDim)];
                double dist_local=LARGE_POSITIVE_FLOAT;
                index_t node_id_local=-1;
#pragma omp for
                for (dim_t n=0; n< mesh-> Nodes->numNodes; n++) {
                    const double D0=coords[INDEX2(0,n,numDim)] - X0;
                    const double d = D0*D0;
                    if (d < dist_local) {
                        dist_local = d;
                        node_id_local = n;
                    }
                }
#pragma omp critical
                {
                    if (dist_local < dist_p[i]) {
                        dist_p[i] = dist_local;
                        node_id_p[i] = node_id_local;
                    }
                }
            }
        } // end parallel section
    } //numDim

#ifdef ESYS_MPI
    // now we need to reduce this across all processors
    {
        MPI_Op op;
        int count = 2*numPoints;
        double *sendbuf=new double[count];
        double *recvbuf=new double[count];

        for (dim_t i=0; i<numPoints; ++i) {
            sendbuf[2*i  ]=dist_p[i];
            sendbuf[2*i+1]=static_cast<double>(mesh->Nodes->Id[node_id_p[i]]);
        }
        MPI_Op_create(Finley_Mesh_MPI_minimizeDistance, TRUE, &op);
        MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, op, mpi_info->comm);
        MPI_Op_free(&op);
        // if the node id has changed we found another node which is closer
        // elsewhere
        for (dim_t i=0; i<numPoints; ++i) {
            const int best_fit_Id = static_cast<const int>(recvbuf[2*i+1]+0.5);
            if (best_fit_Id != mesh->Nodes->Id[node_id_p[i]]) {
                node_id_p[i] = -1;
            }
        }
        delete[] sendbuf;
        delete[] recvbuf;
    }
#endif
    delete[] dist_p;

    // we pick the points to be used on this processor
    dim_t numNewPoints=0;
    const index_t firstDOF=Paso_Distribution_getFirstComponent(mesh->Nodes->degreesOfFreedomDistribution);
    const index_t lastDOF=Paso_Distribution_getLastComponent(mesh->Nodes->degreesOfFreedomDistribution);

    for (dim_t i=0; i<numPoints; ++i) {
        if (node_id_p[i]>-1) {
            // this processor uses a node which is identical to point i
            if (mesh->Nodes->globalReducedDOFIndex[node_id_p[i]] > -1) {
                // the point is also used in the reduced mesh
                const index_t global_id=mesh->Nodes->globalDegreesOfFreedom[node_id_p[i]];
                if ( (firstDOF<=global_id) && (global_id<lastDOF) ) {
                    // is this point actually relevant
                    bool notAnOldPoint=true;
                    if (numOldPoints > 0) {
                        // is this point already in the Point table?
                        for (dim_t k=0; k<numOldPoints; ++k) {
                            if (global_id == oldPoints->Nodes[k]) {
                                notAnOldPoint=false;
                                break;
                            }
                        }
                    }
                    if (notAnOldPoint) {
                        // is this point unique in the new list of points?
                        bool notANewPoint=true;
                        for (dim_t k=0; k<numNewPoints; ++k) {
                            if (global_id == mesh->Nodes->globalDegreesOfFreedom[node_id_p[point_index_p[k]]]) {
                                notANewPoint=false;
                                break;
                            }
                        }
                        if (notANewPoint) {
                            point_index_p[numNewPoints]=i;
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
        for(dim_t n=0; n<numOldPoints; n++) {
            newPoints->Owner[n]=oldPoints->Owner[n];
            newPoints->Id[n]   =oldPoints->Id[n];
            newPoints->Tag[n]  =oldPoints->Tag[n];
            newPoints->Nodes[n]=oldPoints->Nodes[n];
            newPoints->Color[n]=0;
        }
    }
#pragma omp parallel for schedule(static)
    for(dim_t n=0; n<numNewPoints; n++) {
        const index_t idx = point_index_p[n];
        newPoints->Owner[numOldPoints+n]=mpi_info->rank;
        newPoints->Id[numOldPoints+n]   =0;
        newPoints->Tag[numOldPoints+n]  =tags_ptr[idx];
        newPoints->Nodes[numOldPoints+n]=node_id_p[idx];
        newPoints->Color[numOldPoints+n]=0;
    }
    newPoints->minColor=0;
    newPoints->maxColor=0;

    // all done, clean up
    delete[] node_id_p;
    delete[] point_index_p;
    Finley_ReferenceElementSet_dealloc(refPoints);
    Esys_MPIInfo_free(mpi_info);
    if (Finley_noError()) {
        delete oldPoints;
        mesh->Points=newPoints;
    } else {
        delete newPoints;
    }
}

