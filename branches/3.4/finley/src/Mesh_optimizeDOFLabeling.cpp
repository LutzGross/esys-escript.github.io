
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


/**********************************************************************************************/

/*   Finley: Mesh: optimizes the labeling of the DOFs on each processor */

/**********************************************************************************************/

#include "Mesh.h"
#include "IndexList.h"

/************************************************************************************/

void Finley_Mesh_optimizeDOFLabeling(Finley_Mesh* in,dim_t *distribution)
{
    int myFirstVertex,myLastVertex, *newGlobalDOFID=NULL, firstVertex, lastVertex;
    int k;
    int mpiSize, myNumVertices,len, p, i;
    Paso_Pattern *pattern=NULL;
    Esys_MPI_rank myRank,current_rank;
    IndexList* index_list=NULL;
#ifdef ESYS_MPI
    Esys_MPI_rank dest,source;
    MPI_Status status;
#endif

    if (in==NULL) return;
    if (in->Nodes == NULL) return;

    myRank=in->MPIInfo->rank;
    mpiSize=in->MPIInfo->size;
    myFirstVertex=distribution[myRank];
    myLastVertex=distribution[myRank+1];
    myNumVertices=myLastVertex-myFirstVertex;
    len=0;
    for (p=0;p<mpiSize;++p) len=MAX(len,distribution[p+1]-distribution[p]);

    index_list=new IndexList[myNumVertices];
    newGlobalDOFID=new int[len];
    /* create the adjacency structure xadj and adjncy */
    if (! ( Finley_checkPtr(index_list) || Finley_checkPtr(newGlobalDOFID) ) ) {
#pragma omp parallel
        {
            /*  insert contributions from element matrices into columns index index_list: */
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                    myFirstVertex, myLastVertex, in->Elements,
                    in->Nodes->globalDegreesOfFreedom,
                    in->Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                    myFirstVertex, myLastVertex, in->FaceElements,
                    in->Nodes->globalDegreesOfFreedom,
                    in->Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                    myFirstVertex, myLastVertex, in->ContactElements,
                    in->Nodes->globalDegreesOfFreedom,
                    in->Nodes->globalDegreesOfFreedom);
            IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list,
                    myFirstVertex, myLastVertex, in->Points,
                    in->Nodes->globalDegreesOfFreedom,
                    in->Nodes->globalDegreesOfFreedom);
        }
        /* create the local matrix pattern */
        pattern=IndexList_createPattern(0,myNumVertices,index_list,myFirstVertex, myLastVertex,-myFirstVertex);

        if (Finley_noError()) Paso_Pattern_reduceBandwidth(pattern,newGlobalDOFID); 

        Paso_Pattern_free(pattern);
    }
    Esys_MPIInfo_noError(in->MPIInfo);
    if (Finley_noError()) {
        /* shift new labeling to create a global id */
#pragma omp parallel for private(i)
        for (i=0;i<myNumVertices;++i) newGlobalDOFID[i]+=myFirstVertex;


        /* distribute new labeling to other processors */
#ifdef ESYS_MPI
        dest=Esys_MPIInfo_mod(mpiSize, myRank + 1);
        source=Esys_MPIInfo_mod(mpiSize, myRank - 1);
#endif
        current_rank=myRank;
        for (p=0; p< mpiSize; ++p) {
            firstVertex=distribution[current_rank];
            lastVertex=distribution[current_rank+1];
#pragma omp parallel for private(i,k)
            for (i=0;i<in->Nodes->numNodes;++i) {
                k=in->Nodes->globalDegreesOfFreedom[i];
                if ( (firstVertex<=k) && (k<lastVertex)) {
                    in->Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                }
            }
   
            if (p<mpiSize-1) {  /* the final send can be skipped */
#ifdef ESYS_MPI
                MPI_Sendrecv_replace(newGlobalDOFID,len, MPI_INT,
                                          dest, in->MPIInfo->msg_tag_counter,
                                          source, in->MPIInfo->msg_tag_counter,
                                          in->MPIInfo->comm,&status);
#endif
                in->MPIInfo->msg_tag_counter++;
                current_rank=Esys_MPIInfo_mod(mpiSize, current_rank-1);
            }
        }
    }
    delete[] index_list;
    delete[] newGlobalDOFID;
}

