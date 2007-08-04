/*
 ************************************************************
 *          Copyright 2007 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: Mesh: optimizes the distribution of DOFs across processors */
/*   using ParMETIS. On return a new distribution is given and the globalDOF are relabled */
/*   accordingly but the mesh has not been redesitributed yet                             */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "IndexList.h"

/**************************************************************/

void Finley_Mesh_optimizeDOFDistribution(Finley_Mesh* in,dim_t *distribution) {

     dim_t dim, i,j,k, myNumVertices,p, mpiSize, len, globalNumVertices,*partition_count=NULL, *new_distribution=NULL, *loc_partition_count=NULL;
     bool_t *setNewDOFId=NULL;
     index_t myFirstVertex, myLastVertex, firstVertex, lastVertex, *newGlobalDOFID=NULL;
     size_t mpiSize_size;
     index_t* partition=NULL;
     Paso_Pattern *pattern=NULL;
     Paso_MPI_rank myRank,dest,source,current_rank, rank;
     Finley_IndexList* index_list=NULL;
     float *xyz=NULL;
     
     #ifdef PASO_MPI
     MPI_Status status;
     #endif

     if (in==NULL) return;
     if (in->Nodes == NULL) return;

     myRank=in->MPIInfo->rank;
     mpiSize=in->MPIInfo->size;
     mpiSize_size=mpiSize*sizeof(dim_t);
     dim=in->Nodes->numDim;
     /* first step is to distribute the elements according to a global X of DOF */

     myFirstVertex=distribution[myRank];
     myLastVertex=distribution[myRank+1];
     myNumVertices=myLastVertex-myFirstVertex;
     globalNumVertices=distribution[mpiSize];
     len=0;
     for (p=0;p<mpiSize;++p) len=MAX(len,distribution[p+1]-distribution[p]);
     partition=TMPMEMALLOC(len,index_t); /* len is used for the sending around of partition later on */
     xyz=TMPMEMALLOC(myNumVertices*dim,float);
     partition_count=TMPMEMALLOC(mpiSize+1,dim_t);
     new_distribution=TMPMEMALLOC(mpiSize+1,dim_t);
     newGlobalDOFID=TMPMEMALLOC(len,index_t);
     setNewDOFId=TMPMEMALLOC(in->Nodes->numNodes,bool_t);
     if (!(Finley_checkPtr(partition) || Finley_checkPtr(xyz) || Finley_checkPtr(partition_count) || Finley_checkPtr(partition_count) || Finley_checkPtr(newGlobalDOFID) || Finley_checkPtr(setNewDOFId))) {

         /* set the coordinates: *?
         /* it is assumed that at least one node on this processor provides a coordinate */
         #pragma omp parallel for private(i,j,k);
         for (i=0;i<in->Nodes->numNodes;++i) {
             k=in->Nodes->globalDegreesOfFreedom[i]-myFirstVertex;
             if ((k>=0) && (k<myNumVertices)) {
                for (j=0;j<dim;++j) xyz[k*dim+j]=(float)(in->Nodes->Coordinates[INDEX2(j,i,dim)]); 
             }
         }

         index_list=TMPMEMALLOC(myNumVertices,Finley_IndexList);
         /* create the adjacency structure xadj and adjncy */
         if (! Finley_checkPtr(index_list)) {
            #pragma omp parallel private(i)
            {
              #pragma omp for schedule(static)
              for(i=0;i<myNumVertices;++i) {
                   index_list[i].extension=NULL;
                   index_list[i].n=0;
              }
              /*  insert contributions from element matrices into colums index index_list: */
              Finley_IndexList_insertElementsWithRowRange(index_list, myFirstVertex, myLastVertex,
                                                          in->Elements,in->Nodes->globalDegreesOfFreedom,
                                                          in->Nodes->globalDegreesOfFreedom);
              Finley_IndexList_insertElementsWithRowRange(index_list, myFirstVertex, myLastVertex,
                                                          in->FaceElements,in->Nodes->globalDegreesOfFreedom,
                                                          in->Nodes->globalDegreesOfFreedom);
              Finley_IndexList_insertElementsWithRowRange(index_list, myFirstVertex, myLastVertex,
                                                          in->ContactElements,in->Nodes->globalDegreesOfFreedom,
                                                          in->Nodes->globalDegreesOfFreedom);
              Finley_IndexList_insertElementsWithRowRange(index_list, myFirstVertex, myLastVertex,
                                                          in->Points,in->Nodes->globalDegreesOfFreedom,
                                                          in->Nodes->globalDegreesOfFreedom);
           }
           
           /* create the matrix pattern */
           pattern=Finley_IndexList_createPattern(myNumVertices,index_list,0,globalNumVertices);

           /* clean up index list */
           if (index_list!=NULL) {
              #pragma omp parallel for private(i) 
              for(i=0;i<myNumVertices;++i) Finley_IndexList_free(index_list[i].extension);
           }

           if (Finley_noError()) {

/*

        ParMETIS_V3_PartGeomKway(distribution,
                                 pattern->ptr,
                                 pattern->index,
                                 idxtype *vwgt, +
                                 idxtype *adjwgt, +
                                 int *wgtﬂag, +
                                 int *numﬂag, +
                                 dim,
                                 xyz,
                                 int *ncon, +
                                 mpiSize, 
                                 ﬂoat *tpwgts, +
                                 ﬂoat *ubvec, +
                                 int *options, +
                                 int *edgecut, +
                                 partition,
                                 in->MPIInfo->comm);
*/
               for (i=0;i<myNumVertices;++i) partition[i]=myRank; /* remove */
           }

           Paso_Pattern_free(pattern);
           
           /* create a new distributioin and labeling of the DOF */
           memset(new_distribution,0,mpiSize_size);
           #pragma omp parallel private(loc_partition_count)
           {
               loc_partition_count=THREAD_MEMALLOC(mpiSize,dim_t);
               memset(loc_partition_count,0,mpiSize_size);
               #pragma omp for private(i)
               for (i=0;i<myNumVertices;++i) loc_partition_count[partition[i]]++ ;
               #pragma omp critical
               {
                  for (i=0;i<mpiSize;++i) new_distribution[i]+=loc_partition_count[i];
               }
               THREAD_MEMFREE(loc_partition_count);
           }
           #ifdef PASO_MPI
              MPI_Allreduce( new_distribution, partition_count, mpiSize, MPI_INT, MPI_SUM, in->MPIInfo->comm );
           #else
               for (i=0;i<mpiSize;++i) partition_count[i]=new_distribution[i];
           #endif
           new_distribution[0]=0;
           for (i=0;i<mpiSize;++i) {
               new_distribution[i+1]=new_distribution[i]+partition_count[i];
               partition_count[i]=0;
           }
           for (i=0;i<myNumVertices;++i) {
              rank=partition[i];
              newGlobalDOFID[i]=new_distribution[rank]+partition_count[rank];
              partition_count[rank]++;
           }
           /* now the overlap needs to be created by sending the partition around*/

           dest=Paso_MPIInfo_mod(mpiSize, myRank + 1);
           source=Paso_MPIInfo_mod(mpiSize, myRank - 1);
           current_rank=myRank;
           #pragma omp parallel for private(i);
           for (i=0;i<in->Nodes->numNodes;++i) setNewDOFId[i]=TRUE;

           for (p=0; p< mpiSize; ++p) {

               firstVertex=distribution[current_rank];
               lastVertex=distribution[current_rank+1];
               #pragma omp parallel for private(i,j,k);
               for (i=0;i<in->Nodes->numNodes;++i) {
                   k=in->Nodes->globalDegreesOfFreedom[i];
                   if (setNewDOFId[i] && (firstVertex<=k) && (k<lastVertex)) {
                        in->Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                        setNewDOFId[i]=FALSE;
                   }
               }

               if (p<mpiSize-1) {  /* the final send can be skipped */
                  #ifdef PASO_MPI
                  MPI_Sendrecv_replace(newGlobalDOFID,len, MPI_INT,
                                       dest, in->MPIInfo->msg_tag_counter,
                                       source, in->MPIInfo->msg_tag_counter,
                                       in->MPIInfo->comm,&status);
                  #endif
                  in->MPIInfo->msg_tag_counter++;
                  current_rank=Paso_MPIInfo_mod(mpiSize, current_rank-1);
              }
           }
           for (i=0;i<mpiSize+1;++i) distribution[i]=new_distribution[i];

           
         }
         TMPMEMFREE(index_list);
     }
     TMPMEMFREE(newGlobalDOFID);
     TMPMEMFREE(setNewDOFId);
     TMPMEMFREE(new_distribution);
     TMPMEMFREE(partition_count);
     TMPMEMFREE(partition);
     TMPMEMFREE(xyz);
     return;
}
