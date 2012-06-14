
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: Mesh: optimizes the distribution of DOFs across processors */
/*   using ParMETIS. On return a new distribution is given and the globalDOF */
/*   are relabeled accordingly but the mesh is not redistributed yet. */

/**************************************************************/

#include "Mesh.h"
#include "IndexList.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_PARMETIS
#include "parmetis.h"
#endif

/**************************************************************
   Check whether there is any node which has no vertex. In case 
   such a node exists, we don't use parmetis since parmetis requires
   that every node has at least 1 vertex (at line 129 of file
   "xyzpart.c" in parmetis 3.1.1, variable "nvtxs" would be 0 if 
   any node has no vertex).
 **************************************************************/
#ifdef USE_PARMETIS
int Check_Inputs_For_Parmetis(dim_t mpiSize, dim_t rank, dim_t *distribution, MPI_Comm *comm)
{
  dim_t i, len;
  int ret_val = 1;

  if (rank == 0){
    for (i=0; i<mpiSize; i++){
      len = distribution[i+1] - distribution[i];
      if (len == 0){
        ret_val = 0;
        break;
      }
    }
  } 
  MPI_Bcast(&ret_val, 1, MPI_INTEGER, 0, *comm);
  if (ret_val == 0) 
    printf("INFO: Parmetis is not used since some nodes have no vertex!\n");
  return ret_val;
}
#endif



/**************************************************************/

void Finley_Mesh_optimizeDOFDistribution(Finley_Mesh* in,dim_t *distribution) {

     dim_t dim, i,j,k, myNumVertices,p, mpiSize, len, globalNumVertices,*partition_count=NULL, *new_distribution=NULL, *loc_partition_count=NULL;
     bool_t *setNewDOFId=NULL;
     index_t myFirstVertex, myLastVertex, firstVertex, lastVertex, *newGlobalDOFID=NULL;
     size_t mpiSize_size;
     index_t* partition=NULL;
     Paso_Pattern *pattern=NULL;
     Esys_MPI_rank myRank,current_rank, rank;
     Finley_IndexList* index_list=NULL;
     float *xyz=NULL;
     int c;
     
     #ifdef ESYS_MPI
     Esys_MPI_rank dest, source;
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
         dim_t *recvbuf=TMPMEMALLOC(mpiSize*mpiSize,dim_t);

         /* set the coordinates: */
         /* it is assumed that at least one node on this processor provides a coordinate */
         #pragma omp parallel for private(i,j,k)
         for (i=0;i<in->Nodes->numNodes;++i) {
             k=in->Nodes->globalDegreesOfFreedom[i]-myFirstVertex;
             if ((k>=0) && (k<myNumVertices)) {
                for (j=0;j<dim;++j) xyz[k*dim+j]=(float)(in->Nodes->Coordinates[INDEX2(j,i,dim)]); 
             }
         }

         index_list=TMPMEMALLOC(myNumVertices,Finley_IndexList);
	 /* ksteube CSR of DOF IDs */
         /* create the adjacency structure xadj and adjncy */
         if (! Finley_checkPtr(index_list)) {
            #pragma omp parallel private(i)
            {
              #pragma omp for schedule(static)
              for(i=0;i<myNumVertices;++i) {
                   index_list[i].extension=NULL;
                   index_list[i].n=0;
              }
	      /* ksteube build CSR format */
              /*  insert contributions from element matrices into columns index index_list: */
              Finley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
                                                                        in->Elements,in->Nodes->globalDegreesOfFreedom,
                                                                        in->Nodes->globalDegreesOfFreedom);
              Finley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
                                                                        in->FaceElements,in->Nodes->globalDegreesOfFreedom,
                                                                        in->Nodes->globalDegreesOfFreedom);
              Finley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
                                                                        in->ContactElements,in->Nodes->globalDegreesOfFreedom,
                                                                        in->Nodes->globalDegreesOfFreedom);
              Finley_IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list, myFirstVertex, myLastVertex,
                                                                        in->Points,in->Nodes->globalDegreesOfFreedom,
                                                                        in->Nodes->globalDegreesOfFreedom);
           }
           
           /* create the local matrix pattern */
           pattern=Finley_IndexList_createPattern(0,myNumVertices,index_list,0,globalNumVertices,0);

           /* clean up index list */
           if (index_list!=NULL) {
              #pragma omp parallel for private(i) 
              for(i=0;i<myNumVertices;++i) Finley_IndexList_free(index_list[i].extension);
           }

           if (Finley_noError()) {

#ifdef USE_PARMETIS

	      if (mpiSize>1 && 
		  Check_Inputs_For_Parmetis(mpiSize, myRank, distribution, &(in->MPIInfo->comm))>0 ) {
		 int i;
		 int wgtflag = 0;
		 int numflag = 0;	/* pattern->ptr is C style: starting from 0 instead of 1 */
		 int ncon = 1;
		 int edgecut;
		 int options[2];
		 float *tpwgts = TMPMEMALLOC(ncon*mpiSize,float);
		 float *ubvec = TMPMEMALLOC(ncon,float);
		 for (i=0; i<ncon*mpiSize; i++) tpwgts[i] = 1.0/(float)mpiSize;
		 for (i=0; i<ncon; i++) ubvec[i] = 1.05;
		 options[0] = 3;
		 options[1] = 15;
	         ParMETIS_V3_PartGeomKway(distribution,
                                          pattern->ptr,
                                          pattern->index,
                                          NULL,
                                          NULL,
                                          &wgtflag,
                                          &numflag,
                                          &dim,
                                          xyz,
                                          &ncon,
                                          &mpiSize, 
                                          tpwgts,
                                          ubvec,
                                          options,
                                          &edgecut,
                                          partition,				/* new CPU ownership of elements */
                                          &(in->MPIInfo->comm));
		 /* printf("ParMETIS number of edges cut by partitioning per processor: %d\n", edgecut/MAX(in->MPIInfo->size,1)); */
                 TMPMEMFREE(ubvec);
                 TMPMEMFREE(tpwgts);
	      } else {
                 for (i=0;i<myNumVertices;++i) partition[i]=0;		/* CPU 0 owns it */
	      }
#else
              for (i=0;i<myNumVertices;++i) partition[i]=myRank;	/* CPU myRank owns it */
#endif

           }

           Paso_Pattern_free(pattern);
           
           /* create a new distribution and labeling of the DOF */
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
           #ifdef ESYS_MPI
	      /* recvbuf will be the concatenation of each CPU's contribution to new_distribution */
	      MPI_Allgather(new_distribution, mpiSize, MPI_INT, recvbuf, mpiSize, MPI_INT, in->MPIInfo->comm);
           #else
               for (i=0;i<mpiSize;++i) recvbuf[i]=new_distribution[i];
           #endif
           new_distribution[0]=0;
	   for (rank=0; rank<mpiSize;rank++) {
	      c=0;
              for (i=0;i<myRank;++i) c+=recvbuf[rank+mpiSize*i];
              for (i=0;i<myNumVertices;++i) {
                 if (rank==partition[i]) {
                    newGlobalDOFID[i]=new_distribution[rank]+c;
                    c++;
	         }
              }
              for (i=myRank+1;i<mpiSize;++i) c+=recvbuf[rank+mpiSize*i];
              new_distribution[rank+1]=new_distribution[rank]+c;
           }
           TMPMEMFREE(recvbuf);

           /* now the overlap needs to be created by sending the partition around*/
#ifdef ESYS_MPI
           dest=Esys_MPIInfo_mod(mpiSize, myRank + 1);
           source=Esys_MPIInfo_mod(mpiSize, myRank - 1);
#endif
           current_rank=myRank;
           #pragma omp parallel for private(i)
           for (i=0;i<in->Nodes->numNodes;++i) setNewDOFId[i]=TRUE;

           for (p=0; p< mpiSize; ++p) {

               firstVertex=distribution[current_rank];
               lastVertex=distribution[current_rank+1];
               #pragma omp parallel for private(i,j,k)
               for (i=0;i<in->Nodes->numNodes;++i) {
                   k=in->Nodes->globalDegreesOfFreedom[i];
                   if (setNewDOFId[i] && (firstVertex<=k) && (k<lastVertex)) {
                        in->Nodes->globalDegreesOfFreedom[i]=newGlobalDOFID[k-firstVertex];
                        setNewDOFId[i]=FALSE;
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
