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

/************************************************************************/

/*   Finley: Mesh: optimizes the labeling of the DOFs on each processor */

/************************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "IndexList.h"

/**************************************************************/

void Finley_Mesh_optimizeDOFLabeling(Finley_Mesh* in,dim_t *distribution) {

     index_t myFirstVertex,myLastVertex, *newGlobalDOFID=NULL, firstVertex, lastVertex;
     register index_t k;
     dim_t mpiSize, myNumVertices,len, p, i;
     Paso_Pattern *pattern=NULL;
     Paso_MPI_rank myRank,dest,source,current_rank, rank;
     Finley_IndexList* index_list=NULL;
     #ifdef PASO_MPI
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

     index_list=TMPMEMALLOC(myNumVertices,Finley_IndexList);
     newGlobalDOFID=TMPMEMALLOC(len,index_t);
     /* create the adjacency structure xadj and adjncy */
     if (! ( Finley_checkPtr(index_list) || Finley_checkPtr(newGlobalDOFID) ) ) {
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
           /* create the local matrix pattern */
           pattern=Finley_IndexList_createPattern(myNumVertices,index_list,myFirstVertex, myLastVertex,-myFirstVertex);

           /* clean up index list */
           if (index_list!=NULL) {
              #pragma omp parallel for private(i) 
              for(i=0;i<myNumVertices;++i) Finley_IndexList_free(index_list[i].extension);
           }

           if (Finley_noError()) Paso_Pattern_reduceBandwidth(pattern,newGlobalDOFID); 

           Paso_Pattern_free(pattern);
      }
      Paso_MPIInfo_noError(in->MPIInfo);
      if (Finley_noError()) {
              /* shift new labeling to create a global id */
              #pragma omp parallel for private(i)
              for (i=0;i<myNumVertices;++i) newGlobalDOFID[i]+=myFirstVertex;


              /* distribute new labeling to other processors */
              dest=Paso_MPIInfo_mod(mpiSize, myRank + 1);
              source=Paso_MPIInfo_mod(mpiSize, myRank - 1);
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
     }
     TMPMEMFREE(index_list);
     TMPMEMFREE(newGlobalDOFID);
#if 0
for (i=0;i<in->Nodes->numNodes;++i) printf("%d ",in->Nodes->globalDegreesOfFreedom[i]);
printf("\n");
#endif
     return;
}
