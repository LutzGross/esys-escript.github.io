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
/*   using ParMETIS. On return the mpiRankOfDOF (includes overlap for the current distribution) */
/*   giving the new processor rank assigned to a DOF. distribution specifies the current distribution of */
/*   DOFs                                                                                                */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_optimizeDOFDistribution(Finley_Mesh* in,dim_t *distribution,Paso_MPI_rank* mpiRankOfDOF) {
     dim_t dim, i,j,k, myNumVertices,p, mpiSize, len;
     index_t myFirstVertex, myLastVertex, firstVertex, lastVertex;
     index_t* partition=NULL;
     Paso_MPI_rank myRank,dest,source,current_rank;
     float *xyz=NULL;
     #ifdef PASO_MPI
     MPI_Status status;
     #endif

     if (in==NULL) return;
     if (in->Nodes == NULL) return;

     myRank=in->MPIInfo->rank;
     mpiSize=in->MPIInfo->size;
     dim=in->Nodes->numDim;
     /* first step is to distribute the elements according to a global X of DOF */

     myFirstVertex=distribution[myRank];
     myLastVertex=distribution[myRank+1];
     myNumVertices=myLastVertex-myFirstVertex;
     len=0;
     for (p=0;p<mpiSize;++p) len=MAX(len,distribution[p+1]-distribution[p]);
     partition=TMPMEMALLOC(len,index_t); /* len is used for the sending around of partition later on */
     xyz=TMPMEMALLOC(myNumVertices*dim,float);
     if (!(Finley_checkPtr(partition) || Finley_checkPtr(xyz))) {

         /* set the coordinates: *?
         /* it is assumed that at least one node on this processor provides a coordinate */
         #pragma omp parallel for private(i,j,k);
         for (i=0;i<in->Nodes->numNodes;++i) {
             k=in->Nodes->globalDegreesOfFreedom[i]-myFirstVertex;
             if ((k>=0) && (k<myNumVertices)) {
                for (j=0;j<dim;++j) xyz[k*dim+j]=(float)(in->Nodes->Coordinates[INDEX2(j,i,dim)]); 
             }
         }
/*
         xadj
         adjncy
         

        ParMETIS_V3_PartGeomKway(distribution,
                                 xadj,
                                 adjncy,
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

           /* now the overlap needs to be created by sending the partition around*/

           dest=Paso_MPIInfo_mod(mpiSize, myRank + 1);
           source=Paso_MPIInfo_mod(mpiSize, myRank - 1);
           current_rank=myRank;
           for (p=0; p< mpiSize; ++p) {

               firstVertex=distribution[current_rank];
               lastVertex=distribution[current_rank+1];
               #pragma omp parallel for private(i,j,k);
               for (i=0;i<in->Nodes->numNodes;++i) {
                   k=in->Nodes->globalDegreesOfFreedom[i];
                   if ((firstVertex<=k) && (k<lastVertex)) mpiRankOfDOF[i]=partition[k-firstVertex];
               }

               if (p<mpiSize-1) {  /* the final send can be skipped */
                  #ifdef PASO_MPI
                  MPI_Sendrecv_replace(partition,len, MPI_INT,
                                       dest, in->MPIInfo->msg_tag_counter,
                                       source, in->MPIInfo->msg_tag_counter,
                                       in->MPIInfo->comm,&status);
                  #endif
                  in->MPIInfo->msg_tag_counter++;
                  current_rank=Paso_MPIInfo_mod(mpiSize, current_rank-1);
              }
          }
     }
     TMPMEMFREE(partition);
     TMPMEMFREE(xyz);
     return;
}
