
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: generates rectangular meshes */

/*   Generates a numElements[0] x numElements[1] mesh with second order elements (Rec8) in the rectangle */
/*   [0,Length[0]] x [0,Length[1]]. order is the desired accuracy of the integration scheme. */

/**************************************************************/

#include "RectangularMesh.h"


Finley_Mesh* Finley_RectangularMesh_Rec8(dim_t* numElements,
                                          double* Length,
                                          bool_t* periodic,
                                          index_t order, 
                                          index_t reduced_order, 
                                          bool_t useElementsOnFace,
                                          bool_t useFullElementOrder,
                                          bool_t optimize) 
{
  #define N_PER_E 2
  #define DIM 2
  dim_t N0,N1,NE0,NE1,i0,i1,k,Nstride0,Nstride1;
  dim_t totalNECount,faceNECount,NDOF0,NDOF1,NFaceElements,NN, local_NE0, local_NE1, local_N0, local_N1;
  index_t e_offset1, e_offset0, offset0, offset1, global_i0, global_i1;
  index_t node0, myRank;
  Finley_Mesh* out;
  Paso_MPIInfo *mpi_info = NULL;
  char name[50];
  double time0=Finley_timer();

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError()) {
        return NULL;
  }
  myRank=mpi_info->rank;

  /* set up the global dimensions of the mesh */

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=N_PER_E*NE0+1;
  N1=N_PER_E*NE1+1;

  /*  allocate mesh: */  
  sprintf(name,"Rectangular %d x %d mesh",N0,N1);
  out=Finley_Mesh_alloc(name,DIM,order, reduced_order, mpi_info);
  if (! Finley_noError()) { 
      Paso_MPIInfo_free( mpi_info );
      return NULL;
  }

  if (useFullElementOrder) {
     Finley_Mesh_setElements(out,Finley_ElementFile_alloc(Rec9,
                                            out->order,
                                            out->reduced_order,
                                            mpi_info));
     if (useElementsOnFace) {
         Finley_setError(SYSTEM_ERROR,"rich elements for Rec9 elements is not supported yet.");
     } else {
         Finley_Mesh_setFaceElements(out,Finley_ElementFile_alloc(Line3,
                                                    out->order,
                                                    out->reduced_order,
                                                    mpi_info));
         Finley_Mesh_setContactElements(out,Finley_ElementFile_alloc(Line3_Contact,
                                                       out->order,
                                                       out->reduced_order,
                                                       mpi_info));
     }
  } else  {
     Finley_Mesh_setElements(out,Finley_ElementFile_alloc(Rec8,out->order,out->reduced_order,mpi_info));
     if (useElementsOnFace) {
         Finley_Mesh_setFaceElements(out,Finley_ElementFile_alloc(Rec8Face,
                                                                  out->order,
                                                                  out->reduced_order,
                                                                  mpi_info));
         Finley_Mesh_setContactElements(out,Finley_ElementFile_alloc(Rec8Face_Contact,
                                                                    out->order,
                                                                    out->reduced_order,
                                                                    mpi_info));
     } else {
         Finley_Mesh_setFaceElements(out,Finley_ElementFile_alloc(Line3,
                                                                  out->order,
                                                                  out->reduced_order,
                                                                  mpi_info));
         Finley_Mesh_setContactElements(out,Finley_ElementFile_alloc(Line3_Contact,
                                                                     out->order,
                                                                     out->reduced_order,
                                                                     mpi_info));
     }
  }
  Finley_Mesh_setPoints(out,Finley_ElementFile_alloc(Point1,
                                                 out->order,
                                                 out->reduced_order,
                                                 mpi_info));
  if (! Finley_noError()) {
      Paso_MPIInfo_free( mpi_info );
      Finley_Mesh_free(out);
      return NULL;
  }

  /* work out the largest dimension */
  if (N1==MAX(N0,N1)) {
     Nstride0=1;
     Nstride1=N0;
     local_NE0=NE0;
     e_offset0=0;
     Paso_MPIInfo_Split(mpi_info,NE1,&local_NE1,&e_offset1);
  } else {
     Nstride0=N1;
     Nstride1=1;
     Paso_MPIInfo_Split(mpi_info,NE0,&local_NE0,&e_offset0);
     local_NE1=NE1;
     e_offset1=0;
  }
  offset0=e_offset0*N_PER_E;
  offset1=e_offset1*N_PER_E;
  local_N0=local_NE0*N_PER_E+1;
  local_N1=local_NE1*N_PER_E+1;

  /* get the number of surface elements */

  NFaceElements=0;
  if (!periodic[0]) {
     NDOF0=N0;
     if (e_offset0 == 0) NFaceElements+=local_NE1;
     if (local_NE0+e_offset0 == NE0) NFaceElements+=local_NE1;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
     NDOF1=N1;
     if (e_offset1 == 0) NFaceElements+=local_NE0;
     if (local_NE1+e_offset1 == NE1) NFaceElements+=local_NE0;
  } else {
      NDOF1=N1-1;
  }

  /*  allocate tables: */

  Finley_NodeFile_allocTable(out->Nodes,local_N0*local_N1);
  Finley_ElementFile_allocTable(out->Elements,local_NE0*local_NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);

  if (Finley_noError()) {
     /* create nodes */
   
     #pragma omp parallel for private(i0,i1,k,global_i0,global_i1)
     for (i1=0;i1<local_N1;i1++) {
       for (i0=0;i0<local_N0;i0++) {
           k=i0+local_N0*i1;
           global_i0=i0+offset0;
           global_i1=i1+offset1;
           out->Nodes->Coordinates[INDEX2(0,k,DIM)]=DBLE(global_i0)/DBLE(N0-1)*Length[0];
           out->Nodes->Coordinates[INDEX2(1,k,DIM)]=DBLE(global_i1)/DBLE(N1-1)*Length[1];
           out->Nodes->Id[k]=Nstride0*global_i0+Nstride1*global_i1;
           out->Nodes->Tag[k]=0;
           out->Nodes->globalDegreesOfFreedom[k]=Nstride0*(global_i0%NDOF0) 
                                               +Nstride1*(global_i1%NDOF1);
       }
     }
     /*   set the elements: */
     NN=out->Elements->numNodes;
     #pragma omp parallel for private(i0,i1,k,node0) 
     for (i1=0;i1<local_NE1;i1++) {
         for (i0=0;i0<local_NE0;i0++) {
           
           k=i0+local_NE0*i1;        
           node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1);

           out->Elements->Id[k]=(i0+e_offset0)+NE0*(i1+e_offset1);
           out->Elements->Tag[k]=0;
           out->Elements->Owner[k]=myRank;

           out->Elements->Nodes[INDEX2(0,k,NN)]=node0;
           out->Elements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride0;
           out->Elements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride1;
           out->Elements->Nodes[INDEX2(4,k,NN)]=node0+1*Nstride0;
           out->Elements->Nodes[INDEX2(5,k,NN)]=node0+Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride1+1*Nstride0;
           out->Elements->Nodes[INDEX2(7,k,NN)]=node0+Nstride1;
           if (useFullElementOrder) {
              out->Elements->Nodes[INDEX2(8,k,NN)]=node0+1*Nstride1+1*Nstride0;
           }
         }
     }
     /* face elements */
     NN=out->FaceElements->numNodes;
     totalNECount=NE0*NE1;
     faceNECount=0;
     if (!periodic[0]) {
        /* **  elements on boundary 001 (x1=0): */
     
        if (e_offset0 == 0) {
           #pragma omp parallel for private(i1,k,node0) 
           for (i1=0;i1<local_NE1;i1++) {
      
               k=i1+faceNECount;
               node0=Nstride1*N_PER_E*(i1+e_offset1);

               out->FaceElements->Id[k]=i1+e_offset1+totalNECount;
               out->FaceElements->Tag[k]=1;
               out->FaceElements->Owner[k]=myRank;
               if (useElementsOnFace) {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+2*Nstride1;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+Nstride1;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+2*Nstride1+1*Nstride0;
               } else {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+2*Nstride1;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+Nstride1;
               }
           }
           faceNECount+=local_NE1;
        }
        totalNECount+=NE1;
        /* **  elements on boundary 002 (x1=1): */
        if (local_NE0+e_offset0 == NE0) {
           #pragma omp parallel for private(i1,k,node0) 
           for (i1=0;i1<local_NE1;i1++) {
               k=i1+faceNECount;
               node0=Nstride0*N_PER_E*(NE0-1)+Nstride1*N_PER_E*(i1+e_offset1);

               out->FaceElements->Id[k]=(i1+e_offset1)+totalNECount;
               out->FaceElements->Tag[k]=2;
               out->FaceElements->Owner[k]=myRank;

               if (useElementsOnFace) {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride1;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0;
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+Nstride1;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+1*Nstride0;
               } else {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+Nstride1+2*Nstride0;
               }
           }
           faceNECount+=local_NE1;
         }
         totalNECount+=NE1;
     }
     if (!periodic[1]) {
        /* **  elements on boundary 010 (x2=0): */
        if (e_offset1 == 0) {
           #pragma omp parallel for private(i0,k,node0) 
           for (i0=0;i0<local_NE0;i0++) {
               k=i0+faceNECount;
               node0=Nstride0*N_PER_E*(i0+e_offset0);
         
               out->FaceElements->Id[k]=e_offset0+i0+totalNECount;
               out->FaceElements->Tag[k]=10;
               out->FaceElements->Owner[k]=myRank;
       
               if (useElementsOnFace) {
                   out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0;
                   out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride0;
                   out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride1+2*Nstride0;
                   out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride1;
                   out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+1*Nstride0;
                   out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+Nstride1+2*Nstride0;
                   out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride1+1*Nstride0;
                   out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+Nstride1;
               } else {
                   out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0;
                   out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride0;
                   out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+1*Nstride0;
               }
           }
           faceNECount+=local_NE0;
        }
        totalNECount+=NE0;
        /* **  elements on boundary 020 (x2=1): */
        if (local_NE1+e_offset1 == NE1) {
           #pragma omp parallel for private(i0,k,node0) 
           for (i0=0;i0<local_NE0;i0++) {
               k=i0+faceNECount;
               node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(NE1-1);
   
               out->FaceElements->Id[k]=i0+e_offset0+totalNECount;
               out->FaceElements->Tag[k]=20;
               out->FaceElements->Owner[k]=myRank;
               if (useElementsOnFace) {
                    out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+2*Nstride1+2*Nstride0;
                    out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride1;
                    out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0;
                    out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride0;
                    out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+2*Nstride1+1*Nstride0;
                    out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+Nstride1;
                    out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+1*Nstride0;
                    out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+Nstride1+2*Nstride0;
               } else {
                    out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+2*Nstride1+2*Nstride0;
                    out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride1;
                    out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride1+1*Nstride0;
               }
           }
           faceNECount+=local_NE0;
        }
        totalNECount+=NE0;
     }
     /* add tag names */
     Finley_Mesh_addTagMap(out,"top", 20);
     Finley_Mesh_addTagMap(out,"bottom", 10);
     Finley_Mesh_addTagMap(out,"left", 1);
     Finley_Mesh_addTagMap(out,"right", 2);
   
     /* prepare mesh for further calculatuions:*/
     if (Finley_noError()) {
         Finley_Mesh_resolveNodeIds(out);
     }
     if (Finley_noError()) {
         Finley_Mesh_prepare(out, optimize);
     }
  }
  if (!Finley_noError()) {
      Finley_Mesh_free(out);
  }
  /* free up memory */
  Paso_MPIInfo_free( mpi_info );  
  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  return out;
}
