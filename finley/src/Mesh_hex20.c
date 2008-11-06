
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: generates rectangular meshes  */

/*   Generates a numElements[0] x numElements[1] x numElements[2] mesh with second order elements (Hex20) in the brick */
/*   [0,Length[0]] x [0,Length[1]] x [0,Length[2]]. order is the desired accuracy of the */
/*   integration scheme. */


/**************************************************************/

#include "RectangularMesh.h"

Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t* numElements,
                                          double* Length,
                                          bool_t* periodic,
                                          index_t order, 
                                          index_t reduced_order, 
                                          bool_t useElementsOnFace,
                                          bool_t useFullElementOrder,
                                          bool_t optimize) 
{
  #define N_PER_E 2
  #define DIM 3
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,Nstride0,Nstride1,Nstride2, local_NE0, local_NE1, local_NE2;
  dim_t totalNECount,faceNECount,NDOF0,NDOF1,NDOF2,NFaceElements, local_N0, local_N1, local_N2, NN;
  index_t node0, myRank, e_offset0, e_offset1, e_offset2, offset0, offset1, offset2, global_i0, global_i1, global_i2;
  Finley_Mesh* out;
  Paso_MPIInfo *mpi_info = NULL;
  char name[50];
  #ifdef Finley_TRACE
  double time0=Finley_timer();
  #endif

  /* get MPI information */
  mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );
  if (! Finley_noError()) {
        return NULL;
  }
  myRank=mpi_info->rank;

  /* set up the global dimensions of the mesh */

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=N_PER_E*NE0+1;
  N1=N_PER_E*NE1+1;
  N2=N_PER_E*NE2+1;

  /*  allocate mesh: */  
  sprintf(name,"Rectangular %d x %d x %d mesh",N0,N1,N2);
  out=Finley_Mesh_alloc(name,DIM,order, reduced_order, mpi_info);
  if (! Finley_noError()) { 
      Paso_MPIInfo_free( mpi_info );
      return NULL;
  }

  if (useFullElementOrder) {
     /* Finley_setError(SYSTEM_ERROR,"full element order for Hex elements is not supported yet."); */
     Finley_Mesh_setElements(out,Finley_ElementFile_alloc(Hex27,
                                            out->order,
                                            out->reduced_order,
                                            mpi_info));
     if (useElementsOnFace) {
         Finley_setError(SYSTEM_ERROR,"rich elements for Hex27 elements is not supported yet.");
     } else {
         Finley_Mesh_setFaceElements(out,Finley_ElementFile_alloc(Rec9,
                                                    out->order,
                                                    out->reduced_order,
                                                    mpi_info));
         Finley_Mesh_setContactElements(out,Finley_ElementFile_alloc(Rec9_Contact,
                                                       out->order,
                                                       out->reduced_order,
                                                       mpi_info));
     }

  } else  {
     Finley_Mesh_setElements(out,Finley_ElementFile_alloc(Hex20,out->order,out->reduced_order,mpi_info));
     if (useElementsOnFace) {
         Finley_Mesh_setFaceElements(out,Finley_ElementFile_alloc(Hex20Face,
                                                                  out->order,
                                                                  out->reduced_order,
                                                                  mpi_info));
         Finley_Mesh_setContactElements(out,Finley_ElementFile_alloc(Hex20Face_Contact,
                                                                    out->order,
                                                                    out->reduced_order,
                                                                    mpi_info));
     } else {
         Finley_Mesh_setFaceElements(out,Finley_ElementFile_alloc(Rec8,
                                                                  out->order,
                                                                  out->reduced_order,
                                                                  mpi_info));
         Finley_Mesh_setContactElements(out,Finley_ElementFile_alloc(Rec8_Contact,
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
  if (N2==MAX3(N0,N1,N2)) {
     Nstride0=1;
     Nstride1=N0;
     Nstride2=N0*N1;
     local_NE0=NE0;
     e_offset0=0;
     local_NE1=NE1;
     e_offset1=0;
     Paso_MPIInfo_Split(mpi_info,NE2,&local_NE2,&e_offset2);
  } else if (N1==MAX3(N0,N1,N2)) {
     Nstride0=N2;
     Nstride1=N0*N2;
     Nstride2=1;
     local_NE0=NE0;
     e_offset0=0;
     Paso_MPIInfo_Split(mpi_info,NE1,&local_NE1,&e_offset1);
     local_NE2=NE2;
     e_offset2=0;
  } else {
     Nstride0=N1*N2;
     Nstride1=1;
     Nstride2=N1;
     Paso_MPIInfo_Split(mpi_info,NE0,&local_NE0,&e_offset0);
     local_NE1=NE1;
     e_offset1=0;
     local_NE2=NE2;
     e_offset2=0;
  }
  offset0=e_offset0*N_PER_E;
  offset1=e_offset1*N_PER_E;
  offset2=e_offset2*N_PER_E;
  local_N0=local_NE0>0 ? local_NE0*N_PER_E+1 : 0;
  local_N1=local_NE1>0 ? local_NE1*N_PER_E+1 : 0;
  local_N2=local_NE0>0 ? local_NE2*N_PER_E+1 : 0;

  /* get the number of surface elements */

  NFaceElements=0;
  if (!periodic[2] && (local_NE2>0) ) {
    NDOF2=N2;
    if (offset2==0) NFaceElements+=local_NE1*local_NE0;
    if (local_NE2+e_offset2 == NE2) NFaceElements+=local_NE1*local_NE0;
  } else {
      NDOF2=N2-1;
  }
 
  if (!periodic[0] && (local_NE0>0) ) {
     NDOF0=N0;
     if (e_offset0 == 0) NFaceElements+=local_NE1*local_NE2;
     if (local_NE0+e_offset0 == NE0) NFaceElements+=local_NE1*local_NE2;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1] && (local_NE1>0) ) {
     NDOF1=N1;
     if (e_offset1 == 0) NFaceElements+=local_NE0*local_NE2;
     if (local_NE1+e_offset1 == NE1) NFaceElements+=local_NE0*local_NE2;
  } else {
      NDOF1=N1-1;
  }
  /*  allocate tables: */

  Finley_NodeFile_allocTable(out->Nodes,local_N0*local_N1*local_N2);
  Finley_ElementFile_allocTable(out->Elements,local_NE0*local_NE1*local_NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);

  if (Finley_noError()) {
     /* create nodes */
   
     #pragma omp parallel for private(i0,i1,i2,k,global_i0,global_i1,global_i2)
     for (i2=0;i2<local_N2;i2++) {
       for (i1=0;i1<local_N1;i1++) {
         for (i0=0;i0<local_N0;i0++) {
           k=i0+local_N0*i1+local_N0*local_N1*i2;
           global_i0=i0+offset0;
           global_i1=i1+offset1;
           global_i2=i2+offset2;
           out->Nodes->Coordinates[INDEX2(0,k,DIM)]=DBLE(global_i0)/DBLE(N0-1)*Length[0];
           out->Nodes->Coordinates[INDEX2(1,k,DIM)]=DBLE(global_i1)/DBLE(N1-1)*Length[1];
           out->Nodes->Coordinates[INDEX2(2,k,DIM)]=DBLE(global_i2)/DBLE(N2-1)*Length[2];
           out->Nodes->Id[k]=Nstride0*global_i0+Nstride1*global_i1+Nstride2*global_i2;
           out->Nodes->Tag[k]=0;
           out->Nodes->globalDegreesOfFreedom[k]=Nstride0*(global_i0%NDOF0) 
                                               +Nstride1*(global_i1%NDOF1) 
                                               +Nstride2*(global_i2%NDOF2);
         }
       }
     }
     /*   set the elements: */
     NN=out->Elements->numNodes;
     #pragma omp parallel for private(i0,i1,i2,k,node0) 
     for (i2=0;i2<local_NE2;i2++) {
       for (i1=0;i1<local_NE1;i1++) {
         for (i0=0;i0<local_NE0;i0++) {
           
           k=i0+local_NE0*i1+local_NE0*local_NE1*i2;        
           node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(i2+e_offset2);
   
           out->Elements->Id[k]=(i0+e_offset0)+NE0*(i1+e_offset1)+NE0*NE1*(i2+e_offset2);
           out->Elements->Tag[k]=0;
           out->Elements->Owner[k]=myRank;

           out->Elements->Nodes[INDEX2(0,k,NN)] =node0                                 ;
           out->Elements->Nodes[INDEX2(1,k,NN)] =node0+                      2*Nstride0;
           out->Elements->Nodes[INDEX2(2,k,NN)] =node0+           2*Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(3,k,NN)] =node0+           2*Nstride1;
           out->Elements->Nodes[INDEX2(4,k,NN)] =node0+2*Nstride2                      ;
           out->Elements->Nodes[INDEX2(5,k,NN)] =node0+2*Nstride2           +2*Nstride0;
           out->Elements->Nodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(7,k,NN)] =node0+2*Nstride2+2*Nstride1           ;
           out->Elements->Nodes[INDEX2(8,k,NN)] =node0+                      1*Nstride0;
           out->Elements->Nodes[INDEX2(9,k,NN)] =node0+           1*Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(10,k,NN)]=node0+           2*Nstride1+1*Nstride0;
           out->Elements->Nodes[INDEX2(11,k,NN)]=node0+           1*Nstride1           ;
           out->Elements->Nodes[INDEX2(12,k,NN)]=node0+1*Nstride2                      ;
           out->Elements->Nodes[INDEX2(13,k,NN)]=node0+1*Nstride2           +2*Nstride0;
           out->Elements->Nodes[INDEX2(14,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(15,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
           out->Elements->Nodes[INDEX2(16,k,NN)]=node0+2*Nstride2           +1*Nstride0;
           out->Elements->Nodes[INDEX2(17,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
           out->Elements->Nodes[INDEX2(18,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
           out->Elements->Nodes[INDEX2(19,k,NN)]=node0+2*Nstride2+1*Nstride1           ;
           if (useFullElementOrder) {
              out->Elements->Nodes[INDEX2(20,k,NN)]=node0+           1*Nstride1+1*Nstride0;
              out->Elements->Nodes[INDEX2(21,k,NN)]=node0+1*Nstride2           +1*Nstride0;
              out->Elements->Nodes[INDEX2(22,k,NN)]=node0+1*Nstride2+1*Nstride1+2*Nstride0;
              out->Elements->Nodes[INDEX2(23,k,NN)]=node0+1*Nstride2+2*Nstride1+1*Nstride0;
              out->Elements->Nodes[INDEX2(24,k,NN)]=node0+1*Nstride2+1*Nstride1           ;
              out->Elements->Nodes[INDEX2(25,k,NN)]=node0+2*Nstride2+1*Nstride1+1*Nstride0;
              out->Elements->Nodes[INDEX2(26,k,NN)]=node0+1*Nstride2+1*Nstride1+1*Nstride0;        
           }
         }
       }
     }
     /* face elements */
     NN=out->FaceElements->numNodes;
     totalNECount=NE0*NE1*NE2;
     faceNECount=0;
     /*   these are the quadrilateral elements on boundary 1 (x3=0): */
     if (!periodic[2] && (local_NE2>0)) {
       /* **  elements on boundary 100 (x3=0): */
       if (offset2==0) {
          #pragma omp parallel for private(i0,i1,k,node0) 
          for (i1=0;i1<local_NE1;i1++) {
            for (i0=0;i0<local_NE0;i0++) {
           
              k=i0+local_NE0*i1+faceNECount;
              node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1);
     
              out->FaceElements->Id[k]=(i0+e_offset0)+NE0*(i1+e_offset1)+totalNECount;
              out->FaceElements->Tag[k]=100;
              out->FaceElements->Owner[k]=myRank;
        
              if  (useElementsOnFace) {
                 out->FaceElements->Nodes[INDEX2(0,k,NN)] =node0                                 ;
                 out->FaceElements->Nodes[INDEX2(1,k,NN)] =node0           +2*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(2,k,NN)] =node0           +2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(3,k,NN)] =node0+           2*Nstride0           ;
                 out->FaceElements->Nodes[INDEX2(4,k,NN)] =node0+2*Nstride2                      ;
                 out->FaceElements->Nodes[INDEX2(5,k,NN)] =node0+2*Nstride2+2*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(7,k,NN)] =node0+2*Nstride2           +2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+           1*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(9,k,NN)] =node0+           2*Nstride1+1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(10,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(11,k,NN)]=node0+                      1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(12,k,NN)]=node0+1*Nstride2                      ;
                 out->FaceElements->Nodes[INDEX2(13,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(14,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(15,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(16,k,NN)]=node0+2*Nstride2+1*Nstride1;
                 out->FaceElements->Nodes[INDEX2(17,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(18,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(19,k,NN)]=node0+2*Nstride2           +1*Nstride0;
              } else {
                 out->FaceElements->Nodes[INDEX2(0,k,NN)] =node0                                 ;
                 out->FaceElements->Nodes[INDEX2(1,k,NN)] =node0+           2*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(2,k,NN)] =node0+           2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(3,k,NN)] =node0+                      2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(4,k,NN)] =node0+           1*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(5,k,NN)] =node0+           2*Nstride1+1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(6,k,NN)] =node0+           1*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(7,k,NN)] =node0+                      1*Nstride0;
                 if (useFullElementOrder){
                    out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+           1*Nstride1+1*Nstride0;
                 }
              }
            }
          }
          faceNECount+=local_NE1*local_NE0;
       }
       totalNECount+=NE1*NE0;
       /* **  elements on boundary 200 (x3=1) */
       if (local_NE2+e_offset2 == NE2) {
          #pragma omp parallel for private(i0,i1,k,node0) 
          for (i1=0;i1<local_NE1;i1++) {
            for (i0=0;i0<local_NE0;i0++) {
      
              k=i0+local_NE0*i1+faceNECount;
              node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(NE2-1);
        
              out->FaceElements->Id[k]=(i0+e_offset0)+NE0*(i1+e_offset1)+totalNECount;
              out->FaceElements->Tag[k]=200;
              out->FaceElements->Owner[k]=myRank;
              if  (useElementsOnFace) {
                 out->FaceElements->Nodes[INDEX2(0,k,NN)] =node0+2*Nstride2                      ;
                 out->FaceElements->Nodes[INDEX2(1,k,NN)] =node0+2*Nstride2+           2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(3,k,NN)] =node0+2*Nstride2+2*Nstride1           ;
      
                 out->FaceElements->Nodes[INDEX2(4,k,NN)] =node0                                 ;
                 out->FaceElements->Nodes[INDEX2(5,k,NN)] =node0+2*Nstride0                      ;
                 out->FaceElements->Nodes[INDEX2(6,k,NN)] =node0+           2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(7,k,NN)] =node0+           2*Nstride1;
      
                 out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+2*Nstride2+           1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(9,k,NN)] =node0+2*Nstride2+1*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(10,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(11,k,NN)]=node0+2*Nstride2+1*Nstride1           ;
      
                 out->FaceElements->Nodes[INDEX2(12,k,NN)]=node0+1*Nstride2;
                 out->FaceElements->Nodes[INDEX2(13,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(14,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(15,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
      
                 out->FaceElements->Nodes[INDEX2(16,k,NN)]=node0+                      1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(17,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(18,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(19,k,NN)]=node0+           1*Nstride1           ;
      
              } else {
                 out->FaceElements->Nodes[INDEX2(0,k,NN)] =node0+2*Nstride2                      ;
                 out->FaceElements->Nodes[INDEX2(1,k,NN)] =node0+2*Nstride2           +2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(3,k,NN)] =node0+2*Nstride2+2*Nstride1           ;
                 out->FaceElements->Nodes[INDEX2(4,k,NN)] =node0+2*Nstride2           +1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(5,k,NN)] =node0+2*Nstride2+1*Nstride1+2*Nstride0;
                 out->FaceElements->Nodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+1*Nstride0;
                 out->FaceElements->Nodes[INDEX2(7,k,NN)] =node0+2*Nstride2+1*Nstride1           ;
                 if (useFullElementOrder){
                 out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+2*Nstride2+1*Nstride1+1*Nstride0;
                 }
              }
            }
          }
          faceNECount+=local_NE1*local_NE0;
       }
       totalNECount+=NE1*NE0;
     }
     if (!periodic[0] && (local_NE0>0)) {
        /* **  elements on boundary 001 (x1=0): */
     
        if (e_offset0 == 0) {
           #pragma omp parallel for private(i1,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i1=0;i1<local_NE1;i1++) {
      
               k=i1+local_NE1*i2+faceNECount;
               node0=Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(i2+e_offset2);
               out->FaceElements->Id[k]=(i1+e_offset1)+NE1*(i2+e_offset2)+totalNECount;
               out->FaceElements->Tag[k]=1;
               out->FaceElements->Owner[k]=myRank;
      
               if  (useElementsOnFace) {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)] =node0                                 ;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)] =node0+2*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)] =node0+2*Nstride1                      ;
      
                  out->FaceElements->Nodes[INDEX2(4,k,NN)] =node0+2*Nstride0                      ;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)] =node0+2*Nstride2+2*Nstride0           ;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)] =node0+2*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)] =node0+2*Nstride1+2*Nstride0           ;
      
                  out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+1*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(9,k,NN)] =node0+2*Nstride2+1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(10,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(11,k,NN)]=node0+           1*Nstride1           ;
      
                  out->FaceElements->Nodes[INDEX2(12,k,NN)]=node0+                      1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(13,k,NN)]=node0+2*Nstride2           +1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(14,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(15,k,NN)]=node0+2*Nstride1+           1*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(16,k,NN)]=node0+1*Nstride2+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(17,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(18,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(19,k,NN)]=node0+1*Nstride1+           2*Nstride0;
      
               } else {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)] =node0                                 ;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)] =node0+2*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)] =node0+2*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)] =node0+           2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(4,k,NN)] =node0+1*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)] =node0+2*Nstride2+1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)] =node0+1*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)] =node0+           1*Nstride1           ;
                 if (useFullElementOrder){
                    out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+1*Nstride2+1*Nstride1           ;
                 }
               }
             }
           }
           faceNECount+=local_NE1*local_NE2;
        }
        totalNECount+=NE1*NE2;
     
        /* **  elements on boundary 002 (x1=1): */
        if (local_NE0+e_offset0 == NE0) {
           #pragma omp parallel for private(i1,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i1=0;i1<local_NE1;i1++) {
               k=i1+local_NE1*i2+faceNECount;
               node0=Nstride0*N_PER_E*(NE0-1)+Nstride1*N_PER_E*(i1+e_offset1)+Nstride2*N_PER_E*(i2+e_offset2);
               out->FaceElements->Id[k]=(i1+e_offset1)+NE1*(i2+e_offset2)+totalNECount;
               out->FaceElements->Tag[k]=2;
               out->FaceElements->Owner[k]=myRank;
   
               if  (useElementsOnFace) {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+                      2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+           2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride2+           2*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0                                 ;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+           2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+2*Nstride2                      ;
      
                  out->FaceElements->Nodes[INDEX2(8,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(9,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(10,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(11,k,NN)]=node0+1*Nstride2+           2*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(12,k,NN)]=node0+                      1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(13,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(14,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(15,k,NN)]=node0+2*Nstride2+           1*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(16,k,NN)]=node0+           1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(17,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(18,k,NN)]=node0+2*Nstride2+1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(19,k,NN)]=node0+1*Nstride2                      ;
      
               } else {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0                      +2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+           2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride2+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                 if (useFullElementOrder){
                    out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+1*Nstride2+1*Nstride1+2*Nstride0;
                 }
               }
         
             }
           }
           faceNECount+=local_NE1*local_NE2;
         }
         totalNECount+=NE1*NE2;
     }
     if (!periodic[1] && (local_NE1>0)) {
        /* **  elements on boundary 010 (x2=0): */
        if (e_offset1 == 0) {
           #pragma omp parallel for private(i0,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i0=0;i0<local_NE0;i0++) {
               k=i0+local_NE0*i2+faceNECount;
               node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride2*N_PER_E*(i2+e_offset2);
         
               out->FaceElements->Id[k]=(i2+e_offset2)+NE2*(e_offset0+i0)+totalNECount;
               out->FaceElements->Tag[k]=10;
               out->FaceElements->Owner[k]=myRank;
               if  (useElementsOnFace) {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0                                 ;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+                      2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride2           +2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride2                      ;
      
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+           2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+2*Nstride1+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+2*Nstride2+2*Nstride1           ;
      
                  out->FaceElements->Nodes[INDEX2(8,k,NN)]=node0+                      1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(9,k,NN)]=node0+1*Nstride2+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(10,k,NN)]=node0+2*Nstride2+           1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(11,k,NN)]=node0+1*Nstride2                      ;
      
                  out->FaceElements->Nodes[INDEX2(12,k,NN)]=node0+           1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(13,k,NN)]=node0+           1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(14,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(15,k,NN)]=node0+2*Nstride2+1*Nstride1           ;
   
                  out->FaceElements->Nodes[INDEX2(16,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(17,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(18,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(19,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
      
               } else {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0                                 ;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+                      2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride2+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+                      1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+1*Nstride2+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride2+           1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+1*Nstride2                      ;
                 if (useFullElementOrder){
                    out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+1*Nstride2+         1*Nstride0;
                 }
               }
             }
           }
           faceNECount+=local_NE0*local_NE2;
        }
        totalNECount+=NE0*NE2;
        /* **  elements on boundary 020 (x2=1): */
        if (local_NE1+e_offset1 == NE1) {
           #pragma omp parallel for private(i0,i2,k,node0) 
           for (i2=0;i2<local_NE2;i2++) {
             for (i0=0;i0<local_NE0;i0++) {
               k=i0+local_NE0*i2+faceNECount;
               node0=Nstride0*N_PER_E*(i0+e_offset0)+Nstride1*N_PER_E*(NE1-1)+Nstride2*N_PER_E*(i2+e_offset2);
   
               out->FaceElements->Id[k]=(i2+e_offset2)+NE2*(i0+e_offset0)+totalNECount;
               out->FaceElements->Tag[k]=20;
               out->FaceElements->Owner[k]=myRank;
      
               if  (useElementsOnFace) {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+           2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+2*Nstride1+2*Nstride0           ;
      
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0                                 ;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+2*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+2*Nstride2+           2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+                      2*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(8,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(9,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(10,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(11,k,NN)]=node0+           2*Nstride1+1*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(12,k,NN)]=node0+           1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(13,k,NN)]=node0+2*Nstride2+1*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(14,k,NN)]=node0+2*Nstride2+1*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(15,k,NN)]=node0+           1*Nstride1+2*Nstride0;
      
                  out->FaceElements->Nodes[INDEX2(16,k,NN)]=node0+1*Nstride2                      ;
                  out->FaceElements->Nodes[INDEX2(17,k,NN)]=node0+2*Nstride2           +1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(18,k,NN)]=node0+1*Nstride2           +2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(19,k,NN)]=node0+                      1*Nstride0;
               } else {
                  out->FaceElements->Nodes[INDEX2(0,k,NN)]=node0+           2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(1,k,NN)]=node0+2*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(2,k,NN)]=node0+2*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(3,k,NN)]=node0+           2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(4,k,NN)]=node0+1*Nstride2+2*Nstride1           ;
                  out->FaceElements->Nodes[INDEX2(5,k,NN)]=node0+2*Nstride2+2*Nstride1+1*Nstride0;
                  out->FaceElements->Nodes[INDEX2(6,k,NN)]=node0+1*Nstride2+2*Nstride1+2*Nstride0;
                  out->FaceElements->Nodes[INDEX2(7,k,NN)]=node0+           2*Nstride1+1*Nstride0;
                 if (useFullElementOrder){
                    out->FaceElements->Nodes[INDEX2(8,k,NN)] =node0+1*Nstride2+2*Nstride1+1*Nstride0;
                 }
               }
             }
           }
           faceNECount+=local_NE0*local_NE2;
        }
        totalNECount+=NE0*NE2;
     }
     /* add tag names */
     Finley_Mesh_addTagMap(out,"top", 200);
     Finley_Mesh_addTagMap(out,"bottom", 100);
     Finley_Mesh_addTagMap(out,"left", 1);
     Finley_Mesh_addTagMap(out,"right", 2);
     Finley_Mesh_addTagMap(out,"front", 10);
     Finley_Mesh_addTagMap(out,"back", 20);
   
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
