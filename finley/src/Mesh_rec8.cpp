
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

  Finley: generates rectangular meshes

  Generates a numElements[0] x numElements[1] mesh with second order elements
  (Rec8) in the rectangle [0,Length[0]] x [0,Length[1]].
  order is the desired accuracy of the integration scheme.

*****************************************************************************/

#include "RectangularMesh.h"

namespace finley {

Mesh* RectangularMesh_Rec8(const int* numElements, const double* Length,
                           const bool* periodic, int order, int reduced_order,
                           bool useElementsOnFace, bool useFullElementOrder,
                           bool useMacroElements, bool optimize)
{
#define N_PER_E 2
#define DIM 2
  int N0,N1,NE0,NE1,i0,i1,k,Nstride0=0,Nstride1=0;
  int totalNECount,faceNECount,NDOF0=0,NDOF1=0,NFaceElements,NN, local_NE0, local_NE1, local_N0=0, local_N1=0;
  int e_offset1, e_offset0, offset0=0, offset1=0, global_i0, global_i1;
  int node0, myRank;
  const_ReferenceElementSet_ptr refPoints, refContactElements, refFaceElements, refElements;
  Esys_MPIInfo *mpi_info = NULL;
  char name[50];
  bool generateAllNodes = useFullElementOrder || useMacroElements;
#ifdef Finley_TRACE
  double time0=timer();
#endif

  /* get MPI information */
  mpi_info = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
  if (!noError()) {
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
  Mesh* out = new Mesh(name, DIM, mpi_info);
  if (generateAllNodes) {
     /* setError(SYSTEM_ERROR,"full element order for Hex elements is not supported yet."); */
     if (useMacroElements) {
          refElements.reset(new ReferenceElementSet(Rec9Macro, order, reduced_order));
     } else {
          refElements.reset(new ReferenceElementSet(Rec9, order, reduced_order));
     }
     if (useElementsOnFace) {
         setError(SYSTEM_ERROR,"rich elements for Rec9 elements are not supported yet.");
     } else {
         if (useMacroElements) {
             refFaceElements.reset(new ReferenceElementSet(Line3Macro, order, reduced_order));
         } else {
             refFaceElements.reset(new ReferenceElementSet(Line3, order, reduced_order));
         }
         refContactElements.reset(new ReferenceElementSet(Line3_Contact, order, reduced_order));
     }

  } else  {
     refElements.reset(new ReferenceElementSet(Rec8, order, reduced_order));
     if (useElementsOnFace) {
         refFaceElements.reset(new ReferenceElementSet(Rec8Face, order, reduced_order));
         refContactElements.reset(new ReferenceElementSet(Rec8Face_Contact, order, reduced_order));

     } else {
         refFaceElements.reset(new ReferenceElementSet(Line3, order, reduced_order));
         refContactElements.reset(new ReferenceElementSet(Line3_Contact, order, reduced_order));

     }
  }
  refPoints.reset(new ReferenceElementSet(Point1, order, reduced_order));


  if (noError()) {

      out->setPoints(new ElementFile(refPoints, mpi_info));
      out->setContactElements(new ElementFile(refContactElements, mpi_info));
      out->setFaceElements(new ElementFile(refFaceElements, mpi_info));
      out->setElements(new ElementFile(refElements, mpi_info));

      /* work out the largest dimension */
      if (N1==MAX(N0,N1)) {
          Nstride0=1;
          Nstride1=N0;
          local_NE0=NE0;
          e_offset0=0;
          Esys_MPIInfo_Split(mpi_info,NE1,&local_NE1,&e_offset1);
      } else {
          Nstride0=N1;
          Nstride1=1;
          Esys_MPIInfo_Split(mpi_info,NE0,&local_NE0,&e_offset0);
          local_NE1=NE1;
          e_offset1=0;
      }
      offset0=e_offset0*N_PER_E;
      offset1=e_offset1*N_PER_E;
      local_N0=local_NE0>0 ? local_NE0*N_PER_E+1 : 0;
      local_N1=local_NE1>0 ? local_NE1*N_PER_E+1 : 0;

      /* get the number of surface elements */

      NFaceElements=0;
      if (!periodic[0] &&  (local_NE0>0)) {
          NDOF0=N0;
          if (e_offset0 == 0) NFaceElements+=local_NE1;
          if (local_NE0+e_offset0 == NE0) NFaceElements+=local_NE1;
      } else {
          NDOF0=N0-1;
      }
      if (!periodic[1] && (local_NE1>0)) {
          NDOF1=N1;
          if (e_offset1 == 0) NFaceElements+=local_NE0;
          if (local_NE1+e_offset1 == NE1) NFaceElements+=local_NE0;
      } else {
          NDOF1=N1-1;
      }

      /*  allocate tables: */
      out->Nodes->allocTable(local_N0*local_N1);
      out->Elements->allocTable(local_NE0*local_NE1);
      out->FaceElements->allocTable(NFaceElements);
  }

  if (noError()) {
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
           if (generateAllNodes) {
              out->Elements->Nodes[INDEX2(8,k,NN)]=node0+1*Nstride1+1*Nstride0;
           }
         }
     }
     /* face elements */
     NN=out->FaceElements->numNodes;
     totalNECount=NE0*NE1;
     faceNECount=0;
     if (!periodic[0] && (local_NE0>0)) {
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
     if (!periodic[1] && (local_NE1>0)) {
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
  }
  if (noError()) {
     /* add tag names */
     out->addTagMap("top", 20);
     out->addTagMap("bottom", 10);
     out->addTagMap("left", 1);
     out->addTagMap("right", 2);
   }
    // prepare mesh for further calculations
    if (noError()) {
        out->resolveNodeIds();
    }
    if (noError()) {
        out->prepare(optimize);
    }
    if (!noError()) {
        delete out;
        out=NULL;
    }
    Esys_MPIInfo_free(mpi_info);

    return out;
}

} // namespace finley

