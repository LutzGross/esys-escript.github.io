/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: generates rectangular meshes  */

/*   Generates a numElements[0] x numElements[1] x numElements[2] mesh with second order elements (Hex20) in the brick */
/*   [0,Length[0]] x [0,Length[1]] x [0,Length[2]]. order is the desired accuracy of the */
/*   integration scheme. */


/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$

/**************************************************************/

#include "RectangularMesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t* numElements,double* Length,bool_t* periodic,index_t order,bool_t useElementsOnFace) {
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,totalNECount,faceNECount,NDOF0,NDOF1,NDOF2,NFaceElements,NUMNODES,M0,M1,M2;
  index_t node0;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();

  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=2*NE0+1;
  N1=2*NE1+1;
  N2=2*NE2+1;

  if (N0<=MIN(N1,N2)) {
     if (N1 <= N2) {
        M0=1;
        M1=N0;
        M2=N0*N1;
     } else {
        M0=1;
        M2=N0;
        M1=N0*N2;
     }
  } else if (N1<=MIN(N2,N0)) {
     if (N2 <= N0) {
        M1=1;
        M2=N1;
        M0=N2*N1;
     } else {
        M1=1;
        M0=N1;
        M2=N1*N0;
     }
  } else {
     if (N0 <= N1) {
        M2=1;
        M0=N2;
        M1=N2*N0;
     } else {
        M2=1;
        M1=N2;
        M0=N1*N2;
     }
  }

  NFaceElements=0;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements+=2*NE1*NE2;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*NE0*NE2;
  } else {
      NDOF1=N1-1;
  }
  if (!periodic[2]) {
      NDOF2=N2;
      NFaceElements+=2*NE0*NE1;
  } else {
      NDOF2=N2-1;
  }
  
  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d x %d mesh",N0,N1,N2);
  /* TEMPFIX */
#ifndef PASO_MPI
  out=Finley_Mesh_alloc(name,3,order);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Hex20,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex20Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Hex20Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec8,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec8_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
#else
  /* TODO */
  PASO_MPI_TODO;
  out = NULL;
#endif
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  
  /*  allocate tables: */
#ifndef PASO_MPI
  Finley_NodeFile_allocTable(out->Nodes,N0*N1*N2);
#else
  /* TODO */
#endif
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  #pragma omp parallel for private(i0,i1,i2,k)
  for (i2=0;i2<N2;i2++) {
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0;i0++) {
        k=M0*i0+M1*i1+M2*i2;
        out->Nodes->Coordinates[INDEX2(0,k,3)]=DBLE(i0)/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
        out->Nodes->Id[k]=i0+N0*i1+N0*N1*i2;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=M0*(i0%NDOF0) +M1*(i1%NDOF1) +M2*(i2%NDOF2);
      }
    }
  }

  
  /* tags for the faces: */
  if (!periodic[2]) {
     for (i1=0;i1<N1;i1++) {
       for (i0=0;i0<N0;i0++) {
         out->Nodes->Tag[M0*i0+M1*i1+M2*0]+=100;
         out->Nodes->Tag[M0*i0+M1*i1+M2*(N2-1)]+=200;
       }
     }
  }
  if (!periodic[1]) {
    for (i2=0;i2<N2;i2++) {
      for (i0=0;i0<N0;i0++) {
         out->Nodes->Tag[M0*i0+M1*0+M2*i2]+=10;
         out->Nodes->Tag[M0*i0+M1*(N1-1)+M2*i2]+=20;
      }
    }
  }
  if (!periodic[0]) {
    for (i2=0;i2<N2;i2++) {
      for (i1=0;i1<N1;i1++) {
        out->Nodes->Tag[M0*0+M1*i1+M2*i2]+=1;
        out->Nodes->Tag[M0*(N0-1)+M1*i1+M2*i2]+=2;
      }
    }
  }
  
  /*   set the elements: */
  
  #pragma omp parallel for private(i0,i1,i2,k,node0) 
  for (i2=0;i2<NE2;i2++) {
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<NE0;i0++) {
        k=i0+NE0*i1+NE0*NE1*i2;
        node0=2*i0+2*i1*N0+2*N0*N1*i2;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);

        out->Elements->Nodes[INDEX2(0,k,20)]=node0;
        out->Elements->Nodes[INDEX2(1,k,20)]=node0+2;
        out->Elements->Nodes[INDEX2(2,k,20)]=node0+2*N0+2;
        out->Elements->Nodes[INDEX2(3,k,20)]=node0+2*N0;
        out->Elements->Nodes[INDEX2(4,k,20)]=node0+2*N0*N1;
        out->Elements->Nodes[INDEX2(5,k,20)]=node0+2*N0*N1+2;
        out->Elements->Nodes[INDEX2(6,k,20)]=node0+2*N0*N1+2*N0+2;
        out->Elements->Nodes[INDEX2(7,k,20)]=node0+2*N0*N1+2*N0;
        out->Elements->Nodes[INDEX2(8,k,20)]=node0+1;
        out->Elements->Nodes[INDEX2(9,k,20)]=node0+N0+2;
        out->Elements->Nodes[INDEX2(10,k,20)]=node0+2*N0+1;
        out->Elements->Nodes[INDEX2(11,k,20)]=node0+N0;
        out->Elements->Nodes[INDEX2(12,k,20)]=node0+N0*N1;
        out->Elements->Nodes[INDEX2(13,k,20)]=node0+N0*N1+2;
        out->Elements->Nodes[INDEX2(14,k,20)]=node0+N0*N1+2*N0+2;
        out->Elements->Nodes[INDEX2(15,k,20)]=node0+N0*N1+2*N0;
        out->Elements->Nodes[INDEX2(16,k,20)]=node0+2*N0*N1+1;
        out->Elements->Nodes[INDEX2(17,k,20)]=node0+2*N0*N1+N0+2;
        out->Elements->Nodes[INDEX2(18,k,20)]=node0+2*N0*N1+2*N0+1;
        out->Elements->Nodes[INDEX2(19,k,20)]=node0+2*N0*N1+N0;
      }
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0)+9*COLOR_MOD(0);
  
  /*   face elements: */
  
  if  (useElementsOnFace) {
     NUMNODES=20;
  } else {
     NUMNODES=8;
  }
  totalNECount=NE0*NE1*NE2;
  faceNECount=0;
  
  /*   these are the quadrilateral elements on boundary 1 (x3=0): */
  
  if (!periodic[2]) {
    /* **  elements on boundary 100 (x3=0): */
    #pragma omp parallel for private(i0,i1,k,node0) 
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<NE0;i0++) {
        k=i0+NE0*i1+faceNECount;
        node0=2*i0+2*i1*N0;
  
        out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
        out->FaceElements->Tag[k]=100;
        out->FaceElements->Color[k]=(i0%2)+2*(i1%2);
  
        if  (useElementsOnFace) {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2;
           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0*N1;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+2*N0;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0*N1+2;
           out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N0;
           out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N0+1;
           out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+N0+2;
           out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+1;
           out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0*N1;
           out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+N0*N1+2*N0;
           out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+N0*N1+2;
           out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+2*N0*N1+N0;
           out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
           out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N0*N1+N0+2;
           out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+2*N0*N1+1;
        } else {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2;
           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0+1;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0+2;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
        }
      }
    }
    totalNECount+=NE1*NE0;
    faceNECount+=NE1*NE0;
  
    /* **  elements on boundary 200 (x3=1) */
    #pragma omp parallel for private(i0,i1,k,node0) 
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<NE0;i0++) {
        k=i0+NE0*i1+faceNECount;
        node0=2*i0+2*i1*N0+2*N0*N1*(NE2-1);
  
        out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
        out->FaceElements->Tag[k]=200;
        out->FaceElements->Color[k]=(i0%2)+2*(i1%2)+4;

        if  (useElementsOnFace) {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0*N1;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2*N0;

           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0+2;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0;

           out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+2*N0*N1+1;
           out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N0*N1+N0+2;
           out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
           out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+2*N0*N1+N0;

           out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0*N1;
           out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+N0*N1+2;
           out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+N0*N1+2*N0;

           out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+1;
           out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+N0+2;
           out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N0+1;
           out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N0;

        } else {
           out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0*N1;
           out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2;
           out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
           out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2*N0;
           out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0*N1+1;
           out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+N0+2;
           out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
           out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0*N1+N0;
        }
      }
    }
    totalNECount+=NE1*NE0;
    faceNECount+=NE1*NE0;
  }
  if (!periodic[0]) {
     /* **  elements on boundary 001 (x1=0): */
  
     #pragma omp parallel for private(i1,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i1=0;i1<NE1;i1++) {
         k=i1+NE1*i2+faceNECount;
         node0=2*i1*N0+2*N0*N1*i2;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=1;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+8;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0+2;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+N0;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+2*N0*N1+1;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+2*N0+1;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N0+2;

         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0;
         }
       }
     }
     totalNECount+=NE1*NE2;
     faceNECount+=NE1*NE2;
  
     /* **  elements on boundary 002 (x1=1): */
  
     #pragma omp parallel for private(i1,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i1=0;i1<NE1;i1++) {
         k=i1+NE1*i2+faceNECount;
         node0=2*(NE0-1)+2*i1*N0+2*N0*N1*i2 ;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=2;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+12;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0*N1;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N0+2;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+N0*N1+2;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N0*N1+2*N0+1;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+2*N0*N1+1;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N0*N1;

         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0*N1+2;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0+2;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0*N1+2;
         }
   
       }
     }
     totalNECount+=NE1*NE2;
     faceNECount+=NE1*NE2;
  }
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x2=0): */
  
     #pragma omp parallel for private(i0,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<NE0;i0++) {
         k=i0+NE0*i2+faceNECount;
         node0=2*i0+2*N0*N1*i2;
   
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Tag[k]=10;
         out->FaceElements->Color[k]=(i2%2)+2*(i0%2)+16;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N1*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N1*N0;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N1*N0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N1*N0+2*N0;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+2*N1*N0+1;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+N1*N0;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+N0+2;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N1*N0+N0+2;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+2*N1*N0+N0;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+2*N1*N0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+N1*N0+2*N0;

         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N1*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N1*N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N1*N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1*N0;
         }
       }
     }
     totalNECount+=NE0*NE2;
     faceNECount+=NE0*NE2;
  
     /* **  elements on boundary 020 (x2=1): */
     
     #pragma omp parallel for private(i0,i2,k,node0) 
     for (i2=0;i2<NE2;i2++) {
       for (i0=0;i0<NE0;i0++) {
         k=i0+NE0*i2+faceNECount;
         node0=2*i0+2*(NE1-1)*N0+2*N0*N1*i2;
   
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Color[k]=(i2%2)+2*(i0%2)+20;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0+2;

            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0*N1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0*N1+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2;

            out->FaceElements->Nodes[INDEX2(8,k,NUMNODES)]=node0+N1*N0+2*N0;
            out->FaceElements->Nodes[INDEX2(9,k,NUMNODES)]=node0+2*N1*N0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(10,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(11,k,NUMNODES)]=node0+2*N0+1;

            out->FaceElements->Nodes[INDEX2(12,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(13,k,NUMNODES)]=node0+2*N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(14,k,NUMNODES)]=node0+2*N0*N1+N0+2;
            out->FaceElements->Nodes[INDEX2(15,k,NUMNODES)]=node0+N0+2;

            out->FaceElements->Nodes[INDEX2(16,k,NUMNODES)]=node0+N1*N0;
            out->FaceElements->Nodes[INDEX2(17,k,NUMNODES)]=node0+2*N1*N0+1;
            out->FaceElements->Nodes[INDEX2(18,k,NUMNODES)]=node0+N0*N1+2;
            out->FaceElements->Nodes[INDEX2(19,k,NUMNODES)]=node0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0*N1+2*N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0+2;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N1*N0+2*N0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N1*N0+2*N0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+2*N0+2;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0+1;
         }
       }
     }
     totalNECount+=NE0*NE2;
     faceNECount+=NE0*NE2;
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=24;

  /*  face elements done: */
  
  /*   condense the nodes: */
  
  Finley_Mesh_resolveNodeIds(out);

  /* prepare mesh for further calculatuions:*/
  Finley_Mesh_prepare(out) ;

  #ifdef Finley_TRACE
  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);
  #endif

  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  return out;
}

