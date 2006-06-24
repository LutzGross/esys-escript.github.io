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

/*   Finley: generates rectangular meshes */

/*   Generates a numElements[0] x numElements[1] mesh with second order elements (Rec8) in the rectangle */
/*   [0,Length[0]] x [0,Length[1]]. order is the desired accuracy of the integration scheme. */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "RectangularMesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Rec8(int* numElements,double* Length,int* periodic,int order,int useElementsOnFace) {
  dim_t N0,N1,NE0,NE1,i0,i1,totalNECount,faceNECount,NDOF0,NDOF1,NFaceElements,NUMNODES,M0,M1;
  index_t k,node0;
  Finley_Mesh* out=NULL;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=2*NE0+1;
  N1=2*NE1+1;

  if (N0<=N1) {
     M0=1;
     M1=N0;
  } else {
     M0=N1;
     M1=1;
  }

  NFaceElements=0;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements+=2*NE1;
  } else {
      NDOF0=N0-1;
  }
  if (!periodic[1]) {
      NDOF1=N1;
      NFaceElements+=2*NE0;
  } else {
      NDOF1=N1-1;
  }

  /*  allocate mesh: */
  
  sprintf(name,"Rectangular %d x %d mesh",N0,N1);
    /* TEMPFIX */
#ifndef PASO_MPI
  out=Finley_Mesh_alloc(name,2,order);

  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Rec8,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec8Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec8Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line3,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Line3_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
#else
  /* TODO */
#endif
  if (! Finley_noError()) {
       Finley_Mesh_dealloc(out);
       return NULL;
  }
  
  /*  allocate tables: */
#ifndef PASO_MPI
  Finley_NodeFile_allocTable(out->Nodes,N0*N1);
#else
  /* TODO */
#endif
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  /*  set nodes: */
                                                                                                                                                                                                     
  #pragma omp parallel for private(i0,i1,k)
  for (i1=0;i1<N1;i1++) {
    for (i0=0;i0<N0;i0++) {
      k=M0*i0+M1*i1;
      out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(i0)/DBLE(N0-1)*Length[0];
      out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
      out->Nodes->Id[k]=i0+N0*i1;
      out->Nodes->Tag[k]=0;
      out->Nodes->degreeOfFreedom[k]=M0*(i0%NDOF0) +M1*(i1%NDOF1);
    }
  }
  /* tags for the faces: */
  if (!periodic[1]) {
    for (i0=0;i0<N0;i0++) {
      out->Nodes->Tag[M0*i0+M1*0]+=10;
      out->Nodes->Tag[M0*i0+M1*(N1-1)]+=20;
    }
  }
  if (!periodic[0]) {
    for (i1=0;i1<N1;i1++) {
      out->Nodes->Tag[M0*0+M1*i1]+=1;
      out->Nodes->Tag[M0*(N0-1)+M1*i1]+=2;
    }
  }

  /*   set the elements: */
  
  #pragma omp parallel for private(i0,i1,k,node0) 
  for (i1=0;i1<NE1;i1++) {
    for (i0=0;i0<NE0;i0++) {
      k=i0+NE0*i1;
      node0=2*i0+2*i1*N0;

      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);

      out->Elements->Nodes[INDEX2(0,k,8)]=node0;
      out->Elements->Nodes[INDEX2(1,k,8)]=node0+2;
      out->Elements->Nodes[INDEX2(2,k,8)]=node0+2*N0+2;
      out->Elements->Nodes[INDEX2(3,k,8)]=node0+2*N0;
      out->Elements->Nodes[INDEX2(4,k,8)]=node0+1;
      out->Elements->Nodes[INDEX2(5,k,8)]=node0+N0+2;
      out->Elements->Nodes[INDEX2(6,k,8)]=node0+2*N0+1;
      out->Elements->Nodes[INDEX2(7,k,8)]=node0+N0;

    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0);
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=8;
  } else {
     NUMNODES=3;
  }
  
  
  totalNECount=NE0*NE1;
  faceNECount=0;
  
  if (!periodic[1]) {
     /* **  elements on boundary 010 (x2=0): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
       k=i0+faceNECount;
       node0=2*i0;

       out->FaceElements->Id[k]=i0+totalNECount;
       out->FaceElements->Tag[k]=10;
       out->FaceElements->Color[k]=i0%2;

       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0+2;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+2*N0+1;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
       }
     }
     totalNECount+=NE0;
     faceNECount+=NE0;
  
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
       k=i0+faceNECount;
       node0=2*i0+2*(NE1-1)*N0;

       out->FaceElements->Id[k]=i0+totalNECount;
       out->FaceElements->Tag[k]=20;
       out->FaceElements->Color[k]=i0%2+2;
       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+2*N0+1;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+1;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0+2;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0+1;
       }
     }
     totalNECount+=NE0;
     faceNECount+=NE0;
  }
  if (!periodic[0]) {
     /* **  elements on boundary 001 (x1=0): */
  
     #pragma omp parallel for private(i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       k=i1+faceNECount;
       node0=2*i1*N0;

       out->FaceElements->Id[k]=i1+totalNECount;
       out->FaceElements->Tag[k]=1;
       out->FaceElements->Color[k]=(i1%2)+4;

       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+1;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0+2;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+2*N0+1;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0;
       }
   
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  
     /* **  elements on boundary 002 (x1=1): */
     
     #pragma omp parallel for private(i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
       k=i1+faceNECount;
       node0=2*(NE0-1)+2*i1*N0;
   
       out->FaceElements->Id[k]=i1+totalNECount;
       out->FaceElements->Tag[k]=2;
       out->FaceElements->Color[k]=(i1%2)+6;
   
       if (useElementsOnFace) {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+2*N0;
          out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
          out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0+2;
          out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+2*N0+1;
          out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0;
          out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
       } else {
          out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+2;
          out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+2*N0+2;
          out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+2;
       }
     }
     totalNECount+=NE1;
     faceNECount+=NE1;

  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=7;

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
