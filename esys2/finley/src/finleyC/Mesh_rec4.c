/**************************************************************/

/*   Finley: generates rectangular meshes  */

/*   Generates a numElements[0] x numElements[1] mesh with first order elements (Rec4) in the rectangle */
/*   [0,Length[0]] x [0,Length[1]]. order is the desired accuracy of the integration scheme. */


/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#include "Mesh.h"
#include "RectangularMesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Rec4(int* numElements,double* Length,int* periodic, int order,int useElementsOnFace) {
  int N0,N1,NE0,NE1,i0,i1,k,node0,totalNECount,faceNECount,NDOF0,NDOF1,NFaceElements,NUMNODES;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  N0=NE0+1;
  N1=NE1+1;

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
  out=Finley_Mesh_alloc(name,2,order);
  if (Finley_ErrorCode!=NO_ERROR) return NULL;

  out->Elements=Finley_ElementFile_alloc(Rec4,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Rec4Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec4Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Line2,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Line2_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
  if (Finley_ErrorCode!=NO_ERROR) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  allocate tables: */
  
  Finley_NodeFile_allocTable(out->Nodes,N0*N1);
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (Finley_ErrorCode!=NO_ERROR) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  set nodes: */
  
  #pragma omp parallel for private(i0,i1,k) 
  for (i1=0;i1<N1;i1++) {
    for (i0=0;i0<N0;i0++) {
      k=i0+N0*i1;
      out->Nodes->Coordinates[INDEX2(0,k,2)]=DBLE(i0)/DBLE(N0-1)*Length[0];
      out->Nodes->Coordinates[INDEX2(1,k,2)]=DBLE(i1)/DBLE(N1-1)*Length[1];
      out->Nodes->Id[k]=k;
      out->Nodes->Tag[k]=0;
      out->Nodes->degreeOfFreedom[k]=(i0%NDOF0) +N0*(i1%NDOF1);
    }
  }
  /* tags for the faces: */
  if (!periodic[1]) {
    for (i0=0;i0<N0;i0++) {
        out->Nodes->Tag[i0+N0*0]+=10;
        out->Nodes->Tag[i0+N0*(N1-1)]+=20;
    }
  }
  if (!periodic[0]) {
    for (i1=0;i1<N1;i1++) {
      out->Nodes->Tag[0+N0*i1]+=1;
      out->Nodes->Tag[(N0-1)+N0*i1]+=2;
    }
  }
  
  /*   set the elements: */
  
  #pragma omp parallel for private(i0,i1,k,node0) 
  for (i1=0;i1<NE1;i1++) {
    for (i0=0;i0<NE0;i0++) {
      k=i0+NE0*i1;
      node0=i0+i1*N0;

      out->Elements->Id[k]=k;
      out->Elements->Tag[k]=0;
      out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1);

      out->Elements->Nodes[INDEX2(0,k,4)]=node0;
      out->Elements->Nodes[INDEX2(1,k,4)]=node0+1;
      out->Elements->Nodes[INDEX2(2,k,4)]=node0+N0+1;
      out->Elements->Nodes[INDEX2(3,k,4)]=node0+N0;

    }
  }
  out->Elements->numColors=COLOR_MOD(0)+3*COLOR_MOD(0)+1;
  
  /*   face elements: */
  
  if (useElementsOnFace) {
     NUMNODES=4;
  } else {
     NUMNODES=2;
  }
  totalNECount=NE0*NE1;
  faceNECount=0;
  
  /* elements on boundary 010 (x2=0): */
  if (!periodic[1]) {
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
          k=i0+faceNECount;
          node0=i0;
   
          out->FaceElements->Id[k]=i0+totalNECount;
          out->FaceElements->Tag[k]=10;
          out->FaceElements->Color[k]=i0%2;
   
          if (useElementsOnFace) {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
              out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+1;
              out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0;
          } else {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
          }
     }
     totalNECount+=NE0;
     faceNECount+=NE0;
     
     /* **  elements on boundary 020 (x2=1): */
  
     #pragma omp parallel for private(i0,k,node0) 
     for (i0=0;i0<NE0;i0++) {
         k=i0+faceNECount;
         node0=i0+(NE1-1)*N0;
   
         out->FaceElements->Id[k]=i0+totalNECount;
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Color[k]=i0%2+2;
   
         if (useElementsOnFace) {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
             out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0;
             out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
         } else {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
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
          node0=i1*N0;

          out->FaceElements->Id[k]=i1+totalNECount;
          out->FaceElements->Tag[k]=1;
          out->FaceElements->Color[k]=i1%2+4;

          if (useElementsOnFace) {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
              out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+1;
              out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0+1;
          } else {
              out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
              out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0;
   
          }
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  
     /* **  elements on boundary 002 (x1=1): */
  
     #pragma omp parallel for private(i1,k,node0) 
     for (i1=0;i1<NE1;i1++) {
          k=i1+faceNECount;
          node0=(NE0-1)+i1*N0;

          out->FaceElements->Id[k]=i1+totalNECount;
          out->FaceElements->Tag[k]=2;
          out->FaceElements->Color[k]=i1%2+6;
   
          if (useElementsOnFace) {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
             out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0;
             out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0;
          } else {
             out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
             out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
          }
     }
     totalNECount+=NE1;
     faceNECount+=NE1;
  }
  out->FaceElements->numColors=8;

  /*  face elements done: */
  
  /*   condense the nodes: */
  
  Finley_Mesh_resolveNodeIds(out);

  /* prepare mesh for further calculatuions:*/

  Finley_Mesh_prepare(out) ;

  printf("timing: mesh generation: %.4e sec\n",Finley_timer()-time0);

  if (Finley_ErrorCode!=NO_ERROR) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  return out;
}

/*
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

