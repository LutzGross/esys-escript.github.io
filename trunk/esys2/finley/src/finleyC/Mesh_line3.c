/**************************************************************/

/*   Finley: generates rectangular meshes  */

/*   Generates numElements[0] mesh with second order elements (Line3) in the interval */
/*   [0,Length[0]]. order is the desired accuracy of the integration scheme. */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#include "Mesh.h"
#include "RectangularMesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Line3(int* numElements,double* Length,int* periodic,int order,int useElementsOnFace) {
  int N0,NE0,i0,k,node0,NDOF0,NFaceElements,NUMNODES;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  N0=2*NE0+1;
  if (!periodic[0]) {
      NDOF0=N0;
      NFaceElements=2;
  } else {
      NDOF0=N0-1;
      NFaceElements=0;
  }
  
  /*  allocate mesh: */
  
  sprintf(name,"Rectangular mesh with %d nodes",N0);
  out=Finley_Mesh_alloc(name,1,order);
  if (Finley_ErrorCode!=NO_ERROR) return NULL;

  out->Elements=Finley_ElementFile_alloc(Line3,out->order);
  if (useElementsOnFace) {
    out->FaceElements=Finley_ElementFile_alloc(Line3Face,out->order);
    out->ContactElements=Finley_ElementFile_alloc(Line3Face_Contact,out->order);
  } else {
    out->FaceElements=Finley_ElementFile_alloc(Point1,out->order);
    out->ContactElements=Finley_ElementFile_alloc(Point1_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
  if (Finley_ErrorCode!=NO_ERROR) {
       Finley_Mesh_dealloc(out);
       return NULL;
  }
  
  /*  allocate tables: */
  
  Finley_NodeFile_allocTable(out->Nodes,N0);
  Finley_ElementFile_allocTable(out->Elements,NE0);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (Finley_ErrorCode!=NO_ERROR) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  set nodes: */
  
  #pragma omp parallel for private(i0,k) 
  for (i0=0;i0<N0;i0++) {
     k=i0;
     out->Nodes->Coordinates[INDEX2(0,k,1)]=DBLE(i0)/DBLE(N0-1)*Length[0];
     out->Nodes->Id[k]=k;
     out->Nodes->Tag[k]=0;
     out->Nodes->degreeOfFreedom[k]=(i0%NDOF0);
  }
  if (!periodic[0]) {
     out->Nodes->Tag[0]=1;
     out->Nodes->Tag[N0-1]=2;
  }
  
  /*   set the elements: */
  
  #pragma omp parallel for private(i0,k,node0) 
  for (i0=0;i0<NE0;i0++) {
    k=i0;
    node0=2*i0;

    out->Elements->Id[k]=k;
    out->Elements->Tag[k]=0;
    out->Elements->Color[k]=COLOR_MOD(i0);

    out->Elements->Nodes[INDEX2(0,k,3)]=node0;
    out->Elements->Nodes[INDEX2(1,k,3)]=node0+2;
    out->Elements->Nodes[INDEX2(2,k,3)]=node0+1;
  }
  out->Elements->numColors=COLOR_MOD(0)+1;
  
  /*   face elements: */
  if (useElementsOnFace) {
     NUMNODES=3;
  } else {
     NUMNODES=1;
  }
  
  if (!periodic[0]) {
     out->FaceElements->Id[0]=NE0;
     out->FaceElements->Tag[0]=1;
     out->FaceElements->Color[0]=0;
     if (useElementsOnFace) {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
       out->FaceElements->Nodes[INDEX2(1,0,NUMNODES)]=2;
       out->FaceElements->Nodes[INDEX2(2,0,NUMNODES)]=1;
     } else {
       out->FaceElements->Nodes[INDEX2(0,0,NUMNODES)]=0;
     }

     out->FaceElements->Id[1]=NE0+1;
     out->FaceElements->Tag[1]=2;
     out->FaceElements->Color[1]=1;
     if (useElementsOnFace) {
        out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
        out->FaceElements->Nodes[INDEX2(1,1,NUMNODES)]=N0-3;
        out->FaceElements->Nodes[INDEX2(2,1,NUMNODES)]=N0-2;
     } else {
        out->FaceElements->Nodes[INDEX2(0,1,NUMNODES)]=N0-1;
     }
  }
  out->FaceElements->numColors=2;

  /*  face elements done: */
  
  /*   condense the nodes: */
  
  Finley_Mesh_resolveNodeIds(out);

  /* prepare mesh for further calculations:*/

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

