/**************************************************************/

/*   Finley: prints Mesh */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

/*  prints the mesh to the standarts output: */

void Finley_Mesh_print(Finley_Mesh *in) {
  int NN,i,j,numDim;

  /* write header */

  printf("Mesh name: %s\n",in->Name);
  
  /*  write nodes: */
  
  if (in->Nodes!=NULL) {
    numDim=in->Nodes->numDim;
    printf("=== %1dD-Nodes:\nnumber of nodes=%d\ndegrees of freedom=%d\nreduced degrees of freedom=%d\nreduced number of nodes=%d\n",
        numDim,in->Nodes->numNodes,in->Nodes->numDegreesOfFreedom,in->Nodes->reducedNumDegreesOfFreedom,in->Nodes->reducedNumNodes);
    printf("Id,Tag,degreeOfFreedom,reducedDegreeOfFreedom,reducedNumNodes,Coordinates\n");
    for (i=0;i<in->Nodes->numNodes;i++) {
      printf("%d,%d,%d,%d,%d,",
            in->Nodes->Id[i],in->Nodes->Tag[i],in->Nodes->degreeOfFreedom[i],in->Nodes->reducedDegreeOfFreedom[i],in->Nodes->toReduced[i]);
      for (j=0;j<numDim;j++) printf(" %20.15e",in->Nodes->Coordinates[INDEX2(j,i,numDim)]);
      printf("\n");
    }
  }
  
  /*  write elements: */

  if (in->Elements!=NULL) {
    printf( "=== %s:\nnumber of elements=%d\nnumber of colors=%d\n",
                       in->Elements->ReferenceElement->Type->Name,in->Elements->numElements,in->Elements->numColors);
    NN=in->Elements->ReferenceElement->Type->numNodes;
    if (in->Elements->numElements>0) {
       printf("Id,Tag,Color,Nodes\n");
       for (i=0;i<in->Elements->numElements;i++) {
         printf("%d,%d,%d,",in->Elements->Id[i],in->Elements->Tag[i],in->Elements->Color[i]);
         for (j=0;j<NN;j++) printf(" %d",in->Nodes->Id[in->Elements->Nodes[INDEX2(j,i,NN)]]);
         printf("\n");
       }
    }
  }

  /*  write face elements: */


  if (in->FaceElements!=NULL) {
    printf( "=== %s:\nnumber of elements=%d\nnumber of colors=%d\n",
                       in->FaceElements->ReferenceElement->Type->Name,in->FaceElements->numElements,in->FaceElements->numColors);
    NN=in->FaceElements->ReferenceElement->Type->numNodes;
    if (in->FaceElements->numElements>0) {
       printf("Id,Tag,Color,Nodes\n");
       for (i=0;i<in->FaceElements->numElements;i++) {
         printf("%d,%d,%d,",in->FaceElements->Id[i],in->FaceElements->Tag[i],in->FaceElements->Color[i]);
         for (j=0;j<NN;j++) printf(" %d",in->Nodes->Id[in->FaceElements->Nodes[INDEX2(j,i,NN)]]);
         printf("\n");
       }
    }
  }

  /*  write Contact elements : */
  if (in->ContactElements!=NULL) {
    printf( "=== %s:\nnumber of elements=%d\nnumber of colors=%d\n",
                       in->ContactElements->ReferenceElement->Type->Name,in->ContactElements->numElements,in->ContactElements->numColors);
    NN=in->ContactElements->ReferenceElement->Type->numNodes;
    if (in->ContactElements->numElements>0) {
       printf("Id,Tag,Color,Nodes\n");
       for (i=0;i<in->ContactElements->numElements;i++) {
         printf("%d,%d,%d,",in->ContactElements->Id[i],in->ContactElements->Tag[i],in->ContactElements->Color[i]);
         for (j=0;j<NN;j++) printf(" %d",in->Nodes->Id[in->ContactElements->Nodes[INDEX2(j,i,NN)]]);
         printf("\n");
       }
    }
  }
  
  /*  write points: */
  if (in->Points!=NULL) {
    printf( "=== %s:\nnumber of elements=%d\nnumber of colors=%d\n",
                       in->Points->ReferenceElement->Type->Name,in->Points->numElements,in->Points->numColors);
    NN=in->Points->ReferenceElement->Type->numNodes;
    if (in->Points->numElements>0) {
       printf("Id,Tag,Color,Nodes\n");
       for (i=0;i<in->Points->numElements;i++) {
         printf("%d,%d,%d,",in->Points->Id[i],in->Points->Tag[i],in->Points->Color[i]);
         for (j=0;j<NN;j++) printf(" %d",in->Nodes->Id[in->Points->Nodes[INDEX2(j,i,NN)]]);
         printf("\n");
       }
    }
  }
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

