/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

/**************************************************************/

/*   Finley: generates rectangular meshes */

/*   Generates a numElements[0] x numElements[1] x numElements[2] mesh with first order elements (Hex8) in the brick */
/*   [0,Length[0]] x [0,Length[1]] x [0,Length[2]]. order is the desired accuracy of the */
/*   integration scheme. */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "RectangularMesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Hex8(dim_t* numElements,double* Length,bool_t* periodic, index_t order,bool_t useElementsOnFace) {
  dim_t N0,N1,N2,NE0,NE1,NE2,i0,i1,i2,k,totalNECount,faceNECount,NDOF0,NDOF1,NDOF2,NFaceElements,NUMNODES;
  index_t node0;
  Finley_Mesh* out;
  char name[50];
  double time0=Finley_timer();
  NE0=MAX(1,numElements[0]);
  NE1=MAX(1,numElements[1]);
  NE2=MAX(1,numElements[2]);
  N0=NE0+1;
  N1=NE1+1;
  N2=NE2+1;

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
  out=Finley_Mesh_alloc(name,3,order);
  if (! Finley_noError()) return NULL;

  out->Elements=Finley_ElementFile_alloc(Hex8,out->order);
  if (useElementsOnFace) {
     out->FaceElements=Finley_ElementFile_alloc(Hex8Face,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Hex8Face_Contact,out->order);
  } else {
     out->FaceElements=Finley_ElementFile_alloc(Rec4,out->order);
     out->ContactElements=Finley_ElementFile_alloc(Rec4_Contact,out->order);
  }
  out->Points=Finley_ElementFile_alloc(Point1,out->order);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }

  
  /*  allocate tables: */
  
  Finley_NodeFile_allocTable(out->Nodes,N0*N1*N2);
  Finley_ElementFile_allocTable(out->Elements,NE0*NE1*NE2);
  Finley_ElementFile_allocTable(out->FaceElements,NFaceElements);
  if (! Finley_noError()) {
      Finley_Mesh_dealloc(out);
      return NULL;
  }
  
  /*  set nodes: */
  
  #pragma omp parallel for private(i0,i1,i2,k) 
  for (i2=0;i2<N2;i2++) {
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0;i0++) {
        k=i0+N0*i1+N0*N1*i2;
        out->Nodes->Coordinates[INDEX2(0,k,3)]=DBLE(i0)/DBLE(N0-1)*Length[0];
        out->Nodes->Coordinates[INDEX2(1,k,3)]=DBLE(i1)/DBLE(N1-1)*Length[1];
        out->Nodes->Coordinates[INDEX2(2,k,3)]=DBLE(i2)/DBLE(N2-1)*Length[2];
        out->Nodes->Id[k]=k;
        out->Nodes->Tag[k]=0;
        out->Nodes->degreeOfFreedom[k]=(i0%NDOF0) +N0*(i1%NDOF1) +N0*N1*(i2%NDOF2);
      }
    }
  }
  /* tags for the faces: */
  if (!periodic[2]) {
    for (i1=0;i1<N1;i1++) {
      for (i0=0;i0<N0;i0++) {
        out->Nodes->Tag[i0+N0*i1+N0*N1*0]+=100;
        out->Nodes->Tag[i0+N0*i1+N0*N1*(N2-1)]+=200;
      }
    }
  }
  if (!periodic[1]) {
    for (i2=0;i2<N2;i2++) {
      for (i0=0;i0<N0;i0++) {
         out->Nodes->Tag[i0+N0*0+N0*N1*i2]+=10;
         out->Nodes->Tag[i0+N0*(N1-1)+N0*N1*i2]+=20;
      }
    }
  }
  if (!periodic[0]) {
    for (i2=0;i2<N2;i2++) {
      for (i1=0;i1<N1;i1++) {
        out->Nodes->Tag[0+N0*i1+N0*N1*i2]+=1;
        out->Nodes->Tag[(N0-1)+N0*i1+N0*N1*i2]+=2;
      }
    }
  }
  
  /*   set the elements: */
  
  #pragma omp parallel for private(i0,i1,i2,k,node0) 
  for (i2=0;i2<NE2;i2++) {
    for (i1=0;i1<NE1;i1++) {
      for (i0=0;i0<NE0;i0++) {
        k=i0+NE0*i1+NE0*NE1*i2;
        node0=i0+i1*N0+N0*N1*i2;

        out->Elements->Id[k]=k;
        out->Elements->Tag[k]=0;
        out->Elements->Color[k]=COLOR_MOD(i0)+3*COLOR_MOD(i1)+9*COLOR_MOD(i2);;

        out->Elements->Nodes[INDEX2(0,k,8)]=node0;
        out->Elements->Nodes[INDEX2(1,k,8)]=node0+1;
        out->Elements->Nodes[INDEX2(2,k,8)]=node0+N0+1;
        out->Elements->Nodes[INDEX2(3,k,8)]=node0+N0;
        out->Elements->Nodes[INDEX2(4,k,8)]=node0+N0*N1;
        out->Elements->Nodes[INDEX2(5,k,8)]=node0+N0*N1+1;
        out->Elements->Nodes[INDEX2(6,k,8)]=node0+N0*N1+N0+1;
        out->Elements->Nodes[INDEX2(7,k,8)]=node0+N0*N1+N0;

      }
    }
  }
  out->Elements->minColor=0;
  out->Elements->maxColor=COLOR_MOD(0)+3*COLOR_MOD(0)+9*COLOR_MOD(0);
  
  /*   face elements: */
  
  if  (useElementsOnFace) {
     NUMNODES=8;
  } else {
     NUMNODES=4;
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
         node0=i0+i1*N0;
   
         out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
         out->FaceElements->Tag[k]=100;
         out->FaceElements->Color[k]=(i0%2)+2*(i1%2);

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0*N1+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+1;
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
         node0=i0+i1*N0+N0*N1*(NE2-1);
   
         out->FaceElements->Id[k]=i0+NE0*i1+totalNECount;
         out->FaceElements->Tag[k]=200;
         out->FaceElements->Color[k]=(i0%2)+2*(i1%2)+4;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+ N0 * N1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+ N0 * N1+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+ N0 * N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+ N0 * N1+N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+ N0 * N1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+ N0 * N1+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+ N0 * N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+ N0 * N1+N0;
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
         node0=i1*N0+N0*N1*i2;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=1;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+8;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0;
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
         node0=(NE0-1)+i1*N0+N0*N1*i2 ;
   
         out->FaceElements->Id[k]=i1+NE1*i2+totalNECount;
         out->FaceElements->Tag[k]=2;
         out->FaceElements->Color[k]=(i2%2)+2*(i1%2)+12;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0*N1+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N0*N1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0*N1+1;
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
         node0=i0+N0*N1*i2;
   
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Tag[k]=10;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+16;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1*N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*N0;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N1*N0+N0+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+N1*N0+N0;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+1;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N1*N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N1*N0;
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
         node0=i0+(NE1-1)*N0+N0*N1*i2;
   
         out->FaceElements->Tag[k]=20;
         out->FaceElements->Id[k]=i2+NE2*i0+totalNECount;
         out->FaceElements->Color[k]=(i0%2)+2*(i2%2)+20;

         if  (useElementsOnFace) {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0+1;
            out->FaceElements->Nodes[INDEX2(4,k,NUMNODES)]=node0;
            out->FaceElements->Nodes[INDEX2(5,k,NUMNODES)]=node0+N0*N1;
            out->FaceElements->Nodes[INDEX2(6,k,NUMNODES)]=node0+N0*N1+1;
            out->FaceElements->Nodes[INDEX2(7,k,NUMNODES)]=node0+1;
         } else {
            out->FaceElements->Nodes[INDEX2(0,k,NUMNODES)]=node0+N0;
            out->FaceElements->Nodes[INDEX2(1,k,NUMNODES)]=node0+N0*N1+N0;
            out->FaceElements->Nodes[INDEX2(2,k,NUMNODES)]=node0+N0*N1+N0+1;
            out->FaceElements->Nodes[INDEX2(3,k,NUMNODES)]=node0+N0+1;
         }

       }
     }
     totalNECount+=NE0*NE2;
     faceNECount+=NE0*NE2;
  }
  out->FaceElements->minColor=0;
  out->FaceElements->maxColor=23;
   
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

/*
* Revision 1.3  2005/09/01 03:31:35  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-01
*
* Revision 1.2.2.2  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.2.2.1  2005/08/24 02:02:18  gross
* timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
*
* Revision 1.2  2005/07/08 04:07:52  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:52  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

