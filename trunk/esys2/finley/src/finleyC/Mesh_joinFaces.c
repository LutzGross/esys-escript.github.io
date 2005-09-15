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

/*   Finley: Mesh */

/* detects faces in the mesh that match and replaces it by step elements */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_joinFaces(Finley_Mesh* self,double safety_factor,double tolerance) {

   char error_msg[LenErrorMsg_MAX];
   index_t e0,e1,*elem1=NULL,*elem0=NULL,*elem_mask=NULL,*matching_nodes_in_elem1=NULL;
   Finley_ElementFile *newFaceElementsFile=NULL,*newContactElementsFile=NULL;
   dim_t e,i,numPairs;
   if (self->FaceElements==NULL) return;

   if (self->FaceElements->ReferenceElement->Type->numNodesOnFace<=0) {
     sprintf(error_msg,"__FILE__:joining faces cannot be applied to face elements of type %s",self->FaceElements->ReferenceElement->Type->Name);
     Finley_setError(TYPE_ERROR,error_msg);
     return;
   }
   if (self->ContactElements==NULL) {
     Finley_setError(TYPE_ERROR,"__FILE__: no contact element file present.");
     return;
   }

   int NN=self->FaceElements->ReferenceElement->Type->numNodes;
   int NN_Contact=self->ContactElements->ReferenceElement->Type->numNodes;

   if (2*NN!=NN_Contact) {
     sprintf(error_msg,"__FILE__:contact element file for %s cannot hold elements created from face elements %s",
           self->ContactElements->ReferenceElement->Type->Name,self->FaceElements->ReferenceElement->Type->Name);
     Finley_setError(TYPE_ERROR,error_msg);
     return;
   }

   /* allocate work arrays */
   elem1=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   elem0=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   elem_mask=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   matching_nodes_in_elem1=TMPMEMALLOC(self->FaceElements->numElements*NN,index_t);

   if (!(Finley_checkPtr(elem1) || Finley_checkPtr(elem0) || Finley_checkPtr(elem_mask) || Finley_checkPtr(matching_nodes_in_elem1)))  {
      /* find the matching face elements */
      Finley_Mesh_findMatchingFaces(self->Nodes,self->FaceElements,safety_factor,tolerance,&numPairs,elem0,elem1,matching_nodes_in_elem1);
      if (Finley_noError()) {
         /* get a list of the face elements to be kept */
         #pragma omp parallel for private(e) schedule(static)
         for(e=0;e<self->FaceElements->numElements;e++) elem_mask[e]=1;
         for(e=0;e<numPairs;e++) {
             elem_mask[elem0[e]]=0;
             elem_mask[elem1[e]]=0;
         }
         dim_t new_numFaceElements=0;
         /* OMP */
         for(e=0;e<self->FaceElements->numElements;e++) {
             if (elem_mask[e]>0) {
               elem_mask[new_numFaceElements]=e;
               new_numFaceElements++;
             }
         }
         /*  allocate new face element and Contact element files */
         newContactElementsFile=Finley_ElementFile_alloc(self->ContactElements->ReferenceElement->Type->TypeId,self->ContactElements->order);
         newFaceElementsFile=Finley_ElementFile_alloc(self->FaceElements->ReferenceElement->Type->TypeId,self->FaceElements->order);
         if (Finley_noError()) {
               Finley_ElementFile_allocTable(newContactElementsFile,numPairs+self->ContactElements->numElements);
               Finley_ElementFile_allocTable(newFaceElementsFile,new_numFaceElements);
         }
         /* copy the old elements over */
         if (Finley_noError()) {
            /* get the face elements which are still in use:*/
            Finley_ElementFile_gather(elem_mask,self->FaceElements,newFaceElementsFile);
            /* get the Contact elements which are still in use:*/
            Finley_ElementFile_copyTable(0,newContactElementsFile,0,0,self->ContactElements);
            dim_t c=self->ContactElements->numElements;
            /* OMP */
            for (e=0;e<numPairs;e++) {
                 e0=elem0[e];
                 e1=elem1[e];
                 newContactElementsFile->Id[c]=MIN(self->FaceElements->Id[e0],self->FaceElements->Id[e1]);
                 newContactElementsFile->Tag[c]=MIN(self->FaceElements->Tag[e0],self->FaceElements->Tag[e1]);
                 newContactElementsFile->Color[c]=e;
                 for (i=0;i<NN;i++) newContactElementsFile->Nodes[INDEX2(i,c,NN_Contact)]=self->FaceElements->Nodes[INDEX2(i,e0,NN)];
                 for (i=0;i<NN;i++) newContactElementsFile->Nodes[INDEX2(i+NN,c,NN_Contact)]=matching_nodes_in_elem1[INDEX2(i,e,NN)];
                 c++;
            }
            newContactElementsFile->minColor=0;
            newContactElementsFile->maxColor=numPairs-1;
         } 
         /* set new face and Contact elements */
         if (Finley_noError()) {

            Finley_ElementFile_dealloc(self->FaceElements);
            self->FaceElements=newFaceElementsFile;
            Finley_ElementFile_prepare(&(self->FaceElements),self->Nodes->numNodes,self->Nodes->degreeOfFreedom);

            Finley_ElementFile_dealloc(self->ContactElements);
            self->ContactElements=newContactElementsFile;
            Finley_ElementFile_prepare(&(self->ContactElements),self->Nodes->numNodes,self->Nodes->degreeOfFreedom);

            Finley_Mesh_prepareNodes(self);

         } else {
            Finley_ElementFile_dealloc(newFaceElementsFile);
            Finley_ElementFile_dealloc(newContactElementsFile);
         }
      }
   }
   TMPMEMFREE(elem1);
   TMPMEMFREE(elem0);
   TMPMEMFREE(matching_nodes_in_elem1);
   TMPMEMFREE(elem_mask);
}

/*
* $Log$
* Revision 1.6  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:52  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:33  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.2  2005/06/29 02:34:52  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:14  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/

