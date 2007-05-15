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

/*   Finley: Mesh */

/* detects faces in the mesh that match and replaces it by step elements */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_joinFaces(Finley_Mesh* self,double safety_factor,double tolerance, bool_t optimize_labeling) {

   char error_msg[LenErrorMsg_MAX];
   index_t e0,e1,*elem1=NULL,*elem0=NULL,*elem_mask=NULL,*matching_nodes_in_elem1=NULL;
   Finley_ElementFile *newFaceElementsFile=NULL,*newContactElementsFile=NULL;
   dim_t e,i,numPairs, NN, NN_Contact,c, new_numFaceElements;
   if (self->FaceElements==NULL) return;

   if (self->FaceElements->ReferenceElement->Type->numNodesOnFace<=0) {
     sprintf(error_msg,"Finley_Mesh_joinFaces:joining faces cannot be applied to face elements of type %s",self->FaceElements->ReferenceElement->Type->Name);
     Finley_setError(TYPE_ERROR,error_msg);
     return;
   }
   if (self->ContactElements==NULL) {
     Finley_setError(TYPE_ERROR,"Finley_Mesh_joinFaces: no contact element file present.");
     return;
   }

   NN=self->FaceElements->ReferenceElement->Type->numNodes;
   NN_Contact=self->ContactElements->ReferenceElement->Type->numNodes;

   if (2*NN!=NN_Contact) {
     sprintf(error_msg,"Finley_Mesh_joinFaces:contact element file for %s cannot hold elements created from face elements %s",
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
         new_numFaceElements=0;
         /* OMP */
         for(e=0;e<self->FaceElements->numElements;e++) {
             if (elem_mask[e]>0) {
               elem_mask[new_numFaceElements]=e;
               new_numFaceElements++;
             }
         }
         /*  allocate new face element and Contact element files */
#ifndef PASO_MPI
         newContactElementsFile=Finley_ElementFile_alloc(self->ContactElements->ReferenceElement->Type->TypeId,self->ContactElements->order,self->ContactElements->reduced_order);
         newFaceElementsFile=Finley_ElementFile_alloc(self->FaceElements->ReferenceElement->Type->TypeId,self->FaceElements->order,self->FaceElements->reduced_order);
#else
  /* TODO */
  PASO_MPI_TODO;
#endif
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
            c=self->ContactElements->numElements;
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
            newContactElementsFile->isPrepared=self->FaceElements->isPrepared;
         } 
         /* set new face and Contact elements */
         if (Finley_noError()) {

            Finley_ElementFile_dealloc(self->FaceElements);
            self->FaceElements=newFaceElementsFile;
            Finley_ElementFile_prepare(&(self->FaceElements),self->Nodes->numNodes,self->Nodes->degreeOfFreedom);

            Finley_ElementFile_dealloc(self->ContactElements);
            self->ContactElements=newContactElementsFile;
            Finley_ElementFile_prepare(&(self->ContactElements),self->Nodes->numNodes,self->Nodes->degreeOfFreedom);

            Finley_Mesh_prepare(self);

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
   if (Finley_noError()) {
       if ( ! Finley_Mesh_isPrepared(self) ) {
          Finley_setError(SYSTEM_ERROR,"Mesh is not prepared for calculation. Contact the programmers.");
       }
  }
}
