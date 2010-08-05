
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Dudley: Mesh */

/* detects faces in the mesh that match and replaces it by step elements */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Dudley_Mesh_joinFaces(Dudley_Mesh* self,double safety_factor,double tolerance, bool_t optimize) {

   char error_msg[LenErrorMsg_MAX];
   index_t e0,e1,*elem1=NULL,*elem0=NULL,*elem_mask=NULL,*matching_nodes_in_elem1=NULL;
   Dudley_ElementFile *newFaceElementsFile=NULL,*newContactElementsFile=NULL;
   dim_t e,i,numPairs, NN, NN_Contact,c, new_numFaceElements;
   Dudley_ReferenceElement*  faceRefElement=NULL, *contactRefElement=NULL;

   if (self->MPIInfo->size>1) {
     Dudley_setError(TYPE_ERROR,"Dudley_Mesh_joinFaces: MPI is not supported yet.");
     return;
   }
   if (self->ContactElements==NULL) {
     Dudley_setError(TYPE_ERROR,"Dudley_Mesh_joinFaces: no contact element file present.");
     return;
   }
   if (self->FaceElements==NULL) return;
   faceRefElement= Dudley_ReferenceElementSet_borrowReferenceElement(self->FaceElements->referenceElementSet, FALSE);
   contactRefElement= Dudley_ReferenceElementSet_borrowReferenceElement(self->ContactElements->referenceElementSet, FALSE);
   

   NN=self->FaceElements->numNodes;
   NN_Contact=self->ContactElements->numNodes;

   if (faceRefElement->Type->numNodesOnFace<=0) {
     sprintf(error_msg,"Dudley_Mesh_joinFaces:joining faces cannot be applied to face elements of type %s",faceRefElement->Type->Name);
     Dudley_setError(TYPE_ERROR,error_msg);
     return;
   }


   if (contactRefElement->Type->numNodes != 2*faceRefElement->Type->numNodes) {
     sprintf(error_msg,"Dudley_Mesh_joinFaces:contact element file for %s need to hold elements created from face elements %s", contactRefElement->Type->Name,faceRefElement->Type->Name);
     Dudley_setError(TYPE_ERROR,error_msg);
     return;
   }

   /* allocate work arrays */
   elem1=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   elem0=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   elem_mask=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   matching_nodes_in_elem1=TMPMEMALLOC(self->FaceElements->numElements*NN,index_t);

   if (!(Dudley_checkPtr(elem1) || Dudley_checkPtr(elem0) || Dudley_checkPtr(elem_mask) || Dudley_checkPtr(matching_nodes_in_elem1)))  {

      /* find the matching face elements */
      Dudley_Mesh_findMatchingFaces(self->Nodes,self->FaceElements,safety_factor,tolerance,&numPairs,elem0,elem1,matching_nodes_in_elem1);
      if (Dudley_noError()) {
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
         newContactElementsFile=Dudley_ElementFile_alloc(self->ContactElements->referenceElementSet, self->MPIInfo);
         newFaceElementsFile=Dudley_ElementFile_alloc(self->FaceElements->referenceElementSet, self->MPIInfo);
         if (Dudley_noError()) {
               Dudley_ElementFile_allocTable(newContactElementsFile,numPairs+self->ContactElements->numElements);
               Dudley_ElementFile_allocTable(newFaceElementsFile,new_numFaceElements);
         }
         /* copy the old elements over */
         if (Dudley_noError()) {
            /* get the face elements which are still in use:*/
            Dudley_ElementFile_gather(elem_mask,self->FaceElements,newFaceElementsFile);
            /* get the Contact elements which are still in use:*/
            Dudley_ElementFile_copyTable(0,newContactElementsFile,0,0,self->ContactElements);
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
         } 
         /* set new face and Contact elements */
         if (Dudley_noError()) {

            Dudley_ElementFile_free(self->FaceElements);
            self->FaceElements=newFaceElementsFile;
            Dudley_ElementFile_free(self->ContactElements);
            self->ContactElements=newContactElementsFile;
            Dudley_Mesh_prepare(self, optimize);

         } else {
            Dudley_ElementFile_free(newFaceElementsFile);
            Dudley_ElementFile_free(newContactElementsFile);
         }
      }
   }
   TMPMEMFREE(elem1);
   TMPMEMFREE(elem0);
   TMPMEMFREE(matching_nodes_in_elem1);
   TMPMEMFREE(elem_mask);
}
