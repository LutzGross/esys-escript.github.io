/**************************************************************/

/*   Finley: Mesh */

/* detects faces in the mesh that match and replaces it by step elements */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_joinFaces(Finley_Mesh* self,double safety_factor,double tolerance) {

   int numPairs,*elem1=NULL,*elem0=NULL,*elem_mask=NULL,*matching_nodes_in_elem1=NULL;
   Finley_ElementFile *newFaceElementsFile=NULL,*newContactElementsFile=NULL;
   int e,i,e0,e1;
   if (self->FaceElements==NULL) return;

   if (self->FaceElements->ReferenceElement->Type->numNodesOnFace<=0) {
     Finley_ErrorCode=TYPE_ERROR;
     sprintf(Finley_ErrorMsg,"joining faces cannot be applied to face elements of type %s",self->FaceElements->ReferenceElement->Type->Name);
     return;
   }
   if (self->ContactElements==NULL) {
     Finley_ErrorCode=TYPE_ERROR;
     sprintf(Finley_ErrorMsg,"no contact element file present.");
     return;
   }

   int NN=self->FaceElements->ReferenceElement->Type->numNodes;
   int NN_Contact=self->ContactElements->ReferenceElement->Type->numNodes;

   if (2*NN!=NN_Contact) {
     Finley_ErrorCode=TYPE_ERROR;
     sprintf(Finley_ErrorMsg,"contact element file for %s cannot hold elements created from face elements %s",
           self->ContactElements->ReferenceElement->Type->Name,self->FaceElements->ReferenceElement->Type->Name);
     return;
   }

   /* allocate work arrays */
   elem1=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements);
   elem0=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements);
   elem_mask=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements);
   matching_nodes_in_elem1=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements*NN);

   if (!(Finley_checkPtr(elem1) || Finley_checkPtr(elem0) || Finley_checkPtr(elem_mask) || Finley_checkPtr(matching_nodes_in_elem1)))  {
      /* find the matching face elements */
      Finley_Mesh_findMatchingFaces(self->Nodes,self->FaceElements,safety_factor,tolerance,&numPairs,elem0,elem1,matching_nodes_in_elem1);
      if (Finley_ErrorCode==NO_ERROR) {
         /* get a list of the face elements to be kept */
         #pragma omp parallel for private(e) schedule(static)
         for(e=0;e<self->FaceElements->numElements;e++) elem_mask[e]=1;
         for(e=0;e<numPairs;e++) {
             elem_mask[elem0[e]]=0;
             elem_mask[elem1[e]]=0;
         }
         int new_numFaceElements=0;
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
         if (Finley_ErrorCode==NO_ERROR) {
               Finley_ElementFile_allocTable(newContactElementsFile,numPairs+self->ContactElements->numElements);
               Finley_ElementFile_allocTable(newFaceElementsFile,new_numFaceElements);
         }
         /* copy the old elements over */
         if (Finley_ErrorCode==NO_ERROR) {
            /* get the face elements which are still in use:*/
            Finley_ElementFile_gather(elem_mask,self->FaceElements,newFaceElementsFile);
            /* get the Contact elements which are still in use:*/
            Finley_ElementFile_copyTable(0,newContactElementsFile,0,0,self->ContactElements);
            int c=self->ContactElements->numElements;
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
            newContactElementsFile->numColors=numPairs;
         } 
         /* set new face and Contact elements */
         if (Finley_ErrorCode==NO_ERROR) {

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
* Revision 1.3  2004/12/15 03:48:45  jgs
* *** empty log message ***
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

