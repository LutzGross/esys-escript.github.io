/**************************************************************/

/*   Finley: Mesh */

/* removes matching face elements from self */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#include "Mesh.h"

/**************************************************************/


void Finley_Mesh_glueFaces(Finley_Mesh* self,double safety_factor,double tolerance) { 
   Finley_NodeFile *newNodeFile=NULL;
   Finley_ElementFile *newFaceElementsFile=NULL;
   int numPairs,*elem1=NULL,*elem0=NULL,*elem_mask=NULL,*new_node_label=NULL,*new_node_list=NULL,*new_node_mask=NULL,*matching_nodes_in_elem1=NULL;
   int e,i,n,face_node;
   if (self->FaceElements==NULL) return;

   if (self->FaceElements->ReferenceElement->Type->numNodesOnFace<=0) {
     Finley_ErrorCode=TYPE_ERROR;
     sprintf(Finley_ErrorMsg,"glueing faces cannot be applied to face elements pf type %s",self->FaceElements->ReferenceElement->Type->Name);
     return;
   }

   int NNFace=self->FaceElements->ReferenceElement->Type->numNodesOnFace;
   int NN=self->FaceElements->ReferenceElement->Type->numNodes;
   int numDim=self->Nodes->numDim;
   /* allocate work arrays */
   elem1=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements);
   elem0=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements);
   elem_mask=(int*) TMPMEMALLOC(sizeof(int)*self->FaceElements->numElements);
   matching_nodes_in_elem1=(int*) TMPMEMALLOC(sizeof(int)*sizeof(int)*self->FaceElements->numElements*NN);
   new_node_label=(int*) TMPMEMALLOC(sizeof(int)*self->Nodes->numNodes);
   new_node_list=(int*) TMPMEMALLOC(sizeof(int)*self->Nodes->numNodes);
   new_node_mask=(int*) TMPMEMALLOC(sizeof(int)*self->Nodes->numNodes);
   if (!(Finley_checkPtr(elem1) || Finley_checkPtr(elem0) || Finley_checkPtr(elem_mask) || Finley_checkPtr(new_node_label) || Finley_checkPtr(new_node_list) || Finley_checkPtr(new_node_mask) || Finley_checkPtr(matching_nodes_in_elem1)) ) {
      /* find the matching face elements */
      Finley_Mesh_findMatchingFaces(self->Nodes,self->FaceElements,safety_factor,tolerance,&numPairs,elem0,elem1,matching_nodes_in_elem1);
      if (Finley_ErrorCode==NO_ERROR) {
         for(e=0;e<self->FaceElements->numElements;e++) elem_mask[e]=0;
         for(n=0;n<self->Nodes->numNodes;n++) new_node_label[n]=n;
         /* remove mark imatching face elements to be removed */
         for(e=0;e<numPairs;e++) {
             elem_mask[elem0[e]]=1;
             elem_mask[elem1[e]]=1;
             for (i=0;i<NNFace;i++) {
                face_node=self->FaceElements->ReferenceElement->Type->faceNode[i];
                new_node_label[matching_nodes_in_elem1[INDEX2(face_node,e,NN)]]=self->FaceElements->Nodes[INDEX2(face_node,elem0[e],NN)];
             }
         }
         /* create an index of face elements */
         int new_numFaceElements=0;
         for(e=0;e<self->FaceElements->numElements;e++) {
             if (elem_mask[e]<1) {
               elem_mask[new_numFaceElements]=e;
               new_numFaceElements++;
             }
         }
         /* get the new number of nodes */
         int newNumNodes=0;
         for (n=0;n<self->Nodes->numNodes;n++) new_node_mask[n]=-1;
         for (n=0;n<self->Nodes->numNodes;n++) new_node_mask[new_node_label[n]]=1;
         for (n=0;n<self->Nodes->numNodes;n++) {
               if (new_node_mask[n]>0) {
                   new_node_mask[n]=newNumNodes;
                   new_node_list[newNumNodes]=n;
                   newNumNodes++;
               }
         }
         for (n=0;n<self->Nodes->numNodes;n++) new_node_label[n]=new_node_mask[new_node_label[n]];
         /* allocate new node and element files */
         newNodeFile=Finley_NodeFile_alloc(numDim); 
         if (Finley_ErrorCode==NO_ERROR) {
             Finley_NodeFile_allocTable(newNodeFile,newNumNodes);
             if (Finley_ErrorCode==NO_ERROR) {
                newFaceElementsFile=Finley_ElementFile_alloc(self->FaceElements->ReferenceElement->Type->TypeId,self->FaceElements->order);
                if (Finley_ErrorCode==NO_ERROR) {
                   Finley_ElementFile_allocTable(newFaceElementsFile,new_numFaceElements);
                 }
              }
         }
         if (Finley_ErrorCode==NO_ERROR) {
            /* get the new nodes :*/
            Finley_NodeFile_gather(new_node_list,self->Nodes,newNodeFile);
            /* they are the new nodes*/
            Finley_NodeFile_dealloc(self->Nodes);
            self->Nodes=newNodeFile;
            /* get the face elements which are still in use:*/
            Finley_ElementFile_gather(elem_mask,self->FaceElements,newFaceElementsFile);
            /* they are the new face elements */
            Finley_ElementFile_dealloc(self->FaceElements);
            self->FaceElements=newFaceElementsFile;
            
            /* assign new node ids to elements */
            Finley_Mesh_relableElementNodes(new_node_label,0,self);
         } else {
            Finley_NodeFile_dealloc(newNodeFile);
            Finley_ElementFile_dealloc(newFaceElementsFile);
         }
      }
   } 
   TMPMEMFREE(elem1);
   TMPMEMFREE(elem0);
   TMPMEMFREE(elem_mask);
   TMPMEMFREE(new_node_label);
   TMPMEMFREE(new_node_list);
   TMPMEMFREE(new_node_mask);
   TMPMEMFREE(matching_nodes_in_elem1);
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

