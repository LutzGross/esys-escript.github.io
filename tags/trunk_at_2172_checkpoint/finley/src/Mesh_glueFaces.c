
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: Mesh */

/* removes matching face elements from self */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/


void Finley_Mesh_glueFaces(Finley_Mesh* self,double safety_factor,double tolerance,  bool_t optimize) { 
   char error_msg[LenErrorMsg_MAX];
   Finley_NodeFile *newNodeFile=NULL;
   Finley_ElementFile *newFaceElementsFile=NULL;
   dim_t numPairs,e,i,n, NNFace, NN, numDim, new_numFaceElements, newNumNodes;
   index_t face_node, *elem1=NULL,*elem0=NULL,*elem_mask=NULL,*new_node_label=NULL,*new_node_list=NULL,*new_node_mask=NULL,*matching_nodes_in_elem1=NULL;

   if (self->MPIInfo->size>1) {
     Finley_setError(TYPE_ERROR,"Finley_Mesh_glueFaces: MPI is not supported yet.");
     return;
   }
       

   if (self->FaceElements==NULL) return;

   if (self->FaceElements->ReferenceElement->Type->numNodesOnFace<=0) {
     sprintf(error_msg,"Finley_Mesh_glueFaces:glueing faces cannot be applied to face elements of type %s",self->FaceElements->ReferenceElement->Type->Name);
     Finley_setError(TYPE_ERROR,error_msg);
     return;
   }

   NNFace=self->FaceElements->ReferenceElement->Type->numNodesOnFace;
   NN=self->FaceElements->ReferenceElement->Type->numNodes;
   numDim=self->Nodes->numDim;
   /* allocate work arrays */
   elem1=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   elem0=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   elem_mask=TMPMEMALLOC(self->FaceElements->numElements,index_t);
   matching_nodes_in_elem1=TMPMEMALLOC(self->FaceElements->numElements*NN,index_t);
   new_node_label=TMPMEMALLOC(self->Nodes->numNodes,index_t);
   new_node_list=TMPMEMALLOC(self->Nodes->numNodes,index_t);
   new_node_mask=TMPMEMALLOC(self->Nodes->numNodes,index_t);
   if (!(Finley_checkPtr(elem1) || Finley_checkPtr(elem0) || Finley_checkPtr(elem_mask) || Finley_checkPtr(new_node_label) || Finley_checkPtr(new_node_list) || Finley_checkPtr(new_node_mask) || Finley_checkPtr(matching_nodes_in_elem1)) ) {
      /* find the matching face elements */
      Finley_Mesh_findMatchingFaces(self->Nodes,self->FaceElements,safety_factor,tolerance,&numPairs,elem0,elem1,matching_nodes_in_elem1);
      if (Finley_noError()) {
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
         new_numFaceElements=0;
         for(e=0;e<self->FaceElements->numElements;e++) {
             if (elem_mask[e]<1) {
               elem_mask[new_numFaceElements]=e;
               new_numFaceElements++;
             }
         }
         /* get the new number of nodes */
         newNumNodes=0;
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
         newNodeFile=Finley_NodeFile_alloc(numDim, self->MPIInfo); 

         if (Finley_noError()) {
             Finley_NodeFile_allocTable(newNodeFile,newNumNodes);
             if (Finley_noError()) {
                newFaceElementsFile=Finley_ElementFile_alloc(self->FaceElements->ReferenceElement->Type->TypeId,
                                                             self->FaceElements->order, 
                                                             self->FaceElements->reduced_order,
                                                             self->MPIInfo);
                if (Finley_noError()) {
                   Finley_ElementFile_allocTable(newFaceElementsFile,new_numFaceElements);
                 }
              }
         }
         if (Finley_noError()) 
         {
            /* get the new nodes :*/
            Finley_NodeFile_gather(new_node_list,self->Nodes,newNodeFile);
            /* they are the new nodes*/
            Finley_NodeFile_free(self->Nodes);
            self->Nodes=newNodeFile;
            /* get the face elements which are still in use:*/
            Finley_ElementFile_gather(elem_mask,self->FaceElements,newFaceElementsFile);
            /* they are the new face elements */
            Finley_ElementFile_free(self->FaceElements);
            self->FaceElements=newFaceElementsFile;
            
            /* assign new node ids to elements */
            Finley_Mesh_relableElementNodes(new_node_label,0,self);

            Finley_Mesh_prepare(self, optimize);
         } 
         else 
         {
            Finley_NodeFile_free(newNodeFile);
            Finley_ElementFile_free(newFaceElementsFile);
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
