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

/*   at input the element nodes refers to the numbering defined the Id assigned to the nodes in the */
/*   NodeFile. At the output, the numbering of the element nodes is between 0 and numNodes */
/*   degreesOfFreedom are not neccessarily referening to a dense numbering */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/

void  Finley_Mesh_resolveNodeIds(Finley_Mesh* in) {
  char error_msg[LenErrorMsg_MAX];
  dim_t k,len,numDim,newNumNodes,n;
  index_t min_id,max_id,min_id2,max_id2,*maskNodes=NULL,*maskElements=NULL,*index=NULL;
  Finley_NodeFile *newNodeFile=NULL;
  Finley_resetError();
  numDim=Finley_Mesh_getDim(in);
  
  /*   find the minimum and maximum id used: */
  
  min_id=INDEX_T_MAX;
  max_id=-INDEX_T_MAX;
  Finley_NodeFile_setIdRange(&min_id2,&max_id2,in->Nodes);
  if (min_id2==INDEX_T_MAX || max_id2==-INDEX_T_MAX) {
    Finley_setError(VALUE_ERROR,"__FILE__: Mesh has not been defined completely.");
    goto clean;
  }

  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->Elements);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->FaceElements);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->ContactElements);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  Finley_ElementFile_setNodeRange(&min_id2,&max_id2,in->Points);
  max_id=MAX(max_id,max_id2);
  min_id=MIN(min_id,min_id2);
  #ifdef Finley_TRACE
  printf("Node id range is %d:%d\n",min_id,max_id);
  #endif
  
  /*  allocate a new node file used to gather existing node file: */
  
  len=max_id-min_id+1;
  newNodeFile=Finley_NodeFile_alloc(numDim);
  if (! Finley_noError()) goto clean;

  maskNodes=TMPMEMALLOC(len,index_t);
  if (Finley_checkPtr(maskNodes)) goto clean;

  maskElements=TMPMEMALLOC(len,index_t);
  if (Finley_checkPtr(maskElements)) goto clean;

  index=TMPMEMALLOC(in->Nodes->numNodes,index_t);
  if (Finley_checkPtr(maskElements)) goto clean;

  Finley_NodeFile_allocTable(newNodeFile,len);
  if (! Finley_noError()) goto clean;

  #pragma omp parallel for private(n) schedule(static)
  for (n=0;n<in->Nodes->numNodes;n++) index[n]=-1;
  #pragma omp parallel for private(n) schedule(static)
  for (n=0;n<len;n++) {
         maskNodes[n]=-1;
         maskElements[n]=-1;
  }
  /*  mark the nodes referred by elements in maskElements: */
  
  Finley_Mesh_markNodes(maskElements,min_id,in,FALSE);

  /*  mark defined nodes */

  #pragma omp parallel for private(k) schedule(static)
  for (k=0;k<in->Nodes->numNodes;k++) maskNodes[in->Nodes->Id[k]-min_id]=1;

  /*  check if all referenced nodes are actually defined: */

  #pragma omp parallel for private(k) schedule(static)
  for (k=0;k<len;k++) {
     /* if a node is refered by an element is there a node defined ?*/
     if (maskElements[k]>=0 && maskNodes[k]<0) {
       sprintf(error_msg,"__FILE__:Node id %d is referenced by element but is not defined.",k+min_id);
       Finley_setError(VALUE_ERROR,error_msg);
     }
  }

  if (Finley_noError()) {
  
      /*  scatter the nodefile in->nodes into newNodeFile using index; */
      #pragma omp parallel for private(k) schedule(static)
      for (k=0;k<in->Nodes->numNodes;k++) index[k]=in->Nodes->Id[k]-min_id;
      Finley_NodeFile_scatter(index,in->Nodes,newNodeFile);
  
      /*  relable used nodes: */
      /* index maps the new node labeling onto the old one */
      newNumNodes=Finley_Util_packMask(len,maskElements,index);
      #pragma omp parallel for private(k) schedule(static)
      for (k=0;k<newNumNodes;k++) maskElements[index[k]]=k;

      /*  create a new table of nodes: */
      Finley_NodeFile_deallocTable(in->Nodes);
      Finley_NodeFile_allocTable(in->Nodes,newNumNodes);
      if (! Finley_noError()) goto clean;

      /* gather the new nodefile into in->Nodes */
      Finley_NodeFile_gather(index,newNodeFile,in->Nodes);

      /*  relable nodes of the elements: */
      Finley_Mesh_relableElementNodes(maskElements,min_id,in);

  }

  /*  clean-up: */
  
  clean: TMPMEMFREE(maskNodes);
         TMPMEMFREE(maskElements);
         TMPMEMFREE(index);
         Finley_NodeFile_deallocTable(newNodeFile);
         Finley_NodeFile_dealloc(newNodeFile);
}

/*
* $Log$
* Revision 1.6  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:20  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:54  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:33  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.2  2005/06/29 02:34:53  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:14  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/

