/**************************************************************/

/*   Finley: Mesh */

/*   at input the element nodes refers to the numbering defined the Id assigned to the nodes in the */
/*   NodeFile. At the output, the numbering of the element nodes is between 0 and numNodes */
/*   degreesOfFreedom are not neccessarily referening to a dense numbering */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "Mesh.h"

/**************************************************************/

void  Finley_Mesh_resolveNodeIds(Finley_Mesh* in) {
  maybelong k,min_id,max_id,min_id2,max_id2,len,numDim,newNumNodes,n;
  maybelong *maskNodes=NULL,*maskElements=NULL,*index=NULL;
  Finley_NodeFile *newNodeFile=NULL;
  Finley_ErrorCode=NO_ERROR;
  numDim=Finley_Mesh_getDim(in);
  
  /*   find the minimum and maximum id used: */
  
  min_id=MAYBELONG_MAX;
  max_id=-MAYBELONG_MAX;
  Finley_NodeFile_setIdRange(&min_id2,&max_id2,in->Nodes);
  if (min_id2==MAYBELONG_MAX || max_id2==-MAYBELONG_MAX) {
    Finley_ErrorCode=VALUE_ERROR;
    sprintf(Finley_ErrorMsg,"Mesh has not been defined completely.");
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
  if (Finley_ErrorCode!=NO_ERROR) goto clean;

  maskNodes=TMPMEMALLOC(len,maybelong);
  if (Finley_checkPtr(maskNodes)) goto clean;

  maskElements=TMPMEMALLOC(len,maybelong);
  if (Finley_checkPtr(maskElements)) goto clean;

  index=TMPMEMALLOC(in->Nodes->numNodes,maybelong);
  if (Finley_checkPtr(maskElements)) goto clean;

  Finley_NodeFile_allocTable(newNodeFile,len);
  if (Finley_ErrorCode!=NO_ERROR) goto clean;

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
       Finley_ErrorCode=VALUE_ERROR;
       sprintf(Finley_ErrorMsg,"Node id %d is referenced by element but is not defined.",k+min_id);
     }
  }

  if (Finley_ErrorCode==NO_ERROR) {
  
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
      if (Finley_ErrorCode!=NO_ERROR) goto clean;

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
* Revision 1.4  2004/12/15 07:08:33  jgs
* *** empty log message ***
*
*
*
*/

