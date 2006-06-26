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

/*   mark the used nodes with offeset: */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_markNodes(index_t* mask,index_t offset,Finley_Mesh* in,bool_t useLinear) {
          Finley_ElementFile_markNodes(mask,offset,in->Elements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->FaceElements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->ContactElements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->Points,useLinear);
}

#ifdef PASO_MPI

/* assumes that mask has length in->Nodes->numNodes
   hence, it assumes that nodes are dense and offset==0
   offset must ==0
   assumes mask[i]==-1, i=0:in->Nodes->numNodes 
    */
void Finley_Mesh_markOrderedNodesLocation(index_t* mask,index_t offset,Finley_Mesh* in,bool_t useLinear) 
{
  index_t *internalMask=NULL, *boundaryMask=NULL, i;

  if( offset )
  {
    Finley_setError(VALUE_ERROR, "Finley_Mesh_markNodesExternal() : yet to implement offest!=0");
    return;
  }

  internalMask = TMPMEMALLOC( in->Nodes->numNodes, index_t );
  boundaryMask = TMPMEMALLOC( in->Nodes->numNodes, index_t );
  for( i=0; i<in->Nodes->numNodes; i++ )
    internalMask[i] = boundaryMask[i] = 0;

  Finley_ElementFile_markInternalElementNodes(internalMask,offset,in->Elements,useLinear);
  Finley_ElementFile_markInternalElementNodes(internalMask,offset,in->FaceElements,useLinear);
  Finley_ElementFile_markInternalElementNodes(internalMask,offset,in->ContactElements,useLinear);
  Finley_ElementFile_markInternalElementNodes(internalMask,offset,in->Points,useLinear);

  Finley_ElementFile_markBoundaryElementNodes(boundaryMask,offset,in->Elements,useLinear);
  Finley_ElementFile_markBoundaryElementNodes(boundaryMask,offset,in->FaceElements,useLinear);
  Finley_ElementFile_markBoundaryElementNodes(boundaryMask,offset,in->ContactElements,useLinear);
  Finley_ElementFile_markBoundaryElementNodes(boundaryMask,offset,in->Points,useLinear);

  for( i=0; i<in->Nodes->numNodes; i++ )
    if( internalMask[i] && boundaryMask[i] )
      mask[i] = 2;  /* boundary */
    else if( internalMask[i] )
      mask[i] = 1;  /* internal */
    else if( boundaryMask[i] )
      mask[i] = 3;  /* external */
    else
      mask[i] = -1; /* not referenced */
}

void Finley_Mesh_markOrderedDegreesOfFreedomLocation(index_t* mask,index_t offset,Finley_Mesh* in,bool_t useLinear) 
{
  index_t *internalMask=NULL, *boundaryMask=NULL, i;

  if( offset )
  {
    Finley_setError(VALUE_ERROR, "Finley_Mesh_markNodesExternal() : yet to implement offest!=0");
    return;
  }

  internalMask = TMPMEMALLOC( in->Nodes->numDegreesOfFreedom, index_t );
  boundaryMask = TMPMEMALLOC( in->Nodes->numDegreesOfFreedom, index_t );
  for( i=0; i<in->Nodes->numDegreesOfFreedom; i++ )
    internalMask[i] = boundaryMask[i] = 0;

  Finley_ElementFile_markInternalElementDOF(internalMask,offset,in->Nodes->degreeOfFreedom,in->Elements,useLinear);
  Finley_ElementFile_markInternalElementDOF(internalMask,offset,in->Nodes->degreeOfFreedom,in->FaceElements,useLinear);
  Finley_ElementFile_markInternalElementDOF(internalMask,offset,in->Nodes->degreeOfFreedom,in->ContactElements,useLinear);
  Finley_ElementFile_markInternalElementDOF(internalMask,offset,in->Nodes->degreeOfFreedom,in->Points,useLinear);

  Finley_ElementFile_markBoundaryElementDOF(boundaryMask,offset,in->Nodes->degreeOfFreedom,in->Elements,useLinear);
  Finley_ElementFile_markBoundaryElementDOF(boundaryMask,offset,in->Nodes->degreeOfFreedom,in->FaceElements,useLinear);
  Finley_ElementFile_markBoundaryElementDOF(boundaryMask,offset,in->Nodes->degreeOfFreedom,in->ContactElements,useLinear);
  Finley_ElementFile_markBoundaryElementDOF(boundaryMask,offset,in->Nodes->degreeOfFreedom,in->Points,useLinear);

  for( i=0; i<in->Nodes->numDegreesOfFreedom; i++ )
    if( internalMask[i] && boundaryMask[i] )
      mask[i] = 2;  /* boundary */
    else if( internalMask[i] )
      mask[i] = 1;  /* internal */
    else if( boundaryMask[i] )
      mask[i] = 3;  /* external */
    else
      mask[i] = -1; /* not referenced */
}
#endif

/*
* $Log$
* Revision 1.3  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.2.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.2  2005/07/08 04:07:53  jgs
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

