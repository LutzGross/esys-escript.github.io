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

void Finley_Mesh_markDOFsConnectedToRange(index_t* mask, index_t offset, index_t marker, 
                                          index_t firstDOF,index_t lastDOF,Finley_Mesh* in, bool_t useLinear)
{
   index_t *dofIndex;
   if (useLinear) {
       dofIndex=in->Nodes->globalReducedDOFIndex;
   } else {
       dofIndex=in->Nodes->globalDegreesOfFreedom;
   }
   Finley_ElementFile_markDOFsConnectedToRange(mask,offset,marker,firstDOF,lastDOF,dofIndex,in->Elements,useLinear);
   Finley_ElementFile_markDOFsConnectedToRange(mask,offset,marker,firstDOF,lastDOF,dofIndex,in->FaceElements,useLinear);
   Finley_ElementFile_markDOFsConnectedToRange(mask,offset,marker,firstDOF,lastDOF,dofIndex,in->ContactElements,useLinear);
   Finley_ElementFile_markDOFsConnectedToRange(mask,offset,marker,firstDOF,lastDOF,dofIndex,in->Points,useLinear);
}
