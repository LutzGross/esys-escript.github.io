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
/*  Version: $Id: Mesh_markNodes.c 1223 2007-08-03 02:40:39Z gross $ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

Finley_Mesh_markDOFsConnectedToRange(Finley_Mesh* in, index_t* mask,
                                     index_t firstDOF,index_t lastDOF,
                                     index_t marker,index_t mask_offset,bool_t useLinear) {
          Finley_Mesh_markDOFsConnectedToRange(mask,offset,in->Elements,useLinear);
          Finley_Mesh_markDOFsConnectedToRange(mask,offset,in->FaceElements,useLinear);
          Finley_Mesh_markDOFsConnectedToRange(mask,offset,in->ContactElements,useLinear);
          Finley_Mesh_markDOFsConnectedToRange(mask,offset,in->Points,useLinear);
}
