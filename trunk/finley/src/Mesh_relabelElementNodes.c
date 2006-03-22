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

/*                                                                      */
/*   assigns new node reference numbers to elements in element file in. */
/*   if k is the old node, the new node is newNode[k-offset].           */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_relableElementNodes(index_t* newNode,index_t offset,Finley_Mesh* in) {
      Finley_ElementFile_relableNodes(newNode,offset,in->Elements);
      Finley_ElementFile_relableNodes(newNode,offset,in->FaceElements);
      Finley_ElementFile_relableNodes(newNode,offset,in->ContactElements);
      Finley_ElementFile_relableNodes(newNode,offset,in->Points);
}
/*
* $Log$
* Revision 1.3  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.2.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.2  2005/07/08 04:07:54  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:53  gross
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

