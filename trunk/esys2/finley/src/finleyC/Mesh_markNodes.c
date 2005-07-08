/**************************************************************/

/*   Finley: Mesh */

/*   mark the used nodes with offeset: */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_markNodes(index_t* mask,index_t offset,Finley_Mesh* in,bool_t useLinear) {
          Finley_ElementFile_markNodes(mask,offset,in->Elements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->FaceElements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->ContactElements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->Points,useLinear);
}
/*
* $Log$
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

