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

void Finley_Mesh_markNodes(int* mask,int offset,Finley_Mesh* in,int useLinear) {
          Finley_ElementFile_markNodes(mask,offset,in->Elements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->FaceElements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->ContactElements,useLinear);
          Finley_ElementFile_markNodes(mask,offset,in->Points,useLinear);
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

