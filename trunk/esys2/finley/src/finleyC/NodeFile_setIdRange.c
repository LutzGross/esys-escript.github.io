/* $Id$ */
/**************************************************************/

/*   Finley: Mesh: NodeFile */

/*   returns the maximum and minimum node id number of nodes: */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"
#include "Common.h"
#include "Util.h"

/**************************************************************/


void Finley_NodeFile_setIdRange(index_t* min_id,index_t* max_id,Finley_NodeFile* in) {
   *min_id=Finley_Util_getMinInt(1,in->numNodes,in->Id);
   *max_id=Finley_Util_getMaxInt(1,in->numNodes,in->Id);
}
/* 
* $Log$
* Revision 1.2  2005/07/08 04:07:56  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:54  gross
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
