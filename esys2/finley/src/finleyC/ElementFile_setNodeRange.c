/* $Id$ */
/**************************************************************/
/*                                                                                            */
/*   Finley: ElementFile                                                                      */
/*                                                                                            */
/*   returns the maximum and minimum node reference number of nodes describing the elements:; */
/*                                                                                            */
/*                                                                                            */
/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"
#include "Util.h"

/**************************************************************/

void Finley_ElementFile_setNodeRange(maybelong* min_id,maybelong* max_id,Finley_ElementFile* in) {
   if (in!=NULL) {
      *min_id=Finley_Util_getMinInt(in->ReferenceElement->Type->numNodes,in->numElements,in->Nodes);
      *max_id=Finley_Util_getMaxInt(in->ReferenceElement->Type->numNodes,in->numElements,in->Nodes);
   } else {
       *min_id=MAYBELONG_MAX;
       *max_id=-MAYBELONG_MAX;
   }
}
/* 
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.2  2004/07/02 04:21:13  gross
* Finley C code has been included
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
