/* $Id$ */
/**************************************************************/
/*                                                                                                         */
/*   Finley: ElementFile                                                                                   */
/*                                                                                                         */
/*  reorders the elements in the element file such that the elements are stored close to the nodes         */
/*                                                                                                         */
/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "Util.h"
#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_optimizeDistribution(Finley_ElementFile** in) {
     Finley_Util_ValueAndIndex* item_list=NULL;
     Finley_ElementFile* out=NULL;
     maybelong* index=NULL,e,i;

     if (*in != NULL) {
        if ((*in)->numElements<1) return;
        maybelong NN=(*in)->ReferenceElement->Type->numNodes;
        item_list=(Finley_Util_ValueAndIndex*) TMPMEMALLOC((*in)->numElements*sizeof(Finley_Util_ValueAndIndex));
        index=(maybelong*) TMPMEMALLOC((*in)->numElements*sizeof(maybelong));
        if (! (Finley_checkPtr(item_list) || Finley_checkPtr(index)) ) {
           out=Finley_ElementFile_alloc((*in)->ReferenceElement->Type->TypeId,(*in)->order);
           if (Finley_ErrorCode==NO_ERROR) {
               Finley_ElementFile_allocTable(out,(*in)->numElements);
               if (Finley_ErrorCode==NO_ERROR) {
                     #pragma omp parallel for private(e,i) schedule(static)
                     for (e=0;e<(*in)->numElements;e++) {
                          item_list[e].index=e;
                          item_list[e].value=(*in)->Nodes[INDEX2(0,e,NN)];
                          for (i=1;i<NN;i++) item_list[e].value=MIN(item_list[e].value,(*in)->Nodes[INDEX2(i,e,NN)]);
                     }
                     Finley_Util_sortValueAndIndex((*in)->numElements,item_list);
                     #pragma omp parallel for private(e) schedule(static)
                     for (e=0;e<(*in)->numElements;e++) index[e]=item_list[e].index;
                     Finley_ElementFile_gather(index,*in,out);
                     Finley_ElementFile_dealloc(*in);
                     *in=out;
               } else {
                    Finley_ElementFile_dealloc(out);
               }
           }
        }
        TMPMEMFREE(item_list);
        TMPMEMFREE(index);
   }
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
