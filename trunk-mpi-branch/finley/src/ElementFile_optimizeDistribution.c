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
/*                                                                                                         */
/*   Finley: ElementFile                                                                                   */
/*                                                                                                         */
/*  reorders the elements in the element file such that the elements are stored close to the nodes         */
/*                                                                                                         */
/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Util.h"
#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_optimizeDistribution(Finley_ElementFile** in) {
#ifdef PASO_MPI
    /* TODO */
    /* We will need to take care when reordering elements, because the internal/boundary status of an element is implicit in its
       ordering.  */
    return;
#else
     Finley_Util_ValueAndIndex* item_list=NULL;
     Finley_ElementFile* out=NULL;
     dim_t e,i, NN;
     index_t *index=NULL;
     if (*in != NULL) {
        if ((*in)->numElements<1) return;
        NN=(*in)->ReferenceElement->Type->numNodes;
        item_list=TMPMEMALLOC((*in)->numElements,Finley_Util_ValueAndIndex);
        index=TMPMEMALLOC((*in)->numElements,index_t);
        if (! (Finley_checkPtr(item_list) || Finley_checkPtr(index)) ) {

           out=Finley_ElementFile_alloc((*in)->ReferenceElement->Type->TypeId,(*in)->order, (*in)->reduced_order);
           if (Finley_noError()) {
               Finley_ElementFile_allocTable(out,(*in)->numElements);
               if (Finley_noError()) {
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
#endif
}
/* 
* $Log$
* Revision 1.6  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:18  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:50  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:32  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.2  2005/06/29 02:34:50  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:13  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/
