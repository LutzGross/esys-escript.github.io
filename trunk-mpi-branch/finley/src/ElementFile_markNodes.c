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

/*   Finley: ElementFile */

/*   mark the used nodes with offeset: */

/**************************************************************/

/*  Copyrights by ACcESS Australia 2003,2004,2005 */
/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_markNodes(index_t* mask,index_t offset,Finley_ElementFile* in,bool_t useLinear) {
   dim_t i,NN,NN2,e;
   index_t color,*lin_node,*id=NULL;
   if (in!=NULL) {
     id=TMPMEMALLOC(in->ReferenceElement->Type->numNodes, index_t);
     if (! Finley_checkPtr(id) ){
        for (i=0;i<in->ReferenceElement->Type->numNodes;i++) id[i]=i;
        if (useLinear) {
           NN=in->LinearReferenceElement->Type->numNodes;
           lin_node=in->ReferenceElement->Type->linearNodes;
        } else {
           NN=in->ReferenceElement->Type->numNodes;
           lin_node=id;
        }
        NN2=in->numNodes;
        #pragma omp parallel for private(e,i) schedule(static)
        for (e=0;e<in->numElements;e++) {
            for (i=0;i<NN;i++) {
                mask[in->Nodes[INDEX2(lin_node[i],e,NN2)]-offset]=1;
           }
        }
        TMPMEMFREE(id);
     }
   }
}
