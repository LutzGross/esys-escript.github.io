
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: ElementFile */

/*   mark the used nodes with offeset: */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_markNodes(index_t* mask,index_t offset,Finley_ElementFile* in,bool_t useLinear) {
   dim_t i,NN,NN2,e;
   index_t *lin_node,*id=NULL;
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

void Finley_ElementFile_markDOFsConnectedToRange(index_t* mask,index_t offset,index_t marker,index_t firstDOF,index_t lastDOF,index_t *dofIndex,Finley_ElementFile*in ,bool_t useLinear) 
{
   dim_t i,NN,NN2,e,j;
   index_t color,*lin_node,*id=NULL,k;
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
        for (color=in->minColor;color<=in->maxColor;color++) {
            #pragma omp parallel for private(e,i,j,k) schedule(static)
            for (e=0;e<in->numElements;e++) {
               if (in->Color[e]==color) {
                  for (i=0;i<NN;i++) {
                     k=dofIndex[in->Nodes[INDEX2(lin_node[i],e,NN2)]];
                     if ( (firstDOF<=k) && (k<lastDOF) ) {
                        for (j=0;j<NN;j++) mask[dofIndex[in->Nodes[INDEX2(lin_node[j],e,NN2)]]-offset]=marker;
                        break;
                     }
                  }
               }
            }
        }
     }
     TMPMEMFREE(id);
   }
}

