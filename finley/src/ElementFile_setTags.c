
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

/*   Finley: Mesh: ElementFile */

/*  set tags to newTag where mask>0 */

/**************************************************************/

#include "ElementFile.h"
#include "Util.h"
#include "Assemble.h"

/**************************************************************/


void Finley_ElementFile_setTags(Finley_ElementFile* self,const int newTag, escriptDataC* mask) {
    register dim_t n,q;
    dim_t numElements, numQuad;
    register double *mask_array;
    register bool_t check;
    Finley_resetError();
    if (self==NULL) return;
    numElements=self->numElements;
    if (Finley_Assemble_reducedIntegrationOrder(mask)) {
       numQuad=self->ReferenceElementReducedOrder->numQuadNodes;
    } else {
       numQuad=self->ReferenceElement->numQuadNodes;
    }

    if (1!=getDataPointSize(mask)) {
       Finley_setError(TYPE_ERROR,"Finley_ElementFile_setTags: number of components of mask is 1.");
    } else if (!numSamplesEqual(mask,numQuad,numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_ElementFile_setTags: illegal number of samples of mask Data object");
    }

    /* now we can start */

    if (Finley_noError()) {
         if (isExpanded(mask)) {
             #pragma omp parallel for private(n,check,mask_array) schedule(static)
             for (n=0;n<numElements;n++) {
                 mask_array=getSampleData(mask,n);
                 if (mask_array[0]>0) self->Tag[n]=newTag;
             }
         } else {
             #pragma omp parallel for private(q,n,check,mask_array) schedule(static)
             for (n=0;n<numElements;n++) {
                 mask_array=getSampleData(mask,n);
                 check=FALSE;
                 for (q=0;q<numQuad;q++) check=check || mask_array[q];
                 if (check) self->Tag[n]=newTag;
             }
         }
         Finley_ElementFile_setTagsInUse(self);
    }
}
/*
* $Log$
*
*/
