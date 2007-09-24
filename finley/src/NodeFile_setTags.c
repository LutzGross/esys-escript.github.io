
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

/*   Finley: Mesh: NodeFile */

/*  set tags to newTag where mask>0 */

/**************************************************************/

#include "NodeFile.h"
#include "Util.h"

/**************************************************************/


void Finley_NodeFile_setTags(Finley_NodeFile* self,const int newTag, escriptDataC* mask) {
    register dim_t n;
    dim_t numNodes;
    register double *mask_array;
    register bool_t check;
    Finley_resetError();

    if (self==NULL) return;
    numNodes=self->numNodes;
    if (1!=getDataPointSize(mask)) {
       Finley_setError(TYPE_ERROR,"Finley_NodeFile_setTags: number of components of mask is 1.");
    } else if (!numSamplesEqual(mask,1,numNodes)) {
       Finley_setError(TYPE_ERROR,"Finley_NodeFile_setTags: illegal number of samples of mask Data object");
    }

    /* now we can start */

    if (Finley_noError()) {
             #pragma omp parallel for private(n,check,mask_array) schedule(static)
             for (n=0;n<numNodes;n++) {
                 mask_array=getSampleData(mask,n);
                 if (mask_array[0]>0) self->Tag[n]=newTag;
             }
    }
}
/*
* $Log$
*
*/
