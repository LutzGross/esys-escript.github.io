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

/*   Finley: ElementFile                                                      */

/* copies element file in into element file out starting from offset          */
/* the elements offset to in->numElements+offset-1 in out will be overwritten */
                                                                                                                                                   
/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"

/****************************************************************************/


void Finley_ElementFile_copyTable(index_t offset,Finley_ElementFile* out,index_t node_offset, index_t idOffset,Finley_ElementFile* in) {
    dim_t i,n;
    dim_t NN=out->numNodes;
    if (in==NULL) return;
    if (out->ReferenceElement->Type->TypeId!=in->ReferenceElement->Type->TypeId) {
        Finley_setError(TYPE_ERROR,"Finley_ElementFile_copyTable: dimensions of element files don't match");
    }
    if (out->MPIInfo->comm!=in->MPIInfo->comm) {
        Finley_setError(TYPE_ERROR,"Finley_ElementFile_copyTable: MPI communicators of element files don't match");
    }
    if (Finley_noError()) {
       #pragma omp parallel for private(i,n) schedule(static)
       for(n=0;n<in->numElements;n++) {
          out->Owner[offset+n]=out->Owner[n];
          out->Id[offset+n]=in->Id[n]+idOffset;
          out->Tag[offset+n]=in->Tag[n];
          for(i=0;i<NN;i++) out->Nodes[INDEX2(i,offset+n,NN)]=in->Nodes[INDEX2(i,n,in->ReferenceElement->Type->numNodes)]+node_offset;
       }
    }
}
