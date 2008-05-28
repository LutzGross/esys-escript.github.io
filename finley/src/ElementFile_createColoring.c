
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
/*                                                                                                         */
/*   Finley: ElementFile                                                                                   */
/*                                                                                                         */
/*   This routine tries to reduce the number of colors used to color elements in the Finley_ElementFile in */
/*                                                                                                         */
/**************************************************************/

#include "ElementFile.h"
#include "Util.h"

/**************************************************************/

void Finley_ElementFile_createColoring(Finley_ElementFile* in,dim_t numNodes, index_t* degreeOfFreedom) {
    dim_t e,i,numUncoloredElements,n,len,NN;
    index_t *maskDOF,min_id,max_id;
    bool_t independent;

    if (in==NULL) return;
    if (in->numElements<1) return;
    NN=in->numNodes;

    min_id=Finley_Util_getMinInt(1,numNodes,degreeOfFreedom);
    max_id=Finley_Util_getMaxInt(1,numNodes,degreeOfFreedom);
    len=max_id-min_id+1;
    maskDOF=TMPMEMALLOC(len,index_t);
    if (! Finley_checkPtr(maskDOF) ) {
         #pragma omp parallel for private(e) schedule(static)
         for (e=0;e<in->numElements;e++) in->Color[e]=-1;
         numUncoloredElements=in->numElements;
         in->minColor=0;
         in->maxColor=in->minColor-1;
         while (numUncoloredElements>0) {
            /* initialize the mask marking nodes used by a color */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<len;n++) maskDOF[n]=-1;
            numUncoloredElements=0;
            /* OMP ?*/
            for (e=0;e<in->numElements;e++) {
                  if (in->Color[e]<0) {
                     /* find out if element e is independend from the elements already colored: */
                     independent=TRUE;
                     for (i=0;i<NN;i++) {
#ifdef BOUNDS_CHECK
if (in->Nodes[INDEX2(i,e,NN)] < 0 || in->Nodes[INDEX2(i,e,NN)] >= numNodes) { printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d min_id=%d in->Nodes[INDEX2...]=%d\n", __FILE__, __LINE__, i, e, NN, min_id, in->Nodes[INDEX2(i,e,NN)]); exit(1); }
if ((degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id) >= len || (degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id) < 0) { printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d min_id=%d dof=%d\n", __FILE__, __LINE__, i, e, NN, min_id, degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id); exit(1); }
#endif
		        if (maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]>0) {
                             independent=FALSE;
                             break;
                        }
		     }
                     /* if e is independend a new color is assigned and the nodes are marked as being used */
                     if (independent) {
                           for (i=0;i<NN;i++) maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]=1;
                           in->Color[e]=in->maxColor+1;
                     } else {
                             numUncoloredElements++;
                     }
                  }

            }
            in->maxColor++;
         }  /* end of while loop */
    }
    /* all done : */
    TMPMEMFREE(maskDOF);
}
