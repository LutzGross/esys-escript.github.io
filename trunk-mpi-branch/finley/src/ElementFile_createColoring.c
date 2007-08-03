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
/*   This routine tries to reduce the number of colors used to color elements in the Finley_ElementFile in */
/*                                                                                                         */
/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id: ElementFile_createColoring.c 1140 2007-05-15 03:23:17Z ksteube $ */

/**************************************************************/

#include "ElementFile.h"
#include "Util.h"

/**************************************************************/

void Finley_ElementFile_createColoring(Finley_ElementFile* in,dim_t numNodes, index_t* degreeOfFreedom) {
    dim_t e,i,numUncoloredElements,n,len,NN;
    index_t *maskDOF,color,min_id,max_id;
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
                     for (i=0;i<NN;i++) if (maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]>0) independent=FALSE;
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
