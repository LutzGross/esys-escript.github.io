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
    index_t *maskDOF,*old_Color,color,min_id,max_id,old_maxColor,old_minColor;
    bool_t independent;
    Finley_resetError();

    if (in==NULL) return;
    NN=in->ReferenceElement->Type->numNodes;
    if (in->numElements<1) return;

    min_id=Finley_Util_getMinInt(1,numNodes,degreeOfFreedom);
    max_id=Finley_Util_getMaxInt(1,numNodes,degreeOfFreedom);
    len=max_id-min_id+1;
    maskDOF=TMPMEMALLOC(len,index_t);
    old_Color=TMPMEMALLOC(in->numElements,index_t);
    
    if (! (Finley_checkPtr(maskDOF) || Finley_checkPtr(old_Color) ) ) {
         #pragma omp parallel for private(e) schedule(static)
         for (e=0;e<in->numElements;e++) {
               old_Color[e]=in->Color[e];
               in->Color[e]=-1;
         }
         old_maxColor=in->maxColor;
         old_minColor=in->minColor;
         in->maxColor=-1;
         in->minColor=0;
         numUncoloredElements=in->numElements;
         while (numUncoloredElements>0) {
            /* initialize the mask marking nodes used by a color */
            #pragma omp parallel for private(n) schedule(static)
            for (n=0;n<len;n++) maskDOF[n]=-1;
            /* the existing coloring is used to make sure that the new coloring can be done in parallel */
            numUncoloredElements=0;
            if ((old_maxColor-old_minColor+1)*NN<in->numElements) {
               for (color=old_minColor;color<=old_maxColor;color++) {
                  #pragma omp parallel for private(i,e,independent) schedule(static) reduction(+:numUncoloredElements)
                  for (e=0;e<in->numElements;e++) {
                     if (old_Color[e]==color) {
                        /* find out if element e is independend from the elements already colored: */
                        independent=TRUE;
                        for (i=0;i<NN;i++) if (maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]>0) independent=FALSE;
                        /* if e is independend a new color is assigned and the nodes are marked as being used */
                        if (independent) {
                            for (i=0;i<NN;i++) maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]=1;
                            old_Color[e]=-1;
                            in->Color[e]=in->maxColor+1;
                         } else {
                            numUncoloredElements++;
                         }
                     }
                  }
               } /* end of color loop */
             } else {
               for (e=0;e<in->numElements;e++) {
                  if (old_Color[e]!=-1) {
                     /* find out if element e is independend from the elements already colored: */
                     independent=TRUE;
                     for (i=0;i<NN;i++) if (maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]>0) independent=FALSE;
                     /* if e is independend a new color is assigned and the nodes are marked as being used */
                     if (independent) {
                         for (i=0;i<NN;i++) maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]=1;
                         old_Color[e]=-1;
                         in->Color[e]=in->maxColor+1;
                      } else {
                         numUncoloredElements++;
                      }
                  }
               }
            }
            in->maxColor++;
         }  /* end of while loop */
    }

    /* all done : */

    TMPMEMFREE(maskDOF);
    TMPMEMFREE(old_Color);
}
