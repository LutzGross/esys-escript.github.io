/* $Id$ */
/**************************************************************/
/*                                                                                                         */
/*   Finley: ElementFile                                                                                   */
/*                                                                                                         */
/*   This routine tries to reduce the number of colors used to color elements in the Finley_ElementFile in */
/*                                                                                                         */
/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "ElementFile.h"
#include "Util.h"

/**************************************************************/

void Finley_ElementFile_improveColoring(Finley_ElementFile* in,dim_t numNodes, index_t* degreeOfFreedom) {
    if (in==NULL) return;
    dim_t NN=in->ReferenceElement->Type->numNodes;
    dim_t e,i,numUncoloredElements,n,len;
    index_t *maskDOF,*old_Color,color,min_id,max_id,old_maxColor,old_minColor;
    bool_t independent;
    Finley_ErrorCode=NO_ERROR;

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
/* 
* $Log$
* Revision 1.6  2005/07/22 03:53:07  jgs
* Merge of development branch back to main trunk on 2005-07-22
*
* Revision 1.5  2005/07/08 04:07:49  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:32  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.3  2005/07/18 10:34:54  gross
* some informance improvements when reading meshes
*
* Revision 1.1.1.1.2.2  2005/06/29 02:34:49  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:13  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/
