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

void Finley_ElementFile_improveColoring(Finley_ElementFile* in,maybelong numNodes, maybelong* degreeOfFreedom) {
    if (in==NULL) return;
    maybelong NN=in->ReferenceElement->Type->numNodes;
    maybelong old_numColors,*maskDOF,*old_Color,e,i,numUncoloredElements,n,color,independent;
    Finley_ErrorCode=NO_ERROR;
    int min_id,max_id;
    maybelong len;

    if (in->numElements<1) return;

    min_id=Finley_Util_getMinInt(1,numNodes,degreeOfFreedom);
    max_id=Finley_Util_getMaxInt(1,numNodes,degreeOfFreedom);
    len=max_id-min_id+1;
    maskDOF=TMPMEMALLOC(len,maybelong);
    old_Color=TMPMEMALLOC(in->numElements,maybelong);
    
    if (! (Finley_checkPtr(maskDOF) || Finley_checkPtr(old_Color) ) ) {
         #pragma omp parallel for private(e) schedule(static)
         for (e=0;e<in->numElements;e++) {
               old_Color[e]=in->Color[e];
               in->Color[e]=-1;
         }
         old_numColors=in->numColors;
         in->numColors=0;
         numUncoloredElements=in->numElements;
         while (numUncoloredElements>0) {
            #pragma omp parallel private(color)
            {
               /* initialize the mask marking nodes used by a color */
               #pragma omp for private(n) schedule(static)
               for (n=0;n<len;n++) maskDOF[n]=-1;
               /* the existing coloring is used to make sure that the new coloring can be done in parallel */
               #pragma omp master
               numUncoloredElements=0;
               for (color=0;color<old_numColors;color++) {
                  #pragma omp for private(i,e,independent) schedule(static) reduction(+:numUncoloredElements)
                  for (e=0;e<in->numElements;e++) {
                     if (old_Color[e]==color) {
                        /* find out if element e is independend from the elements already colored: */
                        independent=TRUE;
                        for (i=0;i<NN;i++) if (maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]>0) independent=FALSE;
                        /* if e is independend a new color is assigned and the nodes are marked as being used */
                        if (independent) {
                            for (i=0;i<NN;i++) maskDOF[degreeOfFreedom[in->Nodes[INDEX2(i,e,NN)]]-min_id]=1;
                            old_Color[e]=-1;
                            in->Color[e]=in->numColors;
                         } else {
                            numUncoloredElements++;
                         }
                     }
                  }
               } /* end of color loop */
            }
            in->numColors++;
         }  /* end of while loop */
    }

    /* all done : */

    TMPMEMFREE(maskDOF);
    TMPMEMFREE(old_Color);
}
/* 
* $Log$
* Revision 1.2  2004/12/14 05:39:30  jgs
* *** empty log message ***
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:13  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.2  2004/07/02 04:21:13  gross
* Finley C code has been included
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
