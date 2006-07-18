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

/*   gathers the ElementFile out from the  ElementFile in using index[0:out->numElements-1].  */
/*   index has to be between 0 and in->numElements-1. */
/*   a conservative assumtion on the coloring is made */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/

void Finley_ElementFile_gather(index_t* index, Finley_ElementFile* in, Finley_ElementFile* out) {
   index_t k;
   dim_t e,j;
   dim_t NN_in=in->ReferenceElement->Type->numNodes;
   dim_t NN_out=out->ReferenceElement->Type->numNodes;
   if (in!=NULL) {
     /*OMP */
     #pragma omp parallel for private(e,k,j) schedule(static)
     for (e=0;e<out->numElements;e++) {
        k=index[e];
#ifdef PASO_MPI
        out->Dom[e]=in->Dom[k];
#endif
        out->Id[e]=in->Id[k];
        out->Tag[e]=in->Tag[k];
        out->Color[e]=in->Color[k]+out->maxColor+1;
        for(j=0;j<MIN(NN_out,NN_in);j++) out->Nodes[INDEX2(j,e,NN_out)]=in->Nodes[INDEX2(j,k,NN_in)];
     }
     out->minColor=MIN(out->minColor,in->minColor+out->maxColor+1);
     out->maxColor=MAX(out->maxColor,in->maxColor+out->maxColor+1);
   }
}
/* 
* $Log$
* Revision 1.3  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.2.2.1  2005/09/07 06:26:18  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.2  2005/07/08 04:07:49  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.1.1.1.2.2  2005/06/30 01:53:55  gross
* a bug in coloring fixed
*
* Revision 1.1.1.1.2.1  2005/06/29 02:34:49  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
