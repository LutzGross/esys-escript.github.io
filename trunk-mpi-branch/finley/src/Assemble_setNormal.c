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

/*    assemblage routines: calculates the normal vector at quadrature points on face elements */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/


void Finley_Assemble_setNormal(Finley_NodeFile* nodes, Finley_ElementFile* elements, escriptDataC* normal) {
  double *local_X=NULL, *dVdv=NULL,*normal_array;
  index_t sign,node_offset;
  dim_t e,q;
  if (nodes==NULL || elements==NULL) return;
  dim_t NN=elements->ReferenceElement->Type->numNodes;
  dim_t NS=elements->ReferenceElement->Type->numShapes;
  dim_t numDim=nodes->numDim;
  dim_t numQuad=elements->ReferenceElement->numQuadNodes;
  dim_t numDim_local=elements->ReferenceElement->Type->numDim;
  Finley_resetError();

  /* set some parameter */

  if (getFunctionSpaceType(normal)==FINLEY_CONTACT_ELEMENTS_2) {
      node_offset=NN-NS;
      sign=-1;
  } else {
      node_offset=0;
      sign=1;
  }
  /* check the dimensions of normal */
  if (! (numDim==numDim_local || numDim-1==numDim_local)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: Cannot calculate normal vector");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: illegal number of samples of normal Data object");
  } else if (! numSamplesEqual(normal,numQuad,elements->numElements)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: illegal number of samples of normal Data object");
  } else if (! isDataPointShapeEqual(normal,1,&(numDim))) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: illegal point data shape of normal Data object");
  }  else if (!isExpanded(normal)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_setNormal: expanded Data object is expected for normal.");
  }
   
  /* now we can start */
    if (Finley_noError()) {
          #pragma omp parallel private(local_X,dVdv)
          {
             local_X=dVdv=NULL;
             /* allocation of work arrays */
             local_X=THREAD_MEMALLOC(NS*numDim,double); 
             dVdv=THREAD_MEMALLOC(numQuad*numDim*numDim_local,double); 
             if (!(Finley_checkPtr(local_X) || Finley_checkPtr(dVdv) ) ) {
                       /* open the element loop */
                       #pragma omp for private(e,q,normal_array) schedule(static)
                       for(e=0;e<elements->numElements;e++) {
                          /* gather local coordinates of nodes into local_X: */
                          Finley_Util_Gather_double(NS,&(elements->Nodes[INDEX2(node_offset,e,NN)]),numDim,nodes->Coordinates,local_X);
                          /*  calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q) */
                          Finley_Util_SmallMatMult(numDim,numDim_local*numQuad,dVdv,NS,local_X,elements->ReferenceElement->dSdv);
                          /* get normalized vector:  */
                          normal_array=getSampleData(normal,e);
                          Finley_NormalVector(numQuad,numDim,numDim_local,dVdv,normal_array);
			  for (q=0;q<numQuad*numDim;q++) normal_array[q]*=sign;
                       }
                 }
                 THREAD_MEMFREE(local_X);
                 THREAD_MEMFREE(dVdv);
             }
    }
}
/*
 * $Log$
 * Revision 1.6  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.5.2.1  2005/09/07 06:26:18  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.5  2005/07/08 04:07:48  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2004/12/15 07:08:32  jgs
 * *** empty log message ***
 * Revision 1.1.1.1.2.3  2005/06/29 02:34:48  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.2  2004/11/24 01:37:13  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 *
 *
 */
