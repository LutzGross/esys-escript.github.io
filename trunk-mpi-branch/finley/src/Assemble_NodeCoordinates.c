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

/*    assemblage routines: copies node coordinates into an expanded Data Object */

/**************************************************************/

/*  Copyrights by ACcESS Australia 2003,2004,2005 */
/*  Version: $Id$ */

/**************************************************************/

#include "Util.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

void Finley_Assemble_NodeCoordinates(Finley_NodeFile* nodes,escriptDataC* x) {
  char error_msg[LenErrorMsg_MAX];
  dim_t n;
  size_t dim_size;
  Finley_resetError();
  if (nodes==NULL) return;
  if (! numSamplesEqual(x,1,nodes->numNodes)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_NodeCoordinates: illegal number of samples of Data object");
  } else if (getFunctionSpaceType(x)!=FINLEY_NODES) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_NodeCoordinates: Data object is not defined on nodes.");
  } else if (! isExpanded(x)) {
       Finley_setError(TYPE_ERROR,"Finley_Assemble_NodeCoordinates: expanded Data object expected");
  } else if (! isDataPointShapeEqual(x,1, &(nodes->numDim))) {
       sprintf(error_msg,"Finley_Assemble_NodeCoordinates: Data object of shape (%d,) expected",nodes->numDim);
       Finley_setError(TYPE_ERROR,error_msg);
  } else {
       dim_size=nodes->numDim*sizeof(double);
       #pragma omp parallel for private(n)
       for (n=0;n<nodes->numNodes;n++) 
          memcpy(getSampleDataFast(x,n),&(nodes->Coordinates[INDEX2(0,n,nodes->numDim)]),dim_size);
  }
}
