/* $Id$ */

/**************************************************************/

/*    assemblage routines: copies node coordinates into an expanded Data Object */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Util.h"
#include "Finley.h"
#include "Assemble.h"
#include "NodeFile.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

void Finley_Assemble_NodeCoordinates(Finley_NodeFile* nodes,escriptDataC* x) {
  dim_t n;
  if (nodes==NULL) return;
  if (! numSamplesEqual(x,1,nodes->numNodes)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal number of samples of Data object");
  } else if (getFunctionSpaceType(x)!=FINLEY_NODES) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"Data object is not defined on nodes.");
  } else if (! isExpanded(x)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"expanded Data object expected");
  } else if (! isDataPointShapeEqual(x,1, &(nodes->numDim))) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"Data object of shape (%d,) expected",nodes->numDim);
  } else {
       #pragma omp parallel for private(n)
       for (n=0;n<nodes->numNodes;n++) 
          Finley_copyDouble(nodes->numDim,&(nodes->Coordinates[INDEX2(0,n,nodes->numDim)]),getSampleData(x,n));
  }
}
/*
 * $Log$
 * Revision 1.2  2005/07/08 04:07:45  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:46  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.5  2004/08/26 12:03:52  gross
 * Some other bug in Finley_Assemble_gradient fixed.
 *
 * Revision 1.4  2004/08/05 03:58:27  gross
 * Bug in Assemble_NodeCoordinates fixed
 *
 * Revision 1.3  2004/07/30 04:37:06  gross
 * escript and finley are linking now and RecMeshTest.py has been passed
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
