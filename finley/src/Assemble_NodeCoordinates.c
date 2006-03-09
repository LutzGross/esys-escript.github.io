/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
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
       #pragma omp parallel for private(n)
       for (n=0;n<nodes->numNodes;n++) 
          Finley_copyDouble(nodes->numDim,&(nodes->Coordinates[INDEX2(0,n,nodes->numDim)]),getSampleData(x,n));
  }
}
/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
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
