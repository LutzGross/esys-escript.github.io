/* $Id$ */
/**************************************************************/

/*   Finley: Mesh: NodeFile */

/* copies the array newX into self->coordinates */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#include "NodeFile.h"
#include "Util.h"
#include "escript/Data/DataC.h"

/**************************************************************/


void Finley_NodeFile_setCoordinates(Finley_NodeFile* self,escriptDataC* newX) {
   int n;
   if (getDataPointSize(newX)!=self->numDim)  {
      Finley_ErrorCode=VALUE_ERROR;
      sprintf(Finley_ErrorMsg,"dimension of new coordinates has to be %d.",self->numDim);
   } else if (! numSamplesEqual(newX,1,self->numNodes)) {
         Finley_ErrorCode=VALUE_ERROR;
         sprintf(Finley_ErrorMsg,"number of give nodes must to be %d.",self->numNodes);
   } else {
          #pragma omp parallel for private(n) schedule(static)
          for (n=0;n<self->numNodes;n++)
            Finley_copyDouble(self->numDim,getSampleData(newX,n),&(self->Coordinates[INDEX2(0,n,self->numDim)]));
   }
}
/*
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.5  2004/08/26 12:03:52  gross
* Some other bug in Finley_Assemble_gradient fixed.
*
* Revision 1.4  2004/07/21 05:00:54  gross
* name changes in DataC
*
* Revision 1.3  2004/07/02 04:21:13  gross
* Finley C code has been included
*
* Revision 1.2  2004/07/01 23:54:32  gross
* used DataC now
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
