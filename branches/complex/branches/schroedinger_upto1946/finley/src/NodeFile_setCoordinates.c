
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*   Finley: Mesh: NodeFile */

/* copies the array newX into self->coordinates */

/**************************************************************/

#include "NodeFile.h"
#include "Util.h"

/**************************************************************/


void Finley_NodeFile_setCoordinates(Finley_NodeFile* self,escriptDataC* newX) {
  char error_msg[LenErrorMsg_MAX];
  size_t numDim_size;
  int n;
   if (getDataPointSize(newX)!=self->numDim)  {
      sprintf(error_msg,"Finley_NodeFile_setCoordinates: dimension of new coordinates has to be %d.",self->numDim);
      Finley_setError(VALUE_ERROR,error_msg);
   } else if (! numSamplesEqual(newX,1,self->numNodes)) {
         sprintf(error_msg,"Finley_NodeFile_setCoordinates: number of give nodes must to be %d.",self->numNodes);
         Finley_setError(VALUE_ERROR,error_msg);
   } else {
          numDim_size=self->numDim*sizeof(double);
          Finley_increaseStatus(self);
          #pragma omp parallel for private(n) schedule(static)
          for (n=0;n<self->numNodes;n++) {
            memcpy(&(self->Coordinates[INDEX2(0,n,self->numDim)]), getSampleDataFast(newX,n), numDim_size);
          }
   }
}
