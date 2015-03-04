
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/************************************************************************************/

/*   Dudley: Mesh: NodeFile */

/* copies the array newX into self->coordinates */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "NodeFile.h"
#include "Util.h"

/************************************************************************************/

void Dudley_NodeFile_setCoordinates(Dudley_NodeFile * self, const escript::Data* newX)
{
    char error_msg[LenErrorMsg_MAX];
    size_t numDim_size;
    int n;
    if (getDataPointSize(newX) != self->numDim)
    {
	sprintf(error_msg, "Dudley_NodeFile_setCoordinates: dimension of new coordinates has to be %d.", self->numDim);
	Dudley_setError(VALUE_ERROR, error_msg);
    }
    else if (!numSamplesEqual(newX, 1, self->numNodes))
    {
	sprintf(error_msg, "Dudley_NodeFile_setCoordinates: number of given nodes must to be %d.", self->numNodes);
	Dudley_setError(VALUE_ERROR, error_msg);
    }
    else
    {
	numDim_size = self->numDim * sizeof(double);
	Dudley_increaseStatus(self);
#pragma omp parallel private(n)
	{

#pragma omp for schedule(static)
	    for (n = 0; n < self->numNodes; n++)
	    {
		memcpy(&(self->Coordinates[INDEX2(0, n, self->numDim)]), getSampleDataROFast(newX, n), numDim_size);
	    }
	}
    }
}
