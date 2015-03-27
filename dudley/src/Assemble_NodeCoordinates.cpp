
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

/*    assemblage routines: copies node coordinates into an expanded Data Object */

/************************************************************************************/

#include "Util.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

void Dudley_Assemble_NodeCoordinates(Dudley_NodeFile * nodes, escriptDataC * x)
{
    char error_msg[LenErrorMsg_MAX];
    dim_t n;
    size_t dim_size;
    Dudley_resetError();
    if (nodes == NULL)
	return;
    if (!numSamplesEqual(x, 1, nodes->numNodes))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_NodeCoordinates: illegal number of samples of Data object");
    }
    else if (getFunctionSpaceType(x) != DUDLEY_NODES)
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_NodeCoordinates: Data object is not defined on nodes.");
    }
    else if (!isExpanded(x))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_NodeCoordinates: expanded Data object expected");
    }
    else if (!isDataPointShapeEqual(x, 1, &(nodes->numDim)))
    {
	sprintf(error_msg, "Dudley_Assemble_NodeCoordinates: Data object of shape (%d,) expected", nodes->numDim);
	Dudley_setError(TYPE_ERROR, error_msg);
    }
    else
    {
	dim_size = nodes->numDim * sizeof(double);
	requireWrite(x);
#pragma omp parallel for private(n)
	for (n = 0; n < nodes->numNodes; n++)
	    memcpy(getSampleDataRWFast(x, n), &(nodes->Coordinates[INDEX2(0, n, nodes->numDim)]), dim_size);
    }
}
