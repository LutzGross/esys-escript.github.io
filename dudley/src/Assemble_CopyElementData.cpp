
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

/*    assemblage routines: copies data between elements       */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/****************************************************************************************************************************/
#include "ShapeTable.h"

void Dudley_Assemble_CopyElementData(Dudley_ElementFile * elements, escript::Data * out, const escript::Data * in)
{
    dim_t n, q, numElements, numQuad;
    __const double *in_array;
    double *out_array;
    dim_t numComps = getDataPointSize(out);
    size_t len_size;

    Dudley_resetError();
    if (elements == NULL)
    {
	return;
    }

    numElements = elements->numElements;
    if (Dudley_Assemble_reducedIntegrationOrder(in))
    {
	numQuad = QuadNums[elements->numDim][0];
    }
    else
    {
	numQuad = QuadNums[elements->numDim][1];

    }

    /* check out and in */
    if (numComps != getDataPointSize(in))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_CopyElementData: number of components of input and output Data do not match.");
    }
    else if (!numSamplesEqual(in, numQuad, numElements))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_CopyElementData: illegal number of samples of input Data object");
    }
    else if (!numSamplesEqual(out, numQuad, numElements))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_CopyElementData: illegal number of samples of output Data object");
    }
    else if (!isExpanded(out))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_CopyElementData: expanded Data object is expected for output data.");
    }

    /* now we can start */

    if (Dudley_noError())
    {
	if (isExpanded(in))
	{
	    len_size = numComps * numQuad * sizeof(double);
	    requireWrite(out);
#pragma omp parallel private(n)
	    {
# pragma omp for schedule(static)
		for (n = 0; n < numElements; n++)
		    memcpy(getSampleDataRW(out, n), getSampleDataRO(in, n), len_size);
	    }
	}
	else
	{
	    len_size = numComps * sizeof(double);
	    requireWrite(out);
#pragma omp parallel private(q,n,out_array,in_array)
	    {
# pragma omp for schedule(static)
		for (n = 0; n < numElements; n++)
		{
		    in_array = getSampleDataRO(in, n);
		    out_array = getSampleDataRW(out, n);
		    for (q = 0; q < numQuad; q++)
			memcpy(out_array + q * numComps, in_array, len_size);
		}
	    }
	}
    }
    return;
}
