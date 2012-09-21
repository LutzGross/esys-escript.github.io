
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
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

void Dudley_Assemble_AverageElementData(Dudley_ElementFile * elements, escriptDataC * out, escriptDataC * in)
{
    dim_t n, q, numElements, numQuad_in, numQuad_out, i;
    __const double *in_array;
    double *out_array, vol, volinv, wq;
    register double rtmp;
    dim_t numComps = getDataPointSize(out);
    size_t numComps_size;

    Dudley_resetError();
    if (elements == NULL)
    {
	return;
    }

    numElements = elements->numElements;
    if (Dudley_Assemble_reducedIntegrationOrder(in))
    {
	numQuad_in = QuadNums[elements->numDim][0];
	wq = QuadWeight[elements->numDim][0];

    }
    else
    {
	numQuad_in = QuadNums[elements->numDim][1];
	wq = QuadWeight[elements->numDim][1];
    }
    if (Dudley_Assemble_reducedIntegrationOrder(out))
    {
	numQuad_out = QuadNums[elements->numDim][0];
    }
    else
    {
	numQuad_out = QuadNums[elements->numDim][1];

    }

    /* check out and in */
    if (numComps != getDataPointSize(in))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_AverageElementData: number of components of input and output Data do not match.");
    }
    else if (!numSamplesEqual(in, numQuad_in, numElements))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_AverageElementData: illegal number of samples of input Data object");
    }
    else if (!numSamplesEqual(out, numQuad_out, numElements))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_AverageElementData: illegal number of samples of output Data object");
    }
    else if (!isExpanded(out))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_AverageElementData: expanded Data object is expected for output data.");
    }

    /* now we can start */

    if (Dudley_noError())
    {
	if (isExpanded(in))
	{
	    vol = 0;
	    for (q = 0; q < numQuad_in; ++q)
		vol += wq;
	    volinv = 1. / vol;
	    requireWrite(out);
#pragma omp parallel private(n, i, rtmp, q, in_array, out_array)
	    {
# pragma omp for schedule(static)
		for (n = 0; n < numElements; n++)
		{
		    in_array = getSampleDataRO(in, n);
		    out_array = getSampleDataRW(out, n);
		    for (i = 0; i < numComps; ++i)
		    {
			rtmp = 0;
			for (q = 0; q < numQuad_in; ++q)
			    rtmp += in_array[INDEX2(i, q, numComps)] * wq;
			rtmp *= volinv;
			for (q = 0; q < numQuad_out; ++q)
			    out_array[INDEX2(i, q, numComps)] = rtmp;
		    }
		}
	    }
	}
	else
	{
	    numComps_size = numComps * sizeof(double);
	    requireWrite(out);
#pragma omp parallel private(q,n,out_array,in_array)
	    {
# pragma omp for schedule(static)
		for (n = 0; n < numElements; n++)
		{
		    in_array = getSampleDataRO(in, n);
		    out_array = getSampleDataRW(out, n);
		    for (q = 0; q < numQuad_out; q++)
			memcpy(out_array + q * numComps, in_array, numComps_size);
		}
	    }
	}
    }
    return;
}
