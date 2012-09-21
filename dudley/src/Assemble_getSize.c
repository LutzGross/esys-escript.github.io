
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

/*    assemblage routines: */

/*    calculates the minimum distance between two vertices of elements and assigns the value to each  */
/*    quadrature point in element_size                                                                         */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/
void Dudley_Assemble_getSize(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, escriptDataC * element_size)
{

    double *local_X = NULL, *element_size_array;
    dim_t e, n0, n1, q, i, NVertices, NN, NS, numQuad, numDim;
    double d, diff, max_diff;
    Dudley_resetError();

    if (nodes == NULL || elements == NULL)
    {
	return;
    }

    numDim = nodes->numDim;

    /* now we look up what type of elements we need based on the function space of element_size */
    /* if it is DUDLEY_REDUCED_ELEMENTS or DUDLEY_REDUCED_FACE_ELEMENTS then we have single quad point */

    if (Dudley_Assemble_reducedIntegrationOrder(element_size))
    {
	numQuad = 1;
    }
    else
    {
	numQuad = elements->numDim + 1;
    }

    NN = elements->numNodes;
    NS = elements->numDim + 1;
    NVertices = elements->numDim + 1;

    /* check the dimensions of element_size */

    if (!numSamplesEqual(element_size, numQuad, elements->numElements))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getSize: illegal number of samples of element size Data object");
    }
    else if (!isDataPointShapeEqual(element_size, 0, &(numDim)))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getSize: illegal data point shape of element size Data object");
    }
    else if (!isExpanded(element_size))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getSize: expanded Data object is expected for element size.");
    }
    /* now we can start: */

    if (Dudley_noError())
    {
	requireWrite(element_size);
#pragma omp parallel private(local_X)
	{
	    /* allocation of work arrays */
	    local_X = THREAD_MEMALLOC(NN * numDim, double);
	    if (!Dudley_checkPtr(local_X))
	    {
		/* open the element loop */
#pragma omp for private(e,max_diff,diff,n0,n1,d,q,i,element_size_array) schedule(static)
		for (e = 0; e < elements->numElements; e++)
		{
		    /* gather local coordinates of nodes into local_X(numDim,NN): */
		    Dudley_Util_Gather_double(NS, &(elements->Nodes[INDEX2(0, e, NN)]), numDim, nodes->Coordinates,
					      local_X);
		    /* calculate minimal differences */
		    max_diff = 0;
		    for (n0 = 0; n0 < NVertices; n0++)
		    {
			for (n1 = n0 + 1; n1 < NVertices; n1++)
			{
			    diff = 0;
			    for (i = 0; i < numDim; i++)
			    {
				d = local_X[INDEX2(i, n0, numDim)] - local_X[INDEX2(i, n1, numDim)];
				diff += d * d;
			    }

    			    max_diff = MAX(max_diff, diff);
			    
			}
		    }
		    max_diff = sqrt(max_diff);
		    /* set all values to max_diff */
		    element_size_array = getSampleDataRW(element_size, e);
		    for (q = 0; q < numQuad; q++)
			element_size_array[q] = max_diff;
		}
	    }
	    THREAD_MEMFREE(local_X);
	}			/* end of parallel region */
    }
    return;
}
