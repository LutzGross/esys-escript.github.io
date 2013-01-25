
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

/*	  assemblage routines: calculates the normal vector at quadrature points on face elements */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ShapeTable.h"

/************************************************************************************/

void Dudley_Assemble_setNormal(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, escriptDataC * normal)
{
    double *local_X = NULL, *dVdv = NULL, *normal_array;
    index_t sign;
    dim_t e, q, NN, NS, numDim, numQuad, numDim_local;
    bool_t reduced_integration;
    const double *dSdv = 0;
    if (nodes == NULL || elements == NULL)
	return;

    switch (elements->numDim)
    {
    case 2:
	dSdv = &(DTDV_2D[0][0]);
	break;
    case 3:
	dSdv = &(DTDV_3D[0][0]);
	break;
    default:
	dSdv = &(DTDV_1D[0][0]);
	break;
    }
    Dudley_resetError();
    NN = elements->numNodes;
    numDim = nodes->numDim;
    reduced_integration = Dudley_Assemble_reducedIntegrationOrder(normal);
    numQuad = (!reduced_integration) ? (elements->numDim + 1) : 1;
    numDim_local = elements->numLocalDim;
    NS = elements->numDim + 1;

    /* set some parameters */

    sign = 1;
    /* check the dimensions of normal */
    if (!(numDim == numDim_local || numDim - 1 == numDim_local))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_setNormal: Cannot calculate normal vector");
    }
    else if (!isDataPointShapeEqual(normal, 1, &(numDim)))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_setNormal: illegal number of samples of normal Data object");
    }
    else if (!numSamplesEqual(normal, numQuad, elements->numElements))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_setNormal: illegal number of samples of normal Data object");
    }
    else if (!isDataPointShapeEqual(normal, 1, &(numDim)))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_setNormal: illegal point data shape of normal Data object");
    }
    else if (!isExpanded(normal))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_setNormal: expanded Data object is expected for normal.");
    }

    /* now we can start */
    if (Dudley_noError())
    {
	requireWrite(normal);
#pragma omp parallel private(local_X,dVdv)
	{
	    local_X = dVdv = NULL;
	    /* allocation of work arrays */
	    local_X = THREAD_MEMALLOC(NS * numDim, double);
	    dVdv = THREAD_MEMALLOC(numQuad * numDim * numDim_local, double);
	    if (!(Dudley_checkPtr(local_X) || Dudley_checkPtr(dVdv)))
	    {
		/* open the element loop */
#pragma omp for private(e,q,normal_array) schedule(static)
		for (e = 0; e < elements->numElements; e++)
		{
		    /* gather local coordinates of nodes into local_X: */
		    Dudley_Util_Gather_double(NS, &(elements->Nodes[INDEX2(0, e, NN)]), numDim, nodes->Coordinates,
					      local_X);
		    /*  calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q) */
		    Dudley_Util_SmallMatMult(numDim, numDim_local * numQuad, dVdv, NS, local_X, dSdv);
		    /* get normalized vector:      */
		    normal_array = getSampleDataRW(normal, e);
		    Dudley_NormalVector(numQuad, numDim, numDim_local, dVdv, normal_array);
		    for (q = 0; q < numQuad * numDim; q++)
			normal_array[q] *= sign;
		}
	    }
	    THREAD_MEMFREE(local_X);
	    THREAD_MEMFREE(dVdv);
	}
    }
}
