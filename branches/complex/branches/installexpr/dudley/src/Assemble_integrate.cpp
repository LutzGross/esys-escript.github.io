
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

/*    assemblage routines: integrates data on quadrature points   */

/************************************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/************************************************************************************/

void Dudley_Assemble_integrate(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, escriptDataC * data, double *out)
{
/*    type_t data_type=getFunctionSpaceType(data);*/
    dim_t numQuadTotal;
    dim_t numComps = getDataPointSize(data);
    Dudley_ElementFile_Jacobeans *jac = NULL;
    Esys_MPI_rank my_mpi_rank;

    Dudley_resetError();
    if (nodes == NULL || elements == NULL)
	return;
    my_mpi_rank = nodes->MPIInfo->rank;
    /* set some parameter */
    jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, Dudley_Assemble_reducedIntegrationOrder(data));
    if (Dudley_noError())
    {
	numQuadTotal = jac->numQuad;
	/* check the shape of the data  */
	if (!numSamplesEqual(data, numQuadTotal, elements->numElements))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_integrate: illegal number of samples of integrant kernel Data object");
	}
	/* now we can start */

	if (Dudley_noError())
	{
	    dim_t q, e, i;
	    __const double *data_array = NULL;
	    double *out_local = NULL, rtmp;
	    for (q = 0; q < numComps; q++)
		out[q] = 0;
#pragma omp parallel private(q,i,rtmp,data_array,out_local)
	    {
		out_local = new double[numComps];
		if (!Dudley_checkPtr(out_local))
		{
		    /* initialize local result */

		    for (i = 0; i < numComps; i++)
			out_local[i] = 0;

		    /* open the element loop */

		    if (isExpanded(data))
		    {
#pragma omp for private(e) schedule(static)
			for (e = 0; e < elements->numElements; e++)
			{
			    if (elements->Owner[e] == my_mpi_rank)
			    {
				double vol = jac->absD[e] * jac->quadweight;
				data_array = getSampleDataRO(data, e);
				for (q = 0; q < numQuadTotal; q++)
				{
				    for (i = 0; i < numComps; i++)
					out_local[i] += data_array[INDEX2(i, q, numComps)] * vol;
				}
			    }
			}
		    }
		    else
		    {
#pragma omp for private(e) schedule(static)
			for (e = 0; e < elements->numElements; e++)
			{
			    if (elements->Owner[e] == my_mpi_rank)
			    {
				double vol = jac->absD[e] * jac->quadweight;
				data_array = getSampleDataRO(data, e);
				rtmp = 0.;
				for (q = 0; q < numQuadTotal; q++)
				    rtmp += vol;
				for (i = 0; i < numComps; i++)
				    out_local[i] += data_array[i] * rtmp;
			    }
			}
		    }
		    /* add local results to global result */
#pragma omp critical
		    for (i = 0; i < numComps; i++)
			out[i] += out_local[i];
		}
		delete[] out_local;
	    }
	}
    }
}
