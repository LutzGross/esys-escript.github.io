
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

/*    assemblage routines: copies data between different types nodal representation   */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Util.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Dudley_Assemble_CopyNodalData(Dudley_NodeFile * nodes, escript::Data * out, const escript::Data * in)
{
    dim_t n, k, l, mpiSize;
    dim_t numComps = getDataPointSize(out);
    paso::Coupler_ptr coupler;
    type_t in_data_type = getFunctionSpaceType(in);
    type_t out_data_type = getFunctionSpaceType(out);
    index_t upperBound;
    double *recv_buffer;
    size_t numComps_size = 0;
    Dudley_resetError();
    if (nodes == NULL)
	return;
    mpiSize = nodes->MPIInfo->size;

    /* check out and in */
    if (numComps != getDataPointSize(in))
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_CopyNodalData: number of components of input and output Data do not match.");
    }
    else if (!isExpanded(out))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_CopyNodalData: expanded Data object is expected for output data.");
    }

    /* more sophisticated test needed for overlapping node/DOF counts */
    if (in_data_type == DUDLEY_NODES)
    {
	if (!numSamplesEqual(in, 1, Dudley_NodeFile_getNumNodes(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of input Data object");
	}
    }
    else if (in_data_type == DUDLEY_REDUCED_NODES)
    {
	if (!numSamplesEqual(in, 1, Dudley_NodeFile_getNumReducedNodes(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of input Data object");
	}
    }
    else if (in_data_type == DUDLEY_DEGREES_OF_FREEDOM)
    {
	if (!numSamplesEqual(in, 1, Dudley_NodeFile_getNumDegreesOfFreedom(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of input Data object");
	}
	if ((((out_data_type == DUDLEY_NODES) || (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)) && !isExpanded(in)
	     && (mpiSize > 1)))
	{

	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: DUDLEY_DEGREES_OF_FREEDOM to DUDLEY_NODES or DUDLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
	}
    }
    else if (in_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
    {
	if (!numSamplesEqual(in, 1, Dudley_NodeFile_getNumReducedDegreesOfFreedom(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of input Data object");
	}
	if ((out_data_type == DUDLEY_DEGREES_OF_FREEDOM) && !isExpanded(in) && (mpiSize > 1))
	{

	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: DUDLEY_REDUCED_DEGREES_OF_FREEDOM to DUDLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
	}
    }
    else
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_CopyNodalData: illegal function space type for target object");
    }

    if (out_data_type == DUDLEY_NODES)
    {
	if (!numSamplesEqual(out, 1, Dudley_NodeFile_getNumNodes(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of output Data object");
	}
    }
    else if (out_data_type == DUDLEY_REDUCED_NODES)
    {
	if (!numSamplesEqual(out, 1, Dudley_NodeFile_getNumReducedNodes(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of output Data object");
	}
    }
    else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)
    {
	if (!numSamplesEqual(out, 1, Dudley_NodeFile_getNumDegreesOfFreedom(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of output Data object");
	}
    }
    else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
    {
	if (!numSamplesEqual(out, 1, Dudley_NodeFile_getNumReducedDegreesOfFreedom(nodes)))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_CopyNodalData: illegal number of samples of output Data object");
	}
    }
    else
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_CopyNodalData: illegal function space type for source object");
    }

    /* now we can start */

    if (Dudley_noError())
    {
	/*********************** DUDLEY_NODES **************************************************/
	numComps_size = (size_t) numComps *sizeof(double);
	if (in_data_type == DUDLEY_NODES)
	{
	    requireWrite(out);
	    if (out_data_type == DUDLEY_NODES)
	    {
#pragma omp parallel private(n)
		{

#pragma omp parallel for private(n) schedule(static)
		    for (n = 0; n < nodes->nodesMapping->numNodes; n++)
		    {
			memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, n), numComps_size);
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_REDUCED_NODES)
	    {
#pragma omp parallel private(n)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nodes->reducedNodesMapping->numTargets; n++)
		    {
			memcpy(getSampleDataRWFast(out, n),
			       getSampleDataROFast(in, nodes->reducedNodesMapping->map[n]), numComps_size);
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)
	    {
		int nComps = nodes->degreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel private(n)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nComps; n++)
		    {
			memcpy(getSampleDataRWFast(out, n),
			       getSampleDataROFast(in, nodes->degreesOfFreedomMapping->map[n]), numComps_size);
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
	    {
		int nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel private(n)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nComps; n++)
		    {
			memcpy(getSampleDataRWFast(out, n),
			       getSampleDataROFast(in, nodes->reducedDegreesOfFreedomMapping->map[n]), numComps_size);
		    }
		}
	    }
	/*********************** DUDLEY_REDUCED_NODES **************************************************/
	}
	else if (in_data_type == DUDLEY_REDUCED_NODES)
	{
	    requireWrite(out);
	    if (out_data_type == DUDLEY_NODES)
	    {
		Dudley_setError(TYPE_ERROR, "Dudley_Assemble_CopyNodalData: cannot copy from reduced nodes to nodes.");

	    }
	    else if (out_data_type == DUDLEY_REDUCED_NODES)
	    {
#pragma omp parallel private(n)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nodes->reducedNodesMapping->numNodes; n++)
		    {
			memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, n), numComps_size);
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)
	    {
		Dudley_setError(TYPE_ERROR,
				"Dudley_Assemble_CopyNodalData: cannot copy from reduced nodes to degrees of freedom.");
	    }
	    else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
	    {
		int nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel private(n,k)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nComps; n++)
		    {
			k = nodes->reducedDegreesOfFreedomMapping->map[n];
			memcpy(getSampleDataRWFast(out, n),
			       getSampleDataROFast(in, nodes->reducedNodesMapping->target[k]), numComps_size);
		    }
		}
	    }

	/*********************** DUDLEY_DEGREES_OF_FREEDOM **************************************************/
	}
	else if (in_data_type == DUDLEY_DEGREES_OF_FREEDOM)
	{
	    requireWrite(out);
	    if (out_data_type == DUDLEY_NODES)
	    {
		coupler.reset(new paso::Coupler(nodes->degreesOfFreedomConnector, numComps));
		if (Esys_noError())
		{
		    /* safe provided coupler->copyAll is called before the pointer in "in" is invalidated */
		    const_cast<escript::Data*>(in)->resolve();
		    coupler->startCollect(in->getDataRO());  
		    recv_buffer = coupler->finishCollect();
		    upperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
#pragma omp parallel private(n,k)
		    {
#pragma omp for schedule(static)
			for (n = 0; n < nodes->numNodes; n++)
			{
			    k = nodes->degreesOfFreedomMapping->target[n];
			    if (k < upperBound)
			    {
				memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, k), numComps_size);
			    }
			    else
			    {
				memcpy(getSampleDataRWFast(out, n),
				       &recv_buffer[(k - upperBound) * numComps], numComps_size);
			    }
			}
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_REDUCED_NODES)
	    {
		coupler.reset(new paso::Coupler(nodes->degreesOfFreedomConnector, numComps));
		if (Esys_noError())
		{
		    /* safe provided coupler->copyAll is called before the pointer in "in" is invalidated */
		    const_cast<escript::Data*>(in)->resolve();
		    coupler->startCollect(in->getDataRO());  
		    recv_buffer = coupler->finishCollect();
		    upperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
		    requireWrite(out);

#pragma omp parallel private(n,k,l)
		    {
#pragma omp for schedule(static)
			for (n = 0; n < nodes->reducedNodesMapping->numTargets; n++)
			{
			    l = nodes->reducedNodesMapping->map[n];
			    k = nodes->degreesOfFreedomMapping->target[l];
			    if (k < upperBound)
			    {
				memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, k), numComps_size);
			    }
			    else
			    {
				memcpy(getSampleDataRWFast(out, n),
				       &recv_buffer[(k - upperBound) * numComps], numComps_size);
			    }
			}
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)
	    {
		int nComps = nodes->degreesOfFreedomDistribution->getMyNumComponents();
		requireWrite(out);
#pragma omp parallel private(n)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nComps; n++)
		    {
			memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, n), numComps_size);
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
	    {
		int nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
		requireWrite(out);
#pragma omp parallel private(n,k)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nComps; n++)
		    {
			k = nodes->reducedDegreesOfFreedomMapping->map[n];
			memcpy(getSampleDataRWFast(out, n),
			       getSampleDataROFast(in, nodes->degreesOfFreedomMapping->target[k]), numComps_size);
		    }
		}
	    }

	/*********************** DUDLEY_REDUCED_DEGREES_OF_FREEDOM **************************************************/
	}
	else if (in_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
	{

	    if (out_data_type == DUDLEY_NODES)
	    {
		Dudley_setError(TYPE_ERROR,
				"Dudley_Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to nodes.");
	    }
	    else if (out_data_type == DUDLEY_REDUCED_NODES)
	    {
		coupler.reset(new paso::Coupler(nodes->reducedDegreesOfFreedomConnector, numComps));
		if (Esys_noError())
		{
		    upperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
		    /* safe provided coupler->copyAll is called before the pointer in "in" is invalidated */
		    const_cast<escript::Data*>(in)->resolve();
		    coupler->startCollect(in->getDataRO());  
		    recv_buffer = coupler->finishCollect();
		    requireWrite(out);
#pragma omp parallel private(n,k,l)
		    {
#pragma omp for schedule(static)
			for (n = 0; n < nodes->reducedNodesMapping->numTargets; n++)
			{
			    l = nodes->reducedNodesMapping->map[n];
			    k = nodes->reducedDegreesOfFreedomMapping->target[l];
			    if (k < upperBound)
			    {
				memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, k), numComps_size);
			    }
			    else
			    {
				memcpy(getSampleDataRWFast(out, n),
				       &recv_buffer[(k - upperBound) * numComps], numComps_size);
			    }
			}
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_REDUCED_DEGREES_OF_FREEDOM)
	    {
		int nComps = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
		requireWrite(out);
#pragma omp parallel private(n)
		{
#pragma omp for schedule(static)
		    for (n = 0; n < nComps; n++)
		    {
			memcpy(getSampleDataRWFast(out, n), getSampleDataROFast(in, n), numComps_size);
		    }
		}
	    }
	    else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)
	    {
		Dudley_setError(TYPE_ERROR,
				"Dudley_Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to degrees of freedom.");
	    }
	}
    }
    return;
}
