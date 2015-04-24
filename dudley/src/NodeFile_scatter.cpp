
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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

/*   scatters the NodeFile in into NodeFile out using index[0:in->numNodes-1].  */
/*   index has to be between 0 and out->numNodes-1. */
/*   coloring is chosen for the worst case */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "NodeFile.h"

/************************************************************************************/

void Dudley_NodeFile_scatterEntries(dim_t n, index_t * index, index_t min_index, index_t max_index,
				    index_t * Id_out, index_t * Id_in,
				    index_t * Tag_out, index_t * Tag_in,
				    index_t * globalDegreesOfFreedom_out, index_t * globalDegreesOfFreedom_in,
				    dim_t numDim, double *Coordinates_out, double *Coordinates_in)
{
    dim_t i;
    register index_t k;
    const index_t range = max_index - min_index;
    const size_t numDim_size = (size_t) numDim * sizeof(double);

#pragma omp parallel for private(i,k) schedule(static)
    for (i = 0; i < n; i++)
    {
	k = index[i] - min_index;
	if ((k >= 0) && (k < range))
	{
	    Id_out[k] = Id_in[i];
	    Tag_out[k] = Tag_in[i];
	    globalDegreesOfFreedom_out[k] = globalDegreesOfFreedom_in[i];
	    memcpy(&(Coordinates_out[INDEX2(0, k, numDim)]), &(Coordinates_in[INDEX2(0, i, numDim)]), numDim_size);
	}
    }
}

void Dudley_NodeFile_scatter(index_t * index, Dudley_NodeFile * in, Dudley_NodeFile * out)
{
    Dudley_NodeFile_scatterEntries(out->numNodes, index, 0, in->numNodes,
				   out->Id, in->Id,
				   out->Tag, in->Tag,
				   out->globalDegreesOfFreedom, in->globalDegreesOfFreedom,
				   out->numDim, out->Coordinates, in->Coordinates);
}
