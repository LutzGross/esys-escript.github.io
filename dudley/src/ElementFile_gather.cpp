
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (SUCCESS)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/************************************************************************************/

/*   Dudley: Element File */

/*   gathers the Element File out from the  Element File in using index[0:out->elements-1].  */
/*   index has to be between 0 and in->elements-1. */
/*   a conservative assumption on the coloring is made */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "ElementFile.h"

/************************************************************************************/

void Dudley_ElementFile_gather(index_t * index, Dudley_ElementFile * in, Dudley_ElementFile * out)
{
    index_t k;
    dim_t e, j;
    dim_t NN_in = in->numNodes;
    dim_t NN_out = out->numNodes;
    if (in != NULL)
    {
	/*OMP */
#pragma omp parallel for private(e,k,j) schedule(static)
	for (e = 0; e < out->numElements; e++)
	{
	    k = index[e];
	    out->Id[e] = in->Id[k];
	    out->Tag[e] = in->Tag[k];
	    out->Owner[e] = in->Owner[k];
	    out->Color[e] = in->Color[k] + out->maxColor + 1;
	    for (j = 0; j < MIN(NN_out, NN_in); j++)
		out->Nodes[INDEX2(j, e, NN_out)] = in->Nodes[INDEX2(j, k, NN_in)];
	}
	out->minColor = MIN(out->minColor, in->minColor + out->maxColor + 1);
	out->maxColor = MAX(out->maxColor, in->maxColor + out->maxColor + 1);
    }
}
