
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

/*   Dudley: ElementFile */

/*   mark the used nodes with offset: */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "ElementFile.h"

/************************************************************************************/

void Dudley_ElementFile_markNodes(index_t * mask, index_t offset, dim_t numNodes, Dudley_ElementFile * in,
				  bool useLinear)
{
    dim_t i, NN, e;
    if (in != NULL)
    {
	NN = in->numNodes;
#pragma omp parallel for private(e,i) schedule(static)
	for (e = 0; e < in->numElements; e++)
	{
	    for (i = 0; i < NN; i++)
	    {
		mask[in->Nodes[INDEX2(i, e, NN)] - offset] = 1;
	    }
	}
    }
}

void Dudley_ElementFile_markDOFsConnectedToRange(index_t * mask, index_t offset, index_t marker, index_t firstDOF,
						 index_t lastDOF, index_t * dofIndex, Dudley_ElementFile * in,
						 bool useLinear)
{
    dim_t i, NN, e, j;
    index_t color;
    register index_t k;

    if (in != NULL)
    {
	NN = in->numNodes;
	for (color = in->minColor; color <= in->maxColor; color++)
	{
#pragma omp parallel for private(e,i,j,k) schedule(static)
	    for (e = 0; e < in->numElements; e++)
	    {
		if (in->Color[e] == color)
		{
		    for (i = 0; i < NN; i++)
		    {
			k = dofIndex[in->Nodes[INDEX2(i, e, NN)]];
			if ((firstDOF <= k) && (k < lastDOF))
			{
			    for (j = 0; j < NN; j++)
				mask[dofIndex[in->Nodes[INDEX2(j, e, NN)]] - offset] = marker;
			    break;
			}
		    }
		}
	    }
	}
    }
}
