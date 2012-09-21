
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
/*                                                                                                         */
/*   Dudley: ElementFile                                                                                   */
/*                                                                                                         */
/*  reorders the elements in the element file such that the elements are stored close to the nodes         */
/*                                                                                                         */
/************************************************************************************/

#include "Util.h"
#include "ElementFile.h"

/************************************************************************************/

void Dudley_ElementFile_optimizeOrdering(Dudley_ElementFile ** in)
{
    Dudley_Util_ValueAndIndex *item_list = NULL;
    Dudley_ElementFile *out = NULL;
    dim_t e, i, NN;
    index_t *index = NULL;
    if (*in != NULL)
    {
	if ((*in)->numElements < 1)
	    return;
	NN = (*in)-> /*referenceElementSet-> */ numNodes;
	item_list = TMPMEMALLOC((*in)->numElements, Dudley_Util_ValueAndIndex);
	index = TMPMEMALLOC((*in)->numElements, index_t);
	if (!(Dudley_checkPtr(item_list) || Dudley_checkPtr(index)))
	{

	    out = Dudley_ElementFile_alloc((*in)->etype /*referenceElementSet */ , (*in)->MPIInfo);
	    if (Dudley_noError())
	    {
		Dudley_ElementFile_allocTable(out, (*in)->numElements);
		if (Dudley_noError())
		{
#pragma omp parallel for private(e,i) schedule(static)
		    for (e = 0; e < (*in)->numElements; e++)
		    {
			item_list[e].index = e;
			item_list[e].value = (*in)->Nodes[INDEX2(0, e, NN)];
			for (i = 1; i < NN; i++)
			    item_list[e].value = MIN(item_list[e].value, (*in)->Nodes[INDEX2(i, e, NN)]);
		    }
		    Dudley_Util_sortValueAndIndex((*in)->numElements, item_list);
#pragma omp parallel for private(e) schedule(static)
		    for (e = 0; e < (*in)->numElements; e++)
			index[e] = item_list[e].index;
		    Dudley_ElementFile_gather(index, *in, out);
		    Dudley_ElementFile_free(*in);
		    *in = out;
		}
		else
		{
		    Dudley_ElementFile_free(out);
		}
	    }
	}
	TMPMEMFREE(item_list);
	TMPMEMFREE(index);
    }
}
