
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

/*	 Dudley: Mesh: ElementFile */

/*	set tags to newTag where mask>0 */

/************************************************************************************/

#include "ElementFile.h"
#include "Util.h"
#include "Assemble.h"

/************************************************************************************/

void Dudley_ElementFile_setTags(Dudley_ElementFile * self, const int newTag, escriptDataC * mask)
{
    register dim_t n, q;
    dim_t numElements, numQuad;
    register __const double *mask_array;
    register bool_t check;
    Dudley_resetError();
    if (self == NULL)
	return;
    numElements = self->numElements;

    numQuad = Dudley_Assemble_reducedIntegrationOrder(mask) ? 1 : (self->numDim + 1);
    if (1 != getDataPointSize(mask))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_ElementFile_setTags: number of components of mask is 1.");
    }
    else if (!numSamplesEqual(mask, numQuad, numElements))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_ElementFile_setTags: illegal number of samples of mask Data object");
    }

    /* now we can start */

    if (Dudley_noError())
    {
	if (isExpanded(mask))
	{
#pragma omp parallel private(n,check,mask_array)
	    {
#pragma omp for schedule(static)
		for (n = 0; n < numElements; n++)
		{
		    mask_array = getSampleDataRO(mask, n);
		    if (mask_array[0] > 0)
			self->Tag[n] = newTag;
		}
	    }
	}
	else
	{
#pragma omp parallel private(q,n,check,mask_array)
	    {
#pragma omp for schedule(static)
		for (n = 0; n < numElements; n++)
		{
		    mask_array = getSampleDataRO(mask, n);
		    check = FALSE;
		    for (q = 0; q < numQuad; q++)
			check = check || mask_array[q];
		    if (check)
			self->Tag[n] = newTag;
		}
	    }
	}
	Dudley_ElementFile_setTagsInUse(self);
    }
}
