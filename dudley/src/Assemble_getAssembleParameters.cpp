
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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

/*    assemblage routines: prepares the assemble parameter set */

/************************************************************************************/

#include "Assemble.h"
#include "ShapeTable.h"

/************************************************************************************/

void Dudley_Assemble_getAssembleParameters(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, Paso_SystemMatrix * S,
				    escriptDataC * F, bool reducedIntegrationOrder, Dudley_Assemble_Parameters * parm)
{
    Dudley_resetError();
    parm->shapeFns = NULL;
    if (!isEmpty(F) && !isExpanded(F))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getAssembleParameters: Right hand side is not expanded.");
	return;
    }

    if (!getQuadShape(elements->numDim, reducedIntegrationOrder, &(parm->shapeFns)))
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getAssembleParameters: Can not locate shape functions.");
    }
    /*  check the dimensions of S and F */
    if (S != NULL && !isEmpty(F))
    {
	if (!numSamplesEqual
	    (F, 1,
	     (S->row_distribution->getMyNumComponents() * S->row_block_size) /
	     S->logical_row_block_size))
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_getAssembleParameters: number of rows of matrix and length of right hand side don't match.");
	    return;
	}
    }
    /* get the number of equations and components */
    if (S == NULL)
    {
	if (isEmpty(F))
	{
	    parm->numEqu = 1;
	    parm->numComp = 1;
	}
	else
	{
	    parm->numEqu = getDataPointSize(F);
	    parm->numComp = parm->numEqu;
	}
    }
    else
    {
	if (isEmpty(F))
	{
	    parm->numEqu = S->logical_row_block_size;
	    parm->numComp = S->logical_col_block_size;
	}
	else
	{
	    if (getDataPointSize(F) != S->logical_row_block_size)
	    {
		Dudley_setError(TYPE_ERROR,
				"Dudley_Assemble_getAssembleParameters: matrix row block size and number of components of right hand side don't match.");
		return;
	    }
	    parm->numEqu = S->logical_row_block_size;
	    parm->numComp = S->logical_col_block_size;
	}
    }
    parm->col_DOF = nodes->degreesOfFreedomMapping->target;
    parm->row_DOF = nodes->degreesOfFreedomMapping->target;
    /* get the information for the labeling of the degrees of freedom from matrix */
    if (S != NULL)
    {
	/* Make sure # rows in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
	if (S->row_distribution->getMyNumComponents() * S->row_block_size ==
	    parm->numEqu * nodes->degreesOfFreedomDistribution->getMyNumComponents())
	{
	    parm->row_DOF_UpperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
	    parm->row_DOF = nodes->degreesOfFreedomMapping->target;
	    parm->row_jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, reducedIntegrationOrder);
	}
	else if (S->row_distribution->getMyNumComponents() * S->row_block_size ==
		 parm->numEqu * nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents())
	{
	    parm->row_DOF_UpperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
	    parm->row_DOF = nodes->reducedDegreesOfFreedomMapping->target;
	    parm->row_jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, reducedIntegrationOrder);
	}
	else
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_getAssembleParameters: number of rows in matrix does not match the number of degrees of freedom in mesh");
	}
	/* Make sure # cols in matrix == num DOF for one of: full or reduced (use numLocalDOF for MPI) */
	if (S->col_distribution->getMyNumComponents() * S->col_block_size ==
	    parm->numComp * nodes->degreesOfFreedomDistribution->getMyNumComponents())
	{
	    parm->col_DOF_UpperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
	    parm->col_DOF = nodes->degreesOfFreedomMapping->target;
	    parm->row_jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, reducedIntegrationOrder);
	}
	else if (S->col_distribution->getMyNumComponents() * S->col_block_size ==
		 parm->numComp * nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents())
	{
	    parm->col_DOF_UpperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
	    parm->col_DOF = nodes->reducedDegreesOfFreedomMapping->target;
	    parm->row_jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, reducedIntegrationOrder);
	}
	else
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_getAssembleParameters: number of columns in matrix does not match the number of degrees of freedom in mesh");
	}
    }
    if (!Dudley_noError())
	return;
    /* get the information from right hand side */
    if (!isEmpty(F))
    {
	if (numSamplesEqual(F, 1, nodes->degreesOfFreedomDistribution->getMyNumComponents()))
	{
	    parm->row_DOF_UpperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
	    parm->row_DOF = nodes->degreesOfFreedomMapping->target;
	    parm->row_jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, reducedIntegrationOrder);
	}
	else if (numSamplesEqual
		 (F, 1, nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents()))
	{
	    parm->row_DOF_UpperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
	    parm->row_DOF = nodes->reducedDegreesOfFreedomMapping->target;
	    parm->row_jac = Dudley_ElementFile_borrowJacobeans(elements, nodes, reducedIntegrationOrder);
	}
	else
	{
	    Dudley_setError(TYPE_ERROR,
			    "Dudley_Assemble_getAssembleParameters: length of RHS vector does not match the number of degrees of freedom in mesh");
	}
	if (S == NULL)
	{
	    parm->col_DOF_UpperBound = parm->row_DOF_UpperBound;
	    parm->col_DOF = parm->row_DOF;
	    parm->row_jac = parm->row_jac;
	}
    }

    if (parm->row_jac->numDim != parm->row_jac->numDim)
    {
	Dudley_setError(TYPE_ERROR,
			"Dudley_Assemble_getAssembleParameters: spacial dimension for row and column shape function must match.");
    }

    if (elements->numNodes < parm->row_jac->numShapes)
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getAssembleParameters: too many nodes are expected by row.");
    }
    if (parm->row_jac->numElements != elements->numElements)
    {
	Dudley_setError(TYPE_ERROR, "Dudley_Assemble_getAssembleParameters: number of elements for row is wrong.");
    }

    parm->numQuad = parm->row_jac->numQuad;
    parm->NN = elements->numNodes;
    parm->numElements = elements->numElements;
    parm->numDim = parm->row_jac->numDim;
    parm->numShapes = parm->row_jac->numShapes;

}
