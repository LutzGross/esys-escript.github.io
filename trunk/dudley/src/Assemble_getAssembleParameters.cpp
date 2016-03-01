
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

/****************************************************************************

  Assemblage routines: prepares the assemble parameter set

*****************************************************************************/

#include "Assemble.h"
#include "ShapeTable.h"

#include <paso/SystemMatrix.h>

namespace dudley {

void Assemble_getAssembleParameters(const Dudley_NodeFile* nodes,
        const Dudley_ElementFile* elements, escript::ASM_ptr mat,
        const escript::Data& F, bool reducedIntegrationOrder,
        Assemble_Parameters* parm)
{
    paso::SystemMatrix* S = dynamic_cast<paso::SystemMatrix*>(mat.get());
    if (mat && !S) {
        throw DudleyException("Assemble_getAssembleParameters: Unknown matrix type.");
    }
    parm->shapeFns = NULL;
    if (!F.isEmpty() && !F.actsExpanded()) {
        throw DudleyException("Assemble_getAssembleParameters: Right hand side is not expanded.");
    }

    if (!getQuadShape(elements->numDim, reducedIntegrationOrder, &parm->shapeFns)) {
        throw DudleyException("Assemble_getAssembleParameters: Can not locate shape functions.");
    }
    /*  check the dimensions of S and F */
    if (S != NULL && !F.isEmpty()) {
        if (!F.numSamplesEqual(1, (S->row_distribution->getMyNumComponents()*S->row_block_size) /
             S->logical_row_block_size)) {
            throw DudleyException("Assemble_getAssembleParameters: number of rows of matrix and length of right hand side don't match.");
        }
    }
    /* get the number of equations and components */
    if (S == NULL) {
        if (F.isEmpty()) {
            parm->numEqu = 1;
            parm->numComp = 1;
        } else {
            parm->numEqu = F.getDataPointSize();
            parm->numComp = parm->numEqu;
        }
    } else {
        if (F.isEmpty()) {
            parm->numEqu = S->logical_row_block_size;
            parm->numComp = S->logical_col_block_size;
        } else {
            if (F.getDataPointSize() != S->logical_row_block_size) {
                throw DudleyException("Assemble_getAssembleParameters: matrix row block size and number of components of right hand side don't match.");
            }
            parm->numEqu = S->logical_row_block_size;
            parm->numComp = S->logical_col_block_size;
        }
    }
    parm->col_DOF = nodes->degreesOfFreedomMapping->target;
    parm->row_DOF = nodes->degreesOfFreedomMapping->target;
    /* get the information for the labeling of the degrees of freedom from matrix */
    if (S != NULL) {
        // Make sure # rows in matrix == num DOF for one of:
        // full or reduced (use numLocalDOF for MPI)
        if (S->row_distribution->getMyNumComponents() * S->row_block_size ==
            parm->numEqu * nodes->degreesOfFreedomDistribution->getMyNumComponents()) {
            parm->row_DOF_UpperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
            parm->row_DOF = nodes->degreesOfFreedomMapping->target;
            parm->row_jac = Dudley_ElementFile_borrowJacobians(elements, nodes, reducedIntegrationOrder);
        } else if (S->row_distribution->getMyNumComponents() * S->row_block_size ==
                 parm->numEqu * nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents()) {
            parm->row_DOF_UpperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
            parm->row_DOF = nodes->reducedDegreesOfFreedomMapping->target;
            parm->row_jac = Dudley_ElementFile_borrowJacobians(elements, nodes, reducedIntegrationOrder);
        } else {
            throw DudleyException("Assemble_getAssembleParameters: number of rows in matrix does not match the number of degrees of freedom in mesh");
        }
        // Make sure # cols in matrix == num DOF for one of:
        // full or reduced (use numLocalDOF for MPI)
        if (S->col_distribution->getMyNumComponents() * S->col_block_size ==
            parm->numComp * nodes->degreesOfFreedomDistribution->getMyNumComponents()) {
            parm->col_DOF_UpperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
            parm->col_DOF = nodes->degreesOfFreedomMapping->target;
            parm->row_jac = Dudley_ElementFile_borrowJacobians(elements, nodes, reducedIntegrationOrder);
        } else if (S->col_distribution->getMyNumComponents() * S->col_block_size ==
                 parm->numComp * nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents()) {
            parm->col_DOF_UpperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
            parm->col_DOF = nodes->reducedDegreesOfFreedomMapping->target;
            parm->row_jac = Dudley_ElementFile_borrowJacobians(elements, nodes, reducedIntegrationOrder);
        } else {
            throw DudleyException("Assemble_getAssembleParameters: number of columns in matrix does not match the number of degrees of freedom in mesh");
        }
    }

    // get the information from right hand side
    if (!F.isEmpty()) {
        if (F.numSamplesEqual(1, nodes->degreesOfFreedomDistribution->getMyNumComponents())) {
            parm->row_DOF_UpperBound = nodes->degreesOfFreedomDistribution->getMyNumComponents();
            parm->row_DOF = nodes->degreesOfFreedomMapping->target;
            parm->row_jac = Dudley_ElementFile_borrowJacobians(elements, nodes, reducedIntegrationOrder);
        } else if (F.numSamplesEqual(1, nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents())) {
            parm->row_DOF_UpperBound = nodes->reducedDegreesOfFreedomDistribution->getMyNumComponents();
            parm->row_DOF = nodes->reducedDegreesOfFreedomMapping->target;
            parm->row_jac = Dudley_ElementFile_borrowJacobians(elements, nodes, reducedIntegrationOrder);
        } else {
            throw DudleyException("Assemble_getAssembleParameters: length of RHS vector does not match the number of degrees of freedom in mesh");
        }
        if (S == NULL) {
            parm->col_DOF_UpperBound = parm->row_DOF_UpperBound;
            parm->col_DOF = parm->row_DOF;
            parm->row_jac = parm->row_jac;
        }
    }

    if (parm->row_jac->numDim != parm->row_jac->numDim) {
        throw DudleyException("Assemble_getAssembleParameters: spacial dimension for row and column shape function must match.");
    }

    if (elements->numNodes < parm->row_jac->numShapes) {
        throw DudleyException("Assemble_getAssembleParameters: too many nodes are expected by row.");
    }
    if (parm->row_jac->numElements != elements->numElements) {
        throw DudleyException("Assemble_getAssembleParameters: number of elements for row is wrong.");
    }

    parm->numQuad = parm->row_jac->numQuad;
    parm->NN = elements->numNodes;
    parm->numElements = elements->numElements;
    parm->numDim = parm->row_jac->numDim;
    parm->numShapes = parm->row_jac->numShapes;
}

} // namespace dudley

