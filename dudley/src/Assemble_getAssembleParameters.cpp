
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
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

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif

using escript::ValueError;

namespace dudley {

AssembleParameters::AssembleParameters(const NodeFile* nodes,
                                       const ElementFile* ef,
                                       escript::ASM_ptr sm,
                                       escript::Data& rhs,
                                       bool reducedIntegrationOrder) :
    elements(ef),
    S(sm.get()),
    F(rhs),
    shapeFns(NULL)
{
    if (!rhs.isEmpty() && !rhs.actsExpanded()) {
        throw ValueError("AssembleParameters: Right hand side is not expanded.");
    }

    if (!getQuadShape(elements->numDim, reducedIntegrationOrder, &shapeFns)) {
        throw DudleyException("AssembleParameters: Cannot locate shape functions.");
    }

#ifdef ESYS_HAVE_PASO
    paso::SystemMatrix* pasoMat = sm ?
        dynamic_cast<paso::SystemMatrix*>(sm.get()) : NULL;

    // check the dimensions of matrix and rhs
    if (pasoMat != NULL && !rhs.isEmpty()) {
        const dim_t numRows = pasoMat->row_distribution->getMyNumComponents()*pasoMat->row_block_size;
        if (!rhs.numSamplesEqual(1, numRows/pasoMat->logical_row_block_size)) {
            throw ValueError("AssembleParameters: number of rows of matrix "
                             "and length of right hand side don't match.");
        }
    }
#endif

    // get the number of equations and components
    if (sm == NULL) {
        if (rhs.isEmpty()) {
            numEqu = 1;
        } else {
            numEqu = rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize() != sm->getRowBlockSize()) {
            throw ValueError("AssembleParameters: matrix row block size and "
                      "number of components of right hand side don't match.");
        }
        if (sm->getRowBlockSize() != sm->getColumnBlockSize())
            throw DudleyException("Dudley requires number of equations == number of solutions.");
        numEqu = sm->getRowBlockSize();
    }
    DOF = nodes->borrowTargetDegreesOfFreedom();
    DOF_UpperBound = nodes->getNumDegreesOfFreedom();

#ifdef ESYS_HAVE_PASO
    // get the information for the labeling of the degrees of freedom from
    // the matrix
    if (pasoMat) {
        // Make sure # rows in matrix == num local DOF
        const index_t numRows = pasoMat->row_distribution->getMyNumComponents()*pasoMat->row_block_size;
        const index_t numCols = pasoMat->col_distribution->getMyNumComponents()*pasoMat->col_block_size;

        if (numRows != numEqu * nodes->getNumDegreesOfFreedom()) {
            throw DudleyException("AssembleParameters: number of rows in "
                                  "matrix does not match the number of "
                                  "degrees of freedom in mesh");
        }
        // Make sure # cols in matrix == num local DOF
        if (numCols != numRows) {
            throw DudleyException("AssembleParameters: number of columns in "
                                  "matrix does not match the number of "
                                  "degrees of freedom in mesh");
        }
    }
#endif

    // get the information from right hand side
    if (!rhs.isEmpty() &&
            !rhs.numSamplesEqual(1, nodes->getNumDegreesOfFreedom())) {
        throw DudleyException("AssembleParameters: length of RHS vector does "
                              "not match the number of degrees of freedom "
                              "in mesh");
    }

    jac = elements->borrowJacobians(nodes, reducedIntegrationOrder);

    if (elements->numNodes < jac->numShapes) {
        throw DudleyException("AssembleParameters: too many nodes are "
                              "expected by row.");
    }
    if (jac->numElements != elements->numElements) {
        throw DudleyException("AssembleParameters: number of elements for "
                              "row is wrong.");
    }

    NN = elements->numNodes;
    numQuad = jac->numQuad;
    numDim = jac->numDim;
    numShapes = jac->numShapes;
}

} // namespace dudley

