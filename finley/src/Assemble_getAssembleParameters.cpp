
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

/****************************************************************************

  Assemblage routines: prepares the assemble parameter set

*****************************************************************************/

#include "Assemble.h"

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif

using escript::ValueError;

namespace finley {

AssembleParameters::AssembleParameters(const NodeFile* nodes,
                                       const ElementFile* ef,
                                       escript::ASM_ptr sm,
                                       escript::Data& rhs,
                                       bool reducedOrder) :
    elements(ef),
    S(sm),
    F(rhs)
{
    if (!rhs.isEmpty() && !rhs.actsExpanded()) {
        throw ValueError("AssembleParameters: Right hand side is not expanded.");
    }

#ifdef ESYS_HAVE_PASO
    paso::SystemMatrix<double>* pasoMat = sm ?
        dynamic_cast<paso::SystemMatrix<double>*>(sm.get()) : NULL;

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
            numEqu = numComp = 1;
        } else {
            numEqu = numComp = rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize() != sm->getRowBlockSize()) {
            throw ValueError("AssembleParameters: matrix row block size and "
                      "number of components of right hand side don't match.");
        }
        numEqu = sm->getRowBlockSize();
        numComp = sm->getColumnBlockSize();
    }
    // set some defaults
    row_DOF = nodes->borrowTargetDegreesOfFreedom();
    row_DOF_UpperBound = nodes->getNumDegreesOfFreedom();
    row_jac = ef->borrowJacobians(nodes, false, reducedOrder);
    col_DOF = row_DOF;
    col_DOF_UpperBound = row_DOF_UpperBound;
    col_jac = row_jac;

#ifdef ESYS_HAVE_PASO
    // get the information for the labeling of the degrees of freedom from
    // the matrix
    if (pasoMat) {
        // Make sure # rows in matrix == num DOF for one of:
        // full or reduced (use numLocalDOF for MPI)
        const index_t numRows = pasoMat->row_distribution->getMyNumComponents()*pasoMat->row_block_size;
        const index_t numCols = pasoMat->col_distribution->getMyNumComponents()*pasoMat->col_block_size;

        if (numRows == numEqu * nodes->getNumDegreesOfFreedom()) {
        } else if (numRows == numEqu * nodes->getNumReducedDegreesOfFreedom()) {
            row_DOF_UpperBound = nodes->getNumReducedDegreesOfFreedom();
            row_DOF = nodes->borrowTargetReducedDegreesOfFreedom();
            row_jac = ef->borrowJacobians(nodes, true, reducedOrder);
        } else {
            throw ValueError("AssembleParameters: number of rows in matrix "
                   "does not match the number of degrees of freedom in mesh");
        }
        // Make sure # cols in matrix == num DOF for one of:
        // full or reduced (use numLocalDOF for MPI)
        if (numCols == this->numComp * nodes->getNumDegreesOfFreedom()) {
        } else if (numCols == this->numComp * nodes->getNumReducedDegreesOfFreedom()) {
            col_DOF_UpperBound = nodes->getNumReducedDegreesOfFreedom();
            col_DOF = nodes->borrowTargetReducedDegreesOfFreedom();
            col_jac = ef->borrowJacobians(nodes, true, reducedOrder);
        } else {
            throw ValueError("AssembleParameters: number of columns in matrix "
                   "does not match the number of degrees of freedom in mesh");
        }
    }
#endif

    // get the information from right hand side
    if (!rhs.isEmpty()) {
        if (rhs.numSamplesEqual(1, nodes->getNumDegreesOfFreedom())) {
        } else if (rhs.numSamplesEqual(1, nodes->getNumReducedDegreesOfFreedom())) {
            row_DOF_UpperBound = nodes->getNumReducedDegreesOfFreedom();
            row_DOF = nodes->borrowTargetReducedDegreesOfFreedom();
            row_jac = ef->borrowJacobians(nodes, true, reducedOrder);
        } else {
            throw ValueError("AssembleParameters: length of RHS vector does not match the number of degrees of freedom in mesh");
        }
#ifdef ESYS_HAVE_PASO
        if (sm == NULL) {
            col_DOF_UpperBound = this->row_DOF_UpperBound;
            col_DOF = this->row_DOF;
            col_jac = this->row_jac;
        }
#else // trilinos case
        col_DOF_UpperBound = this->row_DOF_UpperBound;
        col_DOF = this->row_DOF;
        col_jac = this->row_jac;
#endif
    }

    numSub = std::min(row_jac->numSub, col_jac->numSub);
    numQuadSub = row_jac->numQuadTotal / numSub;
    if (row_jac->numSides != col_jac->numSides) {
        throw ValueError("AssembleParameters: number of sides for row and "
                         "column shape functions must match.");
    }
    if (row_jac->numDim != col_jac->numDim) {
        throw ValueError("AssembleParameters: spatial dimension for row and "
                         "column shape function must match.");
    }
    if (ef->numNodes < row_jac->numShapesTotal) {
        throw ValueError("AssembleParameters: too many nodes are expected by row.");
    }
    if (ef->numNodes < col_jac->numShapesTotal) {
        throw ValueError("AssembleParameters: too many nodes are expected by col.");
    }
    if (row_jac->numElements != ef->numElements) {
        throw ValueError("AssembleParameters: number of elements for row is wrong.");
    }
    if (col_jac->numElements != ef->numElements) {
        throw ValueError("AssembleParameters: number of elements for column is wrong.");
    }
    if (row_jac->numQuadTotal != col_jac->numQuadTotal) {
        throw ValueError("AssembleParameters: number of quadrature points for "
                         "row and column shape functions must match.");
    }
    // to consider different basis function for rows and columns this will
    // require some work:
    if (numQuadSub * numSub != row_jac->numQuadTotal) {
        throw ValueError("AssembleParameters: number of quadrature points "
                         "for row is not correct.");
    }
    if (numQuadSub != row_jac->BasisFunctions->numQuadNodes) {
        throw ValueError("AssembleParameters: Incorrect number of quadrature "
                         "points for row.");
    }
    if (numQuadSub != col_jac->BasisFunctions->numQuadNodes) {
        throw ValueError("AssembleParameters: Incorrect number of quadrature "
                         "points for column.");
    }

    NN = elements->numNodes;
    numQuadTotal = row_jac->numQuadTotal;
    numElements = elements->numElements;
    numDim = row_jac->numDim;
    col_node = col_jac->node_selection;
    row_node = row_jac->node_selection;
    numSides = row_jac->numSides;
    row_numShapesTotal = row_jac->numShapesTotal;
    row_numShapes = row_jac->BasisFunctions->Type->numShapes;
    col_numShapesTotal = col_jac->numShapesTotal;
    col_numShapes = col_jac->BasisFunctions->Type->numShapes;
}

} // namespace finley

