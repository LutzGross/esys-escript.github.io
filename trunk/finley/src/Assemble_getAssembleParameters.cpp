
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


/****************************************************************************

  Assemblage routines: prepares the assemble parameter set

*****************************************************************************/

#include "Assemble.h"

namespace finley {

AssembleParameters::AssembleParameters(const NodeFile* nodes,
                                       const ElementFile* ef,
                                       Paso_SystemMatrix* sm,
                                       escript::Data& rhs,
                                       bool reducedOrder)
    : elements(ef),
      S(sm),
      F(rhs)
{
    int numSub, numQuadSub;
    resetError();

    if (!rhs.isEmpty() && !rhs.actsExpanded()) {
        setError(TYPE_ERROR, "AssembleParameters: Right hand side is not expanded.");
        return;
    }
    // check the dimensions of S and rhs
    if (sm!=NULL && !rhs.isEmpty()) {
        if (!rhs.numSamplesEqual(1, (sm->row_distribution->getMyNumComponents()
                           * sm->row_block_size)/sm->logical_row_block_size)) {
            setError(TYPE_ERROR, "AssembleParameters: number of rows of matrix and length of right hand side don't match.");
            return;
        }
    }

    // get the number of equations and components
    if (sm==NULL) {
        if (rhs.isEmpty()) {
            this->numEqu=1;
            this->numComp=1;
        } else {
            this->numEqu=rhs.getDataPointSize();
            this->numComp=this->numEqu;
        }
    } else {
        if (rhs.isEmpty()) {
            this->numEqu=sm->logical_row_block_size;
            this->numComp=sm->logical_col_block_size;
        } else {
            if (rhs.getDataPointSize() != sm->logical_row_block_size) {
                setError(TYPE_ERROR,"AssembleParameters: matrix row block size and number of components of right hand side don't match.");
                return;
            }
            this->numEqu=sm->logical_row_block_size;
            this->numComp=sm->logical_col_block_size;
        }
    }
    this->col_DOF=nodes->borrowTargetDegreesOfFreedom();
    this->row_DOF=nodes->borrowTargetDegreesOfFreedom();
    // get the information for the labeling of the degrees of freedom from
    // the matrix
    if (sm!=NULL) {
        // Make sure # rows in matrix == num DOF for one of:
        // full or reduced (use numLocalDOF for MPI)
        if (sm->row_distribution->getMyNumComponents()*sm->row_block_size ==
                this->numEqu*nodes->getNumDegreesOfFreedom()) {
            this->row_DOF_UpperBound = nodes->getNumDegreesOfFreedom();
            this->row_DOF=nodes->borrowTargetDegreesOfFreedom();
            this->row_jac=ef->borrowJacobians(nodes, false, reducedOrder);
        } else if (sm->row_distribution->getMyNumComponents()*sm->row_block_size ==
                this->numEqu*nodes->getNumReducedDegreesOfFreedom()) {
            this->row_DOF_UpperBound = nodes->getNumReducedDegreesOfFreedom();
            this->row_DOF=nodes->borrowTargetReducedDegreesOfFreedom();
            this->row_jac=ef->borrowJacobians(nodes, true, reducedOrder);
        } else {
            setError(TYPE_ERROR, "AssembleParameters: number of rows in matrix does not match the number of degrees of freedom in mesh");
        }
        // Make sure # cols in matrix == num DOF for one of:
        // full or reduced (use numLocalDOF for MPI)
        if (sm->col_distribution->getMyNumComponents()*sm->col_block_size ==
                this->numComp*nodes->getNumDegreesOfFreedom()) {
            this->col_DOF_UpperBound = nodes->getNumDegreesOfFreedom();
            this->col_DOF=nodes->borrowTargetDegreesOfFreedom();
            this->col_jac=ef->borrowJacobians(nodes, false, reducedOrder);
        } else if (sm->col_distribution->getMyNumComponents()*sm->col_block_size==this->numComp*nodes->getNumReducedDegreesOfFreedom()) {
            this->col_DOF_UpperBound = nodes->getNumReducedDegreesOfFreedom();
            this->col_DOF=nodes->borrowTargetReducedDegreesOfFreedom();
            this->col_jac=ef->borrowJacobians(nodes, true, reducedOrder);
        } else {
            setError(TYPE_ERROR, "AssembleParameters: number of columns in matrix does not match the number of degrees of freedom in mesh");
        }
    }

    if (!noError())
        return;

    // get the information from right hand side
    if (!rhs.isEmpty()) {
        if (rhs.numSamplesEqual(1, nodes->getNumDegreesOfFreedom())) {
            this->row_DOF_UpperBound = nodes->getNumDegreesOfFreedom();
            this->row_DOF=nodes->borrowTargetDegreesOfFreedom();
            this->row_jac=ef->borrowJacobians(nodes, false, reducedOrder);
        } else if (rhs.numSamplesEqual(1, nodes->getNumReducedDegreesOfFreedom())) {
            this->row_DOF_UpperBound = nodes->getNumReducedDegreesOfFreedom();
            this->row_DOF=nodes->borrowTargetReducedDegreesOfFreedom();
            this->row_jac=ef->borrowJacobians(nodes, true, reducedOrder);
        } else {
            setError(TYPE_ERROR, "AssembleParameters: length of RHS vector does not match the number of degrees of freedom in mesh");
        }
        if (sm==NULL) {
            this->col_DOF_UpperBound=this->row_DOF_UpperBound;
            this->col_DOF=this->row_DOF;
            this->col_jac=this->row_jac;
        }
    }

    numSub=std::min(this->row_jac->numSub, this->col_jac->numSub);
    numQuadSub=this->row_jac->numQuadTotal/numSub;
    if (this->row_jac->numSides != this->col_jac->numSides) {
        setError(TYPE_ERROR, "AssembleParameters: number of sides for row and column shape functions must match.");
    }
    if (this->row_jac->numDim != this->col_jac->numDim) {
        setError(TYPE_ERROR, "AssembleParameters: spatial dimension for row and column shape function must match.");
    }
    if (ef->numNodes < this->row_jac->numShapesTotal) {
        setError(TYPE_ERROR, "AssembleParameters: too many nodes are expected by row.");
    }
    if (ef->numNodes < this->col_jac->numShapesTotal) {
        setError(TYPE_ERROR, "AssembleParameters: too many nodes are expected by col.");
    }
    if (this->row_jac->numElements != ef->numElements) {
        setError(TYPE_ERROR, "AssembleParameters: number of elements for row is wrong.");
    }
    if (this->col_jac->numElements != ef->numElements) {
        setError(TYPE_ERROR, "AssembleParameters: number of elements for column is wrong.");
    }
    if (this->row_jac->numQuadTotal != this->col_jac->numQuadTotal) {
        setError(TYPE_ERROR, "AssembleParameters: number of quadrature points for row and column shape functions must match.");
    }
    // to consider different basis function for rows and columns this will
    // require some work:
    if (numQuadSub*numSub != this->row_jac->numQuadTotal) {
        setError(TYPE_ERROR, "AssembleParameters: number of quadrature points for row is not correct.");
    }
    if (numQuadSub != this->row_jac->BasisFunctions->numQuadNodes) {
        setError(TYPE_ERROR, "AssembleParameters: Incorrect number of quadrature points for row.");
    }
    if (numQuadSub != this->col_jac->BasisFunctions->numQuadNodes) {
        setError(TYPE_ERROR, "AssembleParameters: Incorrect number of quadrature points for column.");
    }

    this->numQuadSub=numQuadSub;
    this->numSub=numSub;
    this->numQuadTotal=this->row_jac->numQuadTotal;
    this->NN=elements->numNodes;
    this->numElements=elements->numElements;
    this->numDim=this->row_jac->numDim;
    this->col_node=this->col_jac->node_selection;
    this->row_node=this->row_jac->node_selection;
    this->numSides=this->row_jac->numSides;
    this->row_numShapesTotal=this->row_jac->numShapesTotal;
    this->row_numShapes=this->row_jac->BasisFunctions->Type->numShapes;
    this->col_numShapesTotal=this->col_jac->numShapesTotal;
    this->col_numShapes=this->col_jac->BasisFunctions->Type->numShapes;
}

} // namespace finley

