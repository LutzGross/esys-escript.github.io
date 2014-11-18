
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

  Assemblage routines: header file

*****************************************************************************/

#ifndef __FINLEY_ASSEMBLE_H__
#define __FINLEY_ASSEMBLE_H__

#include "Finley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "paso/SystemMatrix.h"

namespace finley {

struct AssembleParameters {
    AssembleParameters(const NodeFile* nodes, const ElementFile* ef,
                       paso::SystemMatrix_ptr sm, escript::Data& rhs,
                       bool reducedOrder);

    /// element file these parameters apply to
    const ElementFile* elements;
    /// system matrix to be updated
    paso::SystemMatrix_ptr S;
    /// right-hand side to be updated
    escript::Data& F;
    /// total number of quadrature nodes = numQuadSub * numQuadSub
    int numQuadTotal;
    /// number of quadrature nodes per subelements
    int numQuadSub;
    /// number of sides
    int numSides;
    /// number of sub-elements
    int numSub;
    /// number of spatial dimensions
    int numDim;
    /// leading dimension of element node table
    int NN;
    /// number of elements
    int numElements;

    int numEqu;
    const int* row_DOF;
    int row_DOF_UpperBound;
    ElementFile_Jacobians* row_jac;
    const int* row_node;
    int row_numShapesTotal;
    int row_numShapes;
    int numComp;
    const int* col_DOF;
    int col_DOF_UpperBound;
    ElementFile_Jacobians* col_jac;
    const int* col_node;
    int col_numShapesTotal;
    int col_numShapes;
};


/// Entry point for PDE assembly. Checks arguments, populates an
/// AssembleParameters structure and calls appropriate method for the actual
/// work.
void Assemble_PDE(const NodeFile* nodes, const ElementFile* elements,
                  paso::SystemMatrix_ptr S, escript::Data& F,
                  const escript::Data& A, const escript::Data& B,
                  const escript::Data& C, const escript::Data& D,
                  const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Points(const AssembleParameters& p, const escript::Data& d_dirac,
                         const escript::Data& y_dirac);

void Assemble_PDE_Single_1D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Single_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Single_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Single_C(const AssembleParameters& p, const escript::Data& D,
                           const escript::Data& Y);

void Assemble_PDE_System_1D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_System_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_System_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_System_C(const AssembleParameters& p, const escript::Data& D,
                           const escript::Data& Y);

void Assemble_addToSystemMatrix(paso::SystemMatrix_ptr S, const int NN_Equa,
        const int* Nodes_Equa, const int num_Equa, const int NN_Sol,
        const int* Nodes_Sol, const int num_Sol, const double* array);

void Assemble_LumpedSystem(const NodeFile* nodes, const ElementFile* elements,
                           escript::Data& lumpedMat, const escript::Data& D,
                           bool useHRZ);

void Assemble_AverageElementData(const ElementFile* elements,
                                 escript::Data& out, const escript::Data& in);

void Assemble_CopyElementData(const ElementFile* elements, escript::Data& out,
                              const escript::Data& in);

void Assemble_CopyNodalData(const NodeFile* nodes, escript::Data& out,
                            const escript::Data& in);

void Assemble_NodeCoordinates(const NodeFile* nodes, escript::Data& out);

void Assemble_getNormal(const NodeFile* nodes, const ElementFile* elements,
                        escript::Data& normals);

void Assemble_getSize(const NodeFile* nodes, const ElementFile* elements,
                      escript::Data& size);

void Assemble_gradient(const NodeFile* nodes, const ElementFile* elements,
                       escript::Data& gradient, const escript::Data& data);

void Assemble_integrate(const NodeFile* nodes, const ElementFile* elements,
                        const escript::Data& data, double* integrals);

void Assemble_interpolate(const NodeFile* nodes, const ElementFile* elements,
                          const escript::Data& data, escript::Data& output);

void Assemble_jacobians_1D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_2D_M1D_E1D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_2D_M1D_E1D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_2D_M1D_E2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_2D_M1D_E2D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_3D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_3D_M2D_E2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_3D_M2D_E2D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_3D_M2D_E3D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);
void Assemble_jacobians_3D_M2D_E3D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           int numElements, int numNodes, const int* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const int* elementId);

} // namespace finley

#endif // __FINLEY_ASSEMBLE_H__
