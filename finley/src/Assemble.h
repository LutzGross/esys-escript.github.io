
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

  Assemblage routines: header file

*****************************************************************************/

#ifndef __FINLEY_ASSEMBLE_H__
#define __FINLEY_ASSEMBLE_H__

#include "Finley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include <escript/AbstractSystemMatrix.h>

namespace finley {

struct AssembleParameters
{
    AssembleParameters(const NodeFile* nodes, const ElementFile* ef,
                       escript::ASM_ptr sm, escript::Data& rhs,
                       bool reducedOrder);

    /// element file these parameters apply to
    const ElementFile* elements;
    /// system matrix to be updated
    escript::ASM_ptr S;
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
    dim_t numElements;

    int numEqu;
    const index_t* row_DOF;
    index_t row_DOF_UpperBound;
    ElementFile_Jacobians* row_jac;
    const int* row_node;
    int row_numShapesTotal;
    int row_numShapes;
    int numComp;
    const index_t* col_DOF;
    index_t col_DOF_UpperBound;
    ElementFile_Jacobians* col_jac;
    const int* col_node;
    int col_numShapesTotal;
    int col_numShapes;
};


/// Entry point for PDE assembly. Checks arguments, populates an
/// AssembleParameters structure and calls appropriate method for the actual
/// work.
void Assemble_PDE(const NodeFile* nodes, const ElementFile* elements,
                  escript::ASM_ptr S, escript::Data& F,
                  const escript::Data& A, const escript::Data& B,
                  const escript::Data& C, const escript::Data& D,
                  const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_Points(const AssembleParameters& p,
                         const escript::Data& d_dirac,
                         const escript::Data& y_dirac);

void Assemble_PDE_Single_1D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_Single_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_Single_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_Single_C(const AssembleParameters& p, const escript::Data& D,
                           const escript::Data& Y);

void Assemble_PDE_System_1D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_System_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_System_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar>
void Assemble_PDE_System_C(const AssembleParameters& p, const escript::Data& D,
                           const escript::Data& Y);

template<typename Scalar = double>
void Assemble_addToSystemMatrix(escript::ASM_ptr S, int NN_Equa,
                  const index_t* Nodes_Equa, int num_Equa, int NN_Sol,
                  const index_t* Nodes_Sol, int num_Sol, const Scalar* array);

void Assemble_LumpedSystem(const NodeFile* nodes, const ElementFile* elements,
                           escript::Data& lumpedMat, const escript::Data& D,
                           bool useHRZ);

/// averages data
template<typename Scalar>
void Assemble_AverageElementData(const ElementFile* elements,
                                 escript::Data& out, const escript::Data& in);

/// copies data between different types of elements
template<typename Scalar>
void Assemble_CopyElementData(const ElementFile* elements, escript::Data& out,
                              const escript::Data& in);

/// copies data between different types of nodal representations
template<typename Scalar>
void Assemble_CopyNodalData(const NodeFile* nodes, escript::Data& out,
                            const escript::Data& in);

/// copies node coordinates into expanded Data object `x`
void Assemble_NodeCoordinates(const NodeFile* nodes, escript::Data& x);

/// calculates the normal vector at quadrature points on face elements
void Assemble_getNormal(const NodeFile* nodes, const ElementFile* elements,
                        escript::Data& normals);

/// calculates the minimum distance between two vertices of elements and
/// assigns the value to each quadrature point in `size`
void Assemble_getSize(const NodeFile* nodes, const ElementFile* elements,
                      escript::Data& size);

/// Assemblage of Jacobians: calculates the gradient of nodal data at
/// quadrature points
template<typename Scalar>
void Assemble_gradient(const NodeFile* nodes, const ElementFile* elements,
                       escript::Data& gradient, const escript::Data& data);

/// integrates data on quadrature points
template<typename Scalar>
void Assemble_integrate(const NodeFile* nodes, const ElementFile* elements,
                        const escript::Data& data, Scalar* integrals);

template<typename Scalar>
void Assemble_integrate_points(const ElementFile* points,
                        const escript::Data& data, Scalar* integrals);


/// interpolates nodal data in a data array onto elements (=integration points)
template<typename Scalar>
void Assemble_interpolate(const NodeFile* nodes, const ElementFile* elements,
                          const escript::Data& data, escript::Data& output);

void Assemble_jacobians_1D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_2D_M1D_E1D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_2D_M1D_E1D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_2D_M1D_E2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_2D_M1D_E2D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_3D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_3D_M2D_E2D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_3D_M2D_E2D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_3D_M2D_E3D(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

void Assemble_jacobians_3D_M2D_E3D_C(const double* coordinates, int numQuad,
                           const double* QuadWeights, int numShape,
                           dim_t numElements, int numNodes, const index_t* nodes,
                           const double* DSDv, int numTest, const double* DTDv,
                           double* dTdX, double* volume, const index_t* elementId);

} // namespace finley

#endif // __FINLEY_ASSEMBLE_H__

