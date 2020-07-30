
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

  Assemblage routines: header file

*****************************************************************************/

#ifndef __DUDLEY_ASSEMBLE_H__
#define __DUDLEY_ASSEMBLE_H__

#include "Dudley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include <escript/AbstractSystemMatrix.h>

namespace dudley {

struct AssembleParameters
{
    AssembleParameters(const NodeFile* nodes, const ElementFile* ef,
                       escript::ASM_ptr sm, escript::Data& rhs,
                       bool reducedOrder);

    /// element file these parameters apply to
    const ElementFile* elements;
    /// system matrix to be updated
    escript::AbstractSystemMatrix* S;
    /// right-hand side to be updated
    escript::Data& F;
    /// number of quadrature nodes
    int numQuad;
    /// number of spatial dimensions
    int numDim;
    /// leading dimension of element node table
    int NN;
    /// number of equations (= matrix row/column block size)
    int numEqu;
    /// row and column degrees of freedom
    const index_t* DOF;
    /// number of local degrees of freedom
    dim_t DOF_UpperBound;
    /// reference to jacobians
    const ElementFile_Jacobians* jac;
    int numShapes;
    const double* shapeFns;
};

void Assemble_PDE(const NodeFile* nodes, const ElementFile* elements,
                  escript::ASM_ptr S, escript::Data& F,
                  const escript::Data& A, const escript::Data& B,
                  const escript::Data& C, const escript::Data& D,
                  const escript::Data& X, const escript::Data& Y);

template<typename Scalar = double>
void Assemble_PDE_Points(const AssembleParameters& p,
                         const escript::Data& d_dirac,
                         const escript::Data& y_dirac);

template<typename Scalar = double>
void Assemble_PDE_Single_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar = double>
void Assemble_PDE_Single_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar = double>
void Assemble_PDE_System_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

template<typename Scalar = double>
void Assemble_PDE_System_3D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);


/// Adds the matrix array[Eq,Eq,NN,NN] onto the matrix S.
/// The rows/columns are given by i_Eq+Eq*Nodes[Nodes[j_Eq]]
/// (i_Eq=0:Eq; j_Eq=0:NN_Eq).
/// The routine has to be called from a parallel region and assumes that
/// array is fully packed.
template<typename Scalar>
void Assemble_addToSystemMatrix(escript::AbstractSystemMatrix* S,
                                const std::vector<index_t>& Nodes, int numEq,
                                const std::vector<Scalar>& array);

/// Assembles the mass matrix in lumped form.
/// The coefficient D has to be defined on the integration points or not
/// present. `lumpedMat` has to be initialized before the routine is called.
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
                   const escript::Data& data, std::vector<Scalar>& integrals);

/// interpolates nodal data in a data array onto elements (=integration points)
template<typename Scalar>
void Assemble_interpolate(const NodeFile* nodes, const ElementFile* elements,
                          const escript::Data& data, escript::Data& output);

void Assemble_jacobians_2D(const double* coordinates, int numQuad,
                           dim_t numElements, int numNodes,
                           const index_t* nodes, double* dTdX, double* absD,
                           double* quadWeight, const index_t* elementId);

void Assemble_jacobians_2D_M1D_E1D(const double* coordinates, int numQuad,
                           dim_t numElements, int numNodes,
                           const index_t* nodes, double* dTdX, double* absD,
                           double* quadWeight, const index_t* elementId);

void Assemble_jacobians_3D(const double* coordinates, int numQuad,
                           dim_t numElements, int numNodes,
                           const index_t* nodes, double* dTdX, double* abs_D,
                           double* quadWeight, const index_t* elementId);

void Assemble_jacobians_3D_M2D_E2D(const double* coordinates, int numQuad,
                           dim_t numElements, int numNodes,
                           const index_t* nodes, double* dTdX, double* absD,
                           double* quadWeight, const index_t* elementId);


} // namespace dudley

#endif // __DUDLEY_ASSEMBLE_H__

