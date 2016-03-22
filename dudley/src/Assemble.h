
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
    escript::ASM_ptr S;
    /// right-hand side to be updated
    escript::Data& F;
    /// number of quadrature nodes
    int numQuad;
    /// number of spatial dimensions
    int numDim;
    /// leading dimension of element node table
    int NN;
    /// number of equations (= matrix row block size)
    int numEqu;
    /// number of components (= matrix column block size)
    int numComp;
    /// row and column degrees of freedom
    const index_t* DOF;
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

void Assemble_PDE_Points(const AssembleParameters& p,
                         const escript::Data& d_dirac,
                         const escript::Data& y_dirac);

void Assemble_PDE_Single_2D(const AssembleParameters& p,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Single_3D(const AssembleParameters& p,
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


/// Adds the matrix array[Equa,Sol,NN,NN] onto the matrix S.
/// The rows/columns are given by i_Equa+Equa*Nodes_Equa[Nodes[j_Equa]]
/// (i_Equa=0:Equa; j_Equa=0:NN_Equa).
/// The routine has to be called from a parallel region and assumes that
/// Equa=Sol=1, i.e. array is fully packed.
void Assemble_addToSystemMatrix(escript::ASM_ptr S, dim_t NN_Equa,
                                const index_t* Nodes_Equa, dim_t num_Equa,
                                dim_t NN_Sol, const index_t* Nodes_Sol,
                                dim_t num_Sol, const double* array);

/// Assembles the mass matrix in lumped form.
/// The coefficient D has to be defined on the integration points or not
/// present. `lumpedMat` has to be initialized before the routine is called.
void Assemble_LumpedSystem(const NodeFile* nodes, const ElementFile* elements,
                           escript::Data& lumpedMat, const escript::Data& D,
                           bool useHRZ);

/// averages data
void Assemble_AverageElementData(const ElementFile* elements,
                                 escript::Data& out, const escript::Data& in);

/// copies data between different types of elements
void Assemble_CopyElementData(const ElementFile* elements, escript::Data& out,
                              const escript::Data& in);

/// copies data between different types of nodal representations
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
void Assemble_gradient(const NodeFile* nodes, const ElementFile* elements,
                       escript::Data& gradient, const escript::Data& data);

/// integrates data on quadrature points
void Assemble_integrate(const NodeFile* nodes, const ElementFile* elements,
                   const escript::Data& data, std::vector<double>& integrals);

/// interpolates nodal data in a data array onto elements (=integration points)
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

