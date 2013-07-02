
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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


/****************************************************************************

  Assemblage routines: header file

*****************************************************************************/

#ifndef __FINLEY_ASSEMBLE_H__
#define __FINLEY_ASSEMBLE_H__

#include "ReferenceElements.h"
#include "Finley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "paso/SystemMatrix.h"

namespace finley {

struct AssembleParameters {
    AssembleParameters(const NodeFile* nodes, const ElementFile* ef,
                       Paso_SystemMatrix* sm, escript::Data& rhs,
                       bool reducedOrder);

    /// element file these parameters apply to
    const ElementFile* elements;
    /// system matrix to be updated
    Paso_SystemMatrix* S;
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
void Assemble_PDE(NodeFile* nodes, ElementFile* elements, Paso_SystemMatrix* S,
        escript::Data& F, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D, const escript::Data& X,
        const escript::Data& Y);

void Assemble_PDE_Points(const AssembleParameters& p, escript::Data& d_dirac,
                         escript::Data& y_dirac);

void Assemble_PDE_Single_1D(const AssembleParameters& p, escript::Data& A,
        escript::Data& B, escript::Data& C, escript::Data& D, escript::Data& X,
        escript::Data& Y);

void Assemble_PDE_Single_2D(const AssembleParameters& p, escript::Data& A,
        escript::Data& B, escript::Data& C, escript::Data& D, escript::Data& X,
        escript::Data& Y);

void Assemble_PDE_Single_3D(const AssembleParameters& p, escript::Data& A,
        escript::Data& B, escript::Data& C, escript::Data& D, escript::Data& X,
        escript::Data& Y);

void Assemble_PDE_Single_C(const AssembleParameters& p, escript::Data& D,
                           escript::Data& Y);

void Assemble_PDE_System_1D(const AssembleParameters& p, escript::Data& A,
        escript::Data& B, escript::Data& C, escript::Data& D, escript::Data& X,
        escript::Data& Y);

void Assemble_PDE_System_2D(const AssembleParameters& p, escript::Data& A,
        escript::Data& B, escript::Data& C, escript::Data& D, escript::Data& X,
        escript::Data& Y);

void Assemble_PDE_System_3D(const AssembleParameters& p, escript::Data& A,
        escript::Data& B, escript::Data& C, escript::Data& D, escript::Data& X,
        escript::Data& Y);

void Assemble_PDE_System_C(const AssembleParameters& p, escript::Data& D,
                           escript::Data& Y);

void Assemble_addToSystemMatrix(Paso_SystemMatrix*, const int NN_Equa,
        const int* Nodes_Equa, const int num_Equa, const int NN_Sol,
        const int* Nodes_Sol, const int num_Sol, const double* array);

void Assemble_LumpedSystem(NodeFile* nodes, ElementFile* elements,
        escript::Data& lumpedMat, const escript::Data& D, bool useHRZ);

void Assemble_AverageElementData(ElementFile*, escript::Data&, const escript::Data&);
void Assemble_CopyElementData(ElementFile*, escript::Data&, const escript::Data&);
void Assemble_CopyNodalData(const NodeFile*, escript::Data&, const escript::Data&);
void Assemble_NodeCoordinates(NodeFile*, escript::Data&);

void Assemble_gradient(const NodeFile*, const ElementFile*, escript::Data&, const escript::Data&);
void Assemble_integrate(NodeFile*, ElementFile*, const escript::Data&, double*);
void Assemble_interpolate(const NodeFile*, const ElementFile*, const escript::Data&, escript::Data&);
void Assemble_setNormal(NodeFile*, ElementFile*, escript::Data&);
void Assemble_getSize(NodeFile*, ElementFile*, escript::Data&);

void Assemble_jacobians_1D(double*, int, double*, int, int, int, int*, double*,
                           int, double*, double*, double*, int*);
void Assemble_jacobians_2D(double*, int, double*, int, int, int, int*, double*,
                           int, double*, double*, double*, int*);
void Assemble_jacobians_2D_M1D_E2D(double*, int, double*, int, int, int, int*,
                                double*, int, double*, double*, double*, int*);
void Assemble_jacobians_2D_M1D_E2D_C(double*, int, double*, int, int, int,
                          int*, double*, int, double*, double*, double*, int*);
void Assemble_jacobians_2D_M1D_E1D(double*, int, double*, int, int, int, int*,
                                double*, int, double*, double*, double*, int*);
void Assemble_jacobians_2D_M1D_E1D_C(double*, int, double*, int, int, int,
                          int*, double*, int, double*, double*, double*, int*);
void Assemble_jacobians_3D(double*, int, double*, int, int, int, int*, double*,
                           int, double*, double*, double*, int*);
void Assemble_jacobians_3D_M2D_E3D(double*, int, double*, int, int, int, int*,
                                double*, int, double*, double*, double*, int*);
void Assemble_jacobians_3D_M2D_E3D_C(double*, int, double*, int, int, int,
                          int*, double*, int, double*, double*, double*, int*);
void Assemble_jacobians_3D_M2D_E2D(double*, int, double*, int, int, int, int*,
                                double*, int, double*, double*, double*, int*);
void Assemble_jacobians_3D_M2D_E2D_C(double*, int, double*, int, int, int,
                          int*, double*, int, double*, double*, double*, int*);

} // namespace finley

#endif // __FINLEY_ASSEMBLE_H__

