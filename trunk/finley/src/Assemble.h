
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

#ifndef INC_FINLEY_ASSEMBLE
#define INC_FINLEY_ASSEMBLE

#include "ReferenceElements.h"
#include "Finley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "escript/DataC.h"
#include "paso/SystemMatrix.h"

using finley::NodeFile;
using finley::ElementFile;
using finley::ElementFile_Jacobians;

struct Finley_Assemble_Parameters {
    /// total number of quadrature nodes = numQuadSub * numQuadSub
    int numQuadTotal;
    /// number of quadrature nodes per subelements
    int numQuadSub;
    /// number of sides
    int numSides;
    /// number of sub-elements
    int numSub;
    /// spatial dimension
    int numDim;
    /// leading dimension of element node table
    int NN;
    /// number of elements
    int numElements;

    int numEqu;
    int* row_DOF;
    int row_DOF_UpperBound;
    ElementFile_Jacobians* row_jac;
    int* row_node;
    int row_numShapesTotal;
    int row_numShapes;
    int numComp;
    int* col_DOF;
    int col_DOF_UpperBound;
    ElementFile_Jacobians* col_jac;
    int* col_node;
    int col_numShapesTotal;
    int col_numShapes;
};

#define Finley_Assemble_reducedIntegrationOrder(__in__) ( (getFunctionSpaceType(__in__) == FINLEY_REDUCED_ELEMENTS) || (getFunctionSpaceType(__in__) == FINLEY_REDUCED_FACE_ELEMENTS) || (getFunctionSpaceType(__in__) == FINLEY_REDUCED_CONTACT_ELEMENTS_1) || (getFunctionSpaceType(__in__) == FINLEY_REDUCED_CONTACT_ELEMENTS_2) )

void Finley_Assemble_PDE(NodeFile*, ElementFile*,
        Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);

void Finley_Assemble_getAssembleParameters(NodeFile*,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*, bool_t,
        Finley_Assemble_Parameters*);

void Finley_Assemble_PDE_System2_3D(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_System2_2D(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_System2_1D(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_System2_C(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_Single2_3D(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_Single2_2D(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_Single2_1D(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_Single2_C(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_PDE_Points(Finley_Assemble_Parameters,
        ElementFile*, Paso_SystemMatrix*, escriptDataC*,
        escriptDataC*, escriptDataC*);

void Finley_Assemble_NodeCoordinates(NodeFile*, escriptDataC*);

void Finley_Assemble_setNormal(NodeFile*, ElementFile*, escriptDataC*);
void Finley_Assemble_interpolate(NodeFile*, ElementFile*, escriptDataC*, escriptDataC*);
void Finley_Assemble_gradient(NodeFile*, ElementFile*, escriptDataC*, escriptDataC*);
void Finley_Assemble_integrate(NodeFile*, ElementFile*, escriptDataC*, double*);
void Finley_Assemble_getSize(NodeFile*, ElementFile*, escriptDataC*);
void Finley_Assemble_CopyNodalData(NodeFile*, escriptDataC*, escriptDataC*);
void Finley_Assemble_CopyElementData(ElementFile*, escriptDataC*, escriptDataC*);
void Finley_Assemble_AverageElementData(ElementFile*, escriptDataC*, escriptDataC*);

void Finley_Assemble_addToSystemMatrix(Paso_SystemMatrix*, const int NN_Equa,
        const int* Nodes_Equa, const int num_Equa, const int NN_Sol,
        const int* Nodes_Sol, const int num_Sol, const double* array);

void Finley_Assemble_LumpedSystem(NodeFile*, ElementFile*,
        escriptDataC* lumpedMat, escriptDataC* D, const bool_t useHRZ);

void Finley_Assemble_jacobians_1D(double*, int, double*, int, int, int,
        int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_2D(double*, int, double*, int, int, int,
        int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_2D_M1D_E2D(double*, int, double*, int, int,
        int, int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_2D_M1D_E2D_C(double*, int, double*, int,
        int, int, int*, double*, int, double*, double*, double*,
        int*);
void Finley_Assemble_jacobians_2D_M1D_E1D(double*, int, double*, int, int,
        int, int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_2D_M1D_E1D_C(double*, int, double*, int,
        int, int, int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_3D(double*, int, double*, int, int, int,
        int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_3D_M2D_E3D(double*, int, double*, int, int,
        int, int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_3D_M2D_E3D_C(double*, int, double*, int,
        int, int, int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_3D_M2D_E2D(double*, int, double*, int, int,
        int, int*, double*, int, double*, double*, double*, int*);
void Finley_Assemble_jacobians_3D_M2D_E2D_C(double*, int, double*, int,
        int, int, int*, double*, int, double*, double*, double*, int*);

#endif // #ifndef INC_FINLEY_ASSEMBLE

