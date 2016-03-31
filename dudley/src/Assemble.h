
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
#include "escript/AbstractSystemMatrix.h"

namespace dudley {

struct Assemble_Parameters {
    // number of quadrature nodes
    int numQuad;
    // number of spatial dimensions
    int numDim;
    // leading dimension of element node table
    int NN;
    // number of elements
    dim_t numElements;

    int numEqu;
    const index_t* row_DOF;
    dim_t row_DOF_UpperBound;
    Dudley_ElementFile_Jacobians* row_jac;
    int numShapes;

    int numComp;
    const index_t* col_DOF;
    dim_t col_DOF_UpperBound;

    const double* shapeFns;
};

inline bool Assemble_reducedIntegrationOrder(const escript::Data* in)
{
    const int fs = in->getFunctionSpace().getTypeCode();
    return (fs == DUDLEY_REDUCED_ELEMENTS || fs == DUDLEY_REDUCED_FACE_ELEMENTS);
}

void Assemble_getAssembleParameters(const Dudley_NodeFile*, const Dudley_ElementFile*,
                                    escript::ASM_ptr, const escript::Data&,
                                    bool, Assemble_Parameters*);

void Assemble_PDE(const Dudley_NodeFile* nodes, const Dudley_ElementFile* elements,
                  escript::ASM_ptr S, escript::Data& F,
                  const escript::Data& A, const escript::Data& B,
                  const escript::Data& C, const escript::Data& D,
                  const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Points(const Assemble_Parameters& p, const Dudley_ElementFile*,
                         escript::ASM_ptr S, escript::Data& F,
                         const escript::Data& d_dirac,
                         const escript::Data& y_dirac);

void Assemble_PDE_Single_2D(const Assemble_Parameters& p, const Dudley_ElementFile*,
                            escript::ASM_ptr S, escript::Data& F,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_Single_3D(const Assemble_Parameters& p, const Dudley_ElementFile*,
                            escript::ASM_ptr S, escript::Data& F,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_System_2D(const Assemble_Parameters& p, const Dudley_ElementFile*,
                            escript::ASM_ptr S, escript::Data& F,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_PDE_System_3D(const Assemble_Parameters& p, const Dudley_ElementFile*,
                            escript::ASM_ptr S, escript::Data& F,
                            const escript::Data& A, const escript::Data& B,
                            const escript::Data& C, const escript::Data& D,
                            const escript::Data& X, const escript::Data& Y);

void Assemble_NodeCoordinates(Dudley_NodeFile *, escript::Data *);
void Assemble_setNormal(Dudley_NodeFile *, Dudley_ElementFile *, escript::Data *);
void Assemble_interpolate(Dudley_NodeFile *, Dudley_ElementFile *, const escript::Data *, escript::Data *);
void Assemble_gradient(Dudley_NodeFile *, Dudley_ElementFile *, escript::Data *, const escript::Data *);
void Assemble_integrate(Dudley_NodeFile *, Dudley_ElementFile *, const escript::Data *, double *);
void Assemble_getSize(Dudley_NodeFile *, Dudley_ElementFile *, escript::Data *);
void Assemble_CopyNodalData(Dudley_NodeFile * nodes, escript::Data * out, const escript::Data * in);
void Assemble_CopyElementData(Dudley_ElementFile * elements, escript::Data * out, const escript::Data * in);
void Assemble_AverageElementData(Dudley_ElementFile * elements, escript::Data * out, const escript::Data * in);
void Assemble_addToSystemMatrix(escript::ASM_ptr in, const dim_t NN_Equa, const index_t * Nodes_Equa, const dim_t num_Equa,
                                       const dim_t NN_Sol, const index_t * Nodes_Sol, const dim_t num_Sol, const double *array);

void Assemble_jacobians_2D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D, double *quadweight,
                           index_t *);
void Assemble_jacobians_2D_M1D_E1D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D,
                                   double *quadweight, index_t *);
void Assemble_jacobians_3D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D, double *quadweight,
                           index_t *);
void Assemble_jacobians_3D_M2D_E2D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D,
                                   double *quadweight, index_t *);

void Assemble_LumpedSystem(Dudley_NodeFile* nodes, Dudley_ElementFile* elements,
                           escript::Data& lumpedMat, const escript::Data& D,
                           bool useHRZ);

} // namespace dudley

#endif // __DUDLEY_ASSEMBLE_H__

