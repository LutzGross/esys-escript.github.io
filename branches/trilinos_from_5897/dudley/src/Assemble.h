
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

/************************************************************************************/

/*    assemblage routines: header file */

/************************************************************************************/

#ifndef INC_DUDLEY_ASSEMBLE
#define INC_DUDLEY_ASSEMBLE

/************************************************************************************/

#include "Dudley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "escript/AbstractSystemMatrix.h"
#include "escript/Data.h"

struct Dudley_Assemble_Parameters {
    dim_t numQuad;		/* number of quadrature nodes */
    dim_t numDim;		/* spatial dimension */
    dim_t NN;			/* leading dimension of element node table */
    dim_t numElements;		/* number of elements */

    dim_t numEqu;
    index_t *row_DOF;
    dim_t row_DOF_UpperBound;
    Dudley_ElementFile_Jacobeans *row_jac;
    dim_t numShapes;

    dim_t numComp;
    index_t *col_DOF;
    dim_t col_DOF_UpperBound;

    const double *shapeFns;
};

typedef struct Dudley_Assemble_Parameters Dudley_Assemble_Parameters;

#define Dudley_Assemble_reducedIntegrationOrder(__in__) ( (getFunctionSpaceType(__in__) == DUDLEY_REDUCED_ELEMENTS) || (getFunctionSpaceType(__in__) == DUDLEY_REDUCED_FACE_ELEMENTS) )

void Dudley_Assemble_PDE(Dudley_NodeFile*, Dudley_ElementFile*, escript::ASM_ptr, escript::Data*,
			 const escript::Data*, const escript::Data*, const escript::Data*, const escript::Data*, const escript::Data*,
			 const escript::Data*);

void Dudley_Assemble_getAssembleParameters(Dudley_NodeFile*, Dudley_ElementFile*, escript::ASM_ptr, const escript::Data*,
				    bool, Dudley_Assemble_Parameters*);
void Dudley_Assemble_PDE_System2_3D(Dudley_Assemble_Parameters, Dudley_ElementFile*, escript::ASM_ptr, escript::Data*,
				    const escript::Data*, const escript::Data*, const escript::Data*, const escript::Data*, const escript::Data*,
				    const escript::Data *);
void Dudley_Assemble_PDE_System2_2D(Dudley_Assemble_Parameters, Dudley_ElementFile*, escript::ASM_ptr, escript::Data*,
				    const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data*,
				    const escript::Data*);
void Dudley_Assemble_PDE_System2_1D(Dudley_Assemble_Parameters, Dudley_ElementFile *, escript::ASM_ptr, const escript::Data *,
				    escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *,
				    const escript::Data *);

void Dudley_Assemble_PDE_Single2_3D(Dudley_Assemble_Parameters, Dudley_ElementFile *, escript::ASM_ptr, escript::Data *,
				    const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *,
				    const escript::Data *);
void Dudley_Assemble_PDE_Single2_2D(Dudley_Assemble_Parameters, Dudley_ElementFile *, escript::ASM_ptr, escript::Data *,
				    const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *,
				    const escript::Data *);
void Dudley_Assemble_PDE_Single2_1D(Dudley_Assemble_Parameters, Dudley_ElementFile *, escript::ASM_ptr, const escript::Data *,
				    const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *, const escript::Data *,
				    const escript::Data *);
void Dudley_Assemble_PDE_Points(Dudley_Assemble_Parameters, Dudley_ElementFile *, escript::ASM_ptr, escript::Data *, const escript::Data *, const escript::Data *);

void Dudley_Assemble_NodeCoordinates(Dudley_NodeFile *, escript::Data *);
void Dudley_Assemble_setNormal(Dudley_NodeFile *, Dudley_ElementFile *, escript::Data *);
void Dudley_Assemble_interpolate(Dudley_NodeFile *, Dudley_ElementFile *, const escript::Data *, escript::Data *);
void Dudley_Assemble_gradient(Dudley_NodeFile *, Dudley_ElementFile *, escript::Data *, const escript::Data *);
void Dudley_Assemble_integrate(Dudley_NodeFile *, Dudley_ElementFile *, const escript::Data *, double *);
void Dudley_Assemble_getSize(Dudley_NodeFile *, Dudley_ElementFile *, escript::Data *);
void Dudley_Assemble_CopyNodalData(Dudley_NodeFile * nodes, escript::Data * out, const escript::Data * in);
void Dudley_Assemble_CopyElementData(Dudley_ElementFile * elements, escript::Data * out, const escript::Data * in);
void Dudley_Assemble_AverageElementData(Dudley_ElementFile * elements, escript::Data * out, const escript::Data * in);
void Dudley_Assemble_addToSystemMatrix(escript::ASM_ptr in, const dim_t NN_Equa, const index_t * Nodes_Equa, const dim_t num_Equa,
				       const dim_t NN_Sol, const index_t * Nodes_Sol, const dim_t num_Sol, const double *array);

void Dudley_Assemble_jacobeans_2D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D, double *quadweight,
			   index_t *);
void Dudley_Assemble_jacobeans_2D_M1D_E1D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D,
				   double *quadweight, index_t *);
void Dudley_Assemble_jacobeans_3D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D, double *quadweight,
			   index_t *);
void Dudley_Assemble_jacobeans_3D_M2D_E2D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D,
				   double *quadweight, index_t *);

void Dudley_Assemble_LumpedSystem(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, escript::Data * lumpedMat,
				  const escript::Data * D, const bool useHRZ);
#endif				/* #ifndef INC_DUDLEY_ASSEMBLE */
