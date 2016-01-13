
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

/************************************************************************************/

/*    assemblage routines: header file */

/************************************************************************************/

#ifndef INC_DUDLEY_ASSEMBLE
#define INC_DUDLEY_ASSEMBLE

/************************************************************************************/

#include "Dudley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "escript/DataC.h"
#include "paso/SystemMatrix.h"

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

void Dudley_Assemble_PDE(Dudley_NodeFile *, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
			 escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
			 escriptDataC *);


void Dudley_Assemble_getAssembleParameters(Dudley_NodeFile *, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    bool_t, Dudley_Assemble_Parameters *);
void Dudley_Assemble_PDE_System2_3D(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
				    escriptDataC *);
void Dudley_Assemble_PDE_System2_2D(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
				    escriptDataC *);
void Dudley_Assemble_PDE_System2_1D(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
				    escriptDataC *);

void Dudley_Assemble_PDE_Single2_3D(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
				    escriptDataC *);
void Dudley_Assemble_PDE_Single2_2D(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
				    escriptDataC *);
void Dudley_Assemble_PDE_Single2_1D(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *,
				    escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *, escriptDataC *,
				    escriptDataC *);
void Dudley_Assemble_PDE_Points(Dudley_Assemble_Parameters, Dudley_ElementFile *, Paso_SystemMatrix *, escriptDataC *, escriptDataC *, escriptDataC *);

void Dudley_Assemble_NodeCoordinates(Dudley_NodeFile *, escriptDataC *);
void Dudley_Assemble_setNormal(Dudley_NodeFile *, Dudley_ElementFile *, escriptDataC *);
void Dudley_Assemble_interpolate(Dudley_NodeFile *, Dudley_ElementFile *, escriptDataC *, escriptDataC *);
void Dudley_Assemble_gradient(Dudley_NodeFile *, Dudley_ElementFile *, escriptDataC *, escriptDataC *);
void Dudley_Assemble_integrate(Dudley_NodeFile *, Dudley_ElementFile *, escriptDataC *, double *);
void Dudley_Assemble_getSize(Dudley_NodeFile *, Dudley_ElementFile *, escriptDataC *);
void Dudley_Assemble_CopyNodalData(Dudley_NodeFile * nodes, escriptDataC * out, escriptDataC * in);
void Dudley_Assemble_CopyElementData(Dudley_ElementFile * elements, escriptDataC * out, escriptDataC * in);
void Dudley_Assemble_AverageElementData(Dudley_ElementFile * elements, escriptDataC * out, escriptDataC * in);
void Dudley_Assemble_addToSystemMatrix(Paso_SystemMatrix * in, const dim_t NN_Equa, const index_t * Nodes_Equa, const dim_t num_Equa,
				       const dim_t NN_Sol, const index_t * Nodes_Sol, const dim_t num_Sol, const double *array);

void Dudley_Assemble_jacobeans_2D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D, double *quadweight,
			   index_t *);
void Dudley_Assemble_jacobeans_2D_M1D_E1D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D,
				   double *quadweight, index_t *);
void Dudley_Assemble_jacobeans_3D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D, double *quadweight,
			   index_t *);
void Dudley_Assemble_jacobeans_3D_M2D_E2D(double *, dim_t, dim_t, dim_t, index_t *, double *, double *abs_D,
				   double *quadweight, index_t *);

void Dudley_Assemble_LumpedSystem(Dudley_NodeFile * nodes, Dudley_ElementFile * elements, escriptDataC * lumpedMat,
				  escriptDataC * D, const bool_t useHRZ);
#endif				/* #ifndef INC_DUDLEY_ASSEMBLE */
