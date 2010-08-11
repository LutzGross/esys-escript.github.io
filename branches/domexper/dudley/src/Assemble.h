
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assemblage routines: header file */

/**************************************************************/

#ifndef INC_DUDLEY_ASSEMBLE
#define INC_DUDLEY_ASSEMBLE

/**************************************************************/

#include "ReferenceElements.h"
#include "Dudley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "escript/DataC.h"
#include "paso/SystemMatrix.h"

struct Assemble_Parameters {
   dim_t numQuadTotal; /* total number of quadrature nodes = numQuadSub * numQuadSub */
   dim_t numQuadSub; /* number of quadrature nodes per subelements */
   dim_t numSides; /* number of sides */
   dim_t numSub;  /* number of subelements */
   dim_t numDim;  /* spatial dimension */ 
   dim_t NN;     /* leading dimension of element node table */
   dim_t numElements; /* number of elements */

   dim_t numEqu;
   index_t* row_DOF;
   dim_t row_DOF_UpperBound;
   Dudley_ElementFile_Jacobeans* row_jac;
   index_t* row_node;
   dim_t row_numShapesTotal;
   dim_t row_numShapes;
 
   dim_t numComp;
   index_t * col_DOF;
   dim_t col_DOF_UpperBound;
   Dudley_ElementFile_Jacobeans* col_jac;
   index_t* col_node;
   dim_t col_numShapesTotal;
   dim_t col_numShapes;
};

typedef struct Assemble_Parameters Assemble_Parameters;


#define Dudley_Assemble_reducedIntegrationOrder(__in__) ( (getFunctionSpaceType(__in__) == DUDLEY_REDUCED_ELEMENTS) || (getFunctionSpaceType(__in__) == DUDLEY_REDUCED_FACE_ELEMENTS) )

void Dudley_Assemble_PDE(Dudley_NodeFile*,Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                    escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*) ;
void Dudley_Assemble_LumpedSystem(Dudley_NodeFile*,Dudley_ElementFile*, escriptDataC*, escriptDataC*);

void Assemble_getAssembleParameters(Dudley_NodeFile*,Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,bool_t, Assemble_Parameters*);
void  Dudley_Assemble_PDE_System2_3D(Assemble_Parameters, Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_System2_2D(Assemble_Parameters, Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_System2_1D(Assemble_Parameters, Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_System2_C(Assemble_Parameters , Dudley_ElementFile*, Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_Single2_3D(Assemble_Parameters, Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_Single2_2D(Assemble_Parameters, Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_Single2_1D(Assemble_Parameters, Dudley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Dudley_Assemble_PDE_Single2_C(Assemble_Parameters p, Dudley_ElementFile*, Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*);


void Dudley_Assemble_NodeCoordinates(Dudley_NodeFile*,escriptDataC*);
void Dudley_Assemble_setNormal(Dudley_NodeFile*, Dudley_ElementFile*, escriptDataC*);
void Dudley_Assemble_interpolate(Dudley_NodeFile*,Dudley_ElementFile*,escriptDataC*, escriptDataC*);
void Dudley_Assemble_gradient(Dudley_NodeFile*, Dudley_ElementFile*,escriptDataC*, escriptDataC*);
void Dudley_Assemble_integrate(Dudley_NodeFile*,Dudley_ElementFile*,escriptDataC*,double*) ;
void Dudley_Assemble_getSize(Dudley_NodeFile*,Dudley_ElementFile*, escriptDataC*);
void Dudley_Assemble_CopyNodalData(Dudley_NodeFile* nodes,escriptDataC* out,escriptDataC* in);
void Dudley_Assemble_CopyElementData(Dudley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Dudley_Assemble_AverageElementData(Dudley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Dudley_Assemble_addToSystemMatrix(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);

void Assemble_jacobeans_1D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E2D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E2D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E1D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E1D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E3D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E3D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E2D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E2D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);


void Dudley_Assemble_LumpedSystem(Dudley_NodeFile* nodes,Dudley_ElementFile* elements, escriptDataC* lumpedMat, escriptDataC* D);
#endif /* #ifndef INC_DUDLEY_ASSEMBLE */

