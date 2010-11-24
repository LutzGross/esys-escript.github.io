
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

#ifndef INC_FINLEY_ASSEMBLE
#define INC_FINLEY_ASSEMBLE

/**************************************************************/

#include "ReferenceElements.h"
#include "Finley.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "escript/DataC.h"
#include "paso/SystemMatrix.h"

struct Finley_Assemble_Parameters {
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
   Finley_ElementFile_Jacobeans* row_jac;
   index_t* row_node;
   dim_t row_numShapesTotal;
   dim_t row_numShapes;
 
   dim_t numComp;
   index_t * col_DOF;
   dim_t col_DOF_UpperBound;
   Finley_ElementFile_Jacobeans* col_jac;
   index_t* col_node;
   dim_t col_numShapesTotal;
   dim_t col_numShapes;
};

typedef struct Finley_Assemble_Parameters Finley_Assemble_Parameters;


#define Finley_Assemble_reducedIntegrationOrder(__in__) ( (getFunctionSpaceType(__in__) == FINLEY_REDUCED_ELEMENTS) || (getFunctionSpaceType(__in__) == FINLEY_REDUCED_FACE_ELEMENTS) || (getFunctionSpaceType(__in__) == FINLEY_REDUCED_CONTACT_ELEMENTS_1) || (getFunctionSpaceType(__in__) == FINLEY_REDUCED_CONTACT_ELEMENTS_2) )

void Finley_Assemble_PDE(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                    escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*) ;

void Finley_Assemble_getAssembleParameters(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,bool_t, Finley_Assemble_Parameters*);
void  Finley_Assemble_PDE_System2_3D(Finley_Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_System2_2D(Finley_Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_System2_1D(Finley_Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_System2_C(Finley_Assemble_Parameters , Finley_ElementFile*, Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_3D(Finley_Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_2D(Finley_Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_1D(Finley_Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_C(Finley_Assemble_Parameters p, Finley_ElementFile*, Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*);


void Finley_Assemble_NodeCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_Assemble_setNormal(Finley_NodeFile*, Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_interpolate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_gradient(Finley_NodeFile*, Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_integrate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*,double*) ;
void Finley_Assemble_getSize(Finley_NodeFile*,Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_AverageElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_addToSystemMatrix(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);

void Finley_Assemble_jacobeans_1D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_2D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_2D_M1D_E2D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_2D_M1D_E2D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_2D_M1D_E1D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_2D_M1D_E1D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_3D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_3D_M2D_E3D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_3D_M2D_E3D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_3D_M2D_E2D(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);
void Finley_Assemble_jacobeans_3D_M2D_E2D_C(double*, dim_t, double*, dim_t, dim_t, dim_t, index_t*, double*, dim_t, double*, double*, double*, index_t*);


void Finley_Assemble_LumpedSystem(Finley_NodeFile* nodes,Finley_ElementFile* elements, escriptDataC* lumpedMat, escriptDataC* D, const bool_t useHRZ);
#endif /* #ifndef INC_FINLEY_ASSEMBLE */

