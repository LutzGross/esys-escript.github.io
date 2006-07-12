/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*    assemblage routines: header file */

/**************************************************************/

/*  Copyrights by ACcESS Australia 2003,2004,2005 */
/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

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

struct Assemble_Parameters {
   dim_t numQuad;
   dim_t numDim;
   dim_t numElementDim;

   dim_t NN;
   dim_t NS;
   Finley_RefElement* referenceElement;

   dim_t numEqu;
   index_t* label_row;
   Finley_RefElement* referenceElement_row;
   index_t* row_node;
   dim_t NN_row;
   dim_t NS_row;

   dim_t numComp;
   index_t * label_col;
   Finley_RefElement* referenceElement_col;
   index_t* col_node;
   dim_t NN_col;
   dim_t NS_col;

   /* added by Ben Cumming for MPI version */
   dim_t degreeOfFreedomUpperBound;

   index_t id[MAX_numNodes]; /* used to hold a reordering vector, referenced by row_node and col_node */
};

typedef struct Assemble_Parameters Assemble_Parameters;


typedef void (Finley_Assemble_handelShapeMissMatch) (dim_t, dim_t,dim_t, double*,dim_t, dim_t);

void Finley_Assemble_PDE(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                    escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*) ;
void Finley_Assemble_PDE_RHS(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*,escriptDataC*,escriptDataC*) ;
void Finley_Assemble_RobinCondition(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                    escriptDataC*,escriptDataC*,Finley_Assemble_handelShapeMissMatch) ;
void Finley_Assemble_RobinCondition_RHS(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*,escriptDataC*,Finley_Assemble_handelShapeMissMatch);
/* void Finley_Assemble_Points(Finley_Mesh*,Paso_SystemMatrix*,escriptDataC*,escriptDataC*,escriptDataC*) ;*/
void Finley_Assemble_NodeCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_Assemble_setNormal(Finley_NodeFile*, Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_interpolate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_gradient(Finley_NodeFile*, Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_integrate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*,double*) ;
void Finley_Assemble_getSize(Finley_NodeFile*,Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_PDEMatrix_System2(dim_t,dim_t,dim_t,dim_t,dim_t,double*,double*, double*,dim_t, double*, double*,dim_t, double*,dim_t,double*,dim_t,double*,dim_t);
void Finley_Assemble_PDEMatrix_Single2(dim_t,dim_t,dim_t,double*,double*, double*,dim_t, double*, double*,dim_t, double*,dim_t,double*,dim_t,double*,dim_t);
void Finley_Assemble_RHSMatrix_System(dim_t,dim_t,dim_t,dim_t,double*,double*,double*,dim_t, double*,double*,dim_t,double*,dim_t);
void Finley_Assemble_RHSMatrix_Single(dim_t,dim_t,dim_t,double*,double*,double*,dim_t, double*,double*,dim_t,double*,dim_t);


void Assemble_getAssembleParameters(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,Assemble_Parameters*);
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Step_out;
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Step_in;
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Mean_out;
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Mean_in;
void Finley_Assemble_addToSystemMatrix(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);

void Assemble_jacobeans_1D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);
void Assemble_jacobeans_2D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E2D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E1D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);
void Assemble_jacobeans_3D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E3D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E2D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t,double*, double*, double*, index_t*);

#endif /* #ifndef INC_FINLEY_ASSEMBLE */

/*
 * $Log$
 * Revision 1.4  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.3  2005/08/12 01:45:42  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.2.2.2  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2.2.1  2005/08/04 22:41:11  gross
 * some extra routines for finley that might speed-up RHS assembling in some cases (not actived right now)
 *
 * Revision 1.2  2005/07/08 04:07:45  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:46  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:56  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
