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
   dim_t NN;

   dim_t numEqu;
   index_t* row_DOF;
   dim_t row_DOF_UpperBound;
   Finley_ElementFile_Jacobeans* row_jac;
   index_t* row_node;
   dim_t row_NN;
   dim_t row_NS;

   dim_t numComp;
   index_t * col_DOF;
   dim_t col_DOF_UpperBound;
   Finley_ElementFile_Jacobeans* col_jac;
   index_t* col_node;
   dim_t col_NN;
   dim_t col_NS;

   /* added by Ben Cumming for MPI version */

   index_t id[MAX_numNodes]; /* used to hold a reordering vector, referenced by row_node and col_node */
};

typedef struct Assemble_Parameters Assemble_Parameters;


void Finley_Assemble_PDE(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                    escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*) ;
void Assemble_getAssembleParameters(Finley_NodeFile*,Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,bool_t, Assemble_Parameters*);
void  Finley_Assemble_PDE_System2_3D(Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_System2_2D(Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_System2_1D(Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_System2_C(Assemble_Parameters p, Finley_ElementFile*, Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_3D(Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_2D(Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_1D(Assemble_Parameters, Finley_ElementFile*,Paso_SystemMatrix*,escriptDataC*,
                                     escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*);
void  Finley_Assemble_PDE_Single2_C(Assemble_Parameters p, Finley_ElementFile*, Paso_SystemMatrix*, escriptDataC*, escriptDataC*, escriptDataC*);


void Finley_Assemble_NodeCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_Assemble_setNormal(Finley_NodeFile*, Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_interpolate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_gradient(Finley_NodeFile*, Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_integrate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*,double*) ;
void Finley_Assemble_getSize(Finley_NodeFile*,Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_addToSystemMatrix(Paso_SystemMatrix*,dim_t,index_t*, dim_t,dim_t,index_t*,dim_t, double*);
void Assemble_jacobeans_1D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E2D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E2D_C(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E1D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_2D_M1D_E1D_C(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E3D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E3D_C(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E2D(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);
void Assemble_jacobeans_3D_M2D_E2D_C(double*, dim_t,double*, dim_t, dim_t, dim_t, index_t* , double*, dim_t, double*, double*, double*, index_t*);

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
