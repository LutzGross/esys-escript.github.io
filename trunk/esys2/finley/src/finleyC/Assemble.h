/* $Id$ */

#ifndef INC_FINLEY_ASSEMBLE
#define INC_FINLEY_ASSEMBLE

/**************************************************************/

/*    assemblage routines: header file */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/
#include "ReferenceElements.h"
#include "System.h"
#include "ElementFile.h"
#include "NodeFile.h"
#include "escript/Data/DataC.h"

struct Assemble_Parameters {
   int numQuad;
   int numDim;
   int numElementDim;

   int NN;
   int NS;
   Finley_RefElement* referenceElement;

   int numEqu;
   maybelong* label_row;
   Finley_RefElement* referenceElement_row;
   maybelong* row_node;
   int NN_row;
   int NS_row;

   int numComp;
   maybelong * label_col;
   Finley_RefElement* referenceElement_col;
   maybelong* col_node;
   int NN_col;
   int NS_col;

   maybelong id[MAX_numNodes]; /* used to hold a reordering vector, referenced by row_node and col_node */
};

typedef struct Assemble_Parameters Assemble_Parameters;


typedef void (Finley_Assemble_handelShapeMissMatch) (int, int,int, double*,int, int);

void Finley_Assemble_PDE(Finley_NodeFile*,Finley_ElementFile*,Finley_SystemMatrix*,escriptDataC*,
                                    escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*, escriptDataC*) ;
void Finley_Assemble_RobinCondition(Finley_NodeFile*,Finley_ElementFile*,Finley_SystemMatrix*,escriptDataC*,
                                    escriptDataC*,escriptDataC*,Finley_Assemble_handelShapeMissMatch) ;
/* void Finley_Assemble_Points(Finley_Mesh*,Finley_SystemMatrix*,escriptDataC*,escriptDataC*,escriptDataC*) ;*/
void Finley_Assemble_NodeCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_Assemble_setNormal(Finley_NodeFile*, Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_interpolate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_gradient(Finley_NodeFile*, Finley_ElementFile*,escriptDataC*, escriptDataC*);
void Finley_Assemble_integrate(Finley_NodeFile*,Finley_ElementFile*,escriptDataC*,double*) ;
void Finley_Assemble_getSize(Finley_NodeFile*,Finley_ElementFile*, escriptDataC*);
void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_CopyElementData(Finley_ElementFile* elements,escriptDataC* out,escriptDataC* in);
void Finley_Assemble_PDEMatrix_System2(int,int,int,int,int,double*,double*, double*,int, double*, double*,int, double*,int,double*,int,double*,int);
void Finley_Assemble_PDEMatrix_Single2(int,int,int,double*,double*, double*,int, double*, double*,int, double*,int,double*,int,double*,int);
void Finley_Assemble_RHSMatrix_System(int,int,int,int,double*,double*,double*,int, double*,double*,int,double*,int);
void Finley_Assemble_RHSMatrix_Single(int,int,int,double*,double*,double*,int, double*,double*,int,double*,int);


void Assemble_getAssembleParameters(Finley_NodeFile*,Finley_ElementFile*,Finley_SystemMatrix*,escriptDataC*,Assemble_Parameters*);
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Step_out;
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Step_in;
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Mean_out;
Finley_Assemble_handelShapeMissMatch Finley_Assemble_handelShapeMissMatch_Mean_in;

#endif /* #ifndef INC_FINLEY_ASSEMBLE */

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:56  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
