/* $Id$ */

#ifndef INC_FINLEY_SHAPEFUNCTIONS
#define INC_FINLEY_SHAPEFUNCTIONS

/**************************************************************/

/*   Finley: Shape functions header file */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"

/**************************************************************/

#define S_INDEX(_J_,_I_,_NUMNODES_)                INDEX2(_J_,_I_,_NUMNODES_)
#define DSDV_INDEX(_J_,_K_,_I_,_NUMNODES_,_DIM_)   INDEX3(_J_,_K_,_I_,_NUMNODES_,_DIM_)

/**************************************************************/

/*   Interfaces: */

typedef void (Finley_Shape_Function) (int,double*,double*,double*);

Finley_Shape_Function Finley_Shape_Point1;
Finley_Shape_Function Finley_Shape_Line2;
Finley_Shape_Function Finley_Shape_Line3;
Finley_Shape_Function Finley_Shape_Line4;
Finley_Shape_Function Finley_Shape_Tri3;
Finley_Shape_Function Finley_Shape_Tri6;
Finley_Shape_Function Finley_Shape_Tri9;
Finley_Shape_Function Finley_Shape_Tri10;
Finley_Shape_Function Finley_Shape_Rec4;
Finley_Shape_Function Finley_Shape_Rec8;
Finley_Shape_Function Finley_Shape_Rec9;
Finley_Shape_Function Finley_Shape_Rec12;
Finley_Shape_Function Finley_Shape_Rec16;
Finley_Shape_Function Finley_Shape_Tet4;
Finley_Shape_Function Finley_Shape_Tet10;
Finley_Shape_Function Finley_Shape_Tet16;
Finley_Shape_Function Finley_Shape_Hex8;
Finley_Shape_Function Finley_Shape_Hex20;
Finley_Shape_Function Finley_Shape_Hex32;

#endif /* #ifndef INC_FINLEY_SHAPEFUNCTIONS */

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
