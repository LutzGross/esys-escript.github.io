/* $Id$ */

#ifndef INC_FINLEY_RECTANGULARMESH
#define INC_FINLEY_RECTANGULARMESH

/**************************************************************/

/*   Finley: header file for generates rectangular meshes for 1D,2D,3D. */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Hex20(int*,double*,int*,int,int);
Finley_Mesh* Finley_RectangularMesh_Hex8(int*,double*,int*,int,int);
Finley_Mesh* Finley_RectangularMesh_Rec8(int*,double*,int*,int,int);
Finley_Mesh* Finley_RectangularMesh_Rec4(int*,double*,int*,int,int);
Finley_Mesh* Finley_RectangularMesh_Line3(int*,double*,int*,int,int);
Finley_Mesh* Finley_RectangularMesh_Line2(int*,double*,int*,int,int);

#define COLOR_MOD(_I_) ((_I_%2)+((_I_==0) ? 2 : 0))

#endif /* #ifndef INC_FINLEY_RECTANGULARMESH */

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
