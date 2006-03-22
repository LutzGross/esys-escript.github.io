/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

/**************************************************************/

/*   Finley: header file for generates rectangular meshes for 1D,2D,3D. */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#ifndef INC_FINLEY_RECTANGULARMESH
#define INC_FINLEY_RECTANGULARMESH

#include "Mesh.h"

/**************************************************************/

Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t*,double*,bool_t*,dim_t,dim_t);
Finley_Mesh* Finley_RectangularMesh_Hex8(dim_t*,double*,bool_t*,dim_t,dim_t);
Finley_Mesh* Finley_RectangularMesh_Rec8(dim_t*,double*,bool_t*,dim_t,dim_t);
Finley_Mesh* Finley_RectangularMesh_Rec4(dim_t*,double*,bool_t*,dim_t,dim_t);
Finley_Mesh* Finley_RectangularMesh_Line3(dim_t*,double*,bool_t*,dim_t,dim_t);
Finley_Mesh* Finley_RectangularMesh_Line2(dim_t*,double*,bool_t*,dim_t,dim_t);

#define COLOR_MOD(_I_) ((_I_%2)+((_I_==0) ? 2 : 0))

#endif /* #ifndef INC_FINLEY_RECTANGULARMESH */

/*
 * $Log$
 * Revision 1.3  2005/09/15 03:44:23  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.2.2.1  2005/09/07 06:26:21  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2  2005/07/08 04:07:56  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:55  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
