
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

/**************************************************************/

/*   Finley: header file for generates rectangular meshes for 1D,2D,3D. */

/**************************************************************/

#ifndef INC_FINLEY_RECTANGULARMESH
#define INC_FINLEY_RECTANGULARMESH

/**************************************************************/

#include "Mesh.h"

Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);
Finley_Mesh* Finley_RectangularMesh_Hex8(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);
Finley_Mesh* Finley_RectangularMesh_Rec8(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);
Finley_Mesh* Finley_RectangularMesh_Rec4(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);

#endif /* #ifndef INC_FINLEY_RECTANGULARMESH */
