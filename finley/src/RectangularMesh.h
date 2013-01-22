
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

/*   Finley: header file for generating rectangular meshes in 1D,2D,3D. */

/************************************************************************************/

#ifndef INC_FINLEY_RECTANGULARMESH
#define INC_FINLEY_RECTANGULARMESH

/************************************************************************************/

#include "Mesh.h"

Finley_Mesh* Finley_RectangularMesh_Hex20(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t, bool_t);
Finley_Mesh* Finley_RectangularMesh_Hex8(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);
Finley_Mesh* Finley_RectangularMesh_Rec8(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t, bool_t);
Finley_Mesh* Finley_RectangularMesh_Rec4(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);

#endif /* #ifndef INC_FINLEY_RECTANGULARMESH */
