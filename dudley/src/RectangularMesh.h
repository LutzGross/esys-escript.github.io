
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

/*   Dudley: header file for generates rectangular meshes for 1D,2D,3D. */

/**************************************************************/

#ifndef INC_DUDLEY_RECTANGULARMESH
#define INC_DUDLEY_RECTANGULARMESH

/**************************************************************/

#include "Mesh.h"

Dudley_Mesh* Dudley_RectangularMesh_Hex20(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t, bool_t);
Dudley_Mesh* Dudley_RectangularMesh_Hex8(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);
Dudley_Mesh* Dudley_RectangularMesh_Rec8(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t, bool_t);
Dudley_Mesh* Dudley_RectangularMesh_Rec4(dim_t*,double*,bool_t*,index_t,index_t,bool_t, bool_t, bool_t);

#endif /* #ifndef INC_DUDLEY_RECTANGULARMESH */
