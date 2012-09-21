
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

/*   Dudley: header file for generates triangular meshes for 1D,2D,3D. */

/************************************************************************************/

#ifndef INC_DUDLEY_TRIANGULARMESH
#define INC_DUDLEY_TRIANGLULARMESH

/************************************************************************************/

#include "Mesh.h"

Dudley_Mesh *Dudley_TriangularMesh_Tri3(dim_t * numElements, double *Length, index_t order, index_t reduced_order,
					bool_t optimize);
Dudley_Mesh *Dudley_TriangularMesh_Tet4(dim_t * numElements, double *Length, index_t order, index_t reduced_order,
					bool_t optimize);

#endif				/* #ifndef INC_DUDLEY_TRIANGULARMESH */
