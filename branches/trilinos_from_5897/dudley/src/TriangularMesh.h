
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __DUDLEY_TRIANGULARMESH_H__
#define __DUDLEY_TRIANGULARMESH_H__

#include "Mesh.h"

namespace dudley {

/// Generates a 2-dimensional numElements[0] x numElements[1] x 2 mesh with
/// Tri3 elements in the rectangle [0,length[0]] x [0,length[1]].
/// The doubling of elements is due to splitting of rectangular elements.
Mesh* TriangularMesh_Tri3(const dim_t* numElements, const double* length,
                          bool optimize, escript::JMPI mpiInfo);

/// Generates a 3-dimensional numElements[0] x numElements[1] x numElements[2]
/// mesh with Tet4 elements in the brick
/// [0,length[0]] x [0,length[1]] x [0,length[2]].
Mesh* TriangularMesh_Tet4(const dim_t* numElements, const double* length,
                          bool optimize, escript::JMPI mpiInfo);

} // namespace dudley

#endif // __DUDLEY_TRIANGULARMESH_H__

