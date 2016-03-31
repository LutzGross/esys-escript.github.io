
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

#ifndef INC_DUDLEY_TRIANGULARMESH
#define INC_DUDLEY_TRIANGULARMESH

#include "Mesh.h"

namespace dudley {

Dudley_Mesh *Dudley_TriangularMesh_Tri3(dim_t * numElements, double *Length, index_t order, index_t reduced_order,
                    bool optimize, escript::JMPI& mpi_info);
Dudley_Mesh *Dudley_TriangularMesh_Tet4(dim_t * numElements, double *Length, index_t order, index_t reduced_order,
                    bool optimize, escript::JMPI& mpi_info);

} // namespace dudley

#endif /* #ifndef INC_DUDLEY_TRIANGULARMESH */

