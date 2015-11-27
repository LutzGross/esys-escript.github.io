
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************

  Finley: declaration of rectangular mesh generators in 1D, 2D, 3D.

*****************************************************************************/

#ifndef __FINLEY_RECTANGULARMESH_H__
#define __FINLEY_RECTANGULARMESH_H__

#include "Mesh.h"

namespace finley {

Mesh* RectangularMesh_Hex20(const dim_t* numElements, const double* length,
                            const bool* periodic, int order, int reducedOrder,
                            bool useElementsOnFace, bool useFullElementOrder,
                            bool useMacroElements, bool optimize,
                            esysUtils::JMPI& mpi_info);

Mesh* RectangularMesh_Hex8(const dim_t* numElements, const double* length,
                           const bool* periodic, int order, int reducedOrder,
                           bool useElementsOnFace, bool useFullElementOrder,
                           bool useMacroElements, esysUtils::JMPI& mpi_info);

Mesh* RectangularMesh_Rec8(const dim_t* numElements, const double* length,
                           const bool* periodic, int order, int reducedOrder,
                           bool useElementsOnFace, bool useFullElementOrder,
                           bool useMacroElements, bool optimize,
                           esysUtils::JMPI& mpi_info);

Mesh* RectangularMesh_Rec4(const dim_t* numElements, const double* length,
                           const bool* periodic, int order, int reducedOrder,
                           bool useElementsOnFace, bool useFullElementOrder,
                           bool useMacroElements, esysUtils::JMPI& mpi_info);

}

#endif // __FINLEY_RECTANGULARMESH_H__

