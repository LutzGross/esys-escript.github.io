
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

  assemblage routines: copies node coordinates into an expanded Data Object.

*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "Util.h"

#include <sstream>

namespace dudley {

void Assemble_NodeCoordinates(Dudley_NodeFile* nodes, escript::Data* x)
{
    Dudley_resetError();
    if (nodes == NULL)
        return;

    const escript::DataTypes::ShapeType expectedShape(1, nodes->numDim);

    if (!x->numSamplesEqual(1, nodes->numNodes)) {
        Dudley_setError(TYPE_ERROR, "Assemble_NodeCoordinates: illegal number of samples of Data object");
    } else if (x->getFunctionSpace().getTypeCode() != DUDLEY_NODES) {
        Dudley_setError(TYPE_ERROR, "Assemble_NodeCoordinates: Data object is not defined on nodes.");
    } else if (!x->actsExpanded()) {
        Dudley_setError(TYPE_ERROR, "Assemble_NodeCoordinates: expanded Data object expected");
    } else if (x->getDataPointShape() != expectedShape) {
        std::stringstream ss;
        ss << "Assemble_NodeCoordinates: Data object of shape ("
            << nodes->numDim << ",) expected.";
        const std::string errorMsg(ss.str());
        Dudley_setError(TYPE_ERROR, errorMsg.c_str());
    } else {
        const size_t dim_size = nodes->numDim * sizeof(double);
        x->requireWrite();
#pragma omp parallel for
        for (dim_t n = 0; n < nodes->numNodes; n++)
            memcpy(x->getSampleDataRW(n),
                    &nodes->Coordinates[INDEX2(0, n, nodes->numDim)], dim_size);
    }
}

} // namespace dudley

