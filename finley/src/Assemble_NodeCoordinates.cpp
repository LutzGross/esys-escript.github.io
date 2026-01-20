
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


/****************************************************************************

  Assemblage routines: copies node coordinates into an expanded Data object.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

#include <sstream>

namespace finley {

void Assemble_NodeCoordinates(const NodeFile* nodes, escript::Data& x)
{
    if (!nodes) return;

    const escript::DataTypes::ShapeType expectedShape(1, nodes->numDim);

    if (!x.numSamplesEqual(1, nodes->getNumNodes())) {
        throw escript::ValueError("Assemble_NodeCoordinates: illegal number of samples of Data object");
    } else if (x.getFunctionSpace().getTypeCode() != FINLEY_NODES) {
        throw escript::ValueError("Assemble_NodeCoordinates: Data object is not defined on nodes.");
    } else if (!x.actsExpanded()) {
        throw escript::ValueError("Assemble_NodeCoordinates: expanded Data object expected");
    } else if (x.getDataPointShape() != expectedShape) {
        std::stringstream ss;
        ss << "Assemble_NodeCoordinates: Data object of shape ("
            << nodes->numDim << ",) expected.";
        throw escript::ValueError(ss.str());
    } else {
        const size_t dim_size = nodes->numDim * sizeof(double);
        x.requireWrite();
#pragma omp parallel for
        for (dim_t n = 0; n < nodes->getNumNodes(); n++)
            memcpy(x.getSampleDataRW(n),
                    &nodes->Coordinates[INDEX2(0, n, nodes->numDim)], dim_size);
    }
}

} // namespace finley

