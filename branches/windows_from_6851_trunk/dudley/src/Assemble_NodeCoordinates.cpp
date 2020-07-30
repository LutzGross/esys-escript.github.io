
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include "Assemble.h"
#include "Util.h"

#include <escript/index.h>

#include <sstream>

namespace dudley {

void Assemble_NodeCoordinates(const NodeFile* nodes, escript::Data& x)
{
    if (nodes == NULL)
        return;

    if (x.isComplex())
    {
        throw escript::ValueError("Assemble_NodeCoordinates: complex arguments not supported");
    }
    
    const escript::DataTypes::ShapeType expectedShape(1, nodes->numDim);

    if (!x.numSamplesEqual(1, nodes->getNumNodes())) {
        throw escript::ValueError("Assemble_NodeCoordinates: illegal number of samples of Data object");
    } else if (x.getFunctionSpace().getTypeCode() != DUDLEY_NODES) {
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
            memcpy(x.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0)),
                    &nodes->Coordinates[INDEX2(0, n, nodes->numDim)], dim_size);
    }
}

} // namespace dudley

