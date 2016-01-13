
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


/****************************************************************************

  Assemblage routines: copies node coordinates into an expanded Data object.

*****************************************************************************/

#include "Util.h"
#include "Assemble.h"

#include <sstream>

namespace finley {

void Assemble_NodeCoordinates(NodeFile* nodes, escript::Data& x)
{
    Finley_resetError();
    if (!nodes) return;

    const escript::DataTypes::ShapeType expectedShape(1, nodes->numDim);

    if (!x.numSamplesEqual(1, nodes->numNodes)) {
        Finley_setError(TYPE_ERROR, "Assemble_NodeCoordinates: illegal number of samples of Data object");
    } else if (x.getFunctionSpace().getTypeCode() != FINLEY_NODES) {
        Finley_setError(TYPE_ERROR, "Assemble_NodeCoordinates: Data object is not defined on nodes.");
    } else if (!x.actsExpanded()) {
        Finley_setError(TYPE_ERROR, "Assemble_NodeCoordinates: expanded Data object expected");
    } else if (x.getDataPointShape() != expectedShape) {
        std::stringstream ss;
        ss << "Assemble_NodeCoordinates: Data object of shape ("
            << nodes->numDim << ",) expected.";
        std::string errorMsg = ss.str();
        Finley_setError(TYPE_ERROR, errorMsg.c_str());
    } else {
        const size_t dim_size = nodes->numDim*sizeof(double);
        x.requireWrite();
#pragma omp parallel for
        for (int n=0; n<nodes->numNodes; n++)
            memcpy(x.getSampleDataRW(n),
                    &(nodes->Coordinates[INDEX2(0,n,nodes->numDim)]), dim_size);
    }
}

} // namespace finley

