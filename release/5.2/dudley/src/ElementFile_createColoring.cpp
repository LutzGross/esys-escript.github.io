
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

#include "ElementFile.h"
#include "Util.h"

#include <escript/index.h>

namespace dudley {

void ElementFile::createColoring(dim_t nNodes, const index_t* dofMap)
{
    if (numElements < 1)
        return;

    //const std::pair<index_t,index_t> idRange(util::getMinMaxInt(
    //                                        1, dofMap.size(), &dofMap[0]));
    const std::pair<index_t,index_t> idRange(util::getMinMaxInt(
                                            1, nNodes, dofMap));

    const int NN = numNodes;
    const dim_t len = idRange.second - idRange.first + 1;

#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++)
        Color[e] = -1;

    dim_t numUncoloredElements = numElements;
    minColor = 0;
    maxColor = -1;
    index_t* maskDOF = new index_t[len];
    while (numUncoloredElements > 0) {
        // initialize the mask marking nodes used by a color
#pragma omp parallel for
        for (index_t n = 0; n < len; n++)
            maskDOF[n] = -1;
        numUncoloredElements = 0;

        for (index_t e = 0; e < numElements; e++) {
            if (Color[e] < 0) {
                // find out if element e is independent from the elements
                // already colored:
                bool independent = true;
                for (int i = 0; i < NN; i++) {
#ifdef BOUNDS_CHECK
                    ESYS_ASSERT(Nodes[INDEX2(i, e, NN)] >= 0, "BOUNDS_CHECK");
                    ESYS_ASSERT(Nodes[INDEX2(i, e, NN)] < nNodes, "BOUNDS_CHECK");
                    ESYS_ASSERT(dofMap[Nodes[INDEX2(i, e, NN)]] - idRange.first < len, "BOUNDS_CHECK");
                    ESYS_ASSERT(dofMap[Nodes[INDEX2(i, e, NN)]] - idRange.first >= 0, "BOUNDS_CHECK");
#endif
                    if (maskDOF[dofMap[Nodes[INDEX2(i, e, NN)]] - idRange.first] > 0)
                    {
                        independent = false;
                        break;
                    }
                }
                // if e is independent a new color is assigned and the nodes
                // are marked as being used
                if (independent) {
                    for (int i = 0; i < NN; i++)
                        maskDOF[dofMap[Nodes[INDEX2(i, e, NN)]] - idRange.first] = 1;
                    Color[e] = maxColor + 1;
                } else {
                    numUncoloredElements++;
                }
            }
        } // for all elements
        maxColor++;
    } // end of while loop
    delete[] maskDOF;
}

} // namespace dudley

