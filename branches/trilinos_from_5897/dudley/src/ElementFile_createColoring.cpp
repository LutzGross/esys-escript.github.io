
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

#include "ElementFile.h"
#include "Util.h"

namespace dudley {

void ElementFile::createColoring(dim_t nNodes, const index_t* degreeOfFreedom)
{
    if (numElements < 1)
        return;

    //const std::pair<index_t,index_t> idRange(util::getMinMaxInt(
    //                                        1, dofMap.size(), &dofMap[0]));
    const std::pair<index_t,index_t> idRange(util::getMinMaxInt(
                                            1, nNodes, degreeOfFreedom));

    const dim_t NN = numNodes;
    index_t min_id = idRange.first;
    index_t max_id = idRange.second;
    dim_t len = max_id - min_id + 1;
    index_t* maskDOF = new index_t[len];
#pragma omp parallel for
    for (index_t e = 0; e < numElements; e++)
        Color[e] = -1;
    dim_t numUncoloredElements = numElements;
    minColor = 0;
    maxColor = minColor - 1;
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
                    if (Nodes[INDEX2(i, e, NN)] < 0 || Nodes[INDEX2(i, e, NN)] >= nNodes)
                    {
                        printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d min_id=%d Nodes[INDEX2...]=%d\n", __FILE__,
                               __LINE__, i, e, NN, min_id, Nodes[INDEX2(i, e, NN)]);
                        exit(1);
                    }
                    if ((degreeOfFreedom[Nodes[INDEX2(i, e, NN)]] - min_id) >= len
                        || (degreeOfFreedom[Nodes[INDEX2(i, e, NN)]] - min_id) < 0)
                    {
                        printf("BOUNDS_CHECK %s %d i=%d e=%d NN=%d min_id=%d dof=%d\n", __FILE__, __LINE__, i, e,
                               NN, min_id, degreeOfFreedom[Nodes[INDEX2(i, e, NN)]] - min_id);
                        exit(1);
                    }
#endif
                    if (maskDOF[degreeOfFreedom[Nodes[INDEX2(i, e, NN)]] - min_id] > 0)
                    {
                        independent = false;
                        break;
                    }
                }
                // if e is independent a new color is assigned and the nodes
                // are marked as being used
                if (independent) {
                    for (int i = 0; i < NN; i++)
                        maskDOF[degreeOfFreedom[Nodes[INDEX2(i, e, NN)]] - min_id] = 1;
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

