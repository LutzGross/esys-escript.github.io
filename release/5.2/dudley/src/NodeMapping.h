
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

#ifndef __DUDLEY_NODEMAPPING_H__
#define __DUDLEY_NODEMAPPING_H__

#include "Util.h"

namespace dudley {

/// NodeMapping provides a mapping from the local nodes typically to the
/// degrees of freedom, the reduced degrees of freedom or the reduced node set
struct NodeMapping
{
    NodeMapping() : numNodes(0), target(NULL), numTargets(0), map(NULL) {}

    /// resets both map and target
    void clear()
    {
        delete[] map;
        delete[] target;
        target = NULL;
        map = NULL;
        numNodes = 0;
        numTargets = 0;
    }

    /// initializes a node mapping. The target array is copied and a reverse
    /// map created.
    /// theTarget[i]=unused means that no target is defined for FEM node i.
    void assign(const index_t* theTarget, dim_t nNodes, index_t unused)
    {
        clear();

        if (nNodes == 0)
            return;

        numNodes = nNodes;

        std::pair<index_t,index_t> range(
            util::getFlaggedMinMaxInt(numNodes, theTarget, unused));
        if (range.first < 0) {
            throw escript::ValueError("NodeMapping: target has negative entry.");
        }
        numTargets = range.first<=range.second ? range.second+1 : 0;

        target = new index_t[numNodes];
        map = new index_t[numTargets];

        bool err = false;
#pragma omp parallel
        {
#pragma omp for
            for (index_t i=0; i<numNodes; ++i) {
                target[i] = theTarget[i];
                if (target[i] != unused)
                    map[target[i]] = i;
            }
            // sanity check
#pragma omp for
            for (index_t i=0; i<numTargets; ++i) {
                if (map[i] == -1) {
#pragma omp critical
                    err = true;
                }
            }
        }
        if (err)
            throw escript::ValueError("NodeMapping: target does not define a continuous labeling.");
    }

    /// returns the number of target nodes (number of items in the map array)
    inline dim_t getNumTargets() const { return numTargets; }

    /// size of `target` (number of FEM nodes)
    dim_t numNodes;

    /// target[i] defines the target of FEM node i=0,...,numNodes
    index_t* target;

    /// size of `map` (number of target nodes, e.g. DOF, reduced DOF, etc.)
    dim_t numTargets;

    /// maps the target nodes back to the FEM nodes: target[map[i]]=i
    index_t* map;
};

} // namespace dudley

#endif // __DUDLEY_NODEMAPPING_H__

