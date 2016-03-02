
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

/*
  NodeMapping provides a mapping from the local nodes typically to the
  degrees of freedom, the reduced degrees of freedom or the reduced node set.
*/

#ifndef __FINLEY_NODEMAPPING_H__
#define __FINLEY_NODEMAPPING_H__

#include "Util.h"

namespace finley {

struct NodeMapping {
    /// resets both map and target.
    void clear()
    {
        target.clear();
        map.clear();
    }

    /// initializes a node mapping. The target array is copied and a reverse
    /// map created.
    /// theTarget[i]=unused means that no target is defined for FEM node i.
    void assign(const std::vector<index_t>& theTarget, index_t unused)
    {
        if (theTarget.empty())
            return;

        std::pair<index_t,index_t> range(
            util::getFlaggedMinMaxInt(theTarget.size(), &theTarget[0], unused));
        if (range.first < 0) {
            throw escript::ValueError("NodeMapping: target has negative entry.");
        }
        // now we assume min(target)=0!
        const dim_t numTargets = range.first<=range.second ? range.second+1 : 0;
        target.assign(theTarget.begin(), theTarget.end());
        const index_t targetSize = target.size();
        map.assign(numTargets, -1);

        bool err = false;
#pragma omp parallel
        {
#pragma omp for
            for (index_t i=0; i<targetSize; ++i) {
                if (target[i] != unused)
                    map[target[i]]=i;
            }
            // sanity check
#pragma omp for
            for (index_t i=0; i<numTargets; ++i) {
                if (map[i]==-1) {
#pragma omp critical
                    err=true;
                }
            }
        }
        if (err)
            throw escript::ValueError("NodeMapping: target does not define a continuous labeling.");
    }

    /// returns the number of target nodes (number of items in the map array)
    dim_t getNumTargets() const { return map.size(); }

    /// target[i] defines the target of FEM node i=0,...,numNodes-1
    std::vector<index_t> target;
    /// maps the target nodes back to the FEM nodes: target[map[i]]=i
    std::vector<index_t> map;
};

} // namespace finley

#endif // __FINLEY_NODEMAPPING_H__

