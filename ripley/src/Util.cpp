
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/Util.h>
extern "C" {
#include <esysUtils/index.h>
}

using namespace std;

namespace ripley {

IndexPair getFlaggedMinMax(const IndexVector &values, index_t ignore)
{
    index_t vmin=INDEX_T_MAX, vmax=-INDEX_T_MAX;
    if (values.size() > 0) {
#pragma omp parallel
        {
            index_t min_local = vmin;
            index_t max_local = vmax;
#pragma omp for schedule(static) nowait
            for (size_t i=0; i<values.size(); i++) {
                if (values[i] != ignore) {
                    min_local = min(min_local, values[i]);
                    max_local = max(max_local, values[i]);
                }
            }
#pragma omp critical
            vmin = min(vmin, min_local);
            vmax = max(vmax, max_local);
        }
    }
    return IndexPair(vmin,vmax);
}

IndexPair getMinMax(const IndexVector &values)
{
    index_t vmin=INDEX_T_MAX, vmax=-INDEX_T_MAX;
    if (values.size() > 0) {
#pragma omp parallel
        {
            index_t min_local = vmin;
            index_t max_local = vmax;
#pragma omp for schedule(static) nowait
            for (size_t i=0; i<values.size(); i++) {
                min_local = min(min_local, values[i]);
                max_local = max(max_local, values[i]);
            }
#pragma omp critical
            vmin = min(vmin, min_local);
            vmax = max(vmax, max_local);
        }
    }
    return IndexPair(vmin,vmax);
}

IndexPair getGlobalMinMax(const IndexVector &values, Esys_MPIInfo *mpiInfo)
{
    IndexPair minMaxLocal = getMinMax(values);
#ifdef ESYS_MPI
    index_t idRange[2] = {-minMaxLocal.first, minMaxLocal.second}, globalIdRange[2];
    MPI_Allreduce(idRange, globalIdRange, 2, MPI_INT, MPI_MAX, mpiInfo->comm);
    minMaxLocal.first = -globalIdRange[0];
    minMaxLocal.second = globalIdRange[1];
#endif
    return minMaxLocal;
}

IndexVector getUniqueValues(const IndexVector &values, Esys_MPIInfo *mpiInfo)
{
    index_t lastFoundValue = INDEX_T_MIN, minFoundValue, local_minFoundValue;
    IndexVector out;

    while (true) {
        // find smallest value bigger than lastFoundValue 
        minFoundValue = INDEX_T_MAX;
#pragma omp parallel private(local_minFoundValue)
        {
            local_minFoundValue = minFoundValue;
#pragma omp for schedule(static)
            for (size_t i = 0; i < values.size(); i++) {
                const index_t v = values[i];
                if ((v > lastFoundValue) && (v < local_minFoundValue))
                    local_minFoundValue = v;
            }
#pragma omp critical
            {
                if (local_minFoundValue < minFoundValue)
                    minFoundValue = local_minFoundValue;
            }
        }
#ifdef ESYS_MPI
        local_minFoundValue = minFoundValue;
        MPI_Allreduce(&local_minFoundValue, &minFoundValue, 1, MPI_INT, MPI_MIN, mpiInfo->comm);
#endif
        // if we found a new value we need to add this to the return vector

        if (minFoundValue < INDEX_T_MAX) {
            out.push_back(minFoundValue);
            lastFoundValue = minFoundValue;
        } else
            break;
    }
    return out;
}

IndexVector packMask(const IndexVector &mask)
{
    IndexVector out;
    for (size_t k=0; k<mask.size(); k++)
        if (mask[k] >= 0)
            out.push_back(k);

    return out;
}


} // end of namespace ripley

