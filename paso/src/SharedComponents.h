
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


/****************************************************************************/

/*   Paso: shared components                                                */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_SHAREDCOMPONENTS_H__
#define __PASO_SHAREDCOMPONENTS_H__

#include "Paso.h"

namespace paso {

struct SharedComponents;
typedef boost::shared_ptr<SharedComponents> SharedComponents_ptr;
typedef boost::shared_ptr<const SharedComponents> const_SharedComponents_ptr;

PASO_DLL_API
struct SharedComponents
{
    SharedComponents(dim_t localLength, dim_t nNeighbours,
            const Esys_MPI_rank* neighbours, const index_t* sharedArray,
            const index_t* offset, index_t m, index_t b, Esys_MPIInfo *mpiInfo)
        : local_length(localLength*m),
          numNeighbors(nNeighbours)
    {
        mpi_info = Esys_MPIInfo_getReference(mpiInfo);
        neighbor = new Esys_MPI_rank[numNeighbors];
        if (!offset) {
            numSharedComponents = 0;
        } else {
            numSharedComponents = offset[nNeighbours] * m;
        }
        shared = new index_t[numSharedComponents];
        offsetInShared = new index_t[numNeighbors+1];
        if (numNeighbors > 0 && offset != NULL) {
#pragma omp parallel
            {
#pragma omp for
                for (dim_t i=0; i < numNeighbors; i++) {
                    neighbor[i] = neighbours[i];
                    offsetInShared[i] = offset[i] * m;
                }
                offsetInShared[numNeighbors] = offset[nNeighbours] * m;
#pragma omp for
                for (dim_t i=0; i<offset[nNeighbours]; i++) {
                    const index_t itmp=m*sharedArray[i]+b;
                    for (dim_t j=0; j < m; ++j)
                        shared[m*i+j]=itmp+j;
                }
            }
        } else {
            offsetInShared[numNeighbors]=0;
        }
    }

    ~SharedComponents()
    {
        delete[] neighbor;
        delete[] shared;
        delete[] offsetInShared;
        Esys_MPIInfo_free(mpi_info);
    }

    /// local array length shared
    dim_t local_length;

    /// number of processors sharing values with this processor
    dim_t numNeighbors;

    /// offsetInSharedInput[i] points to the first input value in array shared
    /// for processor i. Has length numNeighbors+1
    index_t* offsetInShared;

    /// list of the processors sharing values with this processor
    Esys_MPI_rank* neighbor;

    /// list of the (local) components which are shared with other processors.
    /// Has length numSharedComponents
    index_t* shared;

    /// = offsetInShared[numNeighbors]
    dim_t numSharedComponents;

    Esys_MPIInfo *mpi_info;
};

} // namespace paso

#endif // __PASO_SHAREDCOMPONENTS_H__

