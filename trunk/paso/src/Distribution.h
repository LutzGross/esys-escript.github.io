
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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


/****************************************************************************/

/*   Paso: distribution                                                     */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_DISTRIBUTION_H__
#define __PASO_DISTRIBUTION_H__

#include "Paso.h"
#include "PasoUtil.h"

namespace paso {

struct Distribution;
typedef boost::shared_ptr<Distribution> Distribution_ptr;
typedef boost::shared_ptr<const Distribution> const_Distribution_ptr;

/// describes the distribution of a vector stored on the local process
PASO_DLL_API
struct Distribution
{
    Distribution(const escript::JMPI& mpiInfo,
                 const std::vector<index_t>& firstComponent, index_t m,
                 index_t b) :
        mpi_info(mpiInfo)
    {
        first_component.resize(mpi_info->size+1);
        for (int i=0; i < mpi_info->size+1; ++i)
            first_component[i] = m*firstComponent[i]+b;
    }

    inline index_t getFirstComponent() const
    {
        return first_component[mpi_info->rank];
    }

    inline index_t getLastComponent() const
    {
        return first_component[mpi_info->rank+1];
    }


    inline dim_t getGlobalNumComponents() const
    {
        return getMaxGlobalComponents()-getMinGlobalComponents();
    }

    inline dim_t getMyNumComponents() const
    {
        return getLastComponent()-getFirstComponent();
    }

    inline dim_t getMinGlobalComponents() const
    {
        return first_component[0];
    }

    inline dim_t getMaxGlobalComponents() const
    {
        return first_component[mpi_info->size];
    }

    inline double* createRandomVector(dim_t block) const
    {
        const index_t n_0 = getFirstComponent() * block;
        const index_t n_1 = getLastComponent() * block;
        const index_t n = getGlobalNumComponents() * block;
        const dim_t my_n = n_1-n_0;

        double* out = new double[my_n];

#pragma omp parallel for schedule(static)
        for (index_t i=0; i<my_n; ++i) {
            out[i] = fmod(random_seed*(n_0+i+1), 1.);
        }

        random_seed = fmod(random_seed * (n+1.7), 1.);
        return out;
    }

    // process i has nodes with global indices first_component[i] to
    // first_component[i+1].
    std::vector<index_t> first_component;
    escript::JMPI mpi_info;
    static double random_seed;
};

} // namespace paso

#endif // __PASO_DISTRIBUTION_H__

