
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

/*   Paso: distribution                                                     */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_DISTRIBUTION_H__
#define __PASO_DISTRIBUTION_H__

#include "Common.h"
#include "PasoUtil.h"
#include "esysUtils/Esys_MPI.h"

namespace paso {

struct Distribution;
typedef boost::shared_ptr<Distribution> Distribution_ptr;
typedef boost::shared_ptr<const Distribution> const_Distribution_ptr;

/// describes the distribution of a vector stored on the local process
PASO_DLL_API
struct Distribution
{
    Distribution(Esys_MPIInfo* mpiInfo, const index_t* firstComponent,
                 index_t m, index_t b)
    {
        mpi_info = Esys_MPIInfo_getReference(mpiInfo);
        first_component = new index_t[mpi_info->size+1];
        for (dim_t i=0; i < mpi_info->size+1; ++i)
            first_component[i] = m*firstComponent[i]+b;
    }

    ~Distribution()
    {
        Esys_MPIInfo_free(mpi_info);
        delete[] first_component;
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

    inline dim_t numPositives(const double* x, dim_t block) const
    {
        const dim_t my_n = block * getMyNumComponents();
        dim_t my_out = util::numPositives(my_n, x);
        dim_t out;

#ifdef ESYS_MPI
#pragma omp single
        {
            MPI_Allreduce(&my_out, &out, 1, MPI_INT, MPI_SUM, mpi_info->comm);
        }
#else
        out = my_out;
#endif
        return out;
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
            out[i]=fmod(random_seed*(n_0+i+1), 1.);
        }

        random_seed = fmod(random_seed * (n+1.7), 1.);
        return out;
    }

    // process i has nodes with global indices first_component[i] to
    // first_component[i+1].
    index_t* first_component;
    dim_t reference_counter;
    Esys_MPIInfo *mpi_info;
    static double random_seed;
};

} // namespace paso

#endif // __PASO_DISTRIBUTION_H__

