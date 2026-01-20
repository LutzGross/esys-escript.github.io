
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_DISTRIBUTION_H__
#define __ESCRIPT_DISTRIBUTION_H__

#include <escript/DataTypes.h>

namespace escript {

struct Distribution;
typedef boost::shared_ptr<Distribution> Distribution_ptr;
typedef boost::shared_ptr<const Distribution> const_Distribution_ptr;

/// Describes the distribution of a vector across processes.
/// Process i has entries with global indices first_component[i] to
/// first_component[i+1].
struct Distribution
{
    Distribution(JMPI mpiInfo, const DataTypes::IndexVector& firstComponent,
                 DataTypes::index_t m = 1, DataTypes::index_t b = 0) :
        mpi_info(mpiInfo)
    {
        first_component.resize(mpi_info->size + 1);
        for (int i = 0; i < mpi_info->size+1; ++i)
            first_component[i] = m * firstComponent[i] + b;
    }

    inline DataTypes::index_t getFirstComponent() const
    {
        return first_component[mpi_info->rank];
    }

    inline DataTypes::index_t getLastComponent() const
    {
        return first_component[mpi_info->rank+1];
    }

    inline DataTypes::dim_t getGlobalNumComponents() const
    {
        return getMaxGlobalComponents()-getMinGlobalComponents();
    }

    inline DataTypes::dim_t getMyNumComponents() const
    {
        return getLastComponent()-getFirstComponent();
    }

    inline DataTypes::dim_t getMinGlobalComponents() const
    {
        return first_component[0];
    }

    inline DataTypes::dim_t getMaxGlobalComponents() const
    {
        return first_component[mpi_info->size];
    }

    DataTypes::IndexVector first_component;
    JMPI mpi_info;
};

} // namespace escript

#endif // __ESCRIPT_DISTRIBUTION_H__

