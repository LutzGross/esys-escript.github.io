
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
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

#ifndef __ESYS_TRILINOS_ABSTRACTMATRIXWRAPPER_H__
#define __ESYS_TRILINOS_ABSTRACTMATRIXWRAPPER_H__

#include <trilinoswrap/types.h>

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

template<typename ST>
class AbstractMatrixWrapper
{
public:
    virtual ~AbstractMatrixWrapper() {}

    virtual void resetValues() = 0;

    /// notifies the matrix that changes are about to happen.
    virtual void resumeFill() = 0;

    /// notifies the matrix that a set of changes has occured.
    virtual void fillComplete(bool localOnly) = 0;

    /// sets certain entries to zero, and main diagonal to `mdv`
    virtual void nullifyRowsAndCols(
                              const Teuchos::ArrayView<const real_t>& rowMask,
                              const Teuchos::ArrayView<const real_t>& colView,
                              ST mdv) = 0;

    /// adds entries of an element matrix to this matrix
    virtual void add(const std::vector<LO>& rowIndex,
                     const std::vector<ST>& array) = 0;

    /// computes y += Ax
    virtual void ypAx(const Teuchos::ArrayView<ST>& y,
                      const Teuchos::ArrayView<const ST>& x) const = 0;

    /// solves Ax = b
    virtual void solve(const Teuchos::ArrayView<ST>& x,
                       const Teuchos::ArrayView<const ST>& b,
                       escript::SolverBuddy& sb) const = 0;

    /// saves matrix in Matrix Market (MM) format
    virtual void saveMM(const std::string& filename) const = 0;
};

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_ABSTRACTMATRIXWRAPPER_H__

