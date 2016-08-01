
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

#ifndef __ESYS_TRILINOS_BLOCKCRSMATRIXWRAPPER_H__
#define __ESYS_TRILINOS_BLOCKCRSMATRIXWRAPPER_H__

#include <trilinoswrap/AbstractMatrixWrapper.h>
#include <trilinoswrap/BelosWrapper.h>

#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

namespace esys_trilinos {

template<typename ST>
class BlockCrsMatrixWrapper : public AbstractMatrixWrapper<ST>
{
    typedef Tpetra::Experimental::BlockCrsMatrix<ST,LO,GO,NT> Matrix;

public:
    /**
       \brief
       Creates a new Trilinos Block CRS matrix wrapper using a compatible
       fill-complete Trilinos matrix graph and given block size.
    */
    BlockCrsMatrixWrapper(const_TrilinosGraph_ptr graph, int blocksize);

    void resetValues(bool preserveSolverData = false);

    /// notifies the matrix that changes are about to happen.
    inline void resumeFill() {}

    /// notifies the matrix that a set of changes has occured.
    inline void fillComplete(bool /*localOnly*/) {}

    void nullifyRowsAndCols(const Teuchos::ArrayView<const real_t>& rowMask,
                            const Teuchos::ArrayView<const real_t>& colView,
                            ST mdv);

    void add(const std::vector<LO>& rowIndex, const std::vector<ST>& array);

    void ypAx(const Teuchos::ArrayView<ST>& y,
              const Teuchos::ArrayView<const ST>& x) const;

    void solve(const Teuchos::ArrayView<ST>& x,
               const Teuchos::ArrayView<const ST>& b,
               escript::SolverBuddy& sb) const;

    void saveMM(const std::string& filename) const;

private:
    int blockSize;
    Matrix mat;
    mutable Teuchos::RCP<ProblemType<ST> > m_solver;
    MapType colPointMap;
    LO maxLocalRow;
};

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_BLOCKCRSMATRIXWRAPPER_H__

