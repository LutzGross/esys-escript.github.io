
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

#ifndef __ESYS_TRILINOS_BLOCKCRSMATRIXWRAPPER_H__
#define __ESYS_TRILINOS_BLOCKCRSMATRIXWRAPPER_H__

#include <escript/AbstractSystemMatrix.h>

#include <trilinoswrap/AbstractMatrixWrapper.h>
#include <trilinoswrap/BelosWrapper.h>

// Use standard BlockCrsMatrix (experimental version deprecated in Trilinos 16.2+)
#include <Tpetra_BlockCrsMatrix.hpp>

#include <TpetraExt_MatrixMatrix_def.hpp>


namespace esys_trilinos {

template<typename ST>
class BlockCrsMatrixWrapper : public AbstractMatrixWrapper<ST>
{
    // Use standard BlockCrsMatrix (Trilinos 16.2+)
    typedef Tpetra::BlockCrsMatrix<ST,LO,GO,NT> Matrix;

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

    void IztAIz(const Teuchos::RCP<Tpetra::CrsMatrix<ST,LO,GO,NT>> IZ, long n) ;
    int getNumRows();
    int getNumCols();


private:
    int blockSize;
    Matrix mat;
    mutable Teuchos::RCP<ProblemType<ST> > m_solver;
    MapType colPointMap;
    LO maxLocalRow;
};

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_BLOCKCRSMATRIXWRAPPER_H__

