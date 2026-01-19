
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESYS_TRILINOS_CRSMATRIXWRAPPER_H__
#define __ESYS_TRILINOS_CRSMATRIXWRAPPER_H__

#include <escript/AbstractSystemMatrix.h>

#include <trilinoswrap/AbstractMatrixWrapper.h>
#include <trilinoswrap/Amesos2Wrapper.h>
#include <trilinoswrap/BelosWrapper.h>

#include <Teuchos_RCPDecl.hpp>
#include <TpetraExt_MatrixMatrix_def.hpp>

#include <Tpetra_CrsMatrix.hpp>

namespace esys_trilinos {

template<typename ST>
class CrsMatrixWrapper : public AbstractMatrixWrapper<ST>
{
public:
    typedef Tpetra::CrsMatrix<ST,LO,GO,NT> Matrix;

    /**
       \brief
       Creates a new Trilinos CRS matrix adapter using a compatible
       fill-complete Trilinos matrix graph.
    */
    CrsMatrixWrapper(const_TrilinosGraph_ptr graph);

    void resetValues(bool preserveSolverData = false);

    /// notifies the matrix that changes are about to happen.
    inline void resumeFill()
    {
        mat.resumeFill();
    }

    /// notifies the matrix that a set of changes has occured.
    void fillComplete(bool localOnly);

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

    // Used by Oxley
    void IztAIz(const Teuchos::RCP<Tpetra::CrsMatrix<ST,LO,GO,NT>> IZ, long n);
    int getNumRows();
    int getNumCols();
    void add_single(const LO row, const LO col, const ST value);

protected:
    Matrix mat;
    mutable bool m_resetCalled;
    mutable Teuchos::RCP<ProblemType<ST> > m_solver;
    mutable Teuchos::RCP<OpType<ST> > m_preconditioner;
    mutable Teuchos::RCP<DirectSolverType<Matrix,VectorType<ST> > > m_direct;
    LO maxLocalRow;
};

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_CRSMATRIXWRAPPER_H__

