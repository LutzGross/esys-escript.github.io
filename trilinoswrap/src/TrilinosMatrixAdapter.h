
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

#ifndef __ESYS_TRILINOSMATRIXADAPTER_H__
#define __ESYS_TRILINOSMATRIXADAPTER_H__

#include <escript/AbstractSystemMatrix.h>
#include <escript/FunctionSpace.h>

#include <trilinoswrap/AbstractMatrixWrapper.h>

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

class TrilinosMatrixAdapter : public escript::AbstractSystemMatrix
{
friend class OxleyDomain;

public:
    /**
       \brief
       Creates a new Trilinos CRS/block CRS matrix adapter using a compatible
       fill-complete Trilinos matrix graph.
    */
    TrilinosMatrixAdapter(escript::JMPI mpiInfo, int blocksize,
                          const escript::FunctionSpace& fs,
                          const_TrilinosGraph_ptr graph,
                          bool isComplex = false, bool unroll = false);

    virtual ~TrilinosMatrixAdapter() {}

    virtual void nullifyRowsAndCols(escript::Data& row_q, escript::Data& col_q,
                                    double mdv);  

    virtual void saveMM(const std::string& filename) const
    {
        if (m_isComplex)
            cmat->saveMM(filename);
        else
            mat->saveMM(filename);
    }

    virtual void saveHB(const std::string& filename) const;

    virtual void resetValues(bool preserveSolverData = false)
    {
        if (m_isComplex)
            cmat->resetValues(preserveSolverData);
        else
            mat->resetValues(preserveSolverData);
    }

    /// notifies the matrix that changes are about to happen.
    inline void resumeFill()
    {
        if (m_isComplex)
            cmat->resumeFill();
        else
            mat->resumeFill();
    }

    /// notifies the matrix that a set of changes has occured.
    inline void fillComplete(bool localOnly)
    {
        if (m_isComplex)
            cmat->fillComplete(localOnly);
        else
            mat->fillComplete(localOnly);
    }

    template<typename ST>
    void add(const std::vector<LO>& rowIndex, const std::vector<ST>& array);

    template<typename ST>
    void IztAIz(const Teuchos::RCP<Tpetra::CrsMatrix<ST,LO,GO,NT>> IZ, long n);

    int getNumRows();
    int getNumRowsWorkerCplx();
    int getNumRowsWorkerRealx();
    int getNumCols();
    int getNumColsWorkerCplx();
    int getNumColsWorkerRealx();

    inline int getBlockSize() const { return getRowBlockSize(); }

private:
    virtual void setToSolution(escript::Data& out, escript::Data& in,
                               boost::python::object& options) const;

    virtual void ypAx(escript::Data& y, escript::Data& x) const;

    escript::JMPI m_mpiInfo;
    bool m_isComplex;
    Teuchos::RCP<AbstractMatrixWrapper<real_t> > mat;
    Teuchos::RCP<AbstractMatrixWrapper<cplx_t> > cmat;
};

} // namespace esys_trilinos

#endif // __ESYS_TRILINOSMATRIXADAPTER_H__

