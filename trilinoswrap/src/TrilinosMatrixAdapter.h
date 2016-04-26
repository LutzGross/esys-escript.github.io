
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
public:
    /**
       \brief
       Creates a new Trilinos CRS/block CRS matrix adapter using a compatible
       fill-complete Trilinos matrix graph.
    */
    TrilinosMatrixAdapter(escript::JMPI mpiInfo, int blocksize,
                          const escript::FunctionSpace& fs,
                          const_TrilinosGraph_ptr graph, bool isComplex=false);

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

    virtual void resetValues()
    {
        if (m_isComplex)
            cmat->resetValues();
        else
            mat->resetValues();
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

