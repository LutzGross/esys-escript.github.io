
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

#include "TrilinosMatrixAdapter.h" 
#include "BlockCrsMatrixWrapper.h" 
#include "CrsMatrixWrapper.h" 
#include "TrilinosAdapterException.h" 
#include "UnrolledBlockCrsMatrixWrapper.h" 
#include "util.h" 

#include <escript/index.h>
#include <escript/Data.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/SolverOptions.h>

namespace bp = boost::python;
using Teuchos::rcp;

namespace esys_trilinos {

TrilinosMatrixAdapter::TrilinosMatrixAdapter(escript::JMPI mpiInfo,
        int blocksize, const escript::FunctionSpace& fs,
        const_TrilinosGraph_ptr graph, bool isComplex, bool unroll) :
    AbstractSystemMatrix(blocksize, fs, blocksize, fs),
    m_mpiInfo(mpiInfo),
    m_isComplex(isComplex)
{
    if (isComplex) {
        if (blocksize == 1) {
            cmat = rcp(new CrsMatrixWrapper<cplx_t>(graph));
        } else if (unroll) {
            const_TrilinosGraph_ptr newGraph(util::unrollCrsGraph(graph, blocksize));
            cmat = rcp(new UnrolledBlockCrsMatrixWrapper<cplx_t>(newGraph, blocksize));
        } else {
            cmat = rcp(new BlockCrsMatrixWrapper<cplx_t>(graph, blocksize));
        }
    } else {
        if (blocksize == 1) {
            mat = rcp(new CrsMatrixWrapper<real_t>(graph));
        } else if (unroll) {
            const_TrilinosGraph_ptr newGraph(util::unrollCrsGraph(graph, blocksize));
            mat = rcp(new UnrolledBlockCrsMatrixWrapper<real_t>(newGraph, blocksize));
        } else {
            mat = rcp(new BlockCrsMatrixWrapper<real_t>(graph, blocksize));
        }
    }
}

template<>
void TrilinosMatrixAdapter::add<real_t>(const std::vector<LO>& rowIdx,
                                        const std::vector<real_t>& array)
{
    if (m_isComplex) {
        throw escript::ValueError("Please use complex array to add to complex "
                                  "matrix!");
    } else {
        (*mat).add(rowIdx, array);
    }
}

template<>
void TrilinosMatrixAdapter::add<cplx_t>(const std::vector<LO>& rowIdx,
                                        const std::vector<cplx_t>& array)
{
    if (m_isComplex) {
        (*cmat).add(rowIdx, array);
    } else {
        throw escript::ValueError("Please use real-valued array to add to "
                                  "real-valued matrix!");
    }
}

void TrilinosMatrixAdapter::ypAx(escript::Data& y, escript::Data& x) const
{
    if (x.getDataPointSize() != getBlockSize()) {
        throw TrilinosAdapterException("matrix vector product: block size "
                        "does not match the number of components in input.");
    } else if (y.getDataPointSize() != getBlockSize()) {
        throw TrilinosAdapterException("matrix vector product: block size "
                        "does not match the number of components in output.");
    } else if (x.getFunctionSpace() != getColumnFunctionSpace()) {
        throw TrilinosAdapterException("matrix vector product: matrix "
                   "function space and function space of input don't match.");
    } else if (y.getFunctionSpace() != getRowFunctionSpace()) {
        throw TrilinosAdapterException("matrix vector product: matrix "
                  "function space and function space of output don't match.");
    } else if (y.isComplex() != m_isComplex || x.isComplex() != m_isComplex) {
        throw escript::ValueError("matrix vector product: matrix complexity "
                  "must match vector complexity!");
    }

    // expand data objects
    x.expand();
    y.expand();
    y.requireWrite();

    if (m_isComplex) {
        const Teuchos::ArrayView<const cplx_t> xView(x.getSampleDataRO(0,
                        cplx_t(0)), x.getNumDataPoints()*x.getDataPointSize());
        const Teuchos::ArrayView<cplx_t> yView(y.getSampleDataRW(0, cplx_t(0)),
                                    y.getNumDataPoints()*y.getDataPointSize());
        cmat->ypAx(yView, xView);
    } else {
        const Teuchos::ArrayView<const real_t> xView(x.getSampleDataRO(0),
                                    x.getNumDataPoints()*x.getDataPointSize());
        const Teuchos::ArrayView<real_t> yView(y.getSampleDataRW(0),
                                    y.getNumDataPoints()*y.getDataPointSize());
        mat->ypAx(yView, xView);
    }
}

void TrilinosMatrixAdapter::setToSolution(escript::Data& out, escript::Data& in,
                                 bp::object& options) const
{
    if (out.getDataPointSize() != getBlockSize()) {
        throw TrilinosAdapterException("solve: block size does not match the number of components of solution.");
    } else if (in.getDataPointSize() != getBlockSize()) {
        throw TrilinosAdapterException("solve: block size does not match the number of components of right hand side.");
    } else if (out.getFunctionSpace() != getColumnFunctionSpace()) {
        throw TrilinosAdapterException("solve: matrix function space and function space of solution don't match.");
    } else if (in.getFunctionSpace() != getRowFunctionSpace()) {
        throw TrilinosAdapterException("solve: matrix function space and function space of right hand side don't match.");
    } else if (in.isComplex() != m_isComplex || out.isComplex() != m_isComplex) {
        throw escript::ValueError("solve: matrix complexity must match vector "
                                  "complexity!");
    }

    options.attr("resetDiagnostics")();
    escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy&>(options);
    out.expand();
    out.requireWrite();
    in.expand();

    if (m_isComplex) {
        const Teuchos::ArrayView<const cplx_t> bView(in.getSampleDataRO(0,
                      cplx_t(0)), in.getNumDataPoints()*in.getDataPointSize());
        const Teuchos::ArrayView<cplx_t> outView(out.getSampleDataRW(0,
                    cplx_t(0)), out.getNumDataPoints()*out.getDataPointSize());
        cmat->solve(outView, bView, sb);

    } else {
        const Teuchos::ArrayView<const real_t> bView(in.getSampleDataRO(0),
                                  in.getNumDataPoints()*in.getDataPointSize());
        const Teuchos::ArrayView<real_t> outView(out.getSampleDataRW(0),
                                out.getNumDataPoints()*out.getDataPointSize());
        mat->solve(outView, bView, sb);
    }
}

void TrilinosMatrixAdapter::nullifyRowsAndCols(escript::Data& row_q,
                                               escript::Data& col_q,
                                               double mdv)
{
    if (col_q.getDataPointSize() != getColumnBlockSize()) {
        throw TrilinosAdapterException("nullifyRowsAndCols: column block size does not match the number of components of column mask.");
    } else if (row_q.getDataPointSize() != getRowBlockSize()) {
        throw TrilinosAdapterException("nullifyRowsAndCols: row block size does not match the number of components of row mask.");
    } else if (col_q.getFunctionSpace() != getColumnFunctionSpace()) {
        throw TrilinosAdapterException("nullifyRowsAndCols: column function space and function space of column mask don't match.");
    } else if (row_q.getFunctionSpace() != getRowFunctionSpace()) {
        throw TrilinosAdapterException("nullifyRowsAndCols: row function space and function space of row mask don't match.");
    }

    col_q.expand();
    row_q.expand();
    const Teuchos::ArrayView<const real_t> rowMask(row_q.getSampleDataRO(0),
                            row_q.getNumDataPoints()*row_q.getDataPointSize());
    // we need remote values for col_q
    const Teuchos::ArrayView<const real_t> colView(col_q.getSampleDataRO(0),
                            col_q.getNumDataPoints()*col_q.getDataPointSize());

    if (m_isComplex)
        cmat->nullifyRowsAndCols(rowMask, colView, mdv);
    else
        mat->nullifyRowsAndCols(rowMask, colView, mdv);
}

void TrilinosMatrixAdapter::saveHB(const std::string& filename) const
{
    throw escript::NotImplementedError("Harwell-Boeing interface not available.");
}

}  // end of namespace

