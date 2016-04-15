
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
#include "Amesos2Wrapper.h" 
#include "BelosWrapper.h" 
#include "PreconditionerFactory.h" 
#include "TrilinosAdapterException.h" 

#include <escript/index.h>
#include <escript/Data.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/SolverOptions.h>

#include <BelosTpetraAdapter.hpp>
#include <BelosTypes.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>

namespace bp = boost::python;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace esys_trilinos {

template<typename ST>
void ypAxWorker(RCP<MatrixType<ST> > A, const Teuchos::ArrayView<ST>& y,
                const Teuchos::ArrayView<const ST>& x)
{
    RCP<VectorType<ST> > X = rcp(new VectorType<ST>(
                                              A->getRowMap(), x, x.size(), 1));
    RCP<VectorType<ST> > Y = rcp(new VectorType<ST>(
                                              A->getRowMap(), y, y.size(), 1));

    const ST alpha = Teuchos::ScalarTraits<ST>::one();
    const ST beta = Teuchos::ScalarTraits<ST>::one();

    // Y = beta*Y + alpha*A*X
    A->apply(*X, *Y, Teuchos::NO_TRANS, alpha, beta);
    Y->get1dCopy(y, y.size());
}

template<typename ST>
void solveWorker(RCP<MatrixType<ST> > A, const Teuchos::ArrayView<ST>& x,
                 const Teuchos::ArrayView<const ST>& b,
                 const escript::SolverBuddy& sb)
{
    RCP<VectorType<ST> > X = rcp(new VectorType<ST>(A->getDomainMap(), 1));
    RCP<VectorType<ST> > B = rcp(new VectorType<ST>(A->getRangeMap(), b,
                                                    b.size(), 1));

    if (sb.getSolverMethod() == escript::SO_METHOD_DIRECT) {
        RCP<DirectSolverType<ST> > solver = createDirectSolver<ST>(sb, A, X, B);
        if (sb.isVerbose()) {
            std::cout << solver->description() << std::endl;
            std::cout << "Performing symbolic factorization..." << std::flush;
        }
        solver->symbolicFactorization();
        if (sb.isVerbose()) {
            std::cout << "done\nPerforming numeric factorization..." << std::flush;
        }
        solver->numericFactorization();
        if (sb.isVerbose()) {
            std::cout << "done\nSolving system..." << std::flush;
        }
        solver->solve();
        if (sb.isVerbose()) {
            std::cout << "done" << std::endl;
            RCP<Teuchos::FancyOStream> fos(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)));
            solver->printTiming(*fos, Teuchos::VERB_HIGH);
        }

    } else { // iterative solver
        RCP<SolverType<ST> > solver = createSolver<ST>(sb);
        RCP<OpType<ST> > prec = createPreconditioner<ST>(A, sb);
        RCP<ProblemType<ST> > problem = rcp(new ProblemType<ST>(A, X, B));

        if (!prec.is_null()) {
            // Trilinos BiCGStab does not currently support left preconditioners
            if (sb.getSolverMethod() == escript::SO_METHOD_BICGSTAB)
                problem->setRightPrec(prec);
            else
                problem->setLeftPrec(prec);
        }
        problem->setProblem();
        solver->setProblem(problem);
        Belos::ReturnType result = solver->solve();
        if (sb.isVerbose()) {
            const int numIters = solver->getNumIters();
            if (result == Belos::Converged) {
                std::cout << "The solver took " << numIters
                   << " iteration(s) to reach a relative residual tolerance of "
                   << sb.getTolerance() << "." << std::endl;
            } else {
                std::cout << "The solver took " << numIters
                   << " iteration(s), but did not reach a relative residual "
                   "tolerance of " << sb.getTolerance() << "." << std::endl;
            }
        }
    }
    X->get1dCopy(x, x.size());
}

TrilinosMatrixAdapter::TrilinosMatrixAdapter(escript::JMPI mpiInfo,
        int blocksize, const escript::FunctionSpace& fs,
        const_TrilinosGraph_ptr graph, bool isComplex) :
    AbstractSystemMatrix(blocksize, fs, blocksize, fs),
    m_mpiInfo(mpiInfo),
    m_isComplex(isComplex)
{
    if (blocksize != 1) {
        throw escript::ValueError("Trilinos matrices only support blocksize 1 "
                                  "at the moment!");
    }
    if (isComplex) {
        cmat = rcp(new ComplexMatrix(graph));
        cmat->fillComplete();
        std::cout << "Matrix has " << cmat->getGlobalNumEntries()
                  << " entries." << std::endl;
    } else {
        mat = rcp(new RealMatrix(graph));
        mat->fillComplete();
        std::cout << "Matrix has " << mat->getGlobalNumEntries()
                  << " entries." << std::endl;
    }

}

void TrilinosMatrixAdapter::fillComplete(bool localOnly)
{
    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("No Nonlocal Changes", localOnly);
    if (m_isComplex)
        cmat->fillComplete(cmat->getDomainMap(), cmat->getRangeMap(), params);
    else
        mat->fillComplete(mat->getDomainMap(), mat->getRangeMap(), params);
}

template<>
void TrilinosMatrixAdapter::add<real_t>(const std::vector<LO>& rowIdx,
                                        const std::vector<real_t>& array)
{
    if (m_isComplex) {
        throw escript::ValueError("Please use complex array to add to complex "
                                  "matrix!");
    } else {
        addImpl<real_t>(*mat, rowIdx, array);
    }
}

template<>
void TrilinosMatrixAdapter::add<cplx_t>(const std::vector<LO>& rowIdx,
                                        const std::vector<cplx_t>& array)
{
    if (m_isComplex) {
        addImpl<cplx_t>(*cmat, rowIdx, array);
    } else {
        throw escript::ValueError("Please use real-valued array to add to "
                                  "real-valued matrix!");
    }
}

template<typename ST>
void TrilinosMatrixAdapter::addImpl(MatrixType<ST>& A,
                                    const std::vector<LO>& rowIdx,
                                    const std::vector<ST>& array)
{
    // NOTE: the reason this method takes a reference to the matrix and
    // we do the following with the row map is to avoid messing with shared
    // pointer use counters given that this method may be called from
    // parallel sections!
    const MapType& rowMap = *A.getRowMap();
    const int blockSize = getBlockSize();
    const size_t emSize = rowIdx.size();
    const LO myLast = rowMap.getMaxLocalIndex();
    std::vector<LO> cols(emSize*blockSize);
    std::vector<ST> vals(emSize*blockSize);
    for (size_t i = 0; i < emSize; i++) {
        for (int k = 0; k < blockSize; k++) {
            const LO row = rowIdx[i]*blockSize + k;
            if (row <= myLast) {
                cols.clear();
                vals.clear();
                for (int j = 0; j < emSize; j++) {
                    for (int m = 0; m < blockSize; m++) {
                        const LO col = rowIdx[j]*blockSize + m;
                        cols.push_back(col);
                        const size_t srcIdx =
                            INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                        vals.push_back(array[srcIdx]);
                    }
                }
                A.sumIntoLocalValues(row, cols, vals);
            }
        }
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

    // expand data object if necessary to be able to grab the whole data
    x.expand();
    y.expand();
    y.requireWrite();

    if (m_isComplex) {
        const Teuchos::ArrayView<const cplx_t> xView(
                x.getSampleDataRO(0, cplx_t(0)), x.getNumDataPoints());
        const Teuchos::ArrayView<cplx_t> yView(y.getSampleDataRW(0, cplx_t(0)),
                                               y.getNumDataPoints());
        ypAxWorker<cplx_t>(cmat, yView, xView);
    } else {
        const Teuchos::ArrayView<const real_t> xView(x.getSampleDataRO(0),
                                                     x.getNumDataPoints());
        const Teuchos::ArrayView<real_t> yView(y.getSampleDataRW(0),
                                               y.getNumDataPoints());
        ypAxWorker<real_t>(mat, yView, xView);
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
    escript::SolverBuddy sb = bp::extract<escript::SolverBuddy>(options);
    out.expand();
    in.expand();

    if (m_isComplex) {
        throw escript::NotImplementedError("complex solve not implemented!");
        const Teuchos::ArrayView<const cplx_t> bView(in.getSampleDataRO(0,
                                           cplx_t(0)), in.getNumDataPoints());
        const Teuchos::ArrayView<cplx_t> outView(out.getSampleDataRW(0,
                                           cplx_t(0)), out.getNumDataPoints());
        solveWorker<cplx_t>(cmat, outView, bView, sb);

    } else {
        const Teuchos::ArrayView<const real_t> bView(in.getSampleDataRO(0),
                                                     in.getNumDataPoints());
        const Teuchos::ArrayView<real_t> outView(out.getSampleDataRW(0),
                                                 out.getNumDataPoints());
        solveWorker<real_t>(mat, outView, bView, sb);
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
                                                   row_q.getNumDataPoints());
    // we need remote values for col_q
    const Teuchos::ArrayView<const real_t> colView(col_q.getSampleDataRO(0),
                                                   col_q.getNumDataPoints());

    // TODO:
    if (m_isComplex)
        throw escript::NotImplementedError("nullifyRowsAndCols: complex "
                                           "version not implemented");

    RCP<RealVector> lclCol = rcp(new RealVector(mat->getRowMap(), colView, colView.size(), 1));
    RCP<RealVector> gblCol = rcp(new RealVector(mat->getColMap(), 1));

    const ImportType importer(mat->getRowMap(), mat->getColMap());
    gblCol->doImport(*lclCol, importer, Tpetra::INSERT);
    Teuchos::ArrayRCP<const real_t> colMask(gblCol->getData(0));

    resumeFill();
// Can't use OpenMP here as replaceLocalValues() is not thread-safe.
//#pragma omp parallel for
    for (LO lclrow = 0; lclrow < mat->getNodeNumRows(); lclrow++) {
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const real_t> values;
        std::vector<GO> cols;
        std::vector<real_t> vals;
        mat->getLocalRowView(lclrow, indices, values);
        GO row = mat->getRowMap()->getGlobalElement(lclrow);
        for (size_t c = 0; c < indices.size(); c++) {
            const LO lclcol = indices[c];
            const GO col = mat->getColMap()->getGlobalElement(lclcol);
            if (rowMask[lclrow] > 0. || colMask[lclcol] > 0.) {
                cols.push_back(lclcol);
                vals.push_back(row==col ? (real_t)mdv : (real_t)0);
            }
        }
        if (cols.size() > 0)
            mat->replaceLocalValues(lclrow, cols, vals);
    }
    fillComplete(true);
}

void TrilinosMatrixAdapter::saveMM(const std::string& filename) const
{
    if (m_isComplex) {
        Tpetra::MatrixMarket::Writer<RealMatrix>::writeSparseFile(filename, mat);
    } else {
        Tpetra::MatrixMarket::Writer<ComplexMatrix>::writeSparseFile(filename, cmat);
    }
}

void TrilinosMatrixAdapter::saveHB(const std::string& filename) const
{
    throw escript::NotImplementedError("Harwell-Boeing interface not available.");
}

void TrilinosMatrixAdapter::resetValues()
{
    resumeFill();
    if (m_isComplex) {
        cmat->setAllToScalar(static_cast<const cplx_t>(0.));
    } else {
        mat->setAllToScalar(static_cast<const real_t>(0.));
    }
    fillComplete(true);
}

}  // end of namespace

