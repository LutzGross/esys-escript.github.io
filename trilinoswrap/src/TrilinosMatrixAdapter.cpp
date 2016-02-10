
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
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

#include "TrilinosMatrixAdapter.h" 
#include "TrilinosAdapterException.h" 

#include <esysUtils/index.h>
#include <escript/Data.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/SolverOptions.h>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>

namespace bp = boost::python;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace esys_trilinos {

TrilinosMatrixAdapter::TrilinosMatrixAdapter(esysUtils::JMPI mpiInfo,
        int blocksize, const escript::FunctionSpace& fs,
        const_TrilinosGraph_ptr graph) :
    AbstractSystemMatrix(blocksize, fs, blocksize, fs),
    m_mpiInfo(mpiInfo)
{
    mat = rcp(new MatrixType(graph));
    importer = rcp(new ImportType(mat->getRowMap(), mat->getColMap()));
}

void TrilinosMatrixAdapter::fillComplete(bool localOnly)
{
    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("No Nonlocal Changes", localOnly);
    mat->fillComplete(mat->getDomainMap(), mat->getRangeMap(), params);
}

void TrilinosMatrixAdapter::add(const std::vector<LO>& rowIdx,
                                const std::vector<double>& array)
{
    const int blockSize = getBlockSize();
    const size_t emSize = rowIdx.size();
    RCP<const MapType> rowMap(mat->getRowMap());
    const LO myLast = rowMap->getMaxLocalIndex();
    std::vector<LO> cols(emSize*blockSize);
    std::vector<ST> vals(emSize*blockSize);
    for (size_t i=0; i<emSize; i++) {
        for (int k=0; k<blockSize; k++) {
            const LO row = rowIdx[i]*blockSize + k;
            if (row <= myLast) {
                cols.clear();
                vals.clear();
                for (int j=0; j<emSize; j++) {
                    for (int m=0; m<blockSize; m++) {
                        const LO col = rowIdx[j]*blockSize + m;
                        cols.push_back(col);
                        const index_t srcIdx =
                            INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                        vals.push_back(array[srcIdx]);
                        //std::cout << "A[" << row << ","<<col<<"("<<lclcol<<")"
                        //    <<"]+="<<array[srcIdx]<<std::endl;
                    }
                }
                mat->sumIntoLocalValues(row, cols, vals);
            }
        }
    }
}

void TrilinosMatrixAdapter::ypAx(escript::Data& y, escript::Data& x) const
{
    if (x.getDataPointSize() != getBlockSize()) {
        throw TrilinosAdapterException("matrix vector product: block size does not match the number of components in input.");
    } else if (y.getDataPointSize() != getBlockSize()) {
        throw TrilinosAdapterException("matrix vector product: block size does not match the number of components in output.");
    } else if (x.getFunctionSpace() != getColumnFunctionSpace()) {
        throw TrilinosAdapterException("matrix vector product: matrix function space and function space of input don't match.");
    } else if (y.getFunctionSpace() != getRowFunctionSpace()) {
        throw TrilinosAdapterException("matrix vector product: matrix function space and function space of output don't match.");
    }

    // expand data object if necessary to be able to grab the whole data
    x.expand();
    y.expand();
    y.requireWrite();

    // we need remote values for x
    const Teuchos::ArrayView<const ST> xView(x.getSampleDataRO(0), x.getNumDataPoints());
    RCP<VectorType> lclX = rcp(new VectorType(mat->getRowMap(), xView, xView.size(), 1));
    RCP<VectorType> gblX = rcp(new VectorType(mat->getColMap(), 1));

    gblX->doImport(*lclX, *importer, Tpetra::INSERT);

    const ST alpha = Teuchos::ScalarTraits<ST>::one();
    const ST beta = Teuchos::ScalarTraits<ST>::one();
    const Teuchos::ArrayView<ST> yView(y.getSampleDataRW(0), y.getNumDataPoints());
    RCP<VectorType> Y = rcp(new VectorType(mat->getRowMap(), yView, yView.size(), 1));

    // Y = beta*Y + alpha*A*X
    //mat->apply(*X, *Y, Teuchos::NO_TRANS, alpha, beta);
    mat->localMultiply(*gblX, *Y, Teuchos::NO_TRANS, alpha, beta);
    Y->get1dCopy(yView, yView.size());
}

RCP<SolverType> TrilinosMatrixAdapter::createSolver(
                                         const escript::SolverBuddy& sb) const
{
    Belos::SolverFactory<ST, VectorType, OpType> factory;
    RCP<SolverType> solver;
    RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();

    solverParams->set("Convergence Tolerance", sb.getTolerance());
    solverParams->set("Maximum Iterations", sb.getIterMax());
    if (sb.isVerbose()) {
        solverParams->set("Verbosity", Belos::Errors + Belos::Warnings +
                Belos::TimingDetails + Belos::StatusTestDetails);
    }
    switch (sb.getSolverMethod()) {
        case escript::SO_DEFAULT:
        case escript::SO_METHOD_PCG:
            solver = factory.create("CG", solverParams);
            break;
        case escript::SO_METHOD_GMRES:
            solver = factory.create("GMRES", solverParams);
            break;
        case escript::SO_METHOD_LSQR:
            solver = factory.create("LSQR", solverParams);
            break;
        case escript::SO_METHOD_MINRES:
            solver = factory.create("MINRES", solverParams);
            break;
        case escript::SO_METHOD_TFQMR:
            solver = factory.create("TFQMR", solverParams);
            break;
        default:
            throw TrilinosAdapterException("Unsupported solver type requested.");
    }
    return solver;
}

RCP<OpType> TrilinosMatrixAdapter::createPreconditioner(
                                         const escript::SolverBuddy& sb) const
{
    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    Ifpack2::Factory factory;
    RCP<OpType> prec;
    RCP<Ifpack2::Preconditioner<ST,LO,GO,NT> > ifprec;

    switch (sb.getPreconditioner()) {
        case escript::SO_PRECONDITIONER_NONE:
            break;
        case escript::SO_PRECONDITIONER_AMG:
            {
                params->set("max levels", sb.getLevelMax());
                params->set("number of equations", getBlockSize());
                params->set("cycle type", sb.getCycleType()==1 ? "V" : "W");
                params->set("problem: symmetric", sb.isSymmetric());
                params->set("verbosity", sb.isVerbose()? "high":"low");
                RCP<OpType> A(mat);
                prec = MueLu::CreateTpetraPreconditioner(A, *params);
            }
            break;
        case escript::SO_PRECONDITIONER_ILUT:
            ifprec = factory.create<MatrixType>("ILUT", mat);
            params->set("fact: drop tolerance", sb.getDropTolerance());
            break;
        case escript::SO_PRECONDITIONER_GAUSS_SEIDEL:
            ifprec = factory.create<MatrixType>("RELAXATION", mat);
            params->set("relaxation: type", "Gauss-Seidel");
            params->set("relaxation: sweeps", sb.getNumSweeps());
            params->set("relaxation: damping factor", sb.getRelaxationFactor());
            break;
        case escript::SO_PRECONDITIONER_JACOBI:
            ifprec = factory.create<MatrixType>("RELAXATION", mat);
            params->set("relaxation: type", "Jacobi");
            params->set("relaxation: sweeps", sb.getNumSweeps());
            params->set("relaxation: damping factor", sb.getRelaxationFactor());
            break;
        case escript::SO_PRECONDITIONER_RILU:
            ifprec = factory.create<MatrixType>("RILUK", mat);
            params->set("fact: relax value", sb.getRelaxationFactor());
            break;
        default:
            throw TrilinosAdapterException("Unsupported preconditioner requested.");
    }
    if (!ifprec.is_null()) {
        ifprec->setParameters(*params);
        ifprec->initialize();
        ifprec->compute();
        prec = ifprec;
    }
    return prec;
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
    }

    options.attr("resetDiagnostics")();
    escript::SolverBuddy sb = bp::extract<escript::SolverBuddy>(options);
    out.expand();
    in.expand();

    if (sb.isVerbose()) {
        std::cout << "Matrix has " << mat->getGlobalNumEntries() << " entries." << std::endl;
    }

    const Teuchos::ArrayView<const ST> bView(in.getSampleDataRO(0), in.getNumDataPoints());
    RCP<VectorType> B = rcp(new VectorType(mat->getRangeMap(), bView, bView.size(), 1));
    RCP<VectorType> X = rcp(new VectorType(mat->getDomainMap(), 1));

    //solve...
    RCP<SolverType> solver = createSolver(sb);
    RCP<OpType> prec = createPreconditioner(sb);
    RCP<ProblemType> problem = rcp(new ProblemType(mat, X, B));

    if (!prec.is_null()) {
        problem->setLeftPrec(prec);
    }
    problem->setProblem();
    solver->setProblem(problem);
    Belos::ReturnType result = solver->solve();
    if (sb.isVerbose()) {
        const int numIters = solver->getNumIters();
        if (result == Belos::Converged) {
            std::cout << "The Belos solve took " << numIters
                << " iteration(s) to reach a relative residual tolerance of "
                << sb.getTolerance() << "." << std::endl;
        } else {
            std::cout << "The Belos solve took " << numIters
                << " iteration(s), but did not reach a relative residual "
                "tolerance of " << sb.getTolerance() << "." << std::endl;
        }
    }

    const Teuchos::ArrayView<ST> outView(out.getSampleDataRW(0), out.getNumDataPoints());
    X->get1dCopy(outView, outView.size());
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
    const Teuchos::ArrayView<const ST> rowMask(row_q.getSampleDataRO(0), row_q.getNumDataPoints());
    // we need remote values for col_q
    const Teuchos::ArrayView<const ST> colView(col_q.getSampleDataRO(0), col_q.getNumDataPoints());
    RCP<VectorType> lclCol = rcp(new VectorType(mat->getRowMap(), colView, colView.size(), 1));
    RCP<VectorType> gblCol = rcp(new VectorType(mat->getColMap(), 1));

    gblCol->doImport(*lclCol, *importer, Tpetra::INSERT);
    Teuchos::ArrayRCP<const ST> colMask(gblCol->getData(0));

    mat->resumeFill();
// OpenMP here will interact with OpenMP in Trilinos
//#pragma omp parallel for
    for (LO lclrow=0; lclrow < mat->getNodeNumRows(); lclrow++) {
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const ST> values;
        std::vector<GO> cols;
        std::vector<ST> vals;
        mat->getLocalRowView(lclrow, indices, values);
        GO row = mat->getRowMap()->getGlobalElement(lclrow);
        for (size_t c=0; c<indices.size(); c++) {
            const LO lclcol = indices[c];
            const GO col = mat->getColMap()->getGlobalElement(lclcol);
            if (rowMask[lclrow] > 0. || colMask[lclcol] > 0.) {
                cols.push_back(lclcol);
                vals.push_back(row==col ? (ST)mdv : (ST)0);
            }
        }
        if (cols.size() > 0)
            mat->replaceLocalValues(lclrow, cols, vals);
    }
    mat->fillComplete(mat->getDomainMap(), mat->getRangeMap());
}

void TrilinosMatrixAdapter::saveMM(const std::string& filename) const
{
    Tpetra::MatrixMarket::Writer<MatrixType>::writeSparseFile(filename, mat);
}

void TrilinosMatrixAdapter::saveHB(const std::string& filename) const
{
    throw TrilinosAdapterException("Harwell-Boeing interface not available.");
}

void TrilinosMatrixAdapter::resetValues()
{
    mat->resumeFill();
    mat->setAllToScalar(static_cast<const ST>(0.));
    mat->fillComplete();
}

}  // end of namespace

