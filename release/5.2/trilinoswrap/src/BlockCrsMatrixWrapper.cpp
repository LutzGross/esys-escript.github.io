
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

#include "BlockCrsMatrixWrapper.h" 
#include "BelosWrapper.h" 
#include "PreconditionerFactory.h" 
#include "TrilinosAdapterException.h" 
#include "util.h" 

#include <escript/index.h>
#include <escript/SolverOptions.h>

#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix_Helpers.hpp> // for writing
#include <Tpetra_Vector.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace esys_trilinos {

template<typename ST>
BlockCrsMatrixWrapper<ST>::BlockCrsMatrixWrapper(const_TrilinosGraph_ptr graph,
                                                 int blocksize) :
    blockSize(blocksize),
    mat(*graph, blocksize)
{
    // initialize column point map, needed by nullifyRowsAndCols to communicate
    // remote values
    colPointMap = BlockVectorType<ST>::makePointMap(*mat.getColMap(), blockSize);
    maxLocalRow = graph->getRowMap()->getMaxLocalIndex();
}

template<typename ST>
void BlockCrsMatrixWrapper<ST>::add(const std::vector<LO>& rowIdx,
                                    const std::vector<ST>& array)
{
    const size_t emSize = rowIdx.size();
    std::vector<LO> cols(emSize);
    std::vector<ST> vals(emSize*blockSize*blockSize);
    for (size_t i = 0; i < emSize; i++) {
        const LO row = rowIdx[i];
        if (row <= maxLocalRow) {
            for (int j = 0; j < emSize; j++) {
                cols[j] = rowIdx[j];
                for (int k = 0; k < blockSize; k++) {
                    for (int m = 0; m < blockSize; m++) {
                        const size_t srcIdx =
                            INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                        const size_t destIdx =
                            INDEX3(m, k, j, blockSize, blockSize);
                        vals[destIdx] = array[srcIdx];
                    }
                }
            }
            mat.sumIntoLocalValues(row, &cols[0], &vals[0], emSize);
        }
    }
}

template<typename ST>
void BlockCrsMatrixWrapper<ST>::ypAx(const Teuchos::ArrayView<ST>& y,
                                   const Teuchos::ArrayView<const ST>& x) const
{
    typedef VectorType<ST> Vector;
    RCP<Vector> X = rcp(new Vector(mat.getDomainMap(), x, x.size(), 1));
    RCP<Vector> Y = rcp(new Vector(mat.getDomainMap(), y, y.size(), 1));

    const ST alpha = Teuchos::ScalarTraits<ST>::one();
    const ST beta = Teuchos::ScalarTraits<ST>::one();

    // Y = beta*Y + alpha*A*X
    mat.apply(*X, *Y, Teuchos::NO_TRANS, alpha, beta);
    Y->get1dCopy(y, y.size());
}

template<typename ST>
void BlockCrsMatrixWrapper<ST>::solve(const Teuchos::ArrayView<ST>& x,
                                      const Teuchos::ArrayView<const ST>& b,
                                      escript::SolverBuddy& sb) const
{
    typedef VectorType<ST> Vector;

    RCP<Vector> X = rcp(new Vector(mat.getDomainMap(), 1));
    RCP<Vector> B = rcp(new Vector(mat.getRangeMap(), b, b.size(), 1));
    RCP<const Matrix> A = rcpFromRef(mat);

    if (escript::isDirectSolver(sb.getSolverMethod())) {
        throw TrilinosAdapterException("Amesos2 does not currently support "
                                       "block matrices!");
#if 0
        RCP<DirectSolverType<Matrix,Vector> > solver(m_direct);
        if (solver.is_null()) {
            solver = createDirectSolver<Matrix,Vector>(sb, A, X, B);
            m_direct = solver;
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
                std::cout << "done\n" << std::flush;
            }
        } else {
            solver->setX(X);
            solver->setB(B);
        }
        if (sb.isVerbose()) {
            std::cout << "Solving system..." << std::flush;
        }
        solver->solve();
        if (sb.isVerbose()) {
            std::cout << "done" << std::endl;
            RCP<Teuchos::FancyOStream> fos(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)));
            solver->printTiming(*fos, Teuchos::VERB_HIGH);
        }
#endif
    } else { // iterative solver
        double t0 = Teuchos::Time::wallTime();
        RCP<ProblemType<ST> > problem(m_solver);
        if (problem.is_null()) {
            problem = rcp(new ProblemType<ST>(A, X, B));
            m_solver = problem;
            RCP<OpType<ST> > prec = createPreconditioner<ST>(A, sb);
            if (!prec.is_null()) {
                // Trilinos BiCGStab does not support left preconditioners
                if (sb.getSolverMethod() == escript::SO_METHOD_BICGSTAB)
                    problem->setRightPrec(prec);
                else
                    problem->setLeftPrec(prec);
            }
            problem->setHermitian(sb.isSymmetric());
            problem->setProblem();
        } else {
            for (auto t: problem->getTimers()) {
                t->reset();
            }
            problem->setProblem(X, B);
        }

        double t1 = Teuchos::Time::wallTime();
        RCP<SolverType<ST> > solver = createSolver<ST>(sb);
        solver->setProblem(problem);
        Belos::ReturnType result = solver->solve();
        double t2 = Teuchos::Time::wallTime();
        const int numIters = solver->getNumIters();
        double tol = sb.getTolerance();
        try {
            tol = solver->achievedTol();
        } catch (...) {
        }
        if (sb.isVerbose()) {
            if (result == Belos::Converged) {
                sb.updateDiagnostics("converged", true);
                std::cout << "The solver took " << numIters
                   << " iteration(s) to reach a residual tolerance of "
                   << tol << "." << std::endl;
            } else {
                std::cout << "The solver took " << numIters
                   << " iteration(s), but did not reach a relative residual "
                   "tolerance of " << sb.getTolerance() << "." << std::endl;
            }
        }
        double solverTime = 0.;
        for (auto t: problem->getTimers()) {
            solverTime += t->totalElapsedTime();
        }
        sb.updateDiagnostics("set_up_time", t1-t0);
        sb.updateDiagnostics("net_time", solverTime);
        sb.updateDiagnostics("time", t2-t0);
        sb.updateDiagnostics("num_iter", numIters);
        sb.updateDiagnostics("residual_norm", tol);
    }
    X->get1dCopy(x, x.size());
}

template<typename ST>
void BlockCrsMatrixWrapper<ST>::nullifyRowsAndCols(
                               const Teuchos::ArrayView<const real_t>& rowMask,
                               const Teuchos::ArrayView<const real_t>& colView,
                               ST mdv)
{
    RCP<VectorType<real_t> > lclCol = rcp(new VectorType<real_t>(
                               mat.getRangeMap(), colView, colView.size(), 1));
    RCP<MapType> cpm = rcpFromRef(colPointMap);
    RCP<VectorType<real_t> > gblCol = rcp(new VectorType<real_t>(cpm, 1));

    const ImportType importer(mat.getRangeMap(), cpm);
    gblCol->doImport(*lclCol, importer, Tpetra::INSERT);
    Teuchos::ArrayRCP<const real_t> colMask(gblCol->getData(0));
    const real_t eps = escript::DataTypes::real_t_eps();
    const ST zero = Teuchos::ScalarTraits<ST>::zero();

// Can't use OpenMP here as replaceLocalValues() is not thread-safe.
//#pragma omp parallel for
    // loop through local row blocks
    for (LO lrb = 0; lrb < mat.getNodeNumRows(); lrb++) {
        LO numIndices = 0;
        const LO* indices;
        ST* values;
        mat.getLocalRowView(lrb, indices, values, numIndices);
        std::vector<GO> cols(numIndices);
        std::vector<ST> vals(numIndices*blockSize*blockSize);
        const GO rowblk = mat.getRowMap()->getGlobalElement(lrb);
        for (LO c = 0; c < numIndices; c++) {
            // local/global column block
            const LO lcb = indices[c];
            const GO colblk = mat.getColMap()->getGlobalElement(lcb);
            cols[c] = lcb;
            for (LO ri = 0; ri < blockSize; ri++) {
                const LO lclrow = lrb * blockSize + ri;
                const GO row = rowblk * blockSize + ri;
                for (LO ci = 0; ci < blockSize; ci++) {
                    const LO lclcol = lcb * blockSize + ci;
                    const GO col = colblk * blockSize + ci;
                    const size_t idx = INDEX3(ci, ri, c, blockSize, blockSize);
                    if (std::abs(rowMask[lclrow]) > eps || std::abs(colMask[lclcol]) > eps) {
                        vals[idx] = (row==col ? mdv : zero);
                    } else {
                        // we need to add full blocks so add current value
                        vals[idx] = values[idx];
                    }
                }
            }
        }
        mat.replaceLocalValues(lrb, &cols[0], &vals[0], numIndices);
    }
}

template<typename ST>
void BlockCrsMatrixWrapper<ST>::saveMM(const std::string& filename) const
{
    Teuchos::ParameterList params;
    // for compatibility with paso, not strictly required.
    params.set("precision", 15);
    std::ofstream os(filename);
    Tpetra::Experimental::blockCrsMatrixWriter<ST,LO,GO,NT>(mat, os, params);
    os.close();
}

template<typename ST>
void BlockCrsMatrixWrapper<ST>::resetValues(bool preserveSolverData)
{
    mat.setAllToScalar(static_cast<ST>(0.));
    if (!preserveSolverData) {
        m_solver.reset();
    }
}

// instantiate
template class BlockCrsMatrixWrapper<real_t>;
template class BlockCrsMatrixWrapper<cplx_t>;

}  // end of namespace

