
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

#include <fstream>

#include "BlockCrsMatrixWrapper.h"
#include "BelosWrapper.h"
#include "PreconditionerFactory.h"
#include "TrilinosAdapterException.h"
#include "util.h"

#include <escript/index.h>
#include <escript/SolverOptions.h>

//#ifdef ESYS_TRILINOS_14
#include "Teuchos_ArrayRCPDecl.hpp"
//#else
//#include "Tpetra_createDeepCopy_CrsMatrix.hpp"
//#endif

#include "Tpetra_KokkosCompat_DefaultNode.hpp"
//#ifdef ESYS_NO_KOKKOSCOMPAT
//#include "Kokkos_DefaultNode.hpp"
//#else
//#include "KokkosCompat_DefaultNode.hpp"
//#endif

// Use Tpetra_Core.hpp (Tpetra_DefaultPlatform.hpp is deprecated in Trilinos 16.2+)
#include <Tpetra_Core.hpp>
// Use standard BlockCrsMatrix_Helpers (experimental version deprecated in Trilinos 16.2+)
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Tpetra_Vector.hpp>


#include "TpetraExt_MatrixMatrix_decl.hpp"
#include "TpetraExt_MatrixMatrix_def.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_TripleMatrixMultiply_decl.hpp"
#include "TpetraExt_TripleMatrixMultiply_def.hpp"


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
//#ifdef ESYS_TRILINOS_14
    using Kokkos::ALL;
    using Kokkos::subview;
    using local_inds_host_view_type = typename Tpetra::BlockCrsMatrix<ST,LO,GO,NT>::local_inds_host_view_type;
    using values_host_view_type = typename Tpetra::BlockCrsMatrix<ST,LO,GO,NT>::values_host_view_type;
    for (LO lrb = 0; lrb < static_cast<LO>(mat.getLocalNumRows()); lrb++) {
        LO numIndices = 0;
        local_inds_host_view_type indices;
        values_host_view_type values;
        GO row_number = (GO) lrb;
        numIndices=mat.getNumEntriesInGlobalRow(row_number);
        mat.getLocalRowView(lrb, indices, values);
//#else
 //   for (LO lrb = 0; lrb < mat.getNodeNumRows(); lrb++) {
 //       LO numIndices = 0;
 //       const LO * indices;
 //       ST * values;
 //       mat.getLocalRowView(lrb, indices, values, numIndices);
//#endif
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
    // Use standard blockCrsMatrixWriter (experimental version deprecated in Trilinos 16.2+)
    Tpetra::blockCrsMatrixWriter<ST,LO,GO,NT>(mat, os, params);
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

template<typename ST>
void BlockCrsMatrixWrapper<ST>::IztAIz(const Teuchos::RCP<Tpetra::CrsMatrix<ST,LO,GO,NT>> iz, long n) 
{
// #ifdef ESYS_TRILINOS_14
//     // Initialise some variables
//     Tpetra::global_size_t numGblIndices = n;
//     const esys_trilinos::GO indexBase = 0;
//     escript::JMPI m_mpiInfo;
//     auto comm = Teuchos::DefaultComm<int>::getComm();

//     // Create a new map and rcp matrix
//     typedef Tpetra::Map<LO,GO,NT> map_type;
//     Teuchos::RCP<const map_type> map = Teuchos::rcp ( new map_type (numGblIndices, indexBase, comm));
//     Matrix *result = new Matrix(map, n);

//     // Initialise some more variables
//     RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
//     params->set("No Nonlocal Changes", true);
//     const std::string& label = "ans";

//     const Teuchos::ArrayView<const Scalar> matIn(mat.getSampleDataRO(0, zero), mat.getNumDataPoints())
//     const Teuchos::ArrayView<const Scalar> izIn( iz.getSampleDataRO(0, zero),  iz.getNumDataPoints() )

//     const auto tmp_mat1 = new Teuchos::ArrayRCP<ST>::deepCopy(convertToCrsMatrix(matIn));
//     const auto iz_tmp   = new Teuchos::ArrayRCP<ST>::deepCopy(convertToCrsMatrix(izIn);
    
//     // Do the matrix-matrix multiplication
//     Tpetra::TripleMatrixMultiply::MultiplyRAP<ST,LO,GO,NT>(
//                     iz_tmp,true,tmp_mat1,false,*iz,false,*result,false,label,params);

//     Teuchos::ScalarTraits<ST>::magnitudeType threshold=1E-24;
//     Teuchos::removeCrsMatrixZeros(result, threshold);

//     mat.resumeFill();
//     mat=Tpetra::createDeepCopy(convertToBlockCrsMatrix(*result,blockSize));
//     mat.fillComplete(params);

//     delete tmp_mat1;
//     delete iz_tmp;
//     delete result;
// #else
//     // Initialise some variables
//     Tpetra::global_size_t numGblIndices = n;
//     const esys_trilinos::GO indexBase = 0;
//     escript::JMPI m_mpiInfo;
//     auto comm = Teuchos::DefaultComm<int>::getComm();

//     // Create a new map and rcp matrix
//     typedef Tpetra::Map<LO,GO,NT> map_type;
//     Teuchos::RCP<const map_type> map = Teuchos::rcp ( new map_type (numGblIndices, indexBase, comm));
//     Matrix *result = new Matrix(map, n);

//     // Initialise some more variables
//     RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
//     params->set("No Nonlocal Changes", true);
//     const std::string& label = "ans";
//     const auto tmp_mat1 = new Tpetra::createDeepCopy(convertToCrsMatrix(mat));
//     const auto iz_tmp = new Tpetra::createDeepCopy(*iz);
//     // const auto iz_tmp2 = Tpetra::createDeepCopy(*iz);

//     // Do the matrix-matrix multiplication
//     Tpetra::TripleMatrixMultiply::MultiplyRAP<ST,LO,GO,NT>(
//                     iz_tmp,true,tmp_mat1,false,*iz,false,*result,false,label,params);

//     Teuchos::ScalarTraits<ST>::magnitudeType threshold=1E-24;
//     Teuchos::removeCrsMatrixZeros(result, threshold);

//     mat.resumeFill();
//     mat=Tpetra::createDeepCopy(convertToBlockCrsMatrix(*result,blockSize));
//     mat.fillComplete(params);

//     delete tmp_mat1;
//     delete iz_tmp;
//     delete result;
// #endif
}

template<typename ST>
int BlockCrsMatrixWrapper<ST>::getNumRows()
{
    return mat.getGlobalNumRows();
}

template<typename ST>
int BlockCrsMatrixWrapper<ST>::getNumCols()
{
    return mat.getGlobalNumCols();
}


// instantiate
template class BlockCrsMatrixWrapper<real_t>;
template class BlockCrsMatrixWrapper<cplx_t>;

}  // end of namespace


