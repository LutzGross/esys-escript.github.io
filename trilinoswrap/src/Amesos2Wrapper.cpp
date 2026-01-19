
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

#include <trilinoswrap/Amesos2Wrapper.h>
#include <trilinoswrap/TrilinosAdapterException.h>
#include <trilinoswrap/util.h>

#include <escript/SolverOptions.h>

#include <Amesos2.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <boost/python/dict.hpp>

using Teuchos::RCP;

namespace bp = boost::python;

namespace esys_trilinos {

template<class Matrix, class Vector>
RCP<DirectSolverType<Matrix,Vector> > createDirectSolver(
                                               const escript::SolverBuddy& sb,
                                               RCP<const Matrix> A,
                                               RCP<Vector> X,
                                               RCP<const Vector> B)
{
    using util::extractParamIfSet;

    typedef typename Matrix::scalar_type ST;

    RCP<DirectSolverType<Matrix,Vector> > solver;
    RCP<Teuchos::ParameterList> amesosParams = Teuchos::parameterList("Amesos2");
    const bp::dict& pyParams = sb.getTrilinosParameters();

    const escript::SolverOptions method = sb.getSolverMethod();
    // did user request specific direct solver or not?
    const bool dontcare = method == escript::SO_METHOD_DIRECT;

    if ((dontcare || method == escript::SO_METHOD_DIRECT_TRILINOS) &&
            Amesos2::query("klu2")) {
        solver = Amesos2::create<Matrix, Vector>("klu2", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        // the doco says these params exist but clearly they don't :-(
        extractParamIfSet<std::string>("Trans", pyParams, solverParams);
        extractParamIfSet<bool>("Equil", pyParams, solverParams);
        extractParamIfSet<std::string>("IterRefine", pyParams, solverParams);
        extractParamIfSet<bool>("SymmetricMode", pyParams, solverParams);
        extractParamIfSet<ST>("DiagPivotThresh", pyParams, solverParams);
        extractParamIfSet<std::string>("ColPerm", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if ((dontcare || method == escript::SO_METHOD_DIRECT_MUMPS) &&
            Amesos2::query("MUMPS")) {
        solver = Amesos2::create<Matrix, Vector>("MUMPS", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        if (sb.isVerbose()) {
            solverParams.set("ICNTL(4)", 4);
        }
        extractParamIfSet<int>("ICNTL(1)", pyParams, solverParams);
        extractParamIfSet<int>("ICNTL(2)", pyParams, solverParams);
        extractParamIfSet<int>("ICNTL(3)", pyParams, solverParams);
        extractParamIfSet<int>("ICNTL(4)", pyParams, solverParams);
        extractParamIfSet<int>("ICNTL(6)", pyParams, solverParams);
        extractParamIfSet<int>("ICNTL(9)", pyParams, solverParams);
        extractParamIfSet<int>("ICNTL(11)", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if ((dontcare || method == escript::SO_METHOD_DIRECT_TRILINOS) &&
            Amesos2::query("Basker")) {
        solver = Amesos2::create<Matrix, Vector>("Basker", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        solverParams.set("MatrixType", (sb.isSymmetric() || sb.isHermitian()) ? "symmetric" : "general");
    } else if ((dontcare || method == escript::SO_METHOD_DIRECT_SUPERLU) &&
            Amesos2::query("superludist")) {
        solver = Amesos2::create<Matrix, Vector>("superludist", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());  
        solverParams.set("MatrixType", (sb.isSymmetric() || sb.isHermitian()) ? "symmetric" : "general");
        extractParamIfSet<int>("npcol", pyParams, solverParams);
        extractParamIfSet<int>("nprow", pyParams, solverParams);
        extractParamIfSet<std::string>("ColPerm", pyParams, solverParams);
        extractParamIfSet<bool>("ReplaceTinyPivot", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if ((dontcare || method == escript::SO_METHOD_DIRECT_SUPERLU) &&
            Amesos2::query("superlu")) {
        solver = Amesos2::create<Matrix, Vector>("superlu", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        solverParams.set("SymmetricMode", sb.isSymmetric());
        extractParamIfSet<std::string>("Trans", pyParams, solverParams);
        extractParamIfSet<bool>("Equil", pyParams, solverParams);
        extractParamIfSet<std::string>("IterRefine", pyParams, solverParams);
        extractParamIfSet<bool>("SymmetricMode", pyParams, solverParams);
        extractParamIfSet<ST>("DiagPivotThresh", pyParams, solverParams);
        extractParamIfSet<std::string>("ColPerm", pyParams, solverParams);
        extractParamIfSet<bool>("ILU_Flag", pyParams, solverParams);
        extractParamIfSet<ST>("ILU_DropTol", pyParams, solverParams);
        extractParamIfSet<ST>("ILU_FillFactor", pyParams, solverParams);
        extractParamIfSet<std::string>("ILU_Norm", pyParams, solverParams);
        extractParamIfSet<std::string>("ILU_MILU", pyParams, solverParams);
        extractParamIfSet<ST>("ILU_FillTol", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if ((dontcare || method == escript::SO_METHOD_DIRECT_SUPERLU) &&
            Amesos2::query("superlumt")) {
        solver = Amesos2::create<Matrix, Vector>("superlumt", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
#ifdef _OPENMP  
        solverParams.set("nprocs", omp_get_max_threads());
#else
        solverParams.set("nprocs", 1);
#endif
        solverParams.set("SymmetricMode", sb.isSymmetric());
        extractParamIfSet<int>("nprocs", pyParams, solverParams);
        extractParamIfSet<std::string>("trans", pyParams, solverParams);
        extractParamIfSet<int>("panel_size", pyParams, solverParams);
        extractParamIfSet<int>("relax", pyParams, solverParams);
        extractParamIfSet<bool>("Equil", pyParams, solverParams);
        extractParamIfSet<bool>("SymmetricMode", pyParams, solverParams);
        extractParamIfSet<ST>("DiagPivotThresh", pyParams, solverParams);
        extractParamIfSet<std::string>("ColPerm", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if ((dontcare || method == escript::SO_METHOD_DIRECT_PARDISO) &&
            Amesos2::query("pardiso_mkl")) {
        solver = Amesos2::create<Matrix, Vector>("pardiso_mkl", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        extractParamIfSet<int>("IPARM(2)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(4)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(8)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(10)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(18)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(24)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(25)", pyParams, solverParams);
        extractParamIfSet<int>("IPARM(60)", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if (Amesos2::query("amesos2_cholmod")) {
        solver = Amesos2::create<Matrix, Vector>("amesos2_cholmod", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        solverParams.set("SymmetricMode", sb.isSymmetric());
        extractParamIfSet<std::string>("Trans", pyParams, solverParams);
        extractParamIfSet<bool>("Equil", pyParams, solverParams);
        extractParamIfSet<std::string>("IterRefine", pyParams, solverParams);
        extractParamIfSet<bool>("SymmetricMode", pyParams, solverParams);
        extractParamIfSet<ST>("DiagPivotThresh", pyParams, solverParams);
        extractParamIfSet<std::string>("ColPerm", pyParams, solverParams);
        amesosParams->set(solver->name(), solverParams);
    } else if (Amesos2::query("lapack")) {
        solver = Amesos2::create<Matrix, Vector>("lapack", A, X, B);
        Teuchos::ParameterList solverParams(solver->name());
        solverParams.set("MatrixType", (sb.isSymmetric() || sb.isHermitian()) ? "symmetric" : "general");
    } else {
        if (dontcare) {
            throw TrilinosAdapterException("Could not find an Amesos2 direct solver!");
        } else {
            throw TrilinosAdapterException("The requested direct solver is not available!");
        }
    }
    solver->setParameters(amesosParams);
    return solver;
}

typedef Tpetra::CrsMatrix<real_t,LO,GO,NT> RealMatrix;
typedef Tpetra::CrsMatrix<cplx_t,LO,GO,NT> ComplexMatrix;

// instantiate
template
RCP<DirectSolverType<RealMatrix, RealVector> >
createDirectSolver<RealMatrix,RealVector>(const escript::SolverBuddy& sb,
                                          RCP<const RealMatrix> A,
                                          RCP<RealVector> X,
                                          RCP<const RealVector> B);
template
RCP<DirectSolverType<ComplexMatrix, ComplexVector> >
createDirectSolver<ComplexMatrix, ComplexVector>(
                                          const escript::SolverBuddy& sb,
                                          RCP<const ComplexMatrix> A,
                                          RCP<ComplexVector> X,
                                          RCP<const ComplexVector> B);

/* Amesos2 does not currently support block matrices!
template
RCP<DirectSolverType<RealBlockMatrix, RealBlockVector> >
createDirectSolver<RealBlockMatrix,RealBlockVector>(
                                          const escript::SolverBuddy& sb,
                                          RCP<const RealBlockMatrix> A,
                                          RCP<RealBlockVector> X,
                                          RCP<const RealBlockVector> B);
template
RCP<DirectSolverType<ComplexBlockMatrix, ComplexBlockVector> >
createDirectSolver<ComplexBlockMatrix, ComplexBlockVector>(
                                          const escript::SolverBuddy& sb,
                                          RCP<const ComplexBlockMatrix> A,
                                          RCP<ComplexBlockVector> X,
                                          RCP<const ComplexBlockVector> B);
*/

}  // end of namespace

