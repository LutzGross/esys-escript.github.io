
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

#include <trilinoswrap/Amesos2Wrapper.h>
#include <trilinoswrap/TrilinosAdapterException.h>

#include <escript/SolverOptions.h>

#include <Amesos2.hpp>
#include <Tpetra_CrsMatrix.hpp>

using Teuchos::RCP;

namespace esys_trilinos {

template<class Matrix, class Vector>
RCP<DirectSolverType<Matrix,Vector> > createDirectSolver(
                                               const escript::SolverBuddy& sb,
                                               RCP<const Matrix> A,
                                               RCP<Vector> X,
                                               RCP<const Vector> B)
{
    RCP<DirectSolverType<Matrix,Vector> > solver;
    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();

    // TODO: options!
    if (Amesos2::query("MUMPS")) {
        solver = Amesos2::create<Matrix, Vector>("MUMPS", A, X, B);
        if (sb.isVerbose()) {
            params->set("ICNTL(4)", 4);
        }
    } else if (Amesos2::query("klu2")) {
        solver = Amesos2::create<Matrix, Vector>("klu2", A, X, B);
        params->set("DiagPivotThresh", sb.getDiagonalDominanceThreshold());
        params->set("SymmetricMode", sb.isSymmetric());
    } else if (Amesos2::query("superludist")) {
        solver = Amesos2::create<Matrix, Vector>("superludist", A, X, B);
    } else if (Amesos2::query("superlu")) {
        solver = Amesos2::create<Matrix, Vector>("superlu", A, X, B);
        params->set("DiagPivotThresh", sb.getDiagonalDominanceThreshold());
        params->set("ILU_DropTol", sb.getDropTolerance());
        params->set("SymmetricMode", sb.isSymmetric());
    } else if (Amesos2::query("pardiso_mkl")) {
        solver = Amesos2::create<Matrix, Vector>("pardiso_mkl", A, X, B);
    } else if (Amesos2::query("lapack")) {
        solver = Amesos2::create<Matrix, Vector>("lapack", A, X, B);
    } else if (Amesos2::query("amesos2_cholmod")) {
        solver = Amesos2::create<Matrix, Vector>("amesos2_cholmod", A, X, B);
    } else if (Amesos2::query("superlumt")) {
        solver = Amesos2::create<Matrix, Vector>("superlumt", A, X, B);
        params->set("nprocs", omp_get_max_threads());
        params->set("DiagPivotThresh", sb.getDiagonalDominanceThreshold());
        params->set("SymmetricMode", sb.isSymmetric());
    } else {
        throw TrilinosAdapterException("Could not find an Amesos2 direct solver!");
    }
    solver->setParameters(params);
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

