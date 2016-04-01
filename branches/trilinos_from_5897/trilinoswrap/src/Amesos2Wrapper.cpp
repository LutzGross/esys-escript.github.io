
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

using Teuchos::RCP;

namespace esys_trilinos {

template<typename ST>
RCP<DirectSolverType<ST> > createDirectSolver(const escript::SolverBuddy& sb,
                                              RCP<const MatrixType<ST> > A,
                                              RCP<VectorType<ST> > X,
                                              RCP<const VectorType<ST> > B)
{
    typedef MatrixType<ST> MT;
    typedef VectorType<ST> VT;

    RCP<DirectSolverType<ST> > solver;
    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();

    // TODO: options!
    if (Amesos2::query("MUMPS")) {
        solver = Amesos2::create<MT, VT>("MUMPS", A, X, B);
        if (sb.isVerbose()) {
            params->set("ICNTL(4)", 4);
        }
    } else if (Amesos2::query("klu2")) {
        solver = Amesos2::create<MT, VT>("klu2", A, X, B);
        params->set("DiagPivotThresh", sb.getDiagonalDominanceThreshold());
        params->set("SymmetricMode", sb.isSymmetric());
    } else if (Amesos2::query("superludist")) {
        solver = Amesos2::create<MT, VT>("superludist", A, X, B);
    } else if (Amesos2::query("superlu")) {
        solver = Amesos2::create<MT, VT>("superlu", A, X, B);
        params->set("DiagPivotThresh", sb.getDiagonalDominanceThreshold());
        params->set("ILU_DropTol", sb.getDropTolerance());
        params->set("SymmetricMode", sb.isSymmetric());
    } else if (Amesos2::query("pardiso_mkl")) {
        solver = Amesos2::create<MT, VT>("pardiso_mkl", A, X, B);
    } else if (Amesos2::query("lapack")) {
        solver = Amesos2::create<MT, VT>("lapack", A, X, B);
    } else if (Amesos2::query("amesos2_cholmod")) {
        solver = Amesos2::create<MT, VT>("amesos2_cholmod", A, X, B);
    } else if (Amesos2::query("superlumt")) {
        solver = Amesos2::create<MT, VT>("superlumt", A, X, B);
        params->set("nprocs", omp_get_max_threads());
        params->set("DiagPivotThresh", sb.getDiagonalDominanceThreshold());
        params->set("SymmetricMode", sb.isSymmetric());
    } else {
        throw TrilinosAdapterException("Could not find an Amesos2 direct solver!");
    }
    solver->setParameters(params);
    return solver;
}

// instantiate our two supported versions
template
RCP<DirectSolverType<real_t> > createDirectSolver<real_t>(
                      const escript::SolverBuddy& sb, RCP<const RealMatrix> A,
                      RCP<RealVector> X, RCP<const RealVector> B);
template
RCP<DirectSolverType<cplx_t> > createDirectSolver<cplx_t>(
                   const escript::SolverBuddy& sb, RCP<const ComplexMatrix> A,
                   RCP<ComplexVector> X, RCP<const ComplexVector> B);

}  // end of namespace

