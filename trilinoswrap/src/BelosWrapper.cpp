
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <trilinoswrap/BelosWrapper.h>
#include <trilinoswrap/TrilinosAdapterException.h>
#include <trilinoswrap/util.h>

#include <escript/SolverOptions.h>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

#include <boost/python/dict.hpp>

using Teuchos::RCP;

namespace bp = boost::python;

namespace esys_trilinos {

template<typename ST>
RCP<SolverType<ST> > createSolver(const escript::SolverBuddy& sb)
{
    using util::extractParamIfSet;

    Belos::SolverFactory<ST, VectorType<ST>, OpType<ST> > factory;
    RCP<SolverType<ST> > solver;
    RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();

    solverParams->set("Convergence Tolerance", sb.getTolerance());
    solverParams->set("Maximum Iterations", sb.getIterMax());
    if (sb.isVerbose()) {
        solverParams->set("Verbosity", Belos::Errors + Belos::Warnings +
                Belos::TimingDetails + Belos::StatusTestDetails);
    }

    escript::SolverOptions method = sb.getSolverMethod();
    const bp::dict& pyParams = sb.getTrilinosParameters();

    if (method == escript::SO_DEFAULT) {
        if (sb.isSymmetric() || sb.isHermitian()) {
            method = escript::SO_METHOD_PCG;
        } else {
            method = escript::SO_METHOD_GMRES;
        }
    }

    switch (method) {
        case escript::SO_METHOD_BICGSTAB:
            solver = factory.create("BICGSTAB", solverParams);
            break;
        case escript::SO_METHOD_PCG:
            solver = factory.create("CG", solverParams);
            break;
        case escript::SO_METHOD_PRES20:
            //solverParams->set("Num Blocks", 5);
            //solverParams->set("Maximum Restarts", 20);
            solver = factory.create("GMRES", solverParams);
            break;
        case escript::SO_METHOD_ITERATIVE:
        case escript::SO_METHOD_GMRES:
            extractParamIfSet<int>("Num Blocks", pyParams, *solverParams);
            extractParamIfSet<int>("Maximum Restarts", pyParams, *solverParams);
            extractParamIfSet<std::string>("Orthogonalization", pyParams, *solverParams);
            solver = factory.create("GMRES", solverParams);
            break;
        case escript::SO_METHOD_LSQR:
            // Trilinos breaks when using this combination of methods and packages
            // so switch to another method
            //#ifdef ESYS_TRILINOS_14
            if (sb.getSolverMethod() == escript::SO_METHOD_LSQR
                && sb.getPackage() == escript::SO_PACKAGE_TRILINOS
                && sb.getPreconditioner() == escript::SO_PRECONDITIONER_AMG
                )
            {
                extractParamIfSet<int>("Num Blocks", pyParams, *solverParams);
                extractParamIfSet<int>("Maximum Restarts", pyParams, *solverParams);
                extractParamIfSet<std::string>("Orthogonalization", pyParams, *solverParams);
                solver = factory.create("GMRES", solverParams);
                //notification message to user TODO
                break;
            }
            //#endif
            extractParamIfSet<ST>("Condition Limit", pyParams, *solverParams);
            extractParamIfSet<int>("Term Iter Max", pyParams, *solverParams);
            extractParamIfSet<ST>("Lambda", pyParams, *solverParams);
            solverParams->set("Rel Mat Err", sb.getTolerance());
            solver = factory.create("LSQR", solverParams);
            break;
        case escript::SO_METHOD_MINRES:
            extractParamIfSet<int>("Block Size", pyParams, *solverParams);
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

// instantiate our two supported versions
template
RCP<SolverType<real_t> > createSolver(const escript::SolverBuddy& sb);
template
RCP<SolverType<cplx_t> > createSolver(const escript::SolverBuddy& sb);

}  // end of namespace

