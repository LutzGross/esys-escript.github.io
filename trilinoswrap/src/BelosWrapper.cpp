
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

#include <trilinoswrap/BelosWrapper.h>
#include <trilinoswrap/TrilinosAdapterException.h>

#include <escript/SolverOptions.h>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

using Teuchos::RCP;

namespace esys_trilinos {

template<typename ST>
RCP<SolverType<ST> > createSolver(const escript::SolverBuddy& sb)
{
    Belos::SolverFactory<ST, VectorType<ST>, OpType<ST> > factory;
    RCP<SolverType<ST> > solver;
    RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();

    solverParams->set("Convergence Tolerance", sb.getTolerance());
    solverParams->set("Maximum Iterations", sb.getIterMax());
    if (sb.isVerbose()) {
        solverParams->set("Verbosity", Belos::Errors + Belos::Warnings +
                Belos::TimingDetails + Belos::StatusTestDetails);
    }

    // TODO: options!
    escript::SolverOptions method = sb.getSolverMethod();

    if (method == escript::SO_DEFAULT) {
        if (sb.isSymmetric()) {
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
        case escript::SO_METHOD_GMRES:
            //solverParams->set("Num Blocks", sb.getTruncation());
            //solverParams->set("Maximum Restarts", sb.getRestart());
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

// instantiate our two supported versions
template
RCP<SolverType<real_t> > createSolver(const escript::SolverBuddy& sb);
template
RCP<SolverType<cplx_t> > createSolver(const escript::SolverBuddy& sb);

}  // end of namespace

