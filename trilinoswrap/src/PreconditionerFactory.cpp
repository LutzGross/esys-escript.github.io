
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

#include <trilinoswrap/PreconditionerFactory.h>
#include <trilinoswrap/TrilinosAdapterException.h>
#include <trilinoswrap/util.h>

#include <escript/SolverOptions.h>

#include <Ifpack2_Factory.hpp>
#if 1 //ndef ESYS_INDEXTYPE_LONG
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif

#include <boost/python/dict.hpp>

using Teuchos::RCP;

namespace bp = boost::python;

namespace esys_trilinos {

template<typename ST>
RCP<OpType<ST> > createPreconditioner(RCP<const MatrixType<ST> > mat,
                                      const escript::SolverBuddy& sb)
{
    using util::extractParamIfSet;

    typedef MatrixType<ST> Matrix;

    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    Ifpack2::Factory factory;
    RCP<OpType<ST> > prec;
    RCP<Ifpack2::Preconditioner<ST,LO,GO,NT> > ifprec;
    const bp::dict& pyParams = sb.getTrilinosParameters();

    switch (sb.getPreconditioner()) {
        case escript::SO_PRECONDITIONER_NONE:
            break;
        case escript::SO_PRECONDITIONER_AMG:
            {
#if 1 //ndef ESYS_INDEXTYPE_LONG
                params->set("number of equations", 1);
                params->set("problem: symmetric", sb.isSymmetric() || sb.isHermitian());
                params->set("verbosity", sb.isVerbose()? "high":"none");
                // override parameters if set explicitly for trilinos
                // The set of available parameters is documented in the MueLu
                // user's guide (PDF download)
                // NOTE: passing sub parameter lists is not supported via
                // python due to escript's SolverBuddy constraints.
                // Use the 'xml parameter file' option instead.
                extractParamIfSet<std::string>("problem: type", pyParams, *params);
                extractParamIfSet<std::string>("verbosity", pyParams, *params);
                extractParamIfSet<int>("number of equations", pyParams, *params);
                extractParamIfSet<int>("max levels", pyParams, *params);
                extractParamIfSet<std::string>("cycle type", pyParams, *params);
                extractParamIfSet<bool>("problem: symmetric", pyParams, *params);
                extractParamIfSet<std::string>("xml parameter file", pyParams, *params);
                extractParamIfSet<std::string>("smoother: pre or post", pyParams, *params);
                extractParamIfSet<std::string>("smoother: type", pyParams, *params);
                extractParamIfSet<std::string>("smoother: pre type", pyParams, *params);
                extractParamIfSet<std::string>("smoother: post type", pyParams, *params);
                extractParamIfSet<int>("coarse: max size", pyParams, *params);
                extractParamIfSet<std::string>("coarse: type", pyParams, *params);
                extractParamIfSet<std::string>("aggregation: type", pyParams, *params);
                extractParamIfSet<std::string>("aggregation: ordering", pyParams, *params);
                extractParamIfSet<std::string>("aggregation: drop scheme", pyParams, *params);
                extractParamIfSet<ST>("aggregation: drop tol", pyParams, *params);
                extractParamIfSet<int>("aggregation: min agg size", pyParams, *params);
                extractParamIfSet<int>("aggregation: max agg size", pyParams, *params);
                extractParamIfSet<ST>("aggregation: Dirichlet threshold", pyParams, *params);
                extractParamIfSet<bool>("aggregation: export visualization data", pyParams, *params);
                extractParamIfSet<std::string>("aggregation: output filename", pyParams, *params);
                extractParamIfSet<int>("aggregation: output file: time step", pyParams, *params);
                extractParamIfSet<int>("aggregation: output file: iter", pyParams, *params);
                extractParamIfSet<std::string>("aggregation: output file: agg style", pyParams, *params);
                extractParamIfSet<bool>("aggregation: output file: fine graph edges", pyParams, *params);
                extractParamIfSet<bool>("aggregation: output file: coarse graph edges", pyParams, *params);
                extractParamIfSet<bool>("aggregation: output file: build colormap", pyParams, *params);
                extractParamIfSet<bool>("repartition: enable", pyParams, *params);
                extractParamIfSet<std::string>("repartition: partitioner", pyParams, *params);
                extractParamIfSet<int>("repartition: start level", pyParams, *params);
                extractParamIfSet<int>("repartition: min rows per proc", pyParams, *params);
                extractParamIfSet<ST>("repartition: max imbalance", pyParams, *params);
                extractParamIfSet<bool>("repartition: remap parts", pyParams, *params);
                extractParamIfSet<bool>("repartition: rebalance P and R", pyParams, *params);
                extractParamIfSet<std::string>("multigrid algorithm", pyParams, *params);
                extractParamIfSet<int>("semicoarsen: coarsen rate", pyParams, *params);
                extractParamIfSet<ST>("sa: damping factor", pyParams, *params);
                extractParamIfSet<bool>("sa: use filtered matrix", pyParams, *params);
                extractParamIfSet<bool>("filtered matrix: use lumping", pyParams, *params);
                extractParamIfSet<bool>("filtered matrix: reuse eigenvalue", pyParams, *params);
                extractParamIfSet<std::string>("emin: iterative method", pyParams, *params);
                extractParamIfSet<int>("emin: num iterations", pyParams, *params);
                extractParamIfSet<int>("emin: num reuse iterations", pyParams, *params);
                extractParamIfSet<std::string>("emin: pattern", pyParams, *params);
                extractParamIfSet<int>("emin: pattern order", pyParams, *params);
                extractParamIfSet<std::string>("reuse: type", pyParams, *params);
                extractParamIfSet<bool>("print initial parameters", pyParams, *params);
                extractParamIfSet<bool>("print unused parameters", pyParams, *params);
                extractParamIfSet<bool>("transpose: use implicit", pyParams, *params);
                RCP<OpType<ST> > A(Teuchos::rcp_const_cast<Matrix>(mat));
                prec = MueLu::CreateTpetraPreconditioner(A, *params);
#else
                throw escript::ValueError("MueLu (AMG) is incompatible with index type long!");
#endif
            }
            break;
        case escript::SO_PRECONDITIONER_ILUT:
        {
            ifprec = factory.create<const Matrix>("ILUT", mat);
            params->set("fact: drop tolerance", sb.getDropTolerance());
            params->set("fact: relax value", sb.getRelaxationFactor());
            // override if set explicitly for trilinos
            extractParamIfSet<ST>("fact: relax value", pyParams, *params);
            extractParamIfSet<ST>("fact: drop tolerance", pyParams, *params);
            extractParamIfSet<int>("fact: ilut level-of-fill", pyParams, *params);
            extractParamIfSet<ST>("fact: absolute threshold", pyParams, *params);
            extractParamIfSet<ST>("fact: relative threshold", pyParams, *params);
        }
            break;
        case escript::SO_PRECONDITIONER_GAUSS_SEIDEL:
        {
            ifprec = factory.create<const Matrix>("RELAXATION", mat);
            params->set("relaxation: type", ((sb.isSymmetric() || sb.isHermitian())?
                            "Symmetric Gauss-Seidel" : "Gauss-Seidel"));
            params->set("relaxation: sweeps", sb.getNumSweeps());
            const ST fac = static_cast<ST>(sb.getRelaxationFactor());
            params->set("relaxation: damping factor", fac);
            // override if set explicitly for trilinos
            extractParamIfSet<int>("relaxation: sweeps", pyParams, *params);
            extractParamIfSet<ST>("relaxation: damping factor", pyParams, *params);
            extractParamIfSet<bool>("relaxation: backward mode", pyParams, *params);
            extractParamIfSet<bool>("relaxation: use l1", pyParams, *params);
            extractParamIfSet<ST>("relaxation: l1 eta", pyParams, *params);
            extractParamIfSet<bool>("relaxation: zero starting solution", pyParams, *params);
            extractParamIfSet<bool>("relaxation: fix tiny diagonal entries", pyParams, *params);
            extractParamIfSet<ST>("relaxation: min diagonal value", pyParams, *params);
            extractParamIfSet<bool>("relaxation: check diagonal entries", pyParams, *params);
            // extractParamIfSet<ST>("relaxation: local smoothing indices", pyParams, *params);
        }
            break;
        case escript::SO_PRECONDITIONER_JACOBI:
        {
            ifprec = factory.create<const Matrix>("RELAXATION", mat);
            params->set("relaxation: type", "Jacobi");
            params->set("relaxation: sweeps", sb.getNumSweeps());
            const ST fac = static_cast<ST>(sb.getRelaxationFactor());
            params->set("relaxation: damping factor", fac);
            // override if set explicitly for trilinos
            extractParamIfSet<int>("relaxation: sweeps", pyParams, *params);
            extractParamIfSet<ST>("relaxation: damping factor", pyParams, *params);
            extractParamIfSet<bool>("relaxation: backward mode", pyParams, *params);
            extractParamIfSet<bool>("relaxation: use l1", pyParams, *params);
            extractParamIfSet<ST>("relaxation: l1 eta", pyParams, *params);
            extractParamIfSet<bool>("relaxation: zero starting solution", pyParams, *params);
            extractParamIfSet<bool>("relaxation: fix tiny diagonal entries", pyParams, *params);
            extractParamIfSet<ST>("relaxation: min diagonal value", pyParams, *params);
            extractParamIfSet<bool>("relaxation: check diagonal entries", pyParams, *params);
            // extractParamIfSet<ST>("relaxation: local smoothing indices", pyParams, *params);
        }
            break;
        case escript::SO_PRECONDITIONER_ILU0: // to avoid test failures
        case escript::SO_PRECONDITIONER_RILU:
        {
#if defined(ESYS_HAVE_TPETRA_EXPERIMENTAL_BLOCKCRSH)
            if (dynamic_cast<const Tpetra::Experimental::BlockCrsMatrix<ST,LO,GO,NT>* >(mat.get())) {
#else
            if (dynamic_cast<const Tpetra::BlockCrsMatrix<ST,LO,GO,NT>* >(mat.get())) {
#endif
                ifprec = factory.create<const Matrix>("RBILUK", mat);
            } else {
                ifprec = factory.create<const Matrix>("RILUK", mat);
            }
            params->set("fact: relax value", sb.getRelaxationFactor());
            // override if set explicitly for trilinos
            extractParamIfSet<ST>("fact: relax value", pyParams, *params);
            extractParamIfSet<int>("fact: iluk level-of-fill", pyParams, *params);
            extractParamIfSet<int>("fact: iluk level-of-overlap", pyParams, *params);
            extractParamIfSet<ST>("fact: absolute threshold", pyParams, *params);
            extractParamIfSet<ST>("fact: relative threshold", pyParams, *params);
        }
            break;
        default:
            throw escript::ValueError("Unsupported preconditioner requested.");
    }
    if (!ifprec.is_null()) {
        ifprec->setParameters(*params);
        ifprec->initialize();
        ifprec->compute();
        prec = ifprec;
    }
    return prec;
}

// instantiate our two supported versions
typedef MatrixType<real_t> RealMatrix;
typedef MatrixType<cplx_t> ComplexMatrix;

template
RCP<RealOperator> createPreconditioner<real_t>(RCP<const RealMatrix> mat,
                                               const escript::SolverBuddy& sb);
template
RCP<ComplexOperator> createPreconditioner<cplx_t>(RCP<const ComplexMatrix> mat,
                                               const escript::SolverBuddy& sb);

}  // end of namespace
