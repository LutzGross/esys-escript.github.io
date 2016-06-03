
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

#include <trilinoswrap/PreconditionerFactory.h>
#include <trilinoswrap/TrilinosAdapterException.h>

#include <escript/SolverOptions.h>

#include <Ifpack2_Factory.hpp>
#ifndef ESYS_INDEXTYPE_LONG
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif

using Teuchos::RCP;

namespace esys_trilinos {

template<typename ST>
RCP<OpType<ST> > createPreconditioner(RCP<const MatrixType<ST> > mat,
                                      const escript::SolverBuddy& sb)
{
    typedef MatrixType<ST> Matrix;

    RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    Ifpack2::Factory factory;
    RCP<OpType<ST> > prec;
    RCP<Ifpack2::Preconditioner<ST,LO,GO,NT> > ifprec;

    // TODO: options!
    switch (sb.getPreconditioner()) {
        case escript::SO_PRECONDITIONER_NONE:
            break;
        case escript::SO_PRECONDITIONER_AMG:
            {
#ifndef ESYS_INDEXTYPE_LONG
                params->set("max levels", sb.getLevelMax());
                params->set("number of equations", 1);
                params->set("cycle type", sb.getCycleType()==1 ? "V" : "W");
                params->set("problem: symmetric", sb.isSymmetric());
                params->set("verbosity", sb.isVerbose()? "high":"none");
                RCP<OpType<ST> > A(Teuchos::rcp_const_cast<Matrix>(mat));
                prec = MueLu::CreateTpetraPreconditioner(A, *params);
#else
                throw escript::ValueError("MueLu (AMG) is incompatible with index type long!");
#endif
            }
            break;
        case escript::SO_PRECONDITIONER_ILUT:
            ifprec = factory.create<const Matrix>("ILUT", mat);
            params->set("fact: drop tolerance", sb.getDropTolerance());
            break;
        case escript::SO_PRECONDITIONER_GAUSS_SEIDEL:
        case escript::SO_PRECONDITIONER_JACOBI:
            ifprec = factory.create<const Matrix>("RELAXATION", mat);
            if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_JACOBI) {
                params->set("relaxation: type", "Jacobi");
            } else {
                params->set("relaxation: type", (sb.isSymmetric() ?
                            "Symmetric Gauss-Seidel" : "Gauss-Seidel"));
            }
            params->set("relaxation: sweeps", sb.getNumSweeps());
            params->set("relaxation: damping factor", sb.getRelaxationFactor());
            //params->set("relaxation: backward mode", true);
            break;
        case escript::SO_PRECONDITIONER_ILU0: // to avoid test failures
        case escript::SO_PRECONDITIONER_RILU:
            if (dynamic_cast<const Tpetra::Experimental::BlockCrsMatrix<ST,LO,GO,NT>* >(mat.get())) {
                ifprec = factory.create<const Matrix>("RBILUK", mat);
            } else {
                ifprec = factory.create<const Matrix>("RILUK", mat);
            }
            params->set("fact: relax value", sb.getRelaxationFactor());
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

