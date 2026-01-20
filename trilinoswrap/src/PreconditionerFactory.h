
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

#ifndef __ESYS_TRILINOS_PRECONDITIONERFACTORY_H__
#define __ESYS_TRILINOS_PRECONDITIONERFACTORY_H__

#include <trilinoswrap/types.h>

/// Wrapper for Ifpack2 and MueLu

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

/// creates a preconditioner (Operator) for input matrix A using options in sb.
/// ST is the scalar type used by the matrix.
template<typename ST>
Teuchos::RCP<OpType<ST> > createPreconditioner(
                                      Teuchos::RCP<const MatrixType<ST> > A,
                                      const escript::SolverBuddy& sb);

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_PRECONDITIONERFACTORY_H__

