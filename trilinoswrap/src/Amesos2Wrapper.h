
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

#ifndef __ESYS_TRILINOS_AMESOS2WRAPPER_H__
#define __ESYS_TRILINOS_AMESOS2WRAPPER_H__

#include <trilinoswrap/types.h>

#include <Amesos2_Solver_decl.hpp>

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

template<class Matrix, class Vector>
using DirectSolverType = Amesos2::Solver<Matrix, Vector>;

template<class Matrix, class Vector>
Teuchos::RCP<DirectSolverType<Matrix,Vector> > createDirectSolver(
                                  const escript::SolverBuddy& sb,
                                  Teuchos::RCP<const Matrix> A,
                                  Teuchos::RCP<Vector> X,
                                  Teuchos::RCP<const Vector> B);

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_AMESOS2WRAPPER_H__

