
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __ESYS_AMESOS2WRAPPER_H__
#define __ESYS_AMESOS2WRAPPER_H__

#include <trilinoswrap/types.h>

#include <Amesos2_Solver_decl.hpp>

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

template<typename ST>
using DirectSolverType = Amesos2::Solver< MatrixType<ST>, VectorType<ST> >;

template<typename ST>
Teuchos::RCP<DirectSolverType<ST> > createDirectSolver(
                                  const escript::SolverBuddy& sb,
                                  Teuchos::RCP<const MatrixType<ST> > A,
                                  Teuchos::RCP<VectorType<ST> > X,
                                  Teuchos::RCP<const VectorType<ST> > B);

} // namespace esys_trilinos

#endif // __ESYS_AMESOS2WRAPPER_H__

