
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

#ifndef __ESYS_BELOSWRAPPER_H__
#define __ESYS_BELOSWRAPPER_H__

#include <trilinoswrap/types.h>

#include <BelosSolverManager.hpp>

namespace escript {
    class SolverBuddy;
}

namespace esys_trilinos {

template<typename ST>
using ProblemType = Belos::LinearProblem< ST, VectorType<ST>, OpType<ST> >;

template<typename ST>
using SolverType = Belos::SolverManager< ST, VectorType<ST>, OpType<ST> >;

template<typename ST>
Teuchos::RCP<SolverType<ST> > createSolver(const escript::SolverBuddy& sb);

} // namespace esys_trilinos

#endif // __ESYS_AMESOS2WRAPPER_H__

