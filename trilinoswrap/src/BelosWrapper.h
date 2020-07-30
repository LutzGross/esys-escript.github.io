
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

#ifndef __ESYS_TRILINOS_BELOSWRAPPER_H__
#define __ESYS_TRILINOS_BELOSWRAPPER_H__

#include <trilinoswrap/types.h>

#include <BelosSolverManager.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosTypes.hpp>

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

#endif // __ESYS_TRILINOS_AMESOS2WRAPPER_H__

