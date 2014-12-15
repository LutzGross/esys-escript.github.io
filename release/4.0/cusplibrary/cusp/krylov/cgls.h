
/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
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

/*! \file cgls.h
 *  \brief Conjugate Gradient with Least Squares (CGLS) method
 */

#pragma once

#include <cusp/detail/config.h>

namespace cusp
{
namespace krylov
{

/*! \addtogroup iterative_solvers Iterative Solvers
 *  \addtogroup krylov_methods Krylov Methods
 *  \ingroup iterative_solvers
 *  \{
 */

/*! \p cgls : Conjugate Gradient with Least Squares method
 *
 * Solves the linear system A x = b
 * using the default convergence criteria.
 */
template <class LinearOperator,
          typename ValueType,
          class Vector>
void cgls(LinearOperator& A,
          Vector& x,
          Vector& b,
          ValueType shift);

/*! \p cgls : Conjugate Gradient with Least Squares method
 */
template <class LinearOperator,
          typename ValueType,
          class Vector,
          class Monitor>
void cgls(LinearOperator& A,
          Vector& x,
          Vector& b,
          ValueType shift,
          Monitor& monitor);

} // end namespace krylov
} // end namespace cusp

#include <cusp/krylov/detail/cgls.inl>

