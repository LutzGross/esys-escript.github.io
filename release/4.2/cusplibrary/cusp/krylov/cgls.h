/*
 *  Copyright 2014-2015 The University of Queensland
 *  http://www.uq.edu.au
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

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

