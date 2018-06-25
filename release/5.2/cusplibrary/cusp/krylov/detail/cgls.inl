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

#include <cusp/array1d.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>
#include <cusp/monitor.h>
#include <cusp/linear_operator.h>

namespace blas = cusp::blas;

namespace cusp
{
namespace krylov
{

template <class LinearOperator,
          typename ValueType,
          class Vector>
void cgls(LinearOperator& A,
          Vector& x,
          Vector& b,
          ValueType shift)
{
    cusp::default_monitor<ValueType> monitor(b);

    cusp::krylov::cgls(A, x, b, shift, monitor);
}

template <class LinearOperator,
          typename ValueType,
          class Vector,
          class Monitor>
void cgls(LinearOperator& A,
          Vector& x,
          Vector& b,
          ValueType shift,
          Monitor& monitor)
{
    CUSP_PROFILE_SCOPED();

    typedef typename LinearOperator::memory_space MemorySpace;

    const size_t N = A.num_rows;

    // allocate workspace
    cusp::array1d<ValueType,MemorySpace> y(N);
    cusp::array1d<ValueType,MemorySpace> z(N);
    cusp::array1d<ValueType,MemorySpace> r(N);
    cusp::array1d<ValueType,MemorySpace> p(N);

    // y <- Ax
    cusp::multiply(A, x, y);

    // r <- b - A*x
    blas::axpby(b, y, r, ValueType(1), ValueType(-1));

    // z <- A^T*r - shift*x
    cusp::transposed_multiply(A, r, z);
    blas::axpy(x, z, -shift);

    // p <- z
    blas::copy(z, p);

    // gamma = <r, r>
    ValueType gamma = blas::dotc(z, z);

    while (!monitor.finished(z))
    {
        // y <- Ap
        cusp::multiply(A, p, y);

        // delta <- <y,y> + shift*<p,p>
        ValueType delta = blas::dotc(y,y) + shift*blas::dotc(p,p);

        // alpha <- <r,r> / <y,y>
        ValueType alpha = gamma / delta;

        // x <- x + alpha * p
        blas::axpy(p, x, alpha);

        // r <- r - alpha * y
        blas::axpy(y, r, -alpha);

        // z <- A^T*r - shift*x
        cusp::transposed_multiply(A, r, z);
        blas::axpy(x, z, -shift);

        ValueType gamma_old = gamma;

        // gamma = <r, r>
        gamma = blas::dotc(z, z);

        // beta <- <r_{i+1},r_{i+1}> / <r,r>
        ValueType beta = gamma / gamma_old;

        // p <- r + beta*p
        blas::axpby(z, p, p, ValueType(1), beta);

        ++monitor;
    }
}

} // end namespace krylov
} // end namespace cusp

