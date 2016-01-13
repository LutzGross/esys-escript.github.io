
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

