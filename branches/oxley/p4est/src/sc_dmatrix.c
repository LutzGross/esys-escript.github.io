/*
  This file is part of the SC Library.
  The SC Library provides support for parallel scientific applications.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors

  The SC Library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  The SC Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with the SC Library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
*/

#include <sc_dmatrix.h>
#include <sc_lapack.h>

int
sc_darray_is_valid (const double *darray, size_t nelem)
{
  size_t              zz;

  for (zz = 0; zz < nelem; ++zz) {
    if (darray[zz] != darray[zz]) {     /* ignore the comparison warning */
      return 0;
    }
  }

  return 1;
}

int
sc_darray_is_range (const double *darray, size_t nelem,
                    double low, double high)
{
  size_t              zz;

  for (zz = 0; zz < nelem; ++zz) {
    if (!(low <= darray[zz] && darray[zz] <= high)) {
      return 0;
    }
  }

  return 1;
}

size_t
sc_dmatrix_memory_used (sc_dmatrix_t * dm)
{
  size_t              mem = sizeof (sc_dmatrix_t);

  mem += (dm->m + 1) * sizeof (double *);
  if (!dm->view) {
    mem += dm->m * dm->n * sizeof (double);
  }

  return mem;
}

static void
sc_dmatrix_new_e (sc_dmatrix_t * rdm, sc_bint_t m, sc_bint_t n, double *data)
{
  sc_bint_t           i;

  SC_ASSERT (m >= 0 && n >= 0);
  SC_ASSERT (rdm != NULL);

  rdm->e = SC_ALLOC (double *, m + 1);
  rdm->e[0] = data;

  if (m > 0) {
    for (i = 1; i < m; ++i)
      rdm->e[i] = rdm->e[i - 1] + n;

    rdm->e[m] = NULL;           /* safeguard */
  }

  rdm->m = m;
  rdm->n = n;
}

static sc_dmatrix_t *
sc_dmatrix_new_internal (sc_bint_t m, sc_bint_t n, int init_zero)
{
  sc_dmatrix_t       *rdm;
  double             *data;
  size_t              size = (size_t) (m * n);
#ifdef SC_ENABLE_DEBUG
  double              zero = 0.0;       /* no const to avoid warning */
  const double        anan = 0.0 / zero;
  size_t              zz;
#endif

  SC_ASSERT (m >= 0 && n >= 0);

  rdm = SC_ALLOC (sc_dmatrix_t, 1);

  if (init_zero) {
    data = SC_ALLOC_ZERO (double, size);
  }
  else {
    data = SC_ALLOC (double, size);
#ifdef SC_ENABLE_DEBUG
    /* In debug mode initialize the memory to NaN. */
    for (zz = 0; zz < size; ++zz) {
      data[zz] = anan;
    }
#endif
  }

  sc_dmatrix_new_e (rdm, m, n, data);
  rdm->view = 0;

  return rdm;
}

sc_dmatrix_t       *
sc_dmatrix_new (sc_bint_t m, sc_bint_t n)
{
  return sc_dmatrix_new_internal (m, n, 0);
}

sc_dmatrix_t       *
sc_dmatrix_new_zero (sc_bint_t m, sc_bint_t n)
{
  return sc_dmatrix_new_internal (m, n, 1);
}

sc_dmatrix_t       *
sc_dmatrix_new_data (sc_bint_t m, sc_bint_t n, double *data)
{
  sc_dmatrix_t       *rdm;

  SC_ASSERT (m >= 0 && n >= 0);

  rdm = SC_ALLOC (sc_dmatrix_t, 1);
  sc_dmatrix_new_e (rdm, m, n, data);
  rdm->view = 1;

  return rdm;
}

sc_dmatrix_t       *
sc_dmatrix_new_view (sc_bint_t m, sc_bint_t n, sc_dmatrix_t * orig)
{
  return sc_dmatrix_new_view_offset (0, m, n, orig);
}

sc_dmatrix_t       *
sc_dmatrix_new_view_offset (sc_bint_t o, sc_bint_t m, sc_bint_t n,
                            sc_dmatrix_t * orig)
{
  sc_dmatrix_t       *rdm;

  SC_ASSERT (o >= 0 && m >= 0 && n >= 0);
  SC_ASSERT ((o + m) * n <= orig->m * orig->n);

  rdm = SC_ALLOC (sc_dmatrix_t, 1);
  sc_dmatrix_new_e (rdm, m, n, orig->e[0] + o * n);
  rdm->view = 1;

  return rdm;
}

sc_dmatrix_t       *
sc_dmatrix_new_view_column (sc_dmatrix_t * orig, sc_bint_t j)
{
  sc_dmatrix_t       *rdm;

  SC_ASSERT (orig->m >= 0);
  SC_ASSERT (0 <= j && j < orig->n);

  rdm = SC_ALLOC (sc_dmatrix_t, 1);
  sc_dmatrix_new_e (rdm, orig->m, orig->n, orig->e[0] + j);
  rdm->n = 1;
  rdm->view = 1;

  return rdm;
}

void
sc_dmatrix_view_set_column (sc_dmatrix_t * view,
                            sc_dmatrix_t * orig, sc_bint_t j)
{
  const sc_bint_t     m = view->m;
  sc_bint_t           i;

  SC_ASSERT (view->view);
  SC_ASSERT (view->m == orig->m);
  SC_ASSERT (orig->m >= 0);
  SC_ASSERT (0 <= j && j < orig->n);

  view->e[0] = orig->e[0] + j;

  if (m > 0) {
    for (i = 1; i < m; ++i)
      view->e[i] = view->e[i - 1] + orig->n;

    view->e[m] = NULL;          /* safeguard */
  }

  view->n = 1;
}

void
sc_dmatrix_view_set_row (sc_dmatrix_t * view,
                         sc_dmatrix_t * orig, sc_bint_t i)
{
  SC_ASSERT (view->view);
  SC_ASSERT (view->m == 1);
  SC_ASSERT (orig->n >= 0);
  SC_ASSERT (0 <= i && i < orig->m);

  view->e[0] = orig->e[i];
  view->n = orig->n;
}

sc_dmatrix_t       *
sc_dmatrix_clone (const sc_dmatrix_t * X)
{
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Xdata = X->e[0];
  double             *Ydata;
  sc_dmatrix_t       *clone;

  clone = sc_dmatrix_new (X->m, X->n);
  Ydata = clone->e[0];

  memcpy (Ydata, Xdata, totalsize * sizeof (double));

  return clone;
}

void
sc_dmatrix_reshape (sc_dmatrix_t * dmatrix, sc_bint_t m, sc_bint_t n)
{
  double             *data;

  SC_ASSERT (dmatrix->e != NULL);
  SC_ASSERT (dmatrix->m * dmatrix->n == m * n);

  data = dmatrix->e[0];
  SC_FREE (dmatrix->e);
  sc_dmatrix_new_e (dmatrix, m, n, data);
}

void
sc_dmatrix_resize (sc_dmatrix_t * dmatrix, sc_bint_t m, sc_bint_t n)
{
  double             *data;
  sc_bint_t           size, newsize;

  SC_ASSERT (dmatrix->e != NULL);
  SC_ASSERT (m >= 0 && n >= 0);

  size = dmatrix->m * dmatrix->n;
  newsize = m * n;

  if (!dmatrix->view && size != newsize) {
#ifdef SC_ENABLE_USE_REALLOC
    data = SC_REALLOC (dmatrix->e[0], double, newsize);
#else
    data = SC_ALLOC (double, newsize);
    memcpy (data, dmatrix->e[0],
            (size_t) SC_MIN (newsize, size) * sizeof (double));
    SC_FREE (dmatrix->e[0]);
#endif
  }
  else {
    /* for views you must know that data is large enough */
    data = dmatrix->e[0];
  }
  SC_FREE (dmatrix->e);
  sc_dmatrix_new_e (dmatrix, m, n, data);
}

void
sc_dmatrix_resize_in_place (sc_dmatrix_t * dmatrix, sc_bint_t m, sc_bint_t n)
{
  double             *data;
  sc_bint_t           size, newsize;
  sc_bint_t           i;
  sc_bint_t           old_n = dmatrix->n;
  sc_bint_t           min_m = SC_MIN (m, dmatrix->m);

  SC_ASSERT (dmatrix->e != NULL);
  SC_ASSERT (m >= 0 && n >= 0);
  SC_ASSERT (!dmatrix->view);

  size = dmatrix->m * dmatrix->n;
  newsize = m * n;
  data = dmatrix->e[0];
  if (n < old_n) {
    for (i = 1; i < min_m; i++) {
      memmove (data + i * n, data + i * old_n, n * sizeof (double));
    }
  }
  if (newsize != size) {
#ifdef SC_ENABLE_USE_REALLOC
    data = SC_REALLOC (dmatrix->e[0], double, newsize);
#else
    data = SC_ALLOC (double, newsize);
    memcpy (data, dmatrix->e[0],
            (size_t) SC_MIN (newsize, size) * sizeof (double));
    SC_FREE (dmatrix->e[0]);
#endif
  }
  if (n > old_n) {
    for (i = min_m - 1; i > 0; i--) {
      memmove (data + i * n, data + i * old_n, old_n * sizeof (double));
    }
  }
  SC_FREE (dmatrix->e);
  sc_dmatrix_new_e (dmatrix, m, n, data);
}

void
sc_dmatrix_destroy (sc_dmatrix_t * dmatrix)
{
  if (!dmatrix->view) {
    SC_FREE (dmatrix->e[0]);
  }
  SC_FREE (dmatrix->e);

  SC_FREE (dmatrix);
}

int
sc_dmatrix_is_valid (const sc_dmatrix_t * A)
{
  return sc_darray_is_valid (A->e[0], (size_t) A->m * (size_t) A->n);
}

int
sc_dmatrix_is_symmetric (const sc_dmatrix_t * A, double tolerance)
{
  sc_bint_t           i, j;
  double              diff;

  SC_ASSERT (A->m == A->n);

  for (i = 0; i < A->n; ++i) {
    for (j = i + 1; j < A->n; ++j) {
      diff = fabs (A->e[i][j] - A->e[j][i]);
      if (diff > tolerance) {
        SC_LDEBUGF ("sc_dmatrix not symmetric by %g\n", diff);

        return 0;
      }
    }
  }

  return 1;
}

void
sc_dmatrix_set_zero (sc_dmatrix_t * X)
{
  sc_dmatrix_set_value (X, 0.0);
}

void
sc_dmatrix_set_value (sc_dmatrix_t * X, double value)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  double             *data = X->e[0];

  for (i = 0; i < totalsize; ++i)
    data[i] = value;
}

void
sc_dmatrix_scale (double alpha, sc_dmatrix_t * X)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  double             *Xdata = X->e[0];

  for (i = 0; i < totalsize; ++i)
    Xdata[i] *= alpha;
}

void
sc_dmatrix_shift (double alpha, sc_dmatrix_t * X)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  double             *Xdata = X->e[0];

  for (i = 0; i < totalsize; ++i)
    Xdata[i] += alpha;
}

void
sc_dmatrix_scale_shift (double alpha, double beta, sc_dmatrix_t * X)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  double             *Xdata = X->e[0];

  for (i = 0; i < totalsize; ++i)
    Xdata[i] = alpha * Xdata[i] + beta;
}

void
sc_dmatrix_alphadivide (double alpha, sc_dmatrix_t * X)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  double             *Xdata = X->e[0];

  for (i = 0; i < totalsize; ++i)
    Xdata[i] = alpha / Xdata[i];
}

void
sc_dmatrix_pow (double alpha, sc_dmatrix_t * X)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  double             *Xdata = X->e[0];

  for (i = 0; i < totalsize; ++i)
    Xdata[i] = pow (Xdata[i], alpha);
}

void
sc_dmatrix_fabs (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    Ydata[i] = fabs (Xdata[i]);
}

void
sc_dmatrix_sqrt (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    Ydata[i] = sqrt (Xdata[i]);
}

void
sc_dmatrix_getsign (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *indata = X->e[0];
  double             *outdata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    outdata[i] = (indata[i] >= 0. ? 1 : -1);
}

void
sc_dmatrix_greaterequal (const sc_dmatrix_t * X, double bound,
                         sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *indata = X->e[0];
  double             *outdata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    outdata[i] = (indata[i] >= bound ? 1 : 0);
}

void
sc_dmatrix_lessequal (const sc_dmatrix_t * X, double bound, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *indata = X->e[0];
  double             *outdata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    outdata[i] = (indata[i] <= bound ? 1 : 0);
}

void
sc_dmatrix_maximum (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *indata = X->e[0];
  double             *outdata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    outdata[i] = SC_MAX (indata[i], outdata[i]);
}

void
sc_dmatrix_minimum (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *indata = X->e[0];
  double             *outdata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    outdata[i] = SC_MIN (indata[i], outdata[i]);
}

void
sc_dmatrix_dotmultiply (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    Ydata[i] *= Xdata[i];
}

void
sc_dmatrix_dotdivide (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i)
    Ydata[i] /= Xdata[i];
}

void
sc_dmatrix_dotmultiply_add (const sc_dmatrix_t * A, const sc_dmatrix_t * X,
                            sc_dmatrix_t * Y)
{
  sc_bint_t           i;
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Adata = A->e[0];
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];

  SC_ASSERT (X->m == A->m && X->n == A->n);
  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  for (i = 0; i < totalsize; ++i) {
    Ydata[i] += Adata[i] * Xdata[i];
  }
}

void
sc_dmatrix_copy (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  const sc_bint_t     totalsize = X->m * X->n;
  const double       *Xdata = X->e[0];
  double             *Ydata = Y->e[0];

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  memmove (Ydata, Xdata, totalsize * sizeof (double));
}

void
sc_dmatrix_transpose (const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           i, j, Xrows, Xcols, Xstride, Ystride;
  double             *Ydata;
  const double       *Xdata = X->e[0];

  SC_ASSERT (X->m == Y->n && X->n == Y->m);

  Xrows = X->m;
  Xcols = X->n;
  Xstride = X->n;
  Ystride = Y->n;
  Ydata = Y->e[0];

  for (i = 0; i < Xrows; i++) {
    for (j = 0; j < Xcols; j++) {
      Ydata[j * Ystride + i] = Xdata[i * Xstride + j];
    }
  }
}

void
sc_dmatrix_add (double alpha, const sc_dmatrix_t * X, sc_dmatrix_t * Y)
{
  sc_bint_t           totalsize, inc;

  SC_ASSERT (X->m == Y->m && X->n == Y->n);

  totalsize = X->m * X->n;

  inc = 1;
  if (totalsize > 0) {
    SC_BLAS_DAXPY (&totalsize, &alpha, X->e[0], &inc, Y->e[0], &inc);
  }
}

void
sc_dmatrix_vector (sc_trans_t transa, sc_trans_t transx, sc_trans_t transy,
                   double alpha, const sc_dmatrix_t * A,
                   const sc_dmatrix_t * X, double beta, sc_dmatrix_t * Y)
{
  sc_bint_t           inc = 1;

#ifdef SC_ENABLE_DEBUG
  sc_bint_t           dimX = (transx == SC_NO_TRANS) ? X->m : X->n;
  sc_bint_t           dimY = (transy == SC_NO_TRANS) ? Y->m : Y->n;
  sc_bint_t           dimX1 = (transx == SC_NO_TRANS) ? X->n : X->m;
  sc_bint_t           dimY1 = (transy == SC_NO_TRANS) ? Y->n : Y->m;

  sc_bint_t           Arows = (transa == SC_NO_TRANS) ? A->m : A->n;
  sc_bint_t           Acols = (transa == SC_NO_TRANS) ? A->n : A->m;
#endif

  SC_ASSERT (Acols == dimX && Arows == dimY);
  SC_ASSERT (dimX1 == 1 && dimY1 == 1);

  if (A->n > 0 && A->m > 0) {
    SC_BLAS_DGEMV (&sc_antitranschar[transa], &A->n, &A->m, &alpha,
                   A->e[0], &A->n, X->e[0], &inc, &beta, Y->e[0], &inc);
  }
  else if (beta != 1.) {
    sc_dmatrix_scale (beta, Y);
  }
}

void
sc_dmatrix_multiply (sc_trans_t transa, sc_trans_t transb, double alpha,
                     const sc_dmatrix_t * A, const sc_dmatrix_t * B,
                     double beta, sc_dmatrix_t * C)
{
  sc_bint_t           Acols, Crows, Ccols;
#ifdef SC_ENABLE_DEBUG
  sc_bint_t           Arows, Brows, Bcols;

  Arows = (transa == SC_NO_TRANS) ? A->m : A->n;
  Brows = (transb == SC_NO_TRANS) ? B->m : B->n;
  Bcols = (transb == SC_NO_TRANS) ? B->n : B->m;
#endif

  Acols = (transa == SC_NO_TRANS) ? A->n : A->m;
  Crows = C->m;
  Ccols = C->n;

  SC_ASSERT (Acols == Brows && Arows == Crows && Bcols == Ccols);
  SC_ASSERT (transa == SC_NO_TRANS || transa == SC_TRANS);
  SC_ASSERT (transb == SC_NO_TRANS || transb == SC_TRANS);

  if (Crows > 0 && Ccols > 0) {
    if (Acols > 0) {
      SC_BLAS_DGEMM (&sc_transchar[transb], &sc_transchar[transa], &Ccols,
                     &Crows, &Acols, &alpha, B->e[0], &B->n, A->e[0], &A->n,
                     &beta, C->e[0], &C->n);
    }
    else if (beta != 1.0) {     /* ignore comparison warning */
      sc_dmatrix_scale (beta, C);
    }
  }
}

void
sc_dmatrix_ldivide (sc_trans_t transa, const sc_dmatrix_t * A,
                    const sc_dmatrix_t * B, sc_dmatrix_t * C)
{
  sc_dmatrix_t       *BT;
  sc_trans_t          invtransa =
    (transa == SC_NO_TRANS) ? SC_TRANS : SC_NO_TRANS;

#ifdef SC_ENABLE_DEBUG
  sc_bint_t           A_nrows = (transa == SC_NO_TRANS) ? A->m : A->n;
  sc_bint_t           A_ncols = (transa == SC_NO_TRANS) ? A->n : A->m;
  sc_bint_t           B_nrows = B->m;
  sc_bint_t           B_ncols = B->n;
  sc_bint_t           C_nrows = C->m;
  sc_bint_t           C_ncols = C->n;
#endif

  SC_ASSERT ((C_nrows == A_ncols) && (B_nrows == A_nrows)
             && (B_ncols == C_ncols));

  BT = sc_dmatrix_new (B->n, B->m);
  sc_dmatrix_transpose (B, BT);

  sc_dmatrix_rdivide (invtransa, BT, A, BT);

  sc_dmatrix_transpose (BT, C);

  sc_dmatrix_destroy (BT);
}

void
sc_dmatrix_rdivide (sc_trans_t transb, const sc_dmatrix_t * A,
                    const sc_dmatrix_t * B, sc_dmatrix_t * C)
{
  sc_bint_t           A_nrows = A->m;
  sc_bint_t           B_nrows = (transb == SC_NO_TRANS) ? B->m : B->n;
  sc_bint_t           B_ncols = (transb == SC_NO_TRANS) ? B->n : B->m;
#ifdef SC_ENABLE_DEBUG
  sc_bint_t           A_ncols = A->n;
  sc_bint_t           C_nrows = C->m;
  sc_bint_t           C_ncols = C->n;
#endif
  sc_bint_t           M = B_ncols, N = B_nrows, Nrhs = A_nrows, info = 0;

  SC_ASSERT ((C_nrows == A_nrows) && (B_nrows == C_ncols)
             && (B_ncols == A_ncols));
  SC_ASSERT (N > 0 && Nrhs > 0);

  if (M == N) {
    sc_dmatrix_t       *lu = sc_dmatrix_clone (B);
    sc_bint_t          *ipiv = SC_ALLOC (sc_bint_t, N);

    /* Perform an LU factorization of B. */
    SC_LAPACK_DGETRF (&N, &N, lu->e[0], &N, ipiv, &info);

    SC_CHECK_ABORT (info == 0, "Lapack routine DGETRF failed");

    /* Solve the linear system. */
    sc_dmatrix_copy (A, C);
    SC_LAPACK_DGETRS (&sc_transchar[transb], &N, &Nrhs, lu->e[0], &N,
                      ipiv, C->e[0], &N, &info);

    SC_CHECK_ABORT (info == 0, "Lapack routine DGETRS failed");

    SC_FREE (ipiv);
    sc_dmatrix_destroy (lu);
  }
  else {
    SC_CHECK_ABORT (0, "Only square A's work right now\n");
  }
}

void
sc_dmatrix_solve_transpose_inplace (sc_dmatrix_t * A, sc_dmatrix_t * B)
{
  const sc_bint_t     N = A->m;
  const sc_bint_t     nrhs = B->m;
  sc_bint_t          *ipiv, info;

  SC_ASSERT (A->n == N && B->n == N);

  ipiv = SC_ALLOC (sc_bint_t, N);
  SC_LAPACK_DGESV (&N, &nrhs, A->e[0], &N, ipiv, B->e[0], &N, &info);
  SC_FREE (ipiv);

  SC_CHECK_ABORT (info == 0, "Lapack routine DGESV failed");
}

void
sc_dmatrix_write (const sc_dmatrix_t * dmatrix, FILE * fp)
{
  sc_bint_t           i, j, m, n;

  m = dmatrix->m;
  n = dmatrix->n;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      fprintf (fp, " %16.8e", dmatrix->e[i][j]);
    }
    fprintf (fp, "\n");
  }
}

sc_dmatrix_pool_t  *
sc_dmatrix_pool_new (sc_bint_t m, sc_bint_t n)
{
  sc_dmatrix_pool_t  *dmpool;

  SC_ASSERT (m >= 0 && n >= 0);

  dmpool = SC_ALLOC (sc_dmatrix_pool_t, 1);

  dmpool->m = m;
  dmpool->n = n;
  dmpool->elem_count = 0;
  sc_array_init (&dmpool->freed, sizeof (sc_dmatrix_t *));

  return dmpool;
}

void
sc_dmatrix_pool_destroy (sc_dmatrix_pool_t * dmpool)
{
  size_t              zz;
  sc_dmatrix_t      **pdm;

  SC_ASSERT (dmpool->elem_count == 0);

  for (zz = 0; zz < dmpool->freed.elem_count; ++zz) {
    pdm = (sc_dmatrix_t **) sc_array_index (&dmpool->freed, zz);
    sc_dmatrix_destroy (*pdm);
  }
  sc_array_reset (&dmpool->freed);

  SC_FREE (dmpool);
}

sc_dmatrix_t       *
sc_dmatrix_pool_alloc (sc_dmatrix_pool_t * dmpool)
{
  sc_dmatrix_t       *dm;

  ++dmpool->elem_count;

  if (dmpool->freed.elem_count > 0) {
    dm = *(sc_dmatrix_t **) sc_array_pop (&dmpool->freed);
  }
  else {
    dm = sc_dmatrix_new (dmpool->m, dmpool->n);
  }

#ifdef SC_ENABLE_DEBUG
  sc_dmatrix_set_value (dm, -1.);
#endif

  return dm;
}

void
sc_dmatrix_pool_free (sc_dmatrix_pool_t * dmpool, sc_dmatrix_t * dm)
{
  SC_ASSERT (dmpool->elem_count > 0);
  SC_ASSERT (dm->m == dmpool->m && dm->n == dmpool->n);

  --dmpool->elem_count;

  *(sc_dmatrix_t **) sc_array_push (&dmpool->freed) = dm;
}

sc_darray_work_t   *
sc_darray_work_new (const int n_threads, const int n_blocks,
                    const int n_entries, const int alignment_bytes)
{
  const int           align_dbl = alignment_bytes / 8;
  const int           n_entries_aligned = SC_ALIGN_UP (n_entries, align_dbl);
  sc_darray_work_t   *work;

  SC_ASSERT (0 < n_threads);
  SC_ASSERT (0 < n_blocks);
  SC_ASSERT (alignment_bytes <= 0 || (alignment_bytes % 8) == 0);

  work = SC_ALLOC (sc_darray_work_t, 1);

  work->data = SC_ALLOC (double, n_threads * n_blocks * n_entries_aligned);
  work->n_threads = n_threads;
  work->n_blocks = n_blocks;
  work->n_entries = n_entries_aligned;

  return work;
}

void
sc_darray_work_destroy (sc_darray_work_t * work)
{
  SC_FREE (work->data);
  SC_FREE (work);
}

double             *
sc_darray_work_get (sc_darray_work_t * work, const int thread,
                    const int block)
{
  SC_ASSERT (0 <= thread && thread < work->n_threads);
  SC_ASSERT (0 <= block && block < work->n_blocks);

  return work->data + work->n_entries * (work->n_blocks * thread + block);
}

int
sc_darray_work_get_blockcount (sc_darray_work_t * work)
{
  return work->n_blocks;
}

int
sc_darray_work_get_blocksize (sc_darray_work_t * work)
{
  return work->n_entries;
}
