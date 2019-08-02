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

#include <sc_bspline.h>

int
sc_bspline_min_number_points (int n)
{
  return n + 1;
}

int
sc_bspline_min_number_knots (int n)
{
  return 2 * n + 2;
}

sc_dmatrix_t       *
sc_bspline_knots_new (int n, sc_dmatrix_t * points)
{
#ifdef SC_ENABLE_DEBUG
  const int           d = points->n;
#endif
  const int           p = points->m - 1;
  const int           m = n + p + 1;
  const int           l = m - 2 * n;
  int                 i;
  sc_dmatrix_t       *knots;
  double             *knotse;

  SC_ASSERT (n >= 0 && m >= 1 && d >= 1 && l >= 1);

  knots = sc_dmatrix_new (m + 1, 1);
  knotse = knots->e[0];

  for (i = 0; i < n; ++i) {
    knotse[i] = 0.;
    knotse[m - i] = 1.;
  }
  for (i = 0; i <= l; ++i) {
    knotse[n + i] = i / (double) l;
  }

  return knots;
}

sc_dmatrix_t       *
sc_bspline_knots_new_length (int n, sc_dmatrix_t * points)
{
  const int           d = points->n;
  const int           p = points->m - 1;
  const int           m = n + p + 1;
  const int           l = m - 2 * n;
  int                 i, k;
  double              distsqr, distsum, distalln;
  double             *knotse;
  sc_dmatrix_t       *knots;

  SC_ASSERT (n >= 1);
  SC_ASSERT (n >= 0 && m >= 1 && d >= 1 && l >= 1);

  knots = sc_dmatrix_new_zero (m + 1, 1);
  knotse = knots->e[0];

  /* compute cumulative distance from P_0 and hide inside knots */
  distsum = 0.;
  for (i = 0; i < p; ++i) {
    SC_ASSERT (n + i + 2 >= 0 && n + i + 2 <= m);
    distsqr = 0.;
    for (k = 0; k < d; ++k) {
      distsqr += SC_SQR (points->e[i + 1][k] - points->e[i][k]);
    }
    knotse[n + i + 2] = distsum += sqrt (distsqr);
  }
  distalln = distsum * n;

  /* assign average cumulative distance to knot value */
  for (i = 1; i < l; ++i) {
    distsum = 0.;
    for (k = 0; k < n; ++k) {
      SC_ASSERT (n + i + k + 1 <= m);
      distsum += knotse[n + i + k + 1];
    }
    knotse[n + i] = distsum / distalln;
  }

  /* fill in the beginning and end values */
  for (i = 0; i <= n; ++i) {
    knotse[i] = 0.;
    knotse[m - i] = 1.;
  }

  return knots;
}

sc_dmatrix_t       *
sc_bspline_knots_new_periodic (int n, sc_dmatrix_t * points)
{
#ifdef SC_ENABLE_DEBUG
  const int           d = points->n;
#endif
  const int           p = points->m - 1;
  const int           m = n + p + 1;
  const int           l = m - 2 * n;
  int                 i;
  sc_dmatrix_t       *knots;
  double             *knotse;

  SC_ASSERT (n >= 0 && m >= 1 && d >= 1 && l >= 1);

  knots = sc_dmatrix_new (m + 1, 1);
  knotse = knots->e[0];

  for (i = 0; i <= m; ++i) {
    knotse[i] = (i - n) / (double) l;
  }

  return knots;
}

sc_dmatrix_t       *
sc_bspline_knots_new_length_periodic (int n, sc_dmatrix_t * points)
{
  const int           d = points->n;
  const int           p = points->m - 1;
  const int           m = n + p + 1;
  const int           l = m - 2 * n;
  int                 i, k;
  double              distsqr, distsum, distalln;
  double             *knotse;
  sc_dmatrix_t       *knots;

  SC_ASSERT (n >= 1);
  SC_ASSERT (n >= 0 && m >= 1 && d >= 1 && l >= 1);

  knots = sc_dmatrix_new_zero (m + 1, 1);
  knotse = knots->e[0];

  /* compute cumulative distance from P_0 and hide inside knots */
  distsum = 0.;
  for (i = 0; i < p; ++i) {
    SC_ASSERT (n + i + 2 >= 0 && n + i + 2 <= m);
    distsqr = 0.;
    for (k = 0; k < d; ++k) {
      distsqr += SC_SQR (points->e[i + 1][k] - points->e[i][k]);
    }
    knotse[n + i + 2] = sqrt (distsqr);
    if (i < l) {
      distsum += knotse[n + i + 2];
    }
  }
  distalln = distsum * n;

  /* assign average cumulative distance to knot value */
  knotse[n] = 0.;
  for (i = 1; i < l; ++i) {
    distsum = 0.;
    for (k = 0; k < n; ++k) {
      SC_ASSERT (n + i + k + 1 <= m);
      distsum += knotse[n + i + k + 1];
    }
    knotse[n + i] = knotse[n + i - 1] + distsum / distalln;
  }
  knotse[n + l] = 1.;

  /* fill in the beginning and end values */
  for (i = 0; i < n; ++i) {
    knotse[i] = knotse[i + l] - 1.;
    knotse[m - i] = knotse[2 * n - i] + 1.;
  }

  return knots;
}

void
sc_bspline_make_points_periodic (int n, sc_dmatrix_t * points)
{
  const int           d = points->n;
  const int           l = points->m;
  const int           p = n + l - 1;
  const int           shift = (n / 2);
  int                 i, j;

  SC_ASSERT (d > 0 && n >= 0 && l > 0);

  if (!n) {
    /* already periodic */
    return;
  }

  sc_dmatrix_resize (points, p + 1, d);

  /* shift to make room for the starting points */
  for (i = l - 1; i >= 0; i--) {
    for (j = 0; j < d; j++) {
      SC_ASSERT (i + shift <= p);
      points->e[i + shift][j] = points->e[i][j];
    }
  }

  /* copy the starting points */
  for (i = 0; i < shift; i++) {
    for (j = 0; j < d; j++) {
      SC_ASSERT (i + l <= p);
      points->e[i][j] = points->e[i + l][j];
    }
  }

  /* copy the ending points */
  for (i = shift + l; i <= p; i++) {
    for (j = 0; j < d; j++) {
      SC_ASSERT (i - l >= 0);
      points->e[i][j] = points->e[i - l][j];
    }
  }

  SC_ASSERT (sc_dmatrix_is_valid (points));
}

sc_dmatrix_t       *
sc_bspline_workspace_new (int n, int d)
{
  SC_ASSERT (n >= 0 && d >= 1);

  return sc_dmatrix_new ((n + 1) * (n + 1), d);
}

sc_bspline_t       *
sc_bspline_new (int n, sc_dmatrix_t * points,
                sc_dmatrix_t * knots, sc_dmatrix_t * works)
{
  sc_bspline_t       *bs;

  bs = SC_ALLOC_ZERO (sc_bspline_t, 1);
  bs->d = points->n;
  bs->p = points->m - 1;
  bs->n = n;
  bs->m = bs->n + bs->p + 1;
  bs->l = bs->m - 2 * bs->n;
  bs->cacheknot = n;

  SC_ASSERT (n >= 0 && bs->m >= 1 && bs->d >= 1 && bs->l >= 1);

  bs->points = points;
  if (knots == NULL) {
    bs->knots = sc_bspline_knots_new (bs->n, points);
    bs->knots_owned = 1;
  }
  else {
    SC_ASSERT (knots->m == bs->m + 1);
    SC_ASSERT (knots->n == 1);
    bs->knots = knots;
    bs->knots_owned = 0;
  }
  if (works == NULL) {
    bs->works = sc_bspline_workspace_new (bs->n, bs->d);
    bs->works_owned = 1;
  }
  else {
    SC_ASSERT (works->m == (bs->n + 1) * (bs->n + 1));
    SC_ASSERT (works->n == bs->d);
    bs->works = works;
    bs->works_owned = 0;
  }

  return bs;
}

void
sc_bspline_destroy (sc_bspline_t * bs)
{
  if (bs->knots_owned)
    sc_dmatrix_destroy (bs->knots);
  if (bs->works_owned)
    sc_dmatrix_destroy (bs->works);

  SC_FREE (bs);
}

static int
sc_bspline_find_interval (sc_bspline_t * bs, double t)
{
  int                 i, iguess;
  double              t0, t1;
  const double       *knotse = bs->knots->e[0];

  t0 = knotse[bs->n];
  t1 = knotse[bs->n + bs->l];
  SC_ASSERT (t >= t0 && t <= t1);
  SC_ASSERT (bs->cacheknot >= bs->n && bs->cacheknot < bs->n + bs->l);

  if (t >= t1) {
    iguess = bs->cacheknot = bs->n + bs->l - 1;
  }
  else if (knotse[bs->cacheknot] <= t && t < knotse[bs->cacheknot + 1]) {
    iguess = bs->cacheknot;
  }
  else {
    const int           nshift = 1;
    int                 ileft, iright;
    double              tleft, tright;

    ileft = bs->n;
    iright = bs->n + bs->l - 1;
    iguess = bs->n + (int) floor ((t - t0) / (t1 - t0) * bs->l);
    iguess = SC_MAX (iguess, ileft);
    iguess = SC_MIN (iguess, iright);

    for (i = 0;; ++i) {
      tleft = knotse[iguess];
      tright = knotse[iguess + 1];
      if (t < tleft) {
        iright = iguess - 1;
        if (i < nshift) {
          iguess = iright;
        }
        else {
          iguess = (ileft + iright + 1) / 2;
        }
      }
      else if (t >= tright) {
        ileft = iguess + 1;
        if (i < nshift) {
          iguess = ileft;
        }
        else {
          iguess = (ileft + iright) / 2;
        }
      }
      else {
        if (i > 0) {
          SC_LDEBUGF ("For %g needed %d search steps\n", t, i);
        }
        break;
      }
    }
    bs->cacheknot = iguess;
  }
  SC_ASSERT (iguess >= bs->n && iguess < bs->n + bs->l);
  SC_CHECK_ABORT ((knotse[iguess] <= t && t < knotse[iguess + 1]) ||
                  (t >= t1 && iguess == bs->n + bs->l - 1),
                  "Bug in sc_bspline_find_interval");

  return iguess;
}

void
sc_bspline_evaluate (sc_bspline_t * bs, double t, double *result)
{
  int                 i, k, n;
  int                 iguess;
  int                 toffset;
  double             *wfrom, *wto;
  const double       *knotse = bs->knots->e[0];

  iguess = sc_bspline_find_interval (bs, t);

  toffset = 0;
  wfrom = wto = bs->points->e[iguess - bs->n];
  for (n = bs->n; n > 0; --n) {
    wto = bs->works->e[toffset];

    for (i = 0; i < n; ++i) {
      const double        tleft = knotse[iguess + i - n + 1];
      const double        tright = knotse[iguess + i + 1];
      const double        tfactor = 1. / (tright - tleft);

      for (k = 0; k < bs->d; ++k) {
        wto[bs->d * i + k] =
          ((t - tleft) * wfrom[bs->d * (i + 1) + k] +
           (tright - t) * wfrom[bs->d * i + k]) * tfactor;
      }
    }

    wfrom = wto;
    toffset += n;
  }
  SC_ASSERT (toffset == bs->n * (bs->n + 1) / 2);

  memcpy (result, wfrom, sizeof (double) * bs->d);
}

void
sc_bspline_derivative (sc_bspline_t * bs, double t, double *result)
{
  sc_bspline_derivative_n (bs, 1, t, result);
}

void
sc_bspline_derivative_n (sc_bspline_t * bs, int order,
                         double t, double *result)
{
  int                 i, k, n;
  int                 iguess;
  int                 toffset;
  double             *wfrom, *wto;
  const double       *knotse = bs->knots->e[0];

  SC_ASSERT (order >= 0);

  if (bs->n < order) {
    memset (result, 0, sizeof (double) * bs->d);
    return;
  }

  iguess = sc_bspline_find_interval (bs, t);

  toffset = 0;
  wfrom = wto = bs->points->e[iguess - bs->n];
  for (n = bs->n; n > 0; --n) {
    wto = bs->works->e[toffset];

    if (bs->n < n + order) {
      for (i = 0; i < n; ++i) {
        const double        tleft = knotse[iguess + i - n + 1];
        const double        tright = knotse[iguess + i + 1];
        const double        tfactor = n / (tright - tleft);

        for (k = 0; k < bs->d; ++k) {
          wto[bs->d * i + k] =
            (wfrom[bs->d * (i + 1) + k] - wfrom[bs->d * i + k]) * tfactor;
        }
      }
    }
    else {
      for (i = 0; i < n; ++i) {
        const double        tleft = knotse[iguess + i - n + 1];
        const double        tright = knotse[iguess + i + 1];
        const double        tfactor = 1. / (tright - tleft);

        for (k = 0; k < bs->d; ++k) {
          wto[bs->d * i + k] =
            ((t - tleft) * wfrom[bs->d * (i + 1) + k] +
             (tright - t) * wfrom[bs->d * i + k]) * tfactor;
        }
      }
    }

    wfrom = wto;
    toffset += n;
  }
  SC_ASSERT (toffset == bs->n * (bs->n + 1) / 2);

  memcpy (result, wfrom, sizeof (double) * bs->d);
}

void
sc_bspline_derivative2 (sc_bspline_t * bs, double t, double *result)
{
  int                 i, k, n;
  int                 iguess;
  int                 toffset;
  double             *pfrom, *pto;
  double             *qfrom, *qto;
  const double       *knotse = bs->knots->e[0];

  iguess = sc_bspline_find_interval (bs, t);

  toffset = bs->n + 1;
  pfrom = pto = bs->works->e[0];
  qfrom = qto = bs->points->e[iguess - bs->n];
  memset (pfrom, 0, toffset * bs->d * sizeof (double));
  for (n = bs->n; n > 0; --n) {
    pto = bs->works->e[toffset];
    qto = bs->works->e[toffset + n];

    for (i = 0; i < n; ++i) {
      const double        tleft = knotse[iguess + i - n + 1];
      const double        tright = knotse[iguess + i + 1];
      const double        tfactor = 1. / (tright - tleft);

      for (k = 0; k < bs->d; ++k) {
        pto[bs->d * i + k] =
          ((t - tleft) * pfrom[bs->d * (i + 1) + k] +
           (tright - t) * pfrom[bs->d * i + k] +
           qfrom[bs->d * (i + 1) + k] - qfrom[bs->d * i + k]) * tfactor;
        qto[bs->d * i + k] =
          ((t - tleft) * qfrom[bs->d * (i + 1) + k] +
           (tright - t) * qfrom[bs->d * i + k]) * tfactor;
      }
    }

    pfrom = pto;
    qfrom = qto;
    toffset += 2 * n;
  }
  SC_ASSERT (toffset == (bs->n + 1) * (bs->n + 1));

  memcpy (result, pfrom, sizeof (double) * bs->d);
}
