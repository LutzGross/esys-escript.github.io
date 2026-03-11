##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


"""
Provides :class:`InterpolationTable`, a class-based interface for
interpolating values from a regular-grid lookup table onto a mesh.
Supports order-0 (piecewise constant) and order-1 (linear) interpolation
in 1D, 2D, and 3D.  The lookup table may contain real or complex values;
the coordinate data *x* must always be real.
"""

import numpy as np
from . import escriptcpp as escore


class InterpolationTable:
    """
    Interpolates values from a regular-grid lookup table onto a mesh.

    The coordinate data *x* passed to :meth:`__call__` determines the
    interpolation dimension:

    ========================  ===========  =======================
    ``x.getShape()``          Lookup dim   Required table rank
    ========================  ===========  =======================
    ``()``                    1-D          1  — shape ``(nx,)``
    ``(1,)``                  1-D          1  — shape ``(nx,)``
    ``(2,)``                  2-D          2
    ``(3,)``                  3-D          3
    ========================  ===========  =======================

    The returned `Data` object is always **scalar** (shape ``()``).
    If the table contains complex values the result is complex `Data`;
    the coordinate `Data` *x* must always be real.

    **Table indexing convention** (controlled by *table_indexing*):

    * ``"ijk"`` *(default)* — natural numpy order: the first index varies
      along the x-axis, the second along y, the third along z.
      Shapes: ``(nx,)`` / ``(nx, ny)`` / ``(nx, ny, nz)``.

    * ``"kji"`` — legacy / C++ order: the first index varies along the
      z-axis, the second along y, the third along x.
      Shapes: ``(nx,)`` / ``(ny, nx)`` / ``(nz, ny, nx)``.

    Example usage::

        import numpy as np
        from esys.finley import Rectangle
        from esys.escript import Function
        from esys.escriptcore.interpolation import InterpolationTable

        dom = Rectangle(20, 20)
        x = Function(dom).getX()   # shape (2,)

        # 2-D scalar table, "ijk" convention: shape (nx, ny)
        t = np.random.rand(5, 5)
        interp = InterpolationTable(t, origin=(0., 0.), step=(0.25, 0.25))
        result = interp(x)

        # Equivalent using total extent (covers [0,1] x [0,1])
        interp = InterpolationTable(t, origin=(0., 0.), extend=(1., 1.))
        result = interp(x)

        # Order-0 (piecewise constant) interpolation
        interp0 = InterpolationTable(t, origin=(0., 0.), extend=(1., 1.), order=0)
        result0 = interp0(x)

    :param table: lookup table as a numpy array of rank 1, 2, or 3
    :type table: ``numpy.ndarray``
    :param origin: coordinate(s) of the first table entry; a single float
        for 1-D or a tuple for 2-D / 3-D
    :type origin: ``float`` or ``tuple`` of ``float``
    :param step: cell size(s), i.e. the spacing between adjacent table
        entries; all values must be strictly positive.  Exactly one of
        *step* and *extend* must be supplied.
    :type step: ``float`` or ``tuple`` of ``float``
    :param extend: total extent of the table domain along each coordinate
        axis.  For order-1: distance from the first to the last sample
        point along axis *i* (``step[i] = extend[i] / (n_i - 1)``).
        For order-0: total domain width along axis *i*
        (``step[i] = extend[i] / n_i``), where *n_i* is the number of
        cells along that axis.  Exactly one of *step* and *extend* must
        be supplied.
    :type extend: ``float`` or ``tuple`` of ``float``
    :param order: interpolation order — 0 for piecewise constant, 1 for
        linear (default)
    :type order: ``int``
    :param undef: upper threshold; result values above this trigger a
        ``RuntimeError``
    :type undef: ``float``
    :param check_boundaries: if ``True``, a ``RuntimeError`` is raised when
        a coordinate lies outside the table extent; otherwise the nearest
        boundary value is used
    :type check_boundaries: ``bool``
    :param table_indexing: indexing convention for the table array.
        ``"ijk"`` (default) — first index = x-axis, natural numpy order.
        ``"kji"`` — first index = z-axis, legacy C++ order.
    :type table_indexing: ``str``
    """

    def __init__(self, table, origin, step=None, extend=None, order=1,
                 undef=1.e50, check_boundaries=False, table_indexing="ijk"):
        if not isinstance(table, np.ndarray):
            table = np.asarray(table)
        # Promote integer tables to float; preserve complex dtype
        if np.issubdtype(table.dtype, np.integer):
            table = table.astype(float)
        if np.isscalar(origin):
            origin = (origin,)
        ndim = len(origin)
        if ndim < 1 or ndim > 3:
            raise ValueError("ndim (length of origin) must be 1, 2, or 3")
        if table.ndim != ndim:
            raise ValueError(
                "table rank {} does not match coordinate dimension {} "
                "(set by length of origin)".format(table.ndim, ndim))
        if order not in (0, 1):
            raise ValueError("order must be 0 or 1")
        if table_indexing not in ("ijk", "kji"):
            raise ValueError("table_indexing must be 'ijk' or 'kji'")

        # n[i] = number of table entries along coordinate axis i (0=x,1=y,2=z)
        if table_indexing == "ijk":
            self._n = tuple(int(table.shape[i]) for i in range(ndim))
            # Transpose to internal "kji" storage for C++ workers.
            # 1-D: no change; 2-D: axes (0,1)->(1,0); 3-D: (0,1,2)->(2,1,0)
            if ndim == 1:
                self._table = np.ascontiguousarray(table)
            elif ndim == 2:
                self._table = np.ascontiguousarray(table.T)
            else:
                self._table = np.ascontiguousarray(table.transpose(2, 1, 0))
        else:  # "kji"
            self._n = tuple(int(table.shape[ndim - 1 - i]) for i in range(ndim))
            self._table = np.ascontiguousarray(table)
        self._is_complex = np.iscomplexobj(self._table)

        # Resolve step / extend using n[i] (entries along axis i)
        if step is None and extend is None:
            raise ValueError("Exactly one of 'step' or 'extend' must be provided")
        if step is not None and extend is not None:
            raise ValueError("Only one of 'step' or 'extend' may be provided, not both")
        if extend is not None:
            if np.isscalar(extend):
                extend = (extend,)
            if len(extend) != ndim:
                raise ValueError("extend must have the same length as origin")
            if order == 0:
                # n cells each of width step; domain is [origin, origin+extend]
                step = tuple(float(e) / self._n[i]
                             for i, e in enumerate(extend))
            else:
                # n-1 intervals between n sample points
                step = tuple(float(e) / (self._n[i] - 1)
                             for i, e in enumerate(extend))
        else:
            if np.isscalar(step):
                step = (step,)
            if len(step) != ndim:
                raise ValueError("step must have the same length as origin")
            step = tuple(float(s) for s in step)

        if any(s <= 0 for s in step):
            raise ValueError("All step sizes must be strictly positive")
        self._origin = tuple(float(v) for v in origin)
        self._step = step
        self._ndim = ndim
        self._order = order
        self._undef = float(undef)
        self._check_boundaries = check_boundaries
        self._table_indexing = table_indexing

    # ------------------------------------------------------------------
    # Internal C++ dispatch (table must be a real contiguous numpy array)
    # ------------------------------------------------------------------

    def _cpp_1d(self, x0, table):
        if self._order == 0:
            return x0._interpolateTable1dOrder0(
                table, self._origin[0], self._step[0],
                self._undef, self._check_boundaries)
        else:
            return x0.interpolateTable(
                table, self._origin[0], self._step[0],
                self._undef, self._check_boundaries)

    def _cpp_2d(self, x0, x1, table):
        if self._order == 0:
            return x0._interpolateTable2dOrder0(
                table, self._origin[0], self._step[0],
                x1, self._origin[1], self._step[1],
                self._undef, self._check_boundaries)
        else:
            return x0.interpolateTable(
                table, self._origin[0], self._step[0],
                x1, self._origin[1], self._step[1],
                self._undef, self._check_boundaries)

    def _cpp_3d(self, x0, x1, x2, table):
        if self._order == 0:
            return x0._interpolateTable3dOrder0(
                table, self._origin[0], self._step[0],
                x1, self._origin[1], self._step[1],
                x2, self._origin[2], self._step[2],
                self._undef, self._check_boundaries)
        else:
            return x0._interpolateTable3d(
                table, self._origin[0], self._step[0],
                x1, self._origin[1], self._step[1],
                x2, self._origin[2], self._step[2],
                self._undef, self._check_boundaries)

    def _do_call(self, table, x):
        """Dispatch to the correct C++ worker using a real numpy *table*."""
        sh = x.getShape()
        if sh == ():
            return self._cpp_1d(x, table)
        if sh == (1,):
            return self._cpp_1d(x[0], table)
        if sh == (2,):
            return self._cpp_2d(x[0], x[1], table)
        # sh == (3,)
        return self._cpp_3d(x[0], x[1], x[2], table)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def __call__(self, x):
        """
        Return interpolated values at the mesh points described by *x*.

        :param x: coordinate data with shape ``()``, ``(1,)``, ``(2,)``, or
            ``(3,)``; the shape must be consistent with the table rank.
            *x* must be real — complex coordinate data is not supported.
        :type x: `Data`
        :return: interpolated field (complex if the table contains complex
            values, otherwise real)
        :rtype: `Data` with shape ``()``
        :raises TypeError: if *x* is not a real `Data` object
        :raises ValueError: if *x* has an unsupported shape or is inconsistent
            with the table rank
        """
        if not isinstance(x, escore.Data):
            raise TypeError("x must be a Data object")
        if x.isComplex():
            raise TypeError(
                "coordinate data x must be real; complex x is not supported")

        sh = x.getShape()
        expected_ndim = {(): 1, (1,): 1, (2,): 2, (3,): 3}
        if sh not in expected_ndim:
            raise ValueError(
                "x must have shape (), (1,), (2,), or (3,); got {}".format(sh))
        if expected_ndim[sh] != self._ndim:
            raise ValueError(
                "x with shape {} requires ndim={}, got ndim={}".format(
                    sh, expected_ndim[sh], self._ndim))

        if self._is_complex:
            t_re = np.ascontiguousarray(self._table.real)
            t_im = np.ascontiguousarray(self._table.imag)
            return self._do_call(t_re, x) + 1j * self._do_call(t_im, x)
        return self._do_call(self._table, x)

    def interpolate(self, x):
        """Alias for :meth:`__call__`."""
        return self(x)

    # ------------------------------------------------------------------
    # Read-only accessors
    # ------------------------------------------------------------------

    def getOrder(self):
        """Return the interpolation order (0 or 1)."""
        return self._order

    def getNdim(self):
        """Return the number of coordinate dimensions (1, 2, or 3)."""
        return self._ndim

    def getOrigin(self):
        """Return the tuple of starting coordinates."""
        return self._origin

    def getStep(self):
        """Return the tuple of cell sizes."""
        return self._step

    def getTableIndexing(self):
        """Return the table indexing convention (``'ijk'`` or ``'kji'``)."""
        return self._table_indexing

    def getExtend(self):
        """Return the tuple of domain extents along each coordinate axis.

        For order-1: ``extend[i] = step[i] * (n_i - 1)``
        (distance from first to last sample point along axis *i*).

        For order-0: ``extend[i] = step[i] * n_i``
        (total domain width; domain along axis *i* is
        ``[origin[i], origin[i] + extend[i]]``).

        Here *n_i* is the number of table entries along coordinate axis *i*
        (independent of the *table_indexing* convention).
        """
        if self._order == 0:
            return tuple(s * self._n[i] for i, s in enumerate(self._step))
        else:
            return tuple(s * (self._n[i] - 1) for i, s in enumerate(self._step))
