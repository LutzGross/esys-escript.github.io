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


# This example demonstrates table interpolation and saving data in CSV format.
# It corresponds to Section "Interpolating Data" in the User Guide.

from esys.escript import InterpolationTable, saveDataCSV, sup, mkDir
import numpy
try:
    from esys.finley import Rectangle, Brick
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

if not HAVE_FINLEY:
    print("Finley module not available")
else:
    mkDir("output")
    n = 4
    r = Rectangle(n, n)
    x = r.getX()   # Data with shape (2,)
    toobig = 100   # RuntimeError if any interpolated value exceeds this

    # ------------------------------------------------------------------ #
    # 1-D interpolation — order-1 (linear)
    # Approximate a sine curve over [0, 1]
    # ------------------------------------------------------------------ #

    sine_table = numpy.array([0, 0.70710678118654746, 1, 0.70710678118654746, 0,
                               -0.70710678118654746, -1, -0.70710678118654746, 0])
    minval, maxval = 0., 1.
    # n-1 intervals between n sample points
    step = (maxval - minval) / (len(sine_table) - 1)

    interp1d = InterpolationTable(sine_table, origin=minval, step=step,
                                  order=1, undef=toobig)
    result = interp1d(x[0])
    saveDataCSV("output/1d.csv", inp=x[0], out=result)

    # ------------------------------------------------------------------ #
    # 1-D interpolation — order-0 (piecewise constant / cell-centred)
    # Each cell i covers [origin + i*step, origin + (i+1)*step)
    # Use extend= to specify total domain width directly
    # ------------------------------------------------------------------ #

    interp1d_o0 = InterpolationTable(sine_table, origin=minval,
                                     extend=(maxval - minval),
                                     order=0, undef=toobig)
    result0 = interp1d_o0(x[0])

    # ------------------------------------------------------------------ #
    # 2-D interpolation — order-1 using "ijk" convention (default)
    # table[ix, iy] = value at (xmin + ix*xstep, ymin + iy*ystep)
    # Table shape: (nx, ny)
    # ------------------------------------------------------------------ #

    nx = ny = 10
    xstep = ystep = (maxval - minval) / (nx - 1)
    xmin = ymin = minval

    # Sine amplitude decreases from 1 at y=0 to 0 at y=maxval
    st = sine_table  # length 9
    # Build a 9 x 3 table: table[ix, iy], shape (9, 3)
    # y covers [0, 0.55*2] in 3 rows (just as in the original example)
    table2 = numpy.array([st, 0.5 * st, 0. * st]).T  # shape (9, 3)
    x2step = step        # step along x (sine)
    y2step = 0.55        # step along y

    interp2d = InterpolationTable(table2, origin=(minval, 0.),
                                  step=(x2step, y2step),
                                  order=1, undef=toobig)
    result2 = interp2d(x)   # x has shape (2,)
    saveDataCSV("output/2d.csv", inp0=x[0], inp1=x[1], out=result2)

    # ------------------------------------------------------------------ #
    # 3-D interpolation — order-1 using "ijk" convention (default)
    # table[ix, iy, iz] = ix + iy*10 + iz*100
    # Table shape: (nx, ny, nz)
    # ------------------------------------------------------------------ #

    b = Brick(n, n, n)
    x3 = b.getX()   # shape (3,)
    toobig3 = 1000000

    nx3 = ny3 = nz3 = 10
    xstep3 = ystep3 = zstep3 = (maxval - minval) / (nx3 - 1)

    table3 = numpy.array([[[ix + iy * 10 + iz * 100
                             for iz in range(nz3)]
                            for iy in range(ny3)]
                           for ix in range(nx3)], dtype=float)

    interp3d = InterpolationTable(table3,
                                  origin=(minval, minval, minval),
                                  step=(xstep3, ystep3, zstep3),
                                  order=1, undef=toobig3)
    result3 = interp3d(x3)