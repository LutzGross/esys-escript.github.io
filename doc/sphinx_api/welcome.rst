==================
Welcome
==================

**esys-escript** is a module for implementing mathematical models in Python using
the finite element method (FEM). As users do not access the underlying data
structures it is very easy to use and scripts can run on desktop computers as
well as massive parallel supercomputers without changes.

Application areas for esys-escript include:

* Geophysical inversion
* Earthquakes
* Porous media flow
* Reactive transport
* Plate subduction
* Erosion
* Earth mantle convection
* Tsunamis

esys-escript is designed as an easy-to-use environment for implementing mathematical
models based on non-linear, coupled, time-dependent partial differential equations.
It uses the finite element method (FEM) for spatial discretization and data
representation and is used through Python. It is suitable for rapid prototyping
(e.g. for a student project or thesis) as well as for large software projects.

Scripts are executed sequentially, on multi-core platforms via OpenMP and distributed
computing clusters using MPI. The hybrid mode of OpenMP and MPI is supported and
allows for solving problems with over 200 million unknowns on several thousand cores
on a parallel computer.

For geophysical inversion see also the extensions:

* `gambit <https://github.com/AndreaCodd/gambit>`_
* `fingal <https://github.com/LutzGross/fingal>`_


Main Features
=============

* Python-based user interface
* Two- and three-dimensional finite and spectral element simulations
* Specialized geophysical inversion module
* Support for VTK and SILO file formats
* Unstructured meshes from gmsh
* Parallelization with OpenMP and MPI support
* Flux Controlled Transport solver (FEM-FCT)
* Visualization with VisIt, ParaView, Mayavi and others
* Platform is Linux; there is limited support for macOS and Windows


Access
======

Source code is available at:
`https://github.com/LutzGross/esys-escript.github.io <https://github.com/LutzGross/esys-escript.github.io>`_

A version of the documentation for the current master branch is available
`here <https://lutzgross.github.io/esys-escript.github.io/>`_.


Questions and Bug Reports
=========================

To raise a question or to report a bug please start a
`GitHub issue <https://github.com/esys-escript/esys-escript.github.io/issues>`_.


Debian Distribution
===================

Debian packages ``python3-escript-mpi`` and ``python3-escript`` for **version 5** are available.


Anaconda Installation (Version 5)
=================================

To install *esys-escript* for `anaconda <https://www.anaconda.com>`_, first run ``conda`` and then the command::

    conda install esys-escript -c conda-forge

At present, this is the recommended way to run esys-escript on Windows.


Reference
=========

If you publish work that makes use of esys-escript, we would appreciate it if you would cite the following references:

* R Schaa, L Gross, J du Plessis, *PDE-based geophysical modelling using finite elements: examples from 3D resistivity and 2D magnetotellurics*, Journal of Geophysics and Engineering, Volume 13, Issue 2, April 2016, Pages S59-S73. `DOI:10.1088/1742-2132/13/2/S59 <https://doi.org/10.1088/1742-2132/13/2/S59>`_

* L. Gross, L. Bourgouin, A.J. Hale, H.-B. Muhlhaus, *Interface modeling in incompressible media using level sets in Escript*, Physics of the Earth and Planetary Interiors, Volume 163, Issues 1-4, 2007, Pages 23-34. `DOI:10.1016/j.pepi.2007.04.004 <https://doi.org/10.1016/j.pepi.2007.04.004>`_


License
=======

Copyright (c) 2003-2026 by the esys.escript Group.

The majority of the files contained in this distribution are licensed
under the Apache License, version 2.0: http://www.apache.org/licenses/LICENSE-2.0

See `LICENSE <https://github.com/LutzGross/esys-escript.github.io/blob/master/LICENSE>`_ for details.


Contributors
============

esys-escript has started in 2003 as one of the first projects that made systematic
use of the Python language for large-scale scientific computing. It is the product
of years of work by many people with funding from Australian Commonwealth.

See `CREDITS <https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS>`_ for details.
