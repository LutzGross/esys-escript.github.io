======================================
esys-escript User's Guide
======================================

Solving Partial Differential Equations with Escript and Finley

.. |escript| replace:: **esys-escript**
.. |python| replace:: **Python**
.. |mpi| replace:: **MPI**
.. |openmp| replace:: **OpenMP**
.. |cuda| replace:: **CUDA**

Introduction
============

|escript| is a |python|-based environment for implementing mathematical models, in particular those based on coupled, non-linear, time-dependent partial differential equations.

It consists of five major components:

* escript core library
* finite element solvers finley, dudley, ripley, and speckley (which use fast vendor-supplied solvers or the included PASO linear solver library)
* the meshing interface pycad
* a model library
* an inversion module

The current version supports parallelization through |mpi| for distributed memory, |openmp| for shared memory on CPUs, as well as |cuda| for some GPU-based solvers.

Table of Contents
=================

.. toctree::
   :maxdepth: 2
   :numbered:

   tutorial_pde
   execute
   escript
   linear_pde
   finley
   ripley
   speckley
   weipa
   trilinos
   symbolic
   subworlds
   appendix

Documentation Resources
=======================

* :ref:`genindex` - Complete index of all modules, classes and functions
* :ref:`modindex` - Module index
* `← Back to API Documentation <../index.html>`_
* `Documentation Home <../../index.html>`_
* `Installation Guide <../../installation.html>`_

Citation
========

If you use this software in your research, then we would appreciate (but do not require) a citation. Some relevant references can be found in the Appendix.

Copyright and License
=====================

Copyright (c) 2003-2025 by The University of Queensland

Primary Business: Queensland, Australia

Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0

Development History:
--------------------

* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
