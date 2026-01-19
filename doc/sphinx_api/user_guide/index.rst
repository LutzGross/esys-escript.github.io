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

It consists of the following major components:

* escript core library
* finite element solvers finley for unstructured meshes
* finite element solvers ripley for rectangular meshes
* spectral element solver speckley  for rectangular meshes
* interfaces to sparse matrix solvers including Trilinos, MUMPS, UMFPACK and the PASO linear solver library
* some function supporting large scale inversion

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

See CREDITS file for contributors and development history.
