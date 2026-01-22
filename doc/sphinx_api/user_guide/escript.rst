.. _ESCRIPT CHAP:

==================
The escript Module
==================

This chapter describes the core **esys.escript** module, including the fundamental
concepts of function spaces and Data objects that form the foundation of the
escript framework.

.. note::

   This chapter is being converted from LaTeX. For the complete content with all
   mathematical details and figures, please refer to the `User Guide PDF <../../user/user.pdf>`_.

Concepts
========

**esys.escript** is a Python module that allows you to represent the values of
a function at points in a Domain in such a way that the function will
be useful for Finite Element Method (FEM) simulation. It provides what we call
a *function space* that describes how the data is used in the simulation.

Function Spaces
---------------

To understand the term "function space", consider that the solution of a partial
differential equation (PDE) is a function on a domain. When solving a PDE using FEM,
the solution is piecewise-differentiable but, in general, its gradient is discontinuous.
Different function spaces reflect these different degrees of smoothness.

A function space is described by a ``FunctionSpace`` object. The following generators
are commonly used:

**Continuous Functions:**

* ``Solution(mydomain)`` - Solutions of a PDE
* ``ReducedSolution(mydomain)`` - Solutions with reduced smoothness (lower order approximation)
* ``ContinuousFunction(mydomain)`` - Continuous functions (e.g., temperature distribution)

**Discontinuous Functions:**

* ``Function(mydomain)`` - General functions, not necessarily continuous (e.g., stress field)
* ``ReducedFunction(mydomain)`` - Reduced integration version

**Boundary Functions:**

* ``FunctionOnBoundary(mydomain)`` - Functions on domain boundary (e.g., surface pressure)
* ``ReducedFunctionOnBoundary(mydomain)`` - Reduced integration version

**Contact/Discontinuity Functions:**

* ``FunctionOnContact0(mydomain)`` - Functions on side 0 of a discontinuity
* ``FunctionOnContact1(mydomain)`` - Functions on side 1 of a discontinuity

**Point Functions:**

* ``DiracDeltaFunctions(mydomain)`` - Functions defined on a set of points

Example
^^^^^^^

.. code-block:: python

   from esys.escript import *
   from esys.finley import Rectangle

   mydomain = Rectangle(l0=1., l1=1., n0=40, n1=20)
   solution_space = Solution(mydomain)

Data Objects
============

The ``Data`` class stores functions represented through their values at sample points,
chosen according to the function space. Data objects are used to define PDE coefficients
and store solutions.

Properties
----------

Data objects have:

* **Rank** - Number of indices (0 to 4)
* **Shape** - Range of each index (tuple of integers)

For example, a stress field has rank 2 and shape ``(d, d)`` where ``d`` is the spatial dimension.

Creating Data Objects
---------------------

.. code-block:: python

   from esys.escript import *
   from esys.finley import Rectangle

   mydomain = Rectangle(l0=1., l1=1., n0=40, n1=20)

   # Create Data with shape (2,3), rank 2, constant value 1
   mydat = Data(value=1, what=ContinuousFunction(mydomain), shape=(2,3))

   # Create from numpy array
   import numpy as np
   arr = np.ones((2, 3))
   mydat = Data(arr, ContinuousFunction(mydomain))

Operations on Data
------------------

Data objects support standard mathematical operations:

.. code-block:: python

   # Arithmetic
   c = a + b
   c = a * b
   c = a / b

   # Mathematical functions
   s = sin(data)
   e = exp(data)

   # Spatial operations
   g = grad(data)  # Gradient
   d = div(data)   # Divergence

   # Reductions
   m = sup(data)   # Maximum value
   m = inf(data)   # Minimum value
   n = Lsup(data)  # L-infinity norm

Interpolation
-------------

Functions can be interpolated between compatible function spaces automatically
when required. The ``interpolate`` function explicitly performs this:

.. code-block:: python

   # Interpolate to a different function space
   data_on_boundary = interpolate(data, FunctionOnBoundary(mydomain))

Coordinate Access
-----------------

Access the coordinates of sample points:

.. code-block:: python

   x = mydomain.getX()           # Coordinates on nodes
   x = Function(mydomain).getX() # Coordinates at integration points

Tagged Data
-----------

Regions of a domain can be identified by tags, and Data objects can have
different values in different tagged regions:

.. code-block:: python

   # Set different values based on region tags
   k = Scalar(1.0, Function(mydomain))
   insertTaggedValues(k, tag1=2.0, tag2=3.0)

For Complete Reference
======================

For the complete mathematical details, figures showing function space dependencies,
and comprehensive API documentation, please see:

* `User Guide PDF <../../user/user.pdf>`_
* :doc:`/api_index` - Python API Reference

See Also
========

* :ref:`SEC LinearPDE` - Linear PDE module
* :doc:`finley` - Finley domain module
* :doc:`ripley` - Ripley domain module
