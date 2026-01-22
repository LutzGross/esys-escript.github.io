.. _chap:finley:

=================
The finley Module
=================

The **esys.finley** module provides finite element domains for solving linear, steady
partial differential equations (PDEs) using isoparametric finite elements. It supports
unstructured 1D, 2D, and 3D meshes and is parallelized under both OpenMP and MPI.

.. note::

   This chapter is being converted from LaTeX. For the complete content with all
   element type tables, mesh figures, and mathematical formulations, please refer
   to the `User Guide PDF <../../user/user.pdf>`_.

Overview
========

Finley allows the creation of domains for solving PDEs or systems of PDEs using
isoparametric finite elements. The PDEs are represented by the ``LinearPDE`` class
from escript.

Key features:

* Unstructured meshes in 1D, 2D, and 3D
* Triangular, quadrilateral, tetrahedral, and hexahedral elements
* First and second order elements
* Macro elements for incompressible flow problems
* Contact elements for discontinuities
* OpenMP and MPI parallelization

Creating Domains
================

Built-in Mesh Generators
------------------------

Finley provides factory functions for creating simple rectangular domains:

.. code-block:: python

   from esys.finley import Rectangle, Brick

   # 2D rectangular domain
   domain = Rectangle(l0=1.0, l1=1.0, n0=40, n1=20)

   # 3D brick domain
   domain = Brick(l0=1.0, l1=1.0, l2=1.0, n0=10, n1=10, n2=10)

**Parameters:**

* ``l0``, ``l1``, ``l2`` - Domain dimensions
* ``n0``, ``n1``, ``n2`` - Number of elements in each direction
* ``order`` - Element order (1 or 2, default 1)
* ``useElementsOnFace`` - Use rich face elements (default 0)

Example with Second-Order Elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from esys.finley import Rectangle

   # Second-order elements for higher accuracy
   domain = Rectangle(l0=1.0, l1=1.0, n0=20, n1=20, order=2)

Importing External Meshes
-------------------------

Finley can read meshes from external mesh generators:

**Gmsh Format:**

.. code-block:: python

   from esys.finley import ReadGmsh

   domain = ReadGmsh("mesh.msh", numDim=2)

**Finley Native Format:**

.. code-block:: python

   from esys.finley import ReadMesh, LoadMesh

   domain = ReadMesh("mesh.fly")  # Text format
   domain = LoadMesh("mesh.nc")   # NetCDF format

Element Types
=============

Finley supports various element types:

**2D Elements:**

* ``Tri3`` - Linear triangle (3 nodes)
* ``Tri6`` - Quadratic triangle (6 nodes)
* ``Rec4`` - Bilinear quadrilateral (4 nodes)
* ``Rec8`` - Quadratic quadrilateral (8 nodes)
* ``Rec9`` - Biquadratic quadrilateral (9 nodes)

**3D Elements:**

* ``Tet4`` - Linear tetrahedron (4 nodes)
* ``Tet10`` - Quadratic tetrahedron (10 nodes)
* ``Hex8`` - Trilinear hexahedron (8 nodes)
* ``Hex20`` - Serendipity hexahedron (20 nodes)
* ``Hex27`` - Triquadratic hexahedron (27 nodes)

**Face Elements:**

Each interior element type has corresponding face and contact element types
for boundary conditions and discontinuities.

Mesh Structure
==============

A finley mesh consists of:

1. **Nodes** - Points with coordinates, identified by reference numbers
2. **Elements** - Interior elements defined by their nodes
3. **Face Elements** - Boundary elements for natural boundary conditions
4. **Contact Elements** - Elements across discontinuities
5. **Tags** - Integer labels for identifying regions

Accessing Mesh Information
--------------------------

.. code-block:: python

   from esys.finley import Rectangle
   from esys.escript import *

   domain = Rectangle(l0=1.0, l1=1.0, n0=10, n1=10)

   # Get coordinates
   x = domain.getX()

   # Get spatial dimension
   dim = domain.getDim()

   # Get boundary mask
   boundary_mask = whereZero(x[0]) + whereZero(x[1])

Tagging Regions
---------------

Regions can be tagged for setting different material properties:

.. code-block:: python

   from esys.escript import *

   # Get material property with different values in tagged regions
   k = Scalar(1.0, Function(domain))
   insertTaggedValues(k, left=2.0, right=3.0)

Boundary Conditions
===================

Face elements are used for natural (Neumann) boundary conditions. If no face
elements are specified, finley assumes homogeneous natural boundary conditions
(d=0, y=0) on the entire boundary.

For Dirichlet boundary conditions, use the constraint mechanism via the ``q``
and ``r`` coefficients in ``LinearPDE``.

Example: Poisson Equation
=========================

.. code-block:: python

   from esys.escript import *
   from esys.escript.linearPDEs import Poisson
   from esys.finley import Rectangle

   # Create domain
   domain = Rectangle(l0=1.0, l1=1.0, n0=40, n1=40)

   # Set up boundary conditions
   x = domain.getX()
   gammaD = whereZero(x[0]) + whereZero(x[1])  # Dirichlet on left and bottom

   # Solve Poisson equation: -laplacian(u) = 1
   pde = Poisson(domain=domain)
   pde.setValue(f=1, q=gammaD)
   u = pde.getSolution()

For Complete Reference
======================

For the complete list of element types, mesh file formats, and advanced features,
please see:

* `User Guide PDF <../../user/user.pdf>`_
* :doc:`/api_index` - Python API Reference

See Also
========

* :ref:`SEC LinearPDE` - Linear PDE module
* :doc:`ripley` - Regular grid domains (faster for rectangular meshes)
* :doc:`speckley` - Spectral element domains
