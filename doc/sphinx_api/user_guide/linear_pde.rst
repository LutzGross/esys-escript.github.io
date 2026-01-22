.. _SEC LinearPDE:

=====================
The linearPDEs Module
=====================

The **esys.escript.linearPDEs** module provides classes for defining and solving
linear partial differential equations (PDEs) within escript. The module interfaces
with PDE solver libraries defined through the domain (e.g., finley, ripley).

.. note::

   This chapter is being converted from LaTeX. For the complete mathematical formulations
   and detailed coefficient specifications, please refer to the `User Guide PDF <../../user/user.pdf>`_.

Overview
========

The ``LinearPDE`` class is used to define a general linear, steady, second-order PDE
for an unknown function *u* on a given domain. The module supports:

* Single PDEs with scalar solutions
* Systems of PDEs with vector solutions
* Natural (Neumann) boundary conditions
* Essential (Dirichlet) boundary conditions via constraints
* Contact conditions for discontinuities
* Dirac delta function source terms

PDE Formulation
===============

Single PDE
----------

For a single PDE with a scalar solution, the general form is:

.. math::

   -(A_{jl} u_{,l})_{,j} - (B_j u)_{,j} + C_l u_{,l} + D u = -X_{j,j} + Y

where:

* *A* is a rank-2 tensor (diffusion coefficient)
* *B*, *C*, *X* are rank-1 tensors
* *D*, *Y* are scalars
* Subscript ``,j`` denotes partial derivative with respect to the j-th direction

Natural Boundary Conditions
---------------------------

On the boundary, the following natural boundary condition is applied:

.. math::

   n_j (A_{jl} u_{,l} + B_j u) + d u = n_j X_j + y

where *n* is the outward normal and *d*, *y* are boundary coefficients.

Constraints (Dirichlet Conditions)
----------------------------------

Constraints prescribe the solution value at specific locations:

.. math::

   u = r \text{ where } q > 0

where *q* is a characteristic function defining where the constraint applies.

System of PDEs
--------------

For systems with multiple solution components:

.. math::

   -(A_{ijkl} u_{k,l})_{,j} - (B_{ijk} u_k)_{,j} + C_{ikl} u_{k,l} + D_{ik} u_k = -X_{ij,j} + Y_i

where the coefficients have appropriately higher ranks.

LinearPDE Class
===============

Basic Usage
-----------

.. code-block:: python

   from esys.escript import *
   from esys.escript.linearPDEs import LinearPDE
   from esys.finley import Rectangle

   # Create domain
   mydomain = Rectangle(l0=1., l1=1., n0=40, n1=20)

   # Create PDE
   mypde = LinearPDE(mydomain)

   # Set coefficients
   kappa = 1.0  # diffusion coefficient
   mypde.setValue(A=kappa * kronecker(mydomain), D=1, Y=1)

   # Solve
   u = mypde.getSolution()

Setting Coefficients
--------------------

Use ``setValue()`` to specify PDE coefficients:

.. code-block:: python

   mypde.setValue(
       A=...,    # Rank-2: diffusion tensor
       B=...,    # Rank-1: advection
       C=...,    # Rank-1: reaction gradient
       D=...,    # Scalar: reaction
       X=...,    # Rank-1: flux source
       Y=...,    # Scalar: volume source
       d=...,    # Scalar: boundary reaction
       y=...,    # Scalar: boundary source
       q=...,    # Scalar: constraint mask
       r=...,    # Scalar: constraint value
   )

Symmetry
--------

If the PDE is symmetric (A is symmetric and B=C), enable symmetry for
more efficient solving:

.. code-block:: python

   mypde.setSymmetryOn()

Solver Options
==============

Solver options are accessed through ``getSolverOptions()``:

.. code-block:: python

   from esys.escript.linearPDEs import SolverOptions

   # Enable verbose output
   mypde.getSolverOptions().setVerbosityOn()

   # Set solver method
   mypde.getSolverOptions().setSolverMethod(SolverOptions.PCG)

   # Set preconditioner
   mypde.getSolverOptions().setPreconditioner(SolverOptions.AMG)

   # Set tolerance
   mypde.getSolverOptions().setTolerance(1e-8)

   # Set maximum iterations
   mypde.getSolverOptions().setIterMax(1000)

Available Solvers
-----------------

* ``SolverOptions.DEFAULT`` - Default solver
* ``SolverOptions.DIRECT`` - Direct solver (if available)
* ``SolverOptions.PCG`` - Preconditioned Conjugate Gradient
* ``SolverOptions.GMRES`` - Generalized Minimum Residual
* ``SolverOptions.BICGSTAB`` - Biconjugate Gradient Stabilized

Available Preconditioners
-------------------------

* ``SolverOptions.ILU0`` - Incomplete LU factorization
* ``SolverOptions.JACOBI`` - Jacobi (diagonal)
* ``SolverOptions.AMG`` - Algebraic Multigrid
* ``SolverOptions.RILU`` - Relaxed ILU

Diagnostics
-----------

After solving, retrieve diagnostic information:

.. code-block:: python

   u = mypde.getSolution()

   print("Iterations:", mypde.getDiagnostics("num_iter"))
   print("Total time:", mypde.getDiagnostics("time"))
   print("Setup time:", mypde.getDiagnostics("set_up_time"))
   print("Residual:", mypde.getDiagnostics("residual_norm"))

Specialized PDE Classes
=======================

Poisson
-------

For the Poisson equation: :math:`-\nabla^2 u = f`

.. code-block:: python

   from esys.escript.linearPDEs import Poisson

   pde = Poisson(domain=mydomain)
   pde.setValue(f=1, q=boundary_mask)
   u = pde.getSolution()

Helmholtz
---------

For the Helmholtz equation: :math:`-\nabla^2 u + \omega u = f`

.. code-block:: python

   from esys.escript.linearPDEs import Helmholtz

   pde = Helmholtz(domain=mydomain)
   pde.setValue(f=source, omega=frequency**2, q=boundary_mask)
   u = pde.getSolution()

LameEquation
------------

For linear elasticity problems:

.. code-block:: python

   from esys.escript.linearPDEs import LameEquation

   pde = LameEquation(domain=mydomain)
   pde.setValue(lame_lambda=lam, lame_mu=mu, F=body_force, q=fixed_boundary)
   displacement = pde.getSolution()

Example: Heat Conduction
========================

.. code-block:: python

   from esys.escript import *
   from esys.escript.linearPDEs import LinearPDE
   from esys.finley import Rectangle

   # Domain
   domain = Rectangle(l0=1., l1=1., n0=50, n1=50)
   x = domain.getX()

   # Thermal conductivity (anisotropic)
   kappa = kronecker(domain)
   kappa[0, 0] = 2.0  # Higher conductivity in x-direction

   # Heat source
   Q = whereNegative(length(x - [0.5, 0.5]) - 0.1)  # Point source at center

   # Boundary conditions: T=0 on all boundaries
   q = whereZero(x[0]) + whereZero(x[0] - 1) + whereZero(x[1]) + whereZero(x[1] - 1)

   # Solve
   pde = LinearPDE(domain)
   pde.setSymmetryOn()
   pde.setValue(A=kappa, Y=Q, q=q)
   T = pde.getSolution()

For Complete Reference
======================

For the complete mathematical formulations, system of PDEs, contact conditions,
and advanced solver options, please see:

* `User Guide PDF <../../user/user.pdf>`_
* :doc:`trilinos` - Trilinos solver integration
* :doc:`/api_index` - Python API Reference

See Also
========

* :doc:`escript` - escript module (Data objects, function spaces)
* :doc:`finley` - Finley domain module
* :doc:`symbolic` - Symbolic toolbox for nonlinear PDEs
