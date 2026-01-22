.. _chap:tutorial:

==============================
Tutorial: Solving PDEs
==============================

This chapter provides an introduction to solving partial differential equations (PDEs)
using **esys.escript**. It covers the basics of setting up and solving PDEs with examples
ranging from simple Poisson equations to more complex time-dependent problems.

.. note::

   This chapter is being converted from LaTeX. For the complete content with all figures
   and equations, please refer to the `User Guide PDF <../../user/user.pdf>`_.

Topics Covered
==============

Installation
------------

To download *escript*, please visit https://github.com/LutzGross/esys-escript.github.io.
The website offers information about the installation process, as well as a way to ask
questions at https://github.com/LutzGross/esys-escript.github.io/issues.

First Steps
-----------

This section introduces:

* The Poisson equation and boundary value problems
* Finite element method (FEM) basics
* Creating domains with **esys.finley**
* Setting boundary conditions (Dirichlet and Neumann)
* Solving the PDE and visualizing results

Example - Poisson Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple example solving the Poisson equation:

.. code-block:: python

   from esys.escript import *
   from esys.escript.linearPDEs import Poisson
   from esys.finley import Rectangle
   from esys.weipa import saveVTK

   # Generate domain
   mydomain = Rectangle(l0=1., l1=1., n0=40, n1=20)

   # Define characteristic function of Gamma^D
   x = mydomain.getX()
   gammaD = whereZero(x[0]) + whereZero(x[1])

   # Define PDE and get its solution
   mypde = Poisson(domain=mydomain)
   mypde.setValue(f=1, q=gammaD)
   u = mypde.getSolution()

   # Write solution to file
   saveVTK("u.vtu", sol=u)


Diffusion Problems
------------------

This section covers:

* Time-dependent diffusion equations
* Implicit and explicit time stepping
* Anisotropic diffusion
* Temperature distribution examples


Wave Equations
--------------

This section covers:

* Second-order wave equations
* Acoustic wave propagation
* Time integration schemes


Heated Block Example
--------------------

A practical example of coupled thermal-mechanical analysis.


Dirac Delta Functions
---------------------

Using point sources and sinks in PDE solutions.


Unstructured Meshes
-------------------

Working with meshes imported from external mesh generators like Gmsh.


For Complete Content
====================

For the complete chapter with all mathematical derivations, figures, and detailed examples,
please see the `User Guide PDF <../../user/user.pdf>`_.
