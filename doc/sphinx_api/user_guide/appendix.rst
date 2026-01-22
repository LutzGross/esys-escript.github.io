========
Appendix
========

This appendix provides supplementary mathematical background and advanced topics
for working with **esys.escript**.

.. note::

   This appendix is being converted from LaTeX. For the complete mathematical derivations,
   please refer to the `User Guide PDF <../../user/user.pdf>`_.

.. _EINSTEIN NOTATION:

Einstein Notation
=================

Compact notation is used in equations involving continuum mechanics and linear
algebra, known as Einstein notation or the Einstein summation convention.
It makes conventional tensor equations more compact.

Rules
-----

1. **Index Notation**: The rank of a tensor is represented by an index.

   * *a* - scalar (rank 0)
   * *b_i* - vector (rank 1)
   * *c_ij* - matrix (rank 2)

2. **Summation Convention**: Repeated subscripts are summed over all values.

Examples
--------

The expression:

.. math::

   y = a_0 b_0 + a_1 b_1 + \ldots + a_n b_n = \sum_{i=0}^n a_i b_i

In Einstein notation becomes simply:

.. math::

   y = a_i b_i

The gradient:

.. math::

   \nabla p = \frac{\partial p}{\partial x_0}\mathbf{i} + \frac{\partial p}{\partial x_1}\mathbf{j} + \frac{\partial p}{\partial x_2}\mathbf{k}

In Einstein notation (where comma denotes partial derivative):

.. math::

   \nabla p = p_{,i}

Kronecker Delta
---------------

The Kronecker delta symbol is a matrix with ones on the diagonal and zeros elsewhere:

.. math::

   \delta_{ij} = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \neq j \end{cases}

In escript, use ``kronecker(domain)`` to create this tensor.

.. _APP NEWTON:

Non-Linear PDEs
===============

The ``NonlinearPDE`` class solves nonlinear PDEs using the Newton-Raphson method.

Problem Formulation
-------------------

The solution *u* satisfies:

.. math::

   \int_\Omega v_{i,j} X_{ij} + v_i Y_i \, dx + \int_{\partial\Omega} v_i y_i \, ds = 0

for all smooth test functions *v* with *v = 0* where *q > 0*, where *X* and *Y*
are nonlinear functions of the solution and its gradient.

Newton-Raphson Scheme
---------------------

Starting from an initial guess :math:`u^{(0)}`, the iteration:

.. math::

   u^{(\nu)} = u^{(\nu-1)} - \delta^{(\nu-1)}

converges to the solution. Each correction :math:`\delta` is found by solving
a linearized PDE.

Convergence Criteria
--------------------

The iteration stops when:

.. math::

   \| u - u^{(\nu)} \|_\infty \leq \text{rtol} \cdot \| u \|_\infty

where rtol is the relative tolerance.

Usage in escript
----------------

See the :doc:`symbolic` chapter for details on using the ``NonlinearPDE`` class.

Common Tensors
==============

Identity Tensor
---------------

The identity tensor (Kronecker delta) in 2D:

.. math::

   I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}

In escript:

.. code-block:: python

   I = kronecker(domain)

Fourth-Order Identity
---------------------

For elasticity problems, the fourth-order identity tensor:

.. code-block:: python

   I4 = identityTensor4(domain)

Useful Relations
================

Tensor Operations
-----------------

**Inner product:**

.. math::

   a : b = a_{ij} b_{ij}

In escript: ``inner(a, b)``

**Outer product:**

.. math::

   (a \otimes b)_{ijkl} = a_{ij} b_{kl}

In escript: ``outer(a, b)``

**Trace:**

.. math::

   \text{tr}(A) = A_{ii}

In escript: ``trace(A)``

**Symmetric part:**

.. math::

   \text{sym}(A)_{ij} = \frac{1}{2}(A_{ij} + A_{ji})

In escript: ``symmetric(A)``

Differential Operators
----------------------

**Gradient:**

.. math::

   (\nabla u)_i = u_{,i}

In escript: ``grad(u)``

**Divergence:**

.. math::

   \nabla \cdot \mathbf{v} = v_{i,i}

In escript: ``div(v)``

**Laplacian:**

.. math::

   \nabla^2 u = u_{,ii}

In escript (via PDE): Use ``LinearPDE`` with ``A=kronecker(domain)``

For Complete Reference
======================

For the complete mathematical derivations and additional topics, please see the
`User Guide PDF <../../user/user.pdf>`_.
