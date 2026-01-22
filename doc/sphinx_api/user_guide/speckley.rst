===================
The speckley Module
===================

The **esys.speckley** Module
==============================

.. _chap:speckley:

*speckley* is a high-order form of *ripley*, supporting structured, uniform meshes in two and three dimensions. Uniform meshes allow a more regular division of elements among compute nodes. Possible orders range from 2 to 10, inclusive.

*speckley* domains cannot be created by reading from a mesh file.

The family of domain that will result from a “Rectangle“ or “Brick“ call depends on which module is imported in the specific script. The following line is an example of importing **esys.speckley** domains:

.. code:: python

   from esys.speckley import Rectangle, Brick

Formulation
-----------

For a single PDE that has a solution with a single component the linear PDE is defined in the following form:

.. math::

      .. _SPECKLEY.SINGLE.1:

   \begin{array}{cl} &
   \displaystyle{
   \int_{\Omega}
   D \cdot vu \; d\Omega } + \int_{\Gamma} d \cdot vu \; d{\Gamma}

   = \displaystyle{\int_{\Omega}  X_{j} \cdot v_{,j}+ Y \cdot v \; d\Omega }
   + \displaystyle{\int_{\Gamma} y \cdot v \; d{\Gamma}}
   \end{array}

Meshes
------

.. _SPECKLEY MESHES:

**esys.speckley** meshes are formed of regular elements using Gauss-Labatto-Legendre quadrature points. The number of quadrature points in each axis is dependent on the order of the domain. Examples of small Rectangle domains of different orders are shown in Figure `[SPECKLEY:FIG:MESHES] <#SPECKLEY:FIG:MESHES>`__.

Meshfiles cannot be used to generate **esys.speckley** domains.

.. note::

   See the `User Guide PDF <../../user/user.pdf>`_ for mesh figures showing
   different order domains.

.. _SPECKLEY:FIG:MESHES:

Linear Solvers in “SolverOptions“
---------------------------------

While **esys.speckley** has the same defaults as **esys.ripley**, the “SolverOptions.HRZ_LUMPING“ must be set. “PASO“ is not used in **esys.speckley**.

Cross-domain Interpolation
--------------------------

Data on a **esys.speckley** domain can be interpolated to a matching **esys.ripley** domain provided the two domains have identical dimension, length, and, in multi-process situations, domain sub-divisions.

A utility class, “SpeckleyToRipley“ is available to simplify meeting these conditions. To gain access to the class, the following will be required in the script:

.. code:: python

   from esys.escript.domainCouplers import SpeckleyToRipley

Functions
---------

.. container:: funcdesc

   Brickorder,n0,n1,n2,l0=1.,l1=1.,l2=1.,d0=-1,d1=-1,d2=-1, diracPoints=list(), diracTags=list(), comm=None generates a “Domain“ object representing a three-dimensional brick between :math:`(0,0,0)` and :math:`(l0,l1,l2)` with orthogonal faces. All elements will be regular and of order “order“. The brick is filled with “n0“ elements along the :math:`x_0`-axis, “n1“ elements along the :math:`x_1`-axis and “n2“ elements along the :math:`x_2`-axis. If built with *MPI* support, the domain will be subdivided “d0“ times along the :math:`x_0`-axis, “d1“ times along the :math:`x_1`-axis, and “d2“ times along the :math:`x_2`-axis. “d0“, “d1“, and “d2“ must be factors of the number of *MPI* processes requested. If axial subdivisions are not specified, automatic domain subdivision will take place. This may not be the most efficient construction and will likely result in extra elements being added to ensure proper distribution of work. Any extra elements added in this way will change the length of the domain proportionately. “diracPoints“ is a list of coordinate-tuples of points within the mesh, each point tagged with the respective string within “diracTags“. “comm“ is an optional *MPI* communicator (from ``mpi4py``) that allows using a custom communicator instead of ``MPI_COMM_WORLD``. If “comm“ is not “None“, **esys.escript** must be built with ``mpi4py`` support enabled.

.. container:: funcdesc

   Rectangleorder,n0,n1,l0=1.,l1=1.,d0=-1,d1=-1, diracPoints=list(), diracTags=list(), comm=None generates a “Domain“ object representing a two-dimensional rectangle between :math:`(0,0)` and :math:`(l0,l1)` with orthogonal faces. All elements will be regular and of order “order“. The rectangle is filled with “n0“ elements along the :math:`x_0`-axis and “n1“ elements along the :math:`x_1`-axis. If built with *MPI* support, the domain will be subdivided “d0“ times along the :math:`x_0`-axis and “d1“ times along the :math:`x_1`-axis. “d0“ and “d1“ must be factors of the number of *MPI* processes requested. If axial subdivisions are not specified, automatic domain subdivision will take place. This may not be the most efficient construction and will likely result in extra elements being added to ensure proper distribution of work. Any extra elements added in this way will change the length of the domain proportionately. “diracPoints“ is a list of coordinate-tuples of points within the mesh, each point tagged with the respective string within “diracTags“. “comm“ is an optional *MPI* communicator (from ``mpi4py``) that allows using a custom communicator instead of ``MPI_COMM_WORLD``. If “comm“ is not “None“, **esys.escript** must be built with ``mpi4py`` support enabled.
