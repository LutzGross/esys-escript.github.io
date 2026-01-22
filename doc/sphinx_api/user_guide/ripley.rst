=================
The ripley Module
=================

The **esys.ripley** Module
============================

.. _chap:ripley:

**esys.ripley** is an alternative domain library to **esys.finley**; it supports structured, uniform meshes with rectangular elements in 2D and hexahedral elements in 3D. Uniform meshes allow a straightforward division of elements among processes with *MPI* and allow for a number of optimizations when solving PDEs. **esys.ripley** also supports fast assemblers for certain types of PDE (specifically Lamé and Wave PDEs). These assemblers make use of the regular nature of the domain to optimize the stiffness matrix assembly process for these specific problems. Finally, **esys.ripley** is currently the only domain family that supports GPU-based solvers.

As a result, **esys.ripley** domains cannot be created by reading from a mesh file since only one element type is supported and all elements need to be equally sized. For the same reasons, **esys.ripley** does not allow assigning coordinates via “setX()“.

While **esys.ripley** cannot be used with mesh files, it can be used to read in *GOCAD* data. A script with an example of a voxet reader is included in the examples as ``voxet_reader.py``.

Other than use of meshfiles, **esys.ripley** and **esys.finley** are generally interchangeable in a script with both modules having the “Rectangle“ or “Brick“ functions available. Consider the following example which creates a 2D **esys.ripley** domain:

.. code:: python

   from esys.ripley import Rectangle, Brick
    dom = Rectangle(9, 9)

Multi-resolution domains are supported in **esys.ripley** via “MultiBrick“ and “MultiRectangle“. Each level of one of these domains has twice the elements in each axis of the next lower resolution. The “MultiBrick“ is not currently supported when running **esys.escript** with multiple processes using *MPI*. Interpolation between these multi-resolution domains is possible providing they have matching dimensions and subdivisions, along with a compatible number of elements.

To simplify these conditions the use of “MultiResolutionDomain“ is highly recommended. The following example creates two 2D domains of different resolutions and interpolates between them:

.. code:: python

   from esys.ripley import MultiResolutionDomain
    mrd = MultiResolutionDomain(2, n0=10, n1=10)
    ten_by_ten = mrd.getLevel(0)
    data10 = Vector(..., Function(ten_by_ten))
    ...
    forty_by_forty = mrd.getLevel(2)
    data40 = interpolate(data10, Function(forty_by_forty))

Formulation
-----------

For a single PDE that has a solution with a single component the linear PDE is defined in the following form:

.. math::

      .. _eq:ripleysingle:

   \begin{array}{cl} &
   \displaystyle{
   \int_{\Omega}
   A_{jl} \cdot v_{,j}u_{,l}+ B_{j} \cdot v_{,j} u+ C_{l} \cdot v u_{,l}+D \cdot vu \; d\Omega }
   + \displaystyle{\int_{\Gamma} d \cdot vu \; d{\Gamma} }\\
   = & \displaystyle{\int_{\Omega}  X_{j} \cdot v_{,j}+ Y \cdot v \; d\Omega }
   + \displaystyle{\int_{\Gamma} y \cdot v \; d{\Gamma}}
   \end{array}

Meshes
------

.. _sec:ripleymeshes:

An example 2D mesh from **esys.ripley** is shown in Figure `[fig:ripleyrect] <#fig:ripleyrect>`__. Mesh files cannot be used to generate **esys.ripley** domains, i.e. **esys.ripley** does not have “ReadGmsh()“ or “ReadMesh()“ functions. Instead, **esys.ripley** domains are always created using a call to “Brick()“ or “Rectangle()“, see Section `[sec:ripleyfuncs] <#sec:ripleyfuncs>`__.

.. _fig:ripleyrect:

.. note::

   See the `User Guide PDF <../../user/user.pdf>`_ for mesh figures.

Functions
---------

.. _sec:ripleyfuncs:

.. container:: funcdesc

   Brickn0,n1,n2,l0=1.,l1=1.,l2=1.,d0=-1,d1=-1,d2=-1, diracPoints=list(), diracTags=list(), comm=None generates a “Domain“ object representing a three-dimensional brick between :math:`(0,0,0)` and :math:`(l0,l1,l2)` with orthogonal faces. All elements will be regular. The brick is filled with “n0“ elements along the :math:`x_0`-axis, “n1“ elements along the :math:`x_1`-axis and “n2“ elements along the :math:`x_2`-axis. If built with *MPI* support, the domain will be subdivided “d0“ times along the :math:`x_0`-axis, “d1“ times along the :math:`x_1`-axis, and “d2“ times along the :math:`x_2`-axis. “d0“, “d1“, and “d2“ must be factors of the number of *MPI* processes requested. If axial subdivisions are not specified, automatic domain subdivision will take place. This may not be the most efficient construction and will likely result in extra elements being added to ensure proper distribution of work. Any extra elements added in this way will change the length of the domain proportionately. “diracPoints“ is a list of coordinate-tuples of points within the mesh, each point tagged with the respective string within “diracTags“. “comm“ is an optional *MPI* communicator (from ``mpi4py``) that allows using a custom communicator instead of ``MPI_COMM_WORLD``. If “comm“ is not “None“, **esys.escript** must be built with ``mpi4py`` support enabled.

.. container:: funcdesc

   Rectanglen0,n1,l0=1.,l1=1.,d0=-1,d1=-1, diracPoints=list(), diracTags=list(), comm=None generates a “Domain“ object representing a two-dimensional rectangle between :math:`(0,0)` and :math:`(l0,l1)` with orthogonal faces. All elements will be regular. The rectangle is filled with “n0“ elements along the :math:`x_0`-axis and “n1“ elements along the :math:`x_1`-axis. If built with *MPI* support, the domain will be subdivided “d0“ times along the :math:`x_0`-axis and “d1“ times along the :math:`x_1`-axis. “d0“ and “d1“ must be factors of the number of *MPI* processes requested. If axial subdivisions are not specified, automatic domain subdivision will take place. This may not be the most efficient construction and will likely result in extra elements being added to ensure proper distribution of work. Any extra elements added in this way will change the length of the domain proportionately. “diracPoints“ is a list of coordinate-tuples of points within the mesh, each point tagged with the respective string within “diracTags“. “comm“ is an optional *MPI* communicator (from ``mpi4py``) that allows using a custom communicator instead of ``MPI_COMM_WORLD``. If “comm“ is not “None“, **esys.escript** must be built with ``mpi4py`` support enabled.

The arguments “l0“, “l1“ and “l2“ for “Brick()“ and “Rectangle()“ may also be given as tuples “(x0,x1)“ in which case the coordinates will range between “x0“ and “x1“. For example:

.. code:: python

   from esys.ripley import Rectangle
      dom = Rectangle(10, 10, l0=(5.5, 15.5), l1=(9.0, 14.0))

This will create a rectangle with :math:`10` by :math:`10` elements where the bottom-left node is located at :math:`(5.5, 9.0)` and the top-right node has coordinates :math:`(15.5, 14.0)`, see Figure `[fig:ripleyrect] <#fig:ripleyrect>`__.

The “MultiResolutionDomain“ class is available as a wrapper, taking the dimension of the domain followed by the same arguments as “Brick“ (if a two-dimensional domain is requested, any extra arguments over those used by “Rectangle“ are ignored). All of these standard arguments to “MultiResolutionDomain“ must be supplied as keyword arguments (e.g. “d0“=...). The “MultiResolutionDomain“ can then generate compatible domains for interpolation.

Linear Solvers in “SolverOptions“
---------------------------------

Currently direct solvers and GPU-based solvers are not supported under *MPI* when running with more than one rank. By default, **esys.ripley** uses the iterative solvers "SolverOptions.PCG" for symmetric and "SolverOptions.BICGSTAB" for non-symmetric problems. A GPU will not be used unless explicitly requested via the "setSolverTarget()" method of the solver options. These solvers are only available if **esys.ripley** was built with *CUDA* support. If the direct solver is selected, which can be useful when solving very ill-posed equations, **esys.ripley** uses the "MKL" (Intel Math Kernel Library) solver package. If "MKL" is not available "UMFPACK" is used. If "UMFPACK" is not available a suitable iterative solver from "PASO" is used, but if a direct solver was requested via "SolverOptions" an exception will be raised.

For Complete Reference
======================

For complete mesh figures and additional details, please see the `User Guide PDF <../../user/user.pdf>`_.
