
# *esys-escript* 

*esys-escript* is a module for implementing mathematical models in Python using the finite element method (FEM).  As users do not access the underlying data structures it is very easy to use and scripts can run on desktop computers as well as massive parallel supercomputers without changes. Application areas for esys-escript include geophysical inversion, earthquakes, porous media flow, reactive transport, plate subduction, erosion, earth mantle convection, and tsunamis.

esys-escript is designed as an easy-to-use environment for implementing mathematical models based on non-linear, coupled, time-dependent partial differential equations. It uses the finite element method (FEM) for spatial discretization and data representation and is used through Python. It is suitable for rapid prototyping (e.g. for a student project or thesis) as well as for large software projects. Scripts are executed sequentially, on multi-core platforms via OpenMP and distributed computing clusters using MPI. The hybrid mode of OpenMP and MPI is supported and allows for solving problems with over 200 million unknowns on several thousand cores on a parallel computer.

For geophyscial inversion see also the extensions [gambit](https://github.com/AndreaCodd/gambit) and [fingal](https://github.com/LutzGross/fingal).


## Main Features:

- python based user interface
- two- and three-dimensional finite and spectral element simulations
- specialized geophysical inversion module
- support for VTK and SILO file format
- unstructured meshes from gmsh
- parallelization with OpenMP and MPI support
- Flux Controlled Transport solver (FEM-FCT)
- visualization with VisIt
- platform is Linux; there is limited support for MacOS and Windows  

Further documentation including examples and a user guide for the latest release can be found at
https://esys-escript.github.io/


## Questions & Bugs 

To raise a question or to report a bug please start a [github issue](https://github.com/esys-escript/esys-escript.github.io/issues).

## Linux or MacOS Installation from Source

Please inspect the [installation guide](installation.md) for information on a specific Linux distribution.

## Debian Distribution 

Debian packages [python3-escript-mpi] and [python3-escript] for **version 5** are available.

## Anaconda Installation **Version 5** (needs validation)

To install *esys-escript* for [anaconda](https://www.anaconda.com), first run `conda` and then the command

conda install esys-escript -c conda-forge

At present, this is the recommended way to run esys-escript on Windows.

## Using *esys-escript*

To get started using *esys-escript*, you can use the [user guide](./user.pdf). 

Click on this button to try out *esys-escript* in your web browser (needs validation):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/esys-escript/esys-escript.github.io/HEAD)


## The project was funded by the
   - [AuScope National Collaborative Research Infrastructure Strategy (NCRIS)](https://www.auscope.org.au/) (until end of 2023),
   - Australian Geophysical Observing System (AGOS) (ended 2014),
   - [The University of Queensland](https://www.uq.edu.au)

If you publish work that makes use of esys-escript, we would appreciate it if you would cite the following reference:

- [R Schaa, L Gross, J du Plessis, PDE-based geophysical modelling using finite elements: examples from 3D resistivity and 2D magnetotellurics, Journal of Geophysics and Engineering, Volume 13, Issue 2, April 2016, Pages S59–S73](https://doi.org/10.1088/1742-2132/13/2/S59)
- [L. Gross, L. Bourgouin, A.J. Hale, H.-B. Mühlhaus,
Interface modeling in incompressible media using level sets in Escript,
Physics of the Earth and Planetary Interiors,
Volume 163, Issues 1–4, 2007,
Pages 23-34](doi:10.1016/j.pepi.2007.04.004)


### Contributors
        Lutz Gross
        Adam Ellery
        Andrea Codd
        Joel Fenwick
        Cihan Altinay
        Simon Shaw
        Jaco Du Plessis
        Ralf Schaa
        Peter Hornby
        Thomas Poulet
        Lin Gao
        Artak Amirbekyan
        Ken Steube

