# esys-escript

esys-escript is module for implementing mathematical models in python using the finite element method (FEM). 
As users do not access the underlying data structures it is very easy to use and scripts can run on desktop computers as well as massive 
parallel supercomputer without changes. Application areas for esys-escript include geophysical inversion, earthquakes, porous media flow, reactive transport, plate subduction, erosion, earth mantle convection, and tsunamis.

esys-escript is designed as an easy-to-use environment for implementing mathematical models based on non-linear, coupled, time-dependent partial differential equations. It uses the finite element method (FEM) for spatial discretization and data representation and is used through python.
It is suitable for rapid prototyping (e.g for a student project or thesis) as well as for large software projects. Scripts are executed
sequentially, on multi-core platforms via OpenMP and  
distributed computing clusters using MPI. Hybrid mode of OpenMP and MPI is supportied and allow for 
solving problems with over 200 million unknowns on several thousand cores on a parallel computer.

Esys-escript now includes the esys.downunder module for 3D inversion of geophysical data sets. 
The current version supports gravity, magnetic and joint inversion. 


Main Features:

- python based user interface
- two- and three-dimensional finite and spectral element simulations
- specialized geophysical inversion module
- support for VTK and SILO file format
- unstructured meshes from gmsh
- parallelization with OpenMP and MPI support
- Flux Controlled Transport solver (FEM-FCT)
- visualization with VisIt
- support for Linux, Windows and OSX

Further documentation including examples and a user guide for the latest release can be found at
https://esys-escript.readthedocs.io/en/latest/index.html

The project is funded by the
   - AuScope National Collaborative Research Infrastructure Strategy (NCRIS) (until end of 2019)
   - Australian Geophysical Observing System (AGOS) (ended 2014).
   - School of Earth Sciences at the University of Queensland.

If you publish work which makes use of escript, we would appreciate if you would cite the following reference:

[R Schaa, L Gross, J du Plessis, PDE-based geophysical modelling using finite elements: examples from 3D resistivity and 2D magnetotellurics, Journal of Geophysics and Engineering, Volume 13, Issue 2, April 2016, Pages S59–S73]{https://doi.org/10.1088/1742-2132/13/2/S59}
@article{SchaaGrossDuPlessis2016,
    author={Schaa, R. and Gross, L. and Du Plessis, J.},
    year =2016,
    title={PDE-based geophysical modelling using finite elements: examples from 3D resistivity and 2D magnetotellurics},
    journal = {Journal of Geophysics and Engineering},
    volume=13,
    issue=2,
   pages={S59-S73},
  doi = {doi:10.1088/1742-2132/13/2/S59
}

or

@article{GROSS2006,
        author = {L. Gross and L. Bourgouin and A. J. Hale and H.-B Muhlhaus},
        title = {Interface Modeling in Incompressible Media using Level Sets in Escript},
        journal = {Physics of the Earth and Planetary Interiors},
        year = 2007,
        volume = {163},
        pages = {23--34},
        month = {Aug.},
        doi = {doi:10.1016/j.pepi.2007.04.004},
}

Contributors
        Lutz Gross
        Joel Fenwick
        Adam Ellery
        Andrea Codd
        Cihan Altinay
        Simon Shaw
        Jaco Du Plessis
        Ralf Schaa
        Peter Hornby
        Thomas Poulet
        Lin Gao
        Artak Amirbekyan
        Ken Steube

# Windows Installation

A windows build of esys-escript is available in the conda-forge repository. At present this is the recommended way to run esys-escript on Windows. To install, first run conda and then the command

conda install esys-escript -c conda-forge

# Linux Installation

For the impatient:

Install at least g++, python, scons, boost, numpy
READ the file scons/templates/README_FIRST
Copy a suitable template options file from scons/templates/ to scons/hostname_options.py and modify as required.
type: scons to build escript

For information on a specific Linux distrobution, please consult the install guide (install.pdf).

# Using escript

To get started using escript please consult the user guide (user.pdf) and the (cookbook.pdf) cookbook. 
All of these documents are available here and at https://esys-escript.github.io/

# Bugs

To report a bug, please start a github issue.

This project is supported by AuScope
https://github.com/AuScope
