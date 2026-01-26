# esys-escript Installation Guide

This guide is for **Version 6** of esys-escript.

esys-escript is primarily developed on Linux systems. The current version provides some support for macOS but has not been ported to Windows.

For support and questions, please visit the [Issues site](https://github.com/esys-escript/esys-escript.github.io/issues).

For complete documentation including user guide and API reference, see the [main documentation page](https://lutzgross.github.io/esys-escript.github.io/index.html).

## Dependencies

### Essential Dependencies

- C++ compiler (g++, clang++, or similar)
- Python 3
- SCons build system
- CMake
- python3-numpy 2
- python3-matplotlib
- Boost libraries:
  - boost-python
  - boost-random
  - boost-iostreams


### Non-Essential Dependencies
- LAPACKE - Linear algebra library for matrix operations in PASO solver (liblapacke-dev on Debian/Ubuntu)
- HDF5-serial (for HDF5 file I/O) (libhdf5-serial-dev on Debian/Ubuntu)
- SILO (recommended for VisIt visualization) (libsilo-dev on Debian/Ubuntu)
- netCDF4 (for reading gridded data) (libnetcdf-c++4-dev on Debian/Ubuntu)
- 
- boost-numpy (recommended for Jupyter) (libboost-numpy-dev on Debian/Ubuntu)
- UMFPACK - Direct sparse matrix solver (libsuitesparse-dev on Debian/Ubuntu)
- MUMPS - MUltifrontal Massively Parallel sparse direct Solver, sequential version (libmumps-seq-dev on Debian/Ubuntu)
- zlib - Compression library required by p4est for the oxley domain (zlib1g-dev on Debian/Ubuntu)

### Optional Dependencies
- METIS - Graph partitioning library used by Trilinos for matrix reordering (libmetis-dev on Debian/Ubuntu). Can be used with or without MPI.
- MPI implementation (OpenMPI or MPICH) with python3-mpi4py for distributed parallelization
- ParMETIS - Parallel graph partitioning for MPI builds (libparmetis-dev on Debian/Ubuntu, requires MPI and METIS)
- SymPy (python3-sympy) for symbolic mathematics module - **must be present at build time** (requires version 1.2 or later)
- python3-sphinx and python3-markdown (for building API documentation)
- LaTeX (texlive-latex-base, texlive-latex-extra) for building user guide PDF
- CppUnit (for running tests)

The source code is shipped with [Trilinos](https://trilinos.github.io/) and [p4est](https://www.p4est.org/).

### Recommended Tools

These tools are not required at compile time but enhance functionality:

- [gmsh](https://gmsh.info/) - mesh generation
- [VisIt](https://visit-dav.github.io) - visualization via SILO or VTK
- [ParaView](https://www.paraview.org/) - visualization via VTK
- [Mayavi](https://docs.enthought.com/mayavi/mayavi/) - visualization via VTK
- [Matplotlib](https://matplotlib.org/) - 2D plotting

## Parallelization Support

esys-escript supports multiple parallelization strategies:

- **OpenMP** - Shared memory parallelization (multi-core), enabled at compilation
- **MPI** - Distributed memory parallelization (compute clusters), requires mpi4py
- **Hybrid mode** - Both OpenMP and MPI can be used simultaneously

## Installation from Source

### Quick Start

1. Download the source code:

```bash
git clone --single-branch --branch master https://github.com/esys-escript/esys-escript.github.io.git esys6
cd esys6
```

Alternatively, download a tagged release from [releases page](https://github.com/esys-escript/esys-escript.github.io/releases).

2. Install dependencies for your platform (see Platform-Specific Instructions below)

3. Build esys-escript:

```bash
scons -j4 options_file=scons/templates/<OS>_options.py
```

Replace `<OS>` with your operating system (e.g., `debian`, `ubuntu`, `centos8`, etc.).

4. Set up environment variables:

```bash
export PYTHONPATH=$PWD:$PYTHONPATH
export LD_LIBRARY_PATH=$PWD/lib/esys:$PWD/esys.trilinos/lib:$LD_LIBRARY_PATH
```

5. (Optional) Test your build:

```bash
scons -j4 py_tests options_file=scons/templates/<OS>_options.py
```

### Custom Build Configuration

If you need to customize the build (e.g., disable MPI, change compilers, enable/disable dependencies), create your own options file:

```bash
HOST=`uname -n`
cp scons/templates/<closest_match>_options.py scons/${HOST}_options.py
# Edit scons/${HOST}_options.py as needed
scons -j4
```

**For a complete list of all available configuration options, see the [Build Options Reference](./scons/templates/options.md).**

The options file must include `escript_opts_version = 203` to be accepted by the build system.

### Building Documentation

To build the documentation, install additional dependencies:

```bash
sudo apt-get install python3-sphinx python3-markdown
```

For the user guide PDF, also install LaTeX:

```bash
sudo apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

Then build the docs:

```bash
scons -j4 docs
```

This generates:
- Main entry point: `release/doc/index.html`
- User guide (PDF): `release/doc/user/user.pdf`
- Python API docs (Sphinx): `release/doc/sphinx_api/`
- Example scripts: `release/doc/escript_examples.{zip,tar.gz}`

## Platform-Specific Instructions

### Debian 12+

Install required packages:

```bash
sudo apt-get install python3-dev python3-numpy python3-scipy python3-matplotlib
sudo apt-get install g++ gfortran scons cmake
sudo apt-get install libboost-numpy-dev libboost-python-dev libboost-random-dev libboost-iostreams-dev
sudo apt-get install libhdf5-serial-dev libsilo-dev libnetcdf-dev libsuitesparse-dev liblapacke-dev libmumps-seq-dev zlib1g-dev
sudo apt-get install libmetis-dev
```

For MPI support, additionally install:

```bash
sudo apt-get install python3-mpi4py libparmetis-dev
```

For building documentation, additionally install:

```bash
sudo apt-get install python3-sphinx python3-markdown
sudo apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

For symbolic mathematics support, additionally install:

```bash
sudo apt-get install python3-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/debian_options.py
```

### Ubuntu 24.04+

Install required packages:

```bash
sudo apt-get install python3-dev python3-numpy python3-scipy python3-matplotlib
sudo apt-get install g++ gfortran scons cmake
sudo apt-get install libboost-numpy-dev libboost-python-dev libboost-random-dev libboost-iostreams-dev
sudo apt-get install libhdf5-serial-dev libsilo-dev libnetcdf-dev libsuitesparse-dev liblapacke-dev libmumps-seq-dev zlib1g-dev
sudo apt-get install libmetis-dev
```

For MPI support, additionally install:

```bash
sudo apt-get install python3-mpi4py libparmetis-dev
```

For building documentation, additionally install:

```bash
sudo apt-get install python3-sphinx python3-markdown
sudo apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

For symbolic mathematics support, additionally install:

```bash
sudo apt-get install python3-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/ubuntu_options.py
```

### Linux Mint 20.3+

Mint Linux is based on Ubuntu, use the same commands as Ubuntu above.

### Arch Linux

Install required packages:

```bash
sudo pacman -Sy python python-numpy python-scipy python-matplotlib
sudo pacman -Sy gcc scons cmake
sudo pacman -Sy boost boost-libs suitesparse hdf5 netcdf lapack mumps zlib metis
```

For MPI support:

```bash
sudo pacman -Sy python-mpi4py parmetis
```

For building documentation:

```bash
sudo pacman -Sy python-sphinx python-markdown
sudo pacman -Sy texlive-core texlive-latexextra
```

For symbolic mathematics support:

```bash
sudo pacman -Sy python-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/arch_py3_options.py
```

### Fedora

Install required packages:

```bash
sudo dnf install python3-devel python3-numpy python3-scipy python3-matplotlib
sudo dnf install gcc-c++ gcc-gfortran scons cmake
sudo dnf install boost-devel boost-python3-devel boost-python3 boost-numpy3 boost-iostreams boost-random
sudo dnf install hdf5-devel netcdf-devel suitesparse-devel lapack-devel MUMPS-devel zlib-devel metis-devel
```

For MPI support:

```bash
sudo dnf install python3-mpi4py parmetis-devel
```

For building documentation:

```bash
sudo dnf install python3-sphinx python3-markdown
sudo dnf install texlive-scheme-basic texlive-latex texlive-latex-extra
```

For symbolic mathematics support:

```bash
sudo dnf install python3-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/fedora_py3_options.py
```

### OpenSUSE

Install required packages:

```bash
sudo zypper in python3-devel python3-numpy python3-scipy python3-matplotlib
sudo zypper in gcc gcc-c++ gcc-fortran scons cmake
sudo zypper in libboost_python3-devel libboost_numpy3-devel libboost_random-devel libboost_iostreams-devel
sudo zypper in hdf5-devel netcdf-devel suitesparse-devel lapack-devel mumps-devel zlib-devel metis-devel
```

For MPI support:

```bash
sudo zypper in python3-mpi4py parmetis-devel
```

For building documentation:

```bash
sudo zypper in python3-Sphinx python3-Markdown
sudo zypper in texlive-latex texlive-latex-extra
```

For symbolic mathematics support:

```bash
sudo zypper in python3-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/opensuse_py3_options.py
```

### CentOS 8 / Rocky Linux / AlmaLinux

Enable required repositories:

```bash
sudo dnf install epel-release
sudo dnf config-manager --set-enabled powertools
```

Install required packages:

```bash
sudo dnf install python3-devel python3-numpy python3-scipy python3-matplotlib
sudo dnf install gcc gcc-c++ gcc-gfortran scons cmake
sudo dnf install boost-devel boost-python3 boost-python3-devel
sudo dnf install hdf5-devel netcdf-devel suitesparse suitesparse-devel lapack-devel MUMPS-devel zlib-devel metis-devel
```

For MPI support:

```bash
sudo dnf install python3-mpi4py parmetis-devel
```

For building documentation:

```bash
sudo dnf install python3-sphinx python3-markdown
sudo dnf install texlive-scheme-basic texlive-latex
```

For symbolic mathematics support:

```bash
sudo dnf install python3-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/centos8_options.py
```

### macOS with Homebrew

**Note:** macOS support is limited. Requires Xcode Command Line Tools.

Install Homebrew from [https://brew.sh](https://brew.sh), then:

```bash
brew install python3 numpy scipy matplotlib
brew install scons cmake
brew install boost boost-python3
brew install hdf5 netcdf suite-sparse lapack mumps
```

Some Python packages may need pip:

```bash
pip3 install numpy scipy matplotlib
```

For MPI support:

```bash
brew install open-mpi
pip3 install mpi4py
```

For building documentation:

```bash
pip3 install sphinx markdown
brew install --cask mactex  # For LaTeX/PDF documentation
```

For symbolic mathematics support:

```bash
pip3 install sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/homebrew_options.py
```

### macOS with MacPorts

**Note:** macOS support is limited. Requires Xcode Command Line Tools.

Install MacPorts from [https://www.macports.org](https://www.macports.org), then:

```bash
sudo port install python311
sudo port select --set python python311
sudo port select --set python3 python311
sudo port install py311-numpy py311-scipy py311-matplotlib
sudo port install scons cmake boost hdf5 netcdf suitesparse lapack mumps
```

For MPI support:

```bash
sudo port install openmpi py311-mpi4py
```

For building documentation:

```bash
sudo port install py311-sphinx py311-markdown
sudo port install texlive  # For LaTeX/PDF documentation
```

For symbolic mathematics support:

```bash
sudo port install py311-sympy
```

Build esys-escript:

```bash
scons -j4 options_file=scons/templates/macports_options.py
```

### FreeBSD

**Note:** At the time of testing, numpy installation on FreeBSD had issues. Full testing has been limited.

```bash
sudo pkg install python3 py39-numpy py39-scipy py39-matplotlib
sudo pkg install scons cmake boost-python-libs hdf5 netcdf suitesparse lapack mumps
```

For building documentation:

```bash
sudo pkg install py39-sphinx py39-markdown
sudo pkg install texlive-full  # For LaTeX/PDF documentation
```

For symbolic mathematics support:

```bash
sudo pkg install py39-sympy
```

## Advanced Configuration

### Customizing Build Options

Edit your options file (`scons/<hostname>_options.py`) to customize:

- **OpenMP**: Enable/disable with `openmp = True` or `openmp = False`
- **MPI**: Set `mpi = 'OPENMPI'` or `mpi = 'none'`
- **SymPy**: Enable/disable symbolic module with `sympy = True` or `sympy = False`
- **Compiler**: Change `cxx`, `cc` variables
- **Optimization**: Modify `cxx_extra` flags
- **Libraries**: Adjust library paths and names

### MPI Configuration

#### Basic MPI Setup

To enable MPI with auto-detection (recommended):

```python
mpi = 'auto'       # Auto-detect MPI implementation from mpi4py
mpi4py = True
```

When `mpi='auto'` is set:
- If `mpi4py=True`: Automatically detects the MPI implementation (OpenMPI, MPICH, Intel MPI) from your installed mpi4py package
- If `mpi4py=False`: Sets `mpi='none'` (MPI disabled)

To manually specify the MPI implementation:

```python
mpi = 'OPENMPI'  # or 'MPICH', 'MPICH2', 'INTELMPI'
mpi_prefix = '/usr/lib/x86_64-linux-gnu/openmpi'
mpi_libs = ['mpi_cxx', 'mpi']
```

**Important:** When using `mpi4py`, the MPI flavour must match the MPI implementation that mpi4py was compiled against. The build system will verify compatibility and report an error if there's a mismatch.

To disable MPI:

```python
mpi = 'none'
```

#### Custom MPI Communicators with mpi4py

When `mpi4py = True` is enabled, all domain factory functions (e.g., `Rectangle`, `Brick`, `ReadMesh`, `ReadGmsh`) accept an optional `comm` parameter that takes an mpi4py communicator object:

```python
from mpi4py import MPI
from esys.ripley import Rectangle

# Use a custom communicator instead of MPI_COMM_WORLD
custom_comm = MPI.COMM_WORLD.Split(color=..., key=...)
domain = Rectangle(10, 10, comm=custom_comm)
```

If the `comm` parameter is not provided or is `None`, `MPI_COMM_WORLD` is used by default.

### SymPy Symbolic Module

**Note:** SymPy is a **runtime-only dependency** - it is not needed for compiling C++ code. However, **SymPy must be installed at build time** when `sympy = True` to enable the symbolic module. The build system checks for SymPy availability and sets a feature flag. If SymPy is not found at build time, the symbolic support will be disabled even if SymPy is installed later.

To enable symbolic mathematics support (requires SymPy 1.2 or later):

```python
sympy = True
```

To disable symbolic support:

```python
sympy = False  # Default
```

The symbolic module provides support for symbolic expressions, automatic differentiation, and symbolic PDE formulation. See the user guide for details on using the symbolic toolbox.

**Important:** If you build without SymPy and later install it, you must rebuild with `sympy = True` to enable symbolic support.

### Trilinos Support

esys-escript includes Trilinos for advanced solver capabilities. Trilinos is built automatically during the esys-escript build process.

For complex-valued PDEs and advanced preconditioners, ensure Trilinos support is enabled in your build.

## Testing


### Unit Tests

Run the test suite:

```bash
scons -j4 py_tests
```

Or use the test script:

```bash
./utest.sh build -t4  # Run with 4 OpenMP threads
```

With MPI:

```bash
./utest.sh build -t2 -n2  # 2 threads, 2 MPI processes
```

**Note:** Some tests require optional features (netCDF, SILO, etc.) and will report failures if those features are disabled.

## Installation Layout

After building, the directory structure is:

```
esys-escript.github.io/
├── bin/
│   └── run-escript              # Launcher script
├── esys/                        # Python package
│   ├── escript/
│   ├── finley/
│   ├── ripley/
│   ├── oxley/
│   └── speckley/
├── lib/
│   └── esys/                    # Compiled libraries
├── esys.trilinos/               # Trilinos installation
│   └── lib/
├── release/                     # Documentation (if built)
│   └── doc/
│       ├── index.html
│       └── sphinx_api/
└── README.md
```

## Running esys-escript

### Using the Launcher

```bash
./bin/run-escript myscript.py
```

The launcher automatically sets up the environment (PYTHONPATH, LD_LIBRARY_PATH, etc.).

### Parallel Execution

With OpenMP (4 threads):

```bash
./bin/run-escript -t 4 myscript.py
```

With MPI (4 processes):

```bash
./bin/run-escript -n 4 myscript.py
```

Hybrid OpenMP + MPI (2 processes, 4 threads each):

```bash
./bin/run-escript -n 2 -t 4 myscript.py
```

### Direct Python Usage

If not using the launcher, set environment variables:

```bash
export PYTHONPATH=/path/to/esys-escript.github.io:$PYTHONPATH
export LD_LIBRARY_PATH=/path/to/esys-escript.github.io/lib/esys:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/esys-escript.github.io/esys.trilinos/lib:$LD_LIBRARY_PATH

python3 myscript.py
```

## Troubleshooting

### Common Issues

**Import errors:**
- Ensure `PYTHONPATH` and `LD_LIBRARY_PATH` are set correctly
- Check that all dependencies are installed

**Build failures:**
- Verify all required packages are installed
- Check compiler compatibility (C++11 support required)
- Review `config.log` for detailed error messages

**Test failures:**
- Some tests require optional features (HDF5, SILO, MPI)
- Tests will report failures for disabled features - this is expected

**MPI issues:**
- Ensure mpi4py version matches your MPI implementation
- Check that MPI paths in options file are correct
- Verify MPI is properly initialized with `mpirun` or `mpiexec`

### Getting Help

- Report issues: [GitHub Issues](https://github.com/esys-escript/esys-escript.github.io/issues)
- Documentation: [https://esys-escript.github.io/](https://esys-escript.github.io/)
- User guide: See `release/doc/user/user.pdf` after building docs

## Optional Extras

### Direct Solvers

To enable UMFPACK (direct sparse solver):

```bash
sudo apt-get install libsuitesparse-dev
```

Add to your options file:

```python
umfpack = True
```

### Visualization

Install visualization tools:

```bash
# VisIt (via SILO format)
sudo apt-get install visit

# ParaView (via VTK format)
sudo apt-get install paraview

# Mayavi (via VTK format)
pip3 install mayavi
```

### Mesh Generation

Install gmsh for mesh generation:

```bash
sudo apt-get install gmsh  # Linux
brew install gmsh          # macOS
```

## Compiler Requirements

esys-escript requires a C++11 compliant compiler. Tested compilers include:

- GCC g++ 10.2+
- Clang 11.0+
- Intel icpc v17+

Required C++11 features:
- `std::complex<>`
- `auto` type declarations
- `decltype(T)` type declarations
- Explicit template instantiation
- `std::isnan()` in `std` namespace

## License and Credits

See `LICENSE` and `CREDITS` files in the distribution for license information and contributor acknowledgments.