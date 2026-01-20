# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

esys-escript is a Python-based finite element method (FEM) module for implementing mathematical models. It supports sequential, OpenMP (multi-core), and MPI (distributed) execution modes with a hybrid OpenMP+MPI option. The project is primarily C++ with Python bindings, targeting geophysical modeling applications including inversion, earthquakes, porous media flow, and mantle convection.

## Build System

The project uses **SCons** as its build system, not Make or CMake for the main project.

### Building the Project

1. **Create an options file**: Copy a template from `scons/templates/` to `scons/<hostname>_options.py`
   - Get hostname: `uname -n`
   - Available templates: debian_options.py, ubuntu_options.py, bunya_options.py, etc.

2. **Build commands**:
   ```bash
   # Standard build (uses hostname-based options file)
   scons -j4

   # Build with specific options file
   scons -j4 options_file=scons/templates/debian_options.py

   # Build with tests
   scons -j4 build_full

   # Build and run Python tests
   scons -j4 py_tests

   # Install (default target if all domains are built)
   scons install
   ```

3. **Important build options** (configured in options file):
   - `openmp`: Enable OpenMP parallelization
   - `mpi`: MPI flavor ('OPENMPI', 'MPICH', 'INTELMPI', or 'none')
   - `trilinos`: Enable Trilinos solver library
   - `build_trilinos`: Build Trilinos from bundled source ('make', 'always', 'never')
   - `paso`: Build Paso solver library
   - `domains`: List of domain libraries to build ('finley', 'oxley', 'ripley', 'speckley')
   - `debug`: Enable debug flags
   - `prefix`: Installation prefix (default: current directory)

### Running Tests

```bash
# Run all tests using the generated test script
sh utest.sh /path/to/build -t8 -n1 -p1

# Run sanity check after install
scons sanity

# Run Python tests only
scons py_tests
```

### Building Documentation

The project generates multiple types of documentation:

1. **User Guide (PDF from LaTeX)**:
   ```bash
   # Build user guide PDF from doc/user/
   scons user_pdf

   # Build all PDF documentation
   scons pdfdocs
   ```
   - Source: `doc/user/*.tex`
   - Requires: pdflatex, bibtex, makeindex
   - Output: `<prefix>/release/doc/user/user.pdf`

2. **Python API Documentation (Sphinx)** - Integrated Web Documentation:
   ```bash
   # Build Sphinx-based API documentation with integrated resources
   scons sphinxdoc
   ```
   - Uses Sphinx with autodoc to extract Python docstrings
   - Source: Python modules in `esys/` after installation
   - Script: `doc/sphinx_api/genrst.py` generates RST files from installed modules
   - Requires: sphinx, sphinx-build
   - Configuration: `doc/sphinx_api/conf.py`
   - Output: `<prefix>/release/doc/sphinx_api/` (HTML)
   - Optional: Set `mathjax_path` in options file for math rendering

   **Integrated Documentation Features**:
   - The Sphinx documentation now provides a unified landing page integrating:
     * Python API reference (auto-generated from docstrings)
     * User Guide PDF (linked from `../user/user.pdf`)
     * Example scripts (linked as ZIP and TAR.GZ archives)
   - Custom styling via `doc/sphinx_api/_static/custom.css`
   - The build process automatically copies the user guide PDF and example archives
     into the documentation tree for easy access
   - Start page: `<prefix>/release/doc/sphinx_api/index.html`

3. **C++ API Documentation (Doxygen)**:
   ```bash
   # Build Doxygen documentation for C++ classes
   scons api_doxygen
   ```
   - Requires: doxygen
   - Configuration: `doc/doxygen/doxygen_esys`
   - Output: `<prefix>/release/doc/doxygen/` (HTML)

4. **Example Files**:
   ```bash
   # Package examples as archives
   scons examples_tarfile    # Creates .tar.gz
   scons examples_zipfile    # Creates .zip
   ```
   - Source: `doc/examples/` (usersguide and cookbook examples)
   - Output: `<prefix>/release/doc/escript_examples.{tar.gz,zip}`

5. **Build All Documentation**:
   ```bash
   # Build base documentation (PDF, examples, doxygen)
   scons basedocs

   # Build all documentation (base + sphinx)
   scons docs

   # Build everything for release (docs + install)
   scons release_prep
   ```

**Documentation Directory Structure**:
```
doc/
├── user/          # User guide LaTeX source
├── sphinx_api/    # Sphinx configuration and generation script
├── doxygen/       # Doxygen configuration
└── examples/      # Example scripts and tests
    ├── usersguide/
    └── cookbook/
```

**Important Notes**:
- Sphinx documentation requires the project to be installed first (depends on `env['pyinstall']`)
- The documentation build targets use `AlwaysBuild()`, so they rebuild each time
- PDF generation uses custom PDFLaTeX flags to embed the git revision number
- Example files are tested as part of the build process (see `doc/SConscript`)

### Environment Setup

After building, you need to set environment variables:
```bash
export PYTHONPATH=$ESYSESCRIPT
export LD_LIBRARY_PATH=$ESYSESCRIPT/lib/esys:$ESYSESCRIPT/esys.trilinos/lib
```

Where `ESYSESCRIPT` is the project root directory.

## Code Architecture

### Domain Structure

The project is organized into several independent **domain** implementations, each providing different FEM discretization methods:

- **finley**: Traditional unstructured FEM
- **oxley**: Octree-based adaptive mesh refinement (uses p4est)
- **ripley**: Regular rectangular grids
- **speckley**: Spectral element methods

Each domain follows the same structure:
```
<domain>/
├── src/           # C++ implementation
├── py_src/        # Python bindings and interfaces
└── test/          # Unit tests (C++ and Python)
```

### Core Libraries

- **escriptcore**: Core data structures and abstractions
  - `Data` objects: Distributed field data on FEM meshes
  - `Domain` abstract interface
  - Function spaces and domain couplers

- **escript**: High-level PDE solving interface
  - `LinearPDE` and `NonLinearPDE` classes
  - PDE tools and utilities
  - Cost functions and minimizers

- **paso**: Native iterative solver library (C++)
  - Sparse matrix operations
  - Iterative solvers (PCG, BiCGStab, GMRES)
  - Preconditioners

- **trilinoswrap**: Wrapper for Trilinos solvers
  - Provides access to Trilinos linear algebra
  - Only for single PDEs (not systems) due to block matrix limitations

- **weipa**: Data export library
  - Supports SILO and VTK formats
  - Visit simulation interface

- **pythonMPI**: MPI wrapper for Python (when MPI is enabled)

### Build Products

- Python modules installed to: `<prefix>/esys/`
- Shared libraries installed to: `<prefix>/lib/esys/`
- Executables installed to: `<prefix>/bin/`
  - `run-escript`: Main launcher script (handles MPI/OpenMP setup)
- Trilinos libraries (if built): `<prefix>/esys.trilinos/lib/`

### Configuration System

The build uses a two-stage configuration:

1. **Options file** (`scons/*_options.py`): User-specified build options
   - Must set `escript_opts_version = 203` (current version)
   - Configures compilers, libraries, features

2. **Runtime configuration** (`<prefix>/lib/esys/buildvars`): Generated at build time
   - Contains build flags and paths
   - Read by `run-escript` launcher

### Trilinos Integration

The project includes Trilinos source code (version 16) in `trilinos_source16/` and can build it automatically:

- Set `build_trilinos = "make"` in options file to build Trilinos
- Use `scripts/trilinos_mpi.sh` or `scripts/trilinos_nompi.sh` for build configuration
- Build happens before main esys-escript build
- Trilinos is installed to `esys.trilinos/` directory

### Python-C++ Interface

The project uses Boost.Python for Python bindings:
- C++ classes wrapped in `py_src/` directories
- Each domain provides Python factory functions
- Data objects can be manipulated from Python but computation happens in C++

## Common Development Patterns

### Adding a new test

1. C++ tests: Add to `<domain>/test/` and update corresponding `SConscript`
2. Python tests: Add to `<domain>/test/python/` with `SConscript` entry

### Modifying build options

1. Edit your options file in `scons/<hostname>_options.py`
2. Rebuild with `scons -j4`
3. For Trilinos options, modify `scripts/trilinos_*.sh` scripts

### Working with MPI

- MPI must be enabled at build time (`mpi='OPENMPI'` etc.)
- **ALWAYS use `run-escript` launcher for MPI programs** - DO NOT use `mpirun python3` directly
- Launcher flags: `-n` (processes), `-p` (processes per node), `-t` (threads)

**Critical**: The `pythonMPI` wrapper (used by `run-escript`) initializes MPI *before* Python starts, which is essential for proper Trilinos/Tpetra initialization. Using `mpirun python3` directly causes MPI initialization order issues that lead to index errors and crashes.

**Sub-communicators with mpi4py**:
- You CAN use mpi4py to create sub-communicators and pass them to escript domains
- This replaces the old SplitWorld functionality
- Example: `sub_comm = MPI.COMM_WORLD.Split(color, rank); domain = Rectangle(n0=20, n1=10, comm=sub_comm)`
- See `doc/user/mpi4py.tex` for detailed documentation and examples

### Compiler Selection

The build system detects compilers and applies appropriate flags:
- GNU: g++ (default flags in SConstruct lines 384-398)
- Intel: icpc (lines 355-369)
- Clang: clang++ (lines 412-427)
- MSVC: cl/icl for Windows builds

## Key Files

- `SConstruct`: Main build script (reads options, configures environment, builds all modules)
- `scons/site_init.py`: Build system initialization and helper functions
- `scons/dependencies.py`: Dependency checking functions
- `bin/run-escript`: Main launcher (generated during build)
- `scripts/trilinos_mpi.sh` / `scripts/trilinos_nompi.sh`: Trilinos build scripts
- `doc/SConscript`: Documentation build orchestration
- `doc/sphinx_api/genrst.py`: Auto-generates RST files from installed Python modules
- `doc/sphinx_api/conf.py`: Sphinx configuration for API documentation
- `doc/user/user.tex`: Main LaTeX source for user guide

## Important Notes

- The build system assumes you will NOT move Python after building
- Library paths are set at build time in the launcher script
- Each domain is largely independent but shares the core abstractions
- Changes to C++ code require rebuild, Python-only changes do not
- The project uses C++17 standard