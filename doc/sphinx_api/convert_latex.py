#!/usr/bin/env python3
"""
Create stub RST files for user guide chapters that link to the PDF.

This script creates simple RST files that provide a brief description
and link to the full PDF user guide for detailed content.
"""

import os
import sys
from pathlib import Path


# Chapter definitions with descriptions
CHAPTERS = {
    'tutorial_pde': {
        'name': 'Tutorial: Solving PDEs',
        'description': '''This tutorial introduces the basic concepts of solving partial differential
equations (PDEs) with esys-escript. It covers installation, first steps with
the Poisson equation, time-dependent diffusion problems, wave equations,
elastic deformation, and working with unstructured meshes.'''
    },
    'execute': {
        'name': 'Execution of Scripts',
        'description': '''This chapter explains how to run escript simulations using the ``run-escript``
launcher. It covers command-line options, MPI parallelization, OpenMP threading,
environment variables, and output redirection.'''
    },
    'escript': {
        'name': 'The escript Module',
        'description': '''The core escript module provides the fundamental data structures and operations
for finite element computations. This chapter covers Data objects, function spaces,
mathematical operations, and the symbolic toolbox.'''
    },
    'linear_pde': {
        'name': 'The linearPDEs Module',
        'description': '''This chapter describes the linearPDEs module for solving linear partial
differential equations. It covers the LinearPDE class, coefficient specification,
boundary conditions, solver options, and specialized PDE classes like Poisson,
Helmholtz, and Lame equations.'''
    },
    'finley': {
        'name': 'The finley Module',
        'description': '''The finley module provides unstructured finite element meshes. This chapter
covers mesh generation with Rectangle and Brick, reading mesh files (Gmsh, Fly),
element types, and the mathematical formulation of PDEs in finley.'''
    },
    'ripley': {
        'name': 'The ripley Module',
        'description': '''The ripley module provides structured rectangular meshes optimized for
regular grids. It supports fast assemblers for specific PDE types and
GPU-based solvers. This chapter covers mesh generation and the ripley
formulation.'''
    },
    'speckley': {
        'name': 'The speckley Module',
        'description': '''The speckley module implements spectral element methods on structured grids.
This chapter covers the spectral element formulation and mesh generation
with Rectangle and Brick.'''
    },
    'weipa': {
        'name': 'The weipa Module',
        'description': '''The weipa module handles data export and visualization. It supports VTK and
SILO file formats, and can interface directly with VisIt for in-situ
visualization. This chapter covers the EscriptDataset class and export
functions.'''
    },
    'trilinos': {
        'name': 'Using Trilinos',
        'description': '''This chapter explains how to use the Trilinos solver library with escript.
It covers solver configuration, preconditioners, and advanced settings
using MueLu XML parameter files.'''
    },
    'symbolic': {
        'name': 'Symbolic Toolbox',
        'description': '''The symbolic toolbox enables symbolic manipulation of PDEs using SymPy.
This chapter covers the Symbol class, the Evaluator for efficient evaluation,
and the NonlinearPDE class for solving nonlinear problems.'''
    },
    'appendix': {
        'name': 'Appendix',
        'description': '''The appendix contains reference material including installation instructions,
element type specifications, and additional technical details.'''
    },
    'mpi4py': {
        'name': 'Using mpi4py with escript',
        'description': '''This chapter explains how to use mpi4py to create MPI sub-communicators
for multi-domain simulations. It covers the MPIDomainArray and DataCoupler
classes for managing parallel domain decomposition.'''
    },
}


def create_stub_content(chapter_name, chapter_key, description):
    """Create a stub RST file that links to the PDF."""
    underline = '=' * len(chapter_name)

    return f"""{underline}
{chapter_name}
{underline}

{description}

.. note::

   For complete documentation including mathematical derivations, figures,
   and detailed examples, please see the `User Guide PDF <../user/user.pdf>`_.

"""


def convert_chapter(latex_dir, output_dir, chapter_key, chapter_info, verbose=False):
    """Create a stub RST file for a chapter."""
    chapter_name = chapter_info['name']
    description = chapter_info.get('description', f'This chapter covers {chapter_name.lower()}.')

    if verbose:
        print(f"Creating stub for {chapter_key}: {chapter_name}")

    # Create stub content
    rst_content = create_stub_content(chapter_name, chapter_key, description)

    # Write output
    output_file = output_dir / f'{chapter_key}.rst'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(rst_content)

    if verbose:
        print(f"  -> {output_file}")

    return output_file


def main():
    """Main entry point for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(description='Create stub RST files for user guide chapters')
    parser.add_argument('--latex-dir', type=Path, default=Path('doc/user'),
                        help='Directory containing LaTeX source files')
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='Directory for output RST files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress messages')

    args = parser.parse_args()

    # Create output directory if needed
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.verbose:
        print(f"Creating stub RST files in {args.output_dir}")

    # Process each chapter
    for chapter_key, chapter_info in CHAPTERS.items():
        # Skip mpi4py - it has its own manually written RST file
        if chapter_key == 'mpi4py':
            continue
        convert_chapter(args.latex_dir, args.output_dir, chapter_key, chapter_info, args.verbose)

    if args.verbose:
        print(f"\nCreated {len(CHAPTERS) - 1} stub files")


if __name__ == '__main__':
    main()
