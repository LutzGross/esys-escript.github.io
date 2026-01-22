#!/usr/bin/env python3
"""
Convert LaTeX user guide chapters to reStructuredText for Sphinx.

This script:
1. Pre-processes LaTeX to expand escript-specific macros
2. Converts using pandoc
3. Post-processes to fix RST formatting
"""

import os
import re
import subprocess
import sys
from pathlib import Path

# Macro replacements - simple text substitutions
MACRO_REPLACEMENTS = {
    # Software names
    r'\\escript': '**esys.escript**',
    r'\\esys': '**esys**',
    r'\\finley': '**esys.finley**',
    r'\\ripley': '**esys.ripley**',
    r'\\speckley': '**esys.speckley**',
    r'\\oxley': '**esys.oxley**',
    r'\\weipa': '**esys.weipa**',
    r'\\linearPDEs': '**esys.escript.linearPDEs**',
    r'\\pdetools': '**esys.escript.pdetools**',
    r'\\numpy': '**numpy**',
    r'\\numpyNDA': '**numpy.ndarray**',

    # Tools and libraries
    r'\\PYTHON': '*Python*',
    r'\\pythonthree': '*python3*',
    r'\\pythontwo': '*python2*',
    r'\\MPI': '*MPI*',
    r'\\mpifo': '*mpi4py*',
    r'\\OPENMP': '*OpenMP*',
    r'\\CUDA': '*CUDA*',
    r'\\VTK': '*VTK*',
    r'\\SILO': '*SILO*',
    r'\\VisIt': '*VisIt*',
    r'\\mayavi': '*Mayavi2*',
    r'\\MATPLOTLIB': '*matplotlib*',
    r'\\SCIPY': '*SciPy*',
    r'\\netCDF': '*netCDF*',
    r'\\HDF': '*HDF5*',
    r'\\LINUX': '*Linux*',
    r'\\WINDOWS': '*MS Windows*',
    r'\\gnuplot': '*gnuplot*',
    r'\\gmshextern': '*Gmsh*',
    r'\\GOCAD': '*GOCAD*',
    r'\\OpenCASCADE': '*OpenCASCADE*',

    # Classes
    r'\\LinearPDE': '``LinearPDE``',
    r'\\NLPDE': '``NonlinearPDE``',
    r'\\Poisson': '``Poisson``',
    r'\\Helmholtz': '``Helmholtz``',
    r'\\Lame': '``Lame``',
    r'\\Data': '``Data``',
    r'\\Domain': '``Domain``',
    r'\\FunctionSpace': '``FunctionSpace``',
    r'\\Operator': '``Operator``',
    r'\\SolverOptions': '``SolverOptions``',
    r'\\Point': '``Point``',
    r'\\PropertySet': '``PropertySet``',
    r'\\Design': '``Design``',
    r'\\TagMap': '``TagMap``',
    r'\\EVALUATOR': '``Evaluator``',
    r'\\SYMBOL': '``Symbol``',

    # Function spaces
    r'\\SolutionFS': 'solution ``FunctionSpace``',
    r'\\ReducedSolutionFS': 'reduced solution ``FunctionSpace``',
    r'\\FunctionOnBoundary': 'boundary ``FunctionSpace``',
    r'\\Function': 'general ``FunctionSpace``',
    r'\\ContinuousFunction': 'continuous ``FunctionSpace``',
    r'\\DiracDeltaFunctions': 'Dirac delta-function ``FunctionSpace``',

    # Data types
    r'\\Scalar': 'scalar ``Data`` object',
    r'\\Scalars': 'scalar ``Data`` objects',
    r'\\Vector': 'vector ``Data`` object',
    r'\\Tensor': 'tensor ``Data`` object',
    r'\\RankOne': 'rank-1 ``Data`` object',
    r'\\RankTwo': 'rank-2 ``Data`` object',
    r'\\RankThree': 'rank-3 ``Data`` object',
    r'\\RankFour': 'rank-4 ``Data`` object',
    r'\\EmptyData': 'empty ``Data``',
    r'\\DataSample': 'data sample',
    r'\\DataSamplePoints': 'data sample points',

    # Solver options
    r'\\PCG': '``SolverOptions.PCG``',
    r'\\BiCGStab': '``SolverOptions.BICGSTAB``',
    r'\\Direct': '``SolverOptions.DIRECT``',
    r'\\GMRES': '``SolverOptions.GMRES``',
    r'\\AMG': '``SolverOptions.AMG``',
    r'\\JACOBI': '``SolverOptions.JACOBI``',
    r'\\ILU': '``SolverOptions.ILU0``',
    r'\\ILUT': '``SolverOptions.ILUT``',
    r'\\RILU': '``SolverOptions.RILU``',
    r'\\GAUSSSEIDEL': '``SolverOptions.GAUSS_SEIDEL``',
    r'\\HRZLUMPING': '``SolverOptions.HRZ_LUMPING``',
    r'\\ROWSUMLUMPING': '``SolverOptions.ROWSUM_LUMPING``',

    # Packages
    r'\\MKL': '``MKL``',
    r'\\UMFPACK': '``UMFPACK``',
    r'\\PASO': '``PASO``',

    # Boolean constants
    r'\\True': '``True``',
    r'\\False': '``False``',

    # Other
    r'\\Shape': 'shape',
    r'\\Rank': 'rank',
    r'\\ExampleDirectory': 'example directory',
    r'\\etal': '*et al.*',
    r'\\xspace': '',
}

# Patterns that need regex replacement
REGEX_REPLACEMENTS = [
    # Cross-reference patterns
    (r'\\Sec\{([^}]+)\}', r'Section :ref:`\1`'),
    (r'\\Chap\{([^}]+)\}', r'Chapter :ref:`\1`'),
    (r'\\App\{([^}]+)\}', r'Appendix :ref:`\1`'),
    (r'\\fig\{([^}]+)\}', r'Figure :numref:`\1`'),
    (r'\\eqn\{([^}]+)\}', r'Equation :eq:`\1`'),
    (r'\\tab\{([^}]+)\}', r'Table :numref:`\1`'),
    (r'\\Refe\{([^}]+)\}', r':cite:`\1`'),

    # Code and module references
    (r'\\module\{([^}]+)\}', r'``\1``'),
    (r'\\class\{([^}]+)\}', r'``\1``'),
    (r'\\member\{([^}]+)\}', r'``\1``'),
    (r'\\method\{([^}]+)\}', r'``\1()``'),
    (r'\\function\{([^}]+)\}', r'``\1()``'),
    (r'\\constant\{([^}]+)\}', r'``\1``'),
    (r'\\var\{([^}]+)\}', r'``\1``'),
    (r'\\code\{([^}]+)\}', r'``\1``'),
    (r'\\file\{([^}]+)\}', r'``\1``'),
    (r'\\program\{([^}]+)\}', r'``\1``'),
    (r'\\env\{([^}]+)\}', r'``\1``'),

    # Warning
    (r'\\warning\{([^}]+)\}', r'\n.. warning::\n\n   \1\n'),

    # Partial derivative
    (r'\\fracp\{([^}]+)\}\{([^}]+)\}', r'∂\1/∂\2'),

    # Index entries - remove them
    (r'\\index\{[^}]+\}', ''),

    # Labels - convert to RST labels
    (r'\\label\{([^}]+)\}', r'\n.. _\1:\n'),

    # URL handling
    (r'\\url\{([^}]+)\}', r'`\1 <\1>`_'),

    # Footnotes - simplify
    (r'\\footnote\{([^}]+)\}', r' [#]_'),
]


def preprocess_latex(content):
    """Pre-process LaTeX content to expand macros before pandoc conversion."""

    # Apply simple macro replacements
    for macro, replacement in MACRO_REPLACEMENTS.items():
        content = re.sub(macro + r'(?![a-zA-Z])', replacement, content)

    # Apply regex replacements
    for pattern, replacement in REGEX_REPLACEMENTS:
        content = re.sub(pattern, replacement, content)

    # Handle \input{} commands by noting them
    content = re.sub(r'\\input\{([^}]+)\}', r'\n\n.. note:: Content from \1.tex would be included here\n\n', content)

    # Convert python environment to code-block
    content = re.sub(r'\\begin\{python\}', r'\\begin{lstlisting}[language=Python]', content)
    content = re.sub(r'\\end\{python\}', r'\\end{lstlisting}', content)

    return content


def convert_with_pandoc(latex_content, output_file):
    """Convert LaTeX to RST using pandoc."""
    try:
        result = subprocess.run(
            ['pandoc', '-f', 'latex', '-t', 'rst', '--wrap=none'],
            input=latex_content,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            print(f"Pandoc warning: {result.stderr}")
        return result.stdout
    except Exception as e:
        print(f"Error running pandoc: {e}")
        return None


def postprocess_rst(content, chapter_name):
    """Post-process RST content to fix common issues."""

    # Fix escaped asterisks from bold/italic
    content = re.sub(r'\\\*\*([^*]+)\*\\\*', r'**\1**', content)
    content = re.sub(r'\\\*([^*]+)\\\*', r'*\1*', content)
    content = re.sub(r'\\_', '_', content)

    # Fix double backticks that pandoc sometimes creates
    content = re.sub(r'````', '``', content)

    # Fix code blocks - ensure proper indentation
    lines = content.split('\n')
    fixed_lines = []
    in_code_block = False

    for line in lines:
        # Detect code block markers
        if '.. code::' in line or '.. code-block::' in line:
            in_code_block = True
            fixed_lines.append(line)
            continue

        if in_code_block:
            if line.strip() == '' and len(fixed_lines) > 0:
                # Empty line might end code block
                fixed_lines.append(line)
            elif line and not line[0].isspace() and line.strip() != '':
                # Non-indented, non-empty line ends code block
                in_code_block = False
                fixed_lines.append(line)
            else:
                fixed_lines.append(line)
        else:
            fixed_lines.append(line)

    content = '\n'.join(fixed_lines)

    # Add chapter header if not present
    if not content.strip().startswith('='):
        header = f"{'=' * len(chapter_name)}\n{chapter_name}\n{'=' * len(chapter_name)}\n\n"
        content = header + content

    # Fix math blocks - ensure they use proper RST math directive
    content = re.sub(r'\.\. math::\s*\n\s*\n', '.. math::\n\n   ', content)

    return content


def convert_file(tex_file, output_dir, chapter_name=None):
    """Convert a single LaTeX file to RST."""
    print(f"Converting {tex_file}...")

    with open(tex_file, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    # Pre-process
    content = preprocess_latex(content)

    # Convert with pandoc
    rst_content = convert_with_pandoc(content, None)

    if rst_content is None:
        print(f"  Failed to convert {tex_file}")
        return None

    # Determine chapter name
    if chapter_name is None:
        chapter_name = Path(tex_file).stem.replace('_', ' ').title()

    # Post-process
    rst_content = postprocess_rst(rst_content, chapter_name)

    # Write output
    output_file = output_dir / (Path(tex_file).stem + '.rst')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(rst_content)

    print(f"  -> {output_file}")
    return output_file


def main():
    """Main conversion routine."""
    # Paths
    script_dir = Path(__file__).parent
    latex_dir = script_dir.parent / 'user'
    output_dir = script_dir / 'user_guide'

    # Ensure output directory exists
    output_dir.mkdir(exist_ok=True)

    # Define chapters and their source files
    chapters = {
        'tutorial_pde': {
            'name': 'Tutorial: Solving PDEs',
            'files': ['TutorialPDE.tex', 'firststep.tex', 'diffusion.tex',
                     'wave.tex', 'heatedblock.tex', 'dirac.tex', 'unstructured.tex']
        },
        'execute': {
            'name': 'Execution of Scripts',
            'files': ['execute.tex']
        },
        'escript': {
            'name': 'The escript Module',
            'files': ['escript.tex']
        },
        'linear_pde': {
            'name': 'The linearPDEs Module',
            'files': ['linearPDE.tex']
        },
        'finley': {
            'name': 'The finley Module',
            'files': ['finley.tex', 'finleyelements.tex']
        },
        'ripley': {
            'name': 'The ripley Module',
            'files': ['ripley.tex']
        },
        'speckley': {
            'name': 'The speckley Module',
            'files': ['speckley.tex']
        },
        'weipa': {
            'name': 'The weipa Module',
            'files': ['weipa.tex']
        },
        'trilinos': {
            'name': 'Using Trilinos',
            'files': ['trilinos.tex']
        },
        'symbolic': {
            'name': 'Symbolic Toolbox',
            'files': ['symbolic.tex']
        },
        'appendix': {
            'name': 'Appendix',
            'files': ['appendix.tex', 'notation.tex', 'nonlinearPDE.tex', 'lumping.tex']
        }
    }

    # Convert each chapter
    for chapter_key, chapter_info in chapters.items():
        print(f"\n=== Converting chapter: {chapter_info['name']} ===")

        # For simplicity, convert the main file only (sub-files would need \input handling)
        main_file = latex_dir / chapter_info['files'][0]
        if main_file.exists():
            convert_file(main_file, output_dir, chapter_info['name'])
        else:
            print(f"  Warning: {main_file} not found")

    print("\n=== Conversion complete ===")
    print("Note: The converted files may need manual review and fixes.")
    print("Common issues to check:")
    print("  - Math equations formatting")
    print("  - Code block indentation")
    print("  - Cross-references")
    print("  - Figure paths")


if __name__ == '__main__':
    main()
