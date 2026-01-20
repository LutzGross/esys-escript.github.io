##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
##############################################################################

"""
GOCAD Voxet Reader Example

This script demonstrates how to read GOCAD voxet data and use it with
Ripley domains in escript.

GOCAD voxet files store 3D gridded data, typically used in geosciences
for representing subsurface properties like density, velocity, or resistivity.

Note: This is a placeholder example showing the general structure.
The actual voxet reading functionality may require additional libraries
or custom parsing of GOCAD file formats.

Usage:
    run-escript voxet_reader.py <voxet_file>
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

import sys
import numpy as np
from esys.escript import *
from esys.ripley import Brick


def read_gocad_voxet(filename):
    """
    Read a GOCAD voxet file.

    This is a placeholder function. A full implementation would parse
    the GOCAD voxet format which typically consists of:
    - Header file (.vo) with grid geometry and properties
    - Binary data file with the actual property values

    Parameters
    ----------
    filename : str
        Path to the voxet header file

    Returns
    -------
    dict
        Dictionary containing:
        - 'nx', 'ny', 'nz': Number of cells in each direction
        - 'origin': Origin coordinates (x0, y0, z0)
        - 'spacing': Cell spacing (dx, dy, dz)
        - 'data': 3D numpy array of property values
    """
    print(f"Reading GOCAD voxet file: {filename}")
    print("Note: This is a placeholder implementation")

    # Placeholder: In reality, parse the voxet header file
    # and read the binary data file
    voxet_data = {
        'nx': 50,
        'ny': 50,
        'nz': 20,
        'origin': (0.0, 0.0, 0.0),
        'spacing': (100.0, 100.0, 50.0),
        'data': np.random.rand(50, 50, 20)  # Placeholder data
    }

    return voxet_data


def create_ripley_domain_from_voxet(voxet_data):
    """
    Create a Ripley Brick domain matching the voxet geometry.

    Parameters
    ----------
    voxet_data : dict
        Voxet data as returned by read_gocad_voxet()

    Returns
    -------
    Domain
        Ripley Brick domain
    """
    nx = voxet_data['nx']
    ny = voxet_data['ny']
    nz = voxet_data['nz']

    x0, y0, z0 = voxet_data['origin']
    dx, dy, dz = voxet_data['spacing']

    # Calculate domain extents
    l0 = nx * dx
    l1 = ny * dy
    l2 = nz * dz

    # Create Ripley Brick domain
    domain = Brick(n0=nx, n1=ny, n2=nz, l0=l0, l1=l1, l2=l2,
                   d0=1, d1=1, d2=1)  # No MPI subdivision for simplicity

    print(f"Created Ripley domain: {nx}x{ny}x{nz} cells")
    print(f"Domain extent: {l0}x{l1}x{l2}")

    return domain


def map_voxet_data_to_domain(voxet_data, domain, function_space=None):
    """
    Map voxet data values onto a Ripley domain.

    Parameters
    ----------
    voxet_data : dict
        Voxet data as returned by read_gocad_voxet()
    domain : Domain
        Ripley domain
    function_space : FunctionSpace, optional
        Target function space (default: ReducedFunction)

    Returns
    -------
    Data
        escript Data object containing the voxet values
    """
    if function_space is None:
        function_space = ReducedFunction(domain)

    # Flatten the 3D voxet data to 1D array
    # Note: Ordering may need adjustment depending on voxet format
    values_flat = voxet_data['data'].flatten()

    # Create escript Data object
    # This assumes the voxet grid matches the domain grid
    data = Scalar(0.0, function_space)

    # In a real implementation, you would map the voxet values
    # to the appropriate locations in the escript Data object
    # This may require interpolation if grids don't match exactly

    print(f"Mapped voxet data to domain")
    print(f"Data range: [{np.min(values_flat):.3f}, {np.max(values_flat):.3f}]")

    return data


def main():
    """Main function to demonstrate voxet reading workflow."""

    if len(sys.argv) < 2:
        print("Usage: run-escript voxet_reader.py <voxet_file>")
        print("\nNote: This is a placeholder example.")
        print("Using synthetic data for demonstration.")
        voxet_file = "example.vo"
    else:
        voxet_file = sys.argv[1]

    # Read voxet data
    voxet_data = read_gocad_voxet(voxet_file)

    # Create matching Ripley domain
    domain = create_ripley_domain_from_voxet(voxet_data)

    # Map voxet data onto domain
    property_data = map_voxet_data_to_domain(voxet_data, domain)

    # Example: Use the data in a PDE
    # (This would be application-specific)
    print("\nVoxet data successfully loaded into Ripley domain")
    print("The property_data can now be used in PDEs or exported")

    # Example export
    # from esys.weipa import saveVTK
    # saveVTK("voxet_output.vtu", property=property_data)


if __name__ == "__main__":
    main()
