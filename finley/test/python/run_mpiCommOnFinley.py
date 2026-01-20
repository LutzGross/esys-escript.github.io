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
Test suite for MPI communicator retrieval on Finley domains
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_mpiComm import Test_MPI_Comm_Retrieval

from esys.escript import getMPISizeWorld
from esys.finley import Rectangle, Brick, ReadMesh, ReadGmsh, LoadMesh
import os

# Test data paths
try:
    FINLEY_TEST_DATA = os.environ['FINLEY_TEST_DATA']
except KeyError:
    # Default to finley/test/python directory
    FINLEY_TEST_DATA = os.path.dirname(__file__)

FINLEY_TEST_MESH_PATH = os.path.join(FINLEY_TEST_DATA, "data_meshes")

try:
    FINLEY_WORKDIR = os.environ['FINLEY_WORKDIR']
except KeyError:
    FINLEY_WORKDIR = '.'

# Determine optimal domain subdivisions for current MPI configuration
mpiSize = getMPISizeWorld()

# For 2D domains
for x in [int(mpiSize**0.5), 2, 3, 5, 7, 1]:
    NX = x
    NY = mpiSize // x
    if NX * NY == mpiSize:
        break

# For 3D domains
for x in [(int(mpiSize**(1/3.)), int(mpiSize**(1/3.))), (2, 3), (2, 2), (1, 2), (1, 1)]:
    NXb = x[0]
    NYb = x[1]
    NZb = mpiSize // (x[0] * x[1])
    if NXb * NYb * NZb == mpiSize:
        break


class Test_MPI_Comm_Finley2D(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on 2D Finley Rectangle domains"""

    def createDomain(self, comm=None):
        """Create a 2D Rectangle domain"""
        if comm is None:
            return Rectangle(n0=20, n1=20)
        else:
            return Rectangle(n0=20, n1=20, comm=comm)


class Test_MPI_Comm_Finley3D(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on 3D Finley Brick domains"""

    def createDomain(self, comm=None):
        """Create a 3D Brick domain"""
        if comm is None:
            return Brick(n0=12, n1=12, n2=12)
        else:
            return Brick(n0=12, n1=12, n2=12, comm=comm)


class Test_MPI_Comm_ReadMesh(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on domains loaded via ReadMesh"""

    def createDomain(self, comm=None):
        """Create a domain by reading a mesh file"""
        mesh_path = os.path.join(FINLEY_TEST_MESH_PATH, 'rect_4x4.fly')
        if comm is None:
            return ReadMesh(mesh_path)
        else:
            return ReadMesh(mesh_path, comm=comm)


class Test_MPI_Comm_ReadGmsh(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on domains loaded via ReadGmsh"""

    def createDomain(self, comm=None):
        """Create a domain by reading a Gmsh file"""
        mesh_path = os.path.join(FINLEY_TEST_MESH_PATH, 'tagtest.msh')
        if comm is None:
            return ReadGmsh(mesh_path, 2)
        else:
            return ReadGmsh(mesh_path, 2, comm=comm)


class Test_MPI_Comm_LoadMesh(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on domains loaded via LoadMesh"""

    def createDomain(self, comm=None):
        """Create a domain by loading an HDF5 dump file"""
        from esys.escript import getMPIRankWorld, getMPISizeWorld

        # Use unique filename per MPI world rank to avoid file conflicts
        # when split communicators are used
        world_rank = getMPIRankWorld()
        dump_file = os.path.join(FINLEY_WORKDIR, f'test_load_comm_{world_rank}.dump.h5')

        # Step 1: Create a Rectangle domain with the provided communicator
        if comm is None:
            test_domain = Rectangle(n0=8, n1=8)
        else:
            test_domain = Rectangle(n0=8, n1=8, comm=comm)

        # Step 2: Dump the domain to file
        test_domain.dump(dump_file)

        # Step 3: Verify the dump file exists for this rank
        # Determine which rank we are in the communicator
        if comm is not None:
            my_rank = comm.Get_rank()
            comm_size = comm.Get_size()
        else:
            my_rank = getMPIRankWorld()
            comm_size = getMPISizeWorld()

        # Check if the expected file exists
        # dump() appends rank to filename if comm size > 1: filename.0000, filename.0001, etc.
        if comm_size > 1:
            my_dump_file = f"{dump_file}.{my_rank:04d}"
        else:
            my_dump_file = dump_file

        if not os.path.exists(my_dump_file):
            raise RuntimeError(f"Expected dump file {my_dump_file} not found after dump()")

        # Step 4: Load the domain and return it
        if comm is None:
            return LoadMesh(dump_file)
        else:
            return LoadMesh(dump_file, comm=comm)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
