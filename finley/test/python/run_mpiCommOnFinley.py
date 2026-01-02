##############################################################################
#
# Copyright (c) 2003-2025 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
##############################################################################

"""
Test suite for MPI communicator retrieval on Finley domains
"""

__copyright__ = """Copyright (c) 2003-2025 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_mpiComm import Test_MPI_Comm_Retrieval

from esys.escript import getMPISizeWorld
from esys.finley import Rectangle, Brick, ReadMesh, ReadGmsh, LoadMesh
import os

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
            return Rectangle(n0=20, n1=20, d0=NX, d1=NY)
        else:
            # For custom communicators, let finley auto-detect subdivision
            return Rectangle(n0=20, n1=20, comm=comm)


class Test_MPI_Comm_Finley3D(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on 3D Finley Brick domains"""

    def createDomain(self, comm=None):
        """Create a 3D Brick domain"""
        if comm is None:
            return Brick(n0=12, n1=12, n2=12, d0=NXb, d1=NYb, d2=NZb)
        else:
            # For custom communicators, let finley auto-detect subdivision
            return Brick(n0=12, n1=12, n2=12, comm=comm)


class Test_MPI_Comm_ReadMesh(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on domains loaded via ReadMesh"""

    def createDomain(self, comm=None):
        """Create a domain by reading a mesh file"""
        # Use a test mesh file from the finley test data
        mesh_path = os.path.join(os.environ.get('FINLEY_TEST_MESH_PATH', 'test/python/data'),
                                 'rectangle_8x10.fly')
        if comm is None:
            return ReadMesh(mesh_path)
        else:
            return ReadMesh(mesh_path, comm=comm)


class Test_MPI_Comm_ReadGmsh(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on domains loaded via ReadGmsh"""

    def createDomain(self, comm=None):
        """Create a domain by reading a Gmsh file"""
        # Use a test Gmsh file from the finley test data
        mesh_path = os.path.join(os.environ.get('FINLEY_TEST_MESH_PATH', 'test/python/data'),
                                 'tagtest.msh')
        if comm is None:
            return ReadGmsh(mesh_path, 2)
        else:
            return ReadGmsh(mesh_path, 2, comm=comm)


class Test_MPI_Comm_LoadMesh(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on domains loaded via LoadMesh"""

    def createDomain(self, comm=None):
        """Create a domain by loading an HDF5 dump file"""
        # First create and dump a test domain
        import tempfile
        temp_dir = os.environ.get('FINLEY_WORKDIR', tempfile.gettempdir())
        dump_file = os.path.join(temp_dir, 'test_load_comm.dump.h5')

        # Create a small domain and dump it (only do this once per rank)
        if not os.path.exists(dump_file) or getMPIRankWorld() == 0:
            test_domain = Rectangle(n0=8, n1=8)
            test_domain.dump(dump_file)

        # Wait for all ranks to ensure file is written
        import esys.escript
        if hasattr(esys.escript, 'barrier'):
            esys.escript.barrier()

        # Load the domain with the specified communicator
        if comm is None:
            return LoadMesh(dump_file)
        else:
            return LoadMesh(dump_file, comm=comm)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
