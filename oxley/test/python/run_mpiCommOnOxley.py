##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
##############################################################################

"""
Test suite for MPI communicator retrieval on Oxley domains
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from test_mpiComm import Test_MPI_Comm_Retrieval

from esys.escript import getMPISizeWorld, hasFeature

# Skip if Trilinos not available (required for Oxley)
HAVE_TRILINOS = hasFeature('trilinos')

try:
    from esys.oxley import Rectangle, Brick
except ImportError:
    HAVE_TRILINOS = False

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


@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class Test_MPI_Comm_Oxley2D(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on 2D Oxley Rectangle domains"""

    def createDomain(self, comm=None):
        """Create a 2D Rectangle domain"""
        if comm is None:
            return Rectangle(n0=20, n1=20, d0=NX, d1=NY)
        else:
            # For custom communicators, let oxley auto-detect subdivision
            return Rectangle(n0=20, n1=20, comm=comm)


@unittest.skipIf(not HAVE_TRILINOS, "Trilinos not available")
class Test_MPI_Comm_Oxley3D(Test_MPI_Comm_Retrieval):
    """Test MPI communicator retrieval on 3D Oxley Brick domains"""

    def createDomain(self, comm=None):
        """Create a 3D Brick domain"""
        if comm is None:
            return Brick(n0=12, n1=12, n2=12, d0=NXb, d1=NYb, d2=NZb)
        else:
            # For custom communicators, let oxley auto-detect subdivision
            return Brick(n0=12, n1=12, n2=12, comm=comm)


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
