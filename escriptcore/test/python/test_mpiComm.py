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
Test suite for MPI communicator retrieval (getMPIComm) functionality
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *

# Check if mpi4py is available
try:
    from mpi4py import MPI
    HAVE_MPI4PY = True
except ImportError:
    HAVE_MPI4PY = False
    MPI = None


class Test_MPI_Comm_Retrieval(unittest.TestCase):
    """
    Base test class for MPI communicator retrieval.

    Subclasses must implement createDomain() method that returns a domain.
    Tests verify that getMPIComm() correctly returns the communicator for:
    a) No comm argument (default MPI_COMM_WORLD)
    b) Explicit comm=MPI.COMM_WORLD
    c) Custom split communicator (if multiple processes available)
    """

    def createDomain(self, comm=None):
        """
        Create a domain for testing. Must be implemented by subclasses.

        Args:
            comm: Optional MPI communicator to pass to domain constructor

        Returns:
            Domain object
        """
        raise NotImplementedError("Subclasses must implement createDomain()")

    def setUp(self):
        """Set up test fixtures"""
        self.mpiSize = getMPISizeWorld()
        self.mpiRank = getMPIRankWorld()

    def tearDown(self):
        """Clean up after tests"""
        pass

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    def test_getMPIComm_returns_mpi4py_comm(self):
        """Test that getMPIComm() returns an mpi4py.MPI.Comm object"""
        domain = self.createDomain()
        comm = domain.getMPIComm()

        # Should return mpi4py.MPI.Comm object or None
        self.assertTrue(comm is None or isinstance(comm, MPI.Comm),
                        f"Expected MPI.Comm or None, got {type(comm)}")

        # If MPI is enabled, should not be None
        if self.mpiSize > 0:
            self.assertIsNotNone(comm, "getMPIComm() should not return None when MPI is enabled")

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    def test_getMPIComm_no_arg_uses_comm_world(self):
        """Test case (a): No comm argument defaults to MPI_COMM_WORLD"""
        # Create domain without comm argument
        domain = self.createDomain()
        comm = domain.getMPIComm()

        if comm is not None:
            # Should be IDENT or CONGRUENT to MPI_COMM_WORLD
            result = MPI.Comm.Compare(MPI.COMM_WORLD, comm)
            self.assertIn(result, [MPI.IDENT, MPI.CONGRUENT],
                         f"Expected IDENT or CONGRUENT to MPI_COMM_WORLD, got {result}")

            # Size and rank should match
            self.assertEqual(comm.Get_size(), self.mpiSize,
                           "Communicator size mismatch")
            self.assertEqual(comm.Get_rank(), self.mpiRank,
                           "Communicator rank mismatch")

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    def test_getMPIComm_explicit_comm_world(self):
        """Test case (b): Explicit comm=MPI.COMM_WORLD"""
        # Create domain with explicit MPI_COMM_WORLD
        domain = self.createDomain(comm=MPI.COMM_WORLD)
        comm = domain.getMPIComm()

        if comm is not None:
            # Should be IDENT or CONGRUENT to MPI_COMM_WORLD
            result = MPI.Comm.Compare(MPI.COMM_WORLD, comm)
            self.assertIn(result, [MPI.IDENT, MPI.CONGRUENT],
                         f"Expected IDENT or CONGRUENT to MPI_COMM_WORLD, got {result}")

            # Size and rank should match
            self.assertEqual(comm.Get_size(), self.mpiSize,
                           "Communicator size mismatch")
            self.assertEqual(comm.Get_rank(), self.mpiRank,
                           "Communicator rank mismatch")

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    @unittest.skipIf(getMPISizeWorld() < 2, "Requires at least 2 MPI processes")
    def test_getMPIComm_split_comm_size2(self):
        """Test case (c): Custom split communicator with size 2"""
        # Split into groups of size 2 (or as close as possible)
        color = self.mpiRank // 2
        split_comm = MPI.COMM_WORLD.Split(color, self.mpiRank)

        try:
            # Create domain with split communicator
            domain = self.createDomain(comm=split_comm)
            retrieved_comm = domain.getMPIComm()

            if retrieved_comm is not None:
                # Should be IDENT or CONGRUENT to split_comm
                result = MPI.Comm.Compare(split_comm, retrieved_comm)
                self.assertIn(result, [MPI.IDENT, MPI.CONGRUENT],
                             f"Expected IDENT or CONGRUENT to split communicator, got {result}")

                # Size should match (2 or 1 for odd number of processes)
                expected_size = min(2, self.mpiSize - color * 2)
                self.assertEqual(retrieved_comm.Get_size(), expected_size,
                               f"Expected size {expected_size}, got {retrieved_comm.Get_size()}")

                # Rank should be valid within the split communicator
                self.assertLess(retrieved_comm.Get_rank(), retrieved_comm.Get_size(),
                              "Rank should be less than size")
        finally:
            split_comm.Free()

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    @unittest.skipIf(getMPISizeWorld() < 4, "Requires at least 4 MPI processes")
    def test_getMPIComm_split_comm_size4(self):
        """Test case (c): Custom split communicator with size 4"""
        # Split into groups of size 4 (or as close as possible)
        color = self.mpiRank // 4
        split_comm = MPI.COMM_WORLD.Split(color, self.mpiRank)

        try:
            # Create domain with split communicator
            domain = self.createDomain(comm=split_comm)
            retrieved_comm = domain.getMPIComm()

            if retrieved_comm is not None:
                # Should be IDENT or CONGRUENT to split_comm
                result = MPI.Comm.Compare(split_comm, retrieved_comm)
                self.assertIn(result, [MPI.IDENT, MPI.CONGRUENT],
                             f"Expected IDENT or CONGRUENT to split communicator, got {result}")

                # Size should match (4 or less for non-divisible process counts)
                expected_size = min(4, self.mpiSize - color * 4)
                self.assertEqual(retrieved_comm.Get_size(), expected_size,
                               f"Expected size {expected_size}, got {retrieved_comm.Get_size()}")
        finally:
            split_comm.Free()

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    def test_functionspace_getMPIComm(self):
        """Test that FunctionSpace.getMPIComm() returns domain's communicator"""
        domain = self.createDomain()
        domain_comm = domain.getMPIComm()

        # Create function space
        fs = ContinuousFunction(domain)
        fs_comm = fs.getMPIComm()

        if domain_comm is not None and fs_comm is not None:
            # Should be same communicator
            result = MPI.Comm.Compare(domain_comm, fs_comm)
            self.assertIn(result, [MPI.IDENT, MPI.CONGRUENT],
                         "FunctionSpace communicator should match domain communicator")

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    def test_data_getMPIComm(self):
        """Test that Data.getMPIComm() returns domain's communicator"""
        domain = self.createDomain()
        domain_comm = domain.getMPIComm()

        # Create function space and data
        fs = ContinuousFunction(domain)
        data = Scalar(1.0, fs)
        data_comm = data.getMPIComm()

        if domain_comm is not None and data_comm is not None:
            # Should be same communicator
            result = MPI.Comm.Compare(domain_comm, data_comm)
            self.assertIn(result, [MPI.IDENT, MPI.CONGRUENT],
                         "Data communicator should match domain communicator")

    @unittest.skipIf(not HAVE_MPI4PY, "mpi4py not available")
    @unittest.skipIf(getMPISizeWorld() < 2, "Requires at least 2 MPI processes")
    def test_split_comm_consistency_across_objects(self):
        """Test that split communicator is consistent across Domain, FunctionSpace, and Data"""
        # Create split communicator
        color = self.mpiRank % 2
        split_comm = MPI.COMM_WORLD.Split(color, self.mpiRank)

        try:
            # Create domain, function space, and data
            domain = self.createDomain(comm=split_comm)
            fs = ContinuousFunction(domain)
            data = Scalar(1.0, fs)

            # Get communicators
            domain_comm = domain.getMPIComm()
            fs_comm = fs.getMPIComm()
            data_comm = data.getMPIComm()

            if domain_comm is not None and fs_comm is not None and data_comm is not None:
                # All should be the same communicator
                result1 = MPI.Comm.Compare(domain_comm, fs_comm)
                result2 = MPI.Comm.Compare(fs_comm, data_comm)
                result3 = MPI.Comm.Compare(domain_comm, data_comm)

                self.assertIn(result1, [MPI.IDENT, MPI.CONGRUENT],
                             "Domain and FunctionSpace communicators should match")
                self.assertIn(result2, [MPI.IDENT, MPI.CONGRUENT],
                             "FunctionSpace and Data communicators should match")
                self.assertIn(result3, [MPI.IDENT, MPI.CONGRUENT],
                             "Domain and Data communicators should match")

                # All should match the original split communicator
                result_orig = MPI.Comm.Compare(split_comm, domain_comm)
                self.assertIn(result_orig, [MPI.IDENT, MPI.CONGRUENT],
                             "Retrieved communicator should match original split communicator")
        finally:
            split_comm.Free()


# Helper function for running tests
def get_mpi_comm_tests(domain_factory):
    """
    Create a test class for a specific domain factory function.

    Args:
        domain_factory: Function that creates a domain, signature: domain_factory(comm=None)

    Returns:
        Test class that can be used with unittest
    """
    class DomainSpecificTests(Test_MPI_Comm_Retrieval):
        def createDomain(self, comm=None):
            return domain_factory(comm=comm)

    return DomainSpecificTests
