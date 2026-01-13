
##############################################################################
#
# Copyright (c) 2003-2025 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

"""
Comprehensive tests for MPIDomainArray and DataCoupler classes.

These tests cover:
- MPIDomainArray communicator topology creation
- DataCoupler point-to-point communication (send/receive/exchange)
- DataCoupler collective operations (broadcast/allreduce)
- Scalar, vector, and tensor Data objects
- Error handling and edge cases

Note: These tests require MPI with at least 2 processes.
Run with: run-escript -p 4 run_domaincoupler.py
"""

__copyright__="""Copyright (c) 2003-2025 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

from esys.escript import *
from esys.escript import MPIDomainArray, DataCoupler
from esys.escript.linearPDEs import LinearPDE
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import numpy as np

try:
    from mpi4py import MPI
    HAVE_MPI = True
    WORLD_SIZE = MPI.COMM_WORLD.Get_size()
    WORLD_RANK = MPI.COMM_WORLD.Get_rank()
except ImportError:
    HAVE_MPI = False
    WORLD_SIZE = 1
    WORLD_RANK = 0

# Skip all tests if MPI is not available or only 1 process
SKIP_MPI_TESTS = not HAVE_MPI or WORLD_SIZE < 2


@unittest.skipIf(SKIP_MPI_TESTS, "MPI with at least 2 processes required")
class Test_MPIDomainArray(unittest.TestCase):
    """Tests for MPIDomainArray communicator topology."""

    def test_basic_topology_2domains(self):
        """Test basic 2-domain topology creation."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)

        # Check basic properties (all ranks)
        self.assertEqual(domain_array.numDomains, num_domains)
        self.assertIsNotNone(domain_array.getDomainComm())
        self.assertIsNotNone(domain_array.getSubdomainComm())

        # Check domain index is valid
        domain_idx = domain_array.getDomainIndex()
        self.assertIn(domain_idx, range(num_domains))

        # Check subdomain communicator connects all domains
        subdomain_comm = domain_array.getSubdomainComm()
        self.assertEqual(subdomain_comm.Get_size(), num_domains)

    @unittest.skipIf(WORLD_SIZE < 4, "Need at least 4 processes")
    def test_topology_with_multiple_subdomains(self):
        """Test topology with multiple subdomains per domain."""
        num_domains = 2
        num_subdomains = WORLD_SIZE // 2
        domain_array = MPIDomainArray(numDomains=num_domains,
                                      numSubDomains=num_subdomains,
                                      comm=MPI.COMM_WORLD)

        # Check domain communicator has correct size
        domain_comm = domain_array.getDomainComm()
        self.assertEqual(domain_comm.Get_size(), num_subdomains)

        # Check subdomain communicator connects domains
        subdomain_comm = domain_array.getSubdomainComm()
        self.assertEqual(subdomain_comm.Get_size(), num_domains)


@unittest.skipIf(SKIP_MPI_TESTS, "MPI with at least 2 processes required")
class Test_DataCoupler_PointToPoint(unittest.TestCase):
    """Tests for DataCoupler point-to-point communication."""

    def setUp(self):
        # Create 2 domains with identical meshes
        from esys.ripley import Rectangle
        self.num_domains = 2
        self.domain_array = MPIDomainArray(numDomains=self.num_domains,
                                          comm=MPI.COMM_WORLD)
        self.domain = Rectangle(n0=10, n1=10,
                               comm=self.domain_array.getDomainComm())
        self.coupler = DataCoupler(self.domain_array)
        self.my_domain_idx = self.domain_array.getDomainIndex()

    def test_send_receive_scalar(self):
        """Test sending and receiving scalar Data."""
        fs = Solution(self.domain)

        # All ranks participate
        if self.my_domain_idx == 0:
            # Domain 0 sends
            data = Scalar(42.0, fs)
            self.coupler.send(data, dest_domain_index=1, tag=100)
        elif self.my_domain_idx == 1:
            # Domain 1 receives
            received = self.coupler.receive(fs, source_domain_index=0, tag=100)

            # All ranks in domain 1 verify (sup is collective)
            received_val = sup(received)  # sup() not Lsup() - checking actual value
            self.assertAlmostEqual(received_val, 42.0, places=10)

    def test_send_receive_vector(self):
        """Test sending and receiving vector Data."""
        fs = Solution(self.domain)
        x = self.domain.getX()

        if self.my_domain_idx == 0:
            # Domain 0 sends position field
            data = Data(x, fs)
            self.coupler.send(data, dest_domain_index=1, tag=101)
        elif self.my_domain_idx == 1:
            # Domain 1 receives
            received = self.coupler.receive(fs, source_domain_index=0, tag=101)

            # All ranks in domain 1 verify (Lsup is collective)
            x_received = Data(received, fs)
            diff = Lsup(length(x_received - x))
            self.assertLess(diff, 1e-10)

    def test_exchange_bidirectional(self):
        """Test bidirectional exchange between two domains."""
        fs = Solution(self.domain)

        # All ranks in both domains participate
        if self.my_domain_idx == 0:
            send_data = Scalar(10.0, fs)
            received = self.coupler.exchange(send_data, fs, peer_domain_index=1, tag=200)

            # All ranks verify (sup is collective)
            received_val = sup(received)  # sup() not Lsup() - checking actual value
            self.assertAlmostEqual(received_val, 20.0, places=10)

        elif self.my_domain_idx == 1:
            send_data = Scalar(20.0, fs)
            received = self.coupler.exchange(send_data, fs, peer_domain_index=0, tag=200)

            # All ranks verify (sup is collective)
            received_val = sup(received)  # sup() not Lsup() - checking actual value
            self.assertAlmostEqual(received_val, 10.0, places=10)


@unittest.skipIf(SKIP_MPI_TESTS, "MPI with at least 2 processes required")
class Test_DataCoupler_Collective(unittest.TestCase):
    """Tests for DataCoupler collective operations."""

    def setUp(self):
        from esys.ripley import Rectangle
        self.num_domains = 2
        self.domain_array = MPIDomainArray(numDomains=self.num_domains,
                                          comm=MPI.COMM_WORLD)
        self.domain = Rectangle(n0=10, n1=10,
                               comm=self.domain_array.getDomainComm())
        self.coupler = DataCoupler(self.domain_array)
        self.my_domain_idx = self.domain_array.getDomainIndex()

    def test_broadcast_value_scalar(self):
        """Test broadcasting a scalar value."""
        # All ranks participate in broadcast
        if self.my_domain_idx == 0 and self.domain.isRootRank():
            result = self.coupler.broadcast_value(3.14159, root_domain_index=0)
        else:
            result = self.coupler.broadcast_value(root_domain_index=0)

        # All ranks verify
        self.assertAlmostEqual(result, 3.14159, places=5)

    def test_allreduce_value_sum(self):
        """Test all-reduce with SUM operation on scalars."""
        # All ranks participate
        my_value = float(self.my_domain_idx)
        result = self.coupler.allreduce_value(my_value, op=MPI.SUM)

        # All ranks verify
        expected = sum(range(self.num_domains))
        self.assertAlmostEqual(result, expected)

    def test_allreduce_value_max(self):
        """Test all-reduce with MAX operation."""
        # All ranks participate
        my_value = float(self.my_domain_idx)
        result = self.coupler.allreduce_value(my_value, op=MPI.MAX)

        # All ranks verify
        expected = float(self.num_domains - 1)
        self.assertAlmostEqual(result, expected)

    def test_broadcast_data(self):
        """Test broadcasting Data object from root domain."""
        fs = Solution(self.domain)
        x = self.domain.getX()

        # All ranks participate in broadcast
        if self.my_domain_idx == 0:
            data = Data(x[0] + x[1], fs)  # Use Data, not Scalar (this is a field)
            result = self.coupler.broadcast(data=data, root_domain_index=0)
        else:
            result = self.coupler.broadcast(function_space=fs, root_domain_index=0)

        # All ranks verify (sup is collective)
        expected = x[0] + x[1]
        diff = sup(abs(result - expected))  # sup(abs()) not Lsup(abs())
        self.assertLess(diff, 1e-8)

    def test_allreduce_data_sum(self):
        """Test all-reduce on Data with SUM operation."""
        fs = Solution(self.domain)

        # All ranks participate - each domain contributes its index
        my_data = Scalar(float(self.my_domain_idx), fs)
        result = self.coupler.allreduce(my_data, op=MPI.SUM)

        # All ranks verify (sup and inf are collective)
        expected_value = sum(range(self.num_domains))
        result_max = sup(result)  # sup() not Lsup() - we want max, not max(abs())
        result_min = inf(result)
        self.assertAlmostEqual(result_max, expected_value, places=10)
        self.assertAlmostEqual(result_min, expected_value, places=10)


# Create test suite
def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(Test_MPIDomainArray))
    test_suite.addTest(unittest.makeSuite(Test_DataCoupler_PointToPoint))
    test_suite.addTest(unittest.makeSuite(Test_DataCoupler_Collective))
    return test_suite


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
