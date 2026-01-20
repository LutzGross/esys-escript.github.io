
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""
Comprehensive test suite for MPIDomainArray and DataCoupler classes.

This module provides extensive testing for MPI-based domain decomposition and
inter-domain communication in escript. The test suite validates both the
communicator topology infrastructure (MPIDomainArray) and the data coupling
operations (DataCoupler) across multiple domain types and function spaces.

Test Coverage
-------------
The test suite includes:

**MPIDomainArray Tests:**
    - Basic 2-domain topology creation
    - Multiple subdomains per domain
    - Communicator hierarchy validation

**DataCoupler Point-to-Point Communication:**
    - send() and receive() operations for scalar and vector Data
    - Bidirectional exchange() operations
    - Communication across all supported function spaces

**DataCoupler Collective Operations on Data Objects:**
    - broadcast() operations from root domain to all domains
    - allreduce() operations with SUM operator

**DataCoupler Collective Operations on Scalar Values:**
    - broadcast_value() for Python scalar distribution
    - allreduce_value() with SUM and MAX operators

**Domain Types Tested:**
    - esys.ripley (Rectangle and Brick)
    - esys.finley (Rectangle and Brick)
    - esys.speckley (Rectangle and Brick with order=2)
    - Note: esys.oxley excluded due to deadlock issues (see GitHub issue #123)

**Function Spaces Tested:**
    - Solution
    - ContinuousFunction
    - ReducedSolution
    - Function
    - ReducedFunction
    - FunctionOnBoundary
    - ReducedFunctionOnBoundary
    - DiracDeltaFunctions

**Data Types Tested:**
    - Scalar Data objects
    - Vector Data objects (2D and 3D)
    - Python scalar values

Requirements
------------
These tests require MPI with at least 2 processes. Some tests require exactly
4 processes for multiple subdomain validation.

Usage
-----
Run the complete test suite with::

    run-escript -p 4 run_domaincoupler.py

Or run individual test classes with::

    run-escript -p 4 -m unittest run_domaincoupler.Test_MPIDomainArray

Notes
-----
- Speckley domains do not support ReducedSolution, FunctionOnBoundary, or
  ReducedFunctionOnBoundary. These combinations are automatically skipped.
- All tests use MPI collective operations, requiring all ranks to participate.
- Tests validate results on all participating ranks using collective reduction
  operations (sup, inf).
- Dirac points are automatically added to all domains to enable testing of
  DiracDeltaFunctions function space.

See Also
--------
MPIDomainArray : MPI communicator topology for domain decomposition
DataCoupler : Inter-domain communication interface
"""

__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
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


# Test configurations: (domain_module, domain_class, dimensions, resolution)
DOMAIN_CONFIGS = [
    ('esys.ripley', 'Rectangle', 2, {'n0': 10, 'n1': 10}),
    ('esys.ripley', 'Brick', 3, {'n0': 8, 'n1': 8, 'n2': 8}),
    ('esys.finley', 'Rectangle', 2, {'n0': 10, 'n1': 10}),
    ('esys.finley', 'Brick', 3, {'n0': 6, 'n1': 6, 'n2': 6}),
    ('esys.speckley', 'Rectangle', 2, {'order': 2, 'n0': 10, 'n1': 10}),
    ('esys.speckley', 'Brick', 3, {'order': 2, 'n0': 6, 'n1': 6, 'n2': 6}),
]

# Oxley commented out - may cause deadlocks with DataCoupler operations
# try:
#     import esys.oxley
#     DOMAIN_CONFIGS.extend([
#         ('esys.oxley', 'Rectangle', 2, {'n0': 10, 'n1': 10}),
#         ('esys.oxley', 'Brick', 3, {'n0': 6, 'n1': 6, 'n2': 6}),
#     ])
# except ImportError:
#     pass


def create_domain(domain_module_name, domain_class_name, params, comm):
    """Helper to create a domain instance with Dirac points."""
    import importlib
    module = importlib.import_module(domain_module_name)
    domain_class = getattr(module, domain_class_name)

    # Add Dirac points for testing DiracDeltaFunctions
    dirac_params = params.copy()
    if 'n2' in params:
        # 3D domain - add two Dirac points
        dirac_params['diracPoints'] = [(0.5, 0.5, 0.5), (0.8, 0.8, 0.8)]
        dirac_params['diracTags'] = ['point1', 'point2']
    elif 'n1' in params:
        # 2D domain - add two Dirac points
        dirac_params['diracPoints'] = [(0.5, 0.5), (0.8, 0.8)]
        dirac_params['diracTags'] = ['point1', 'point2']

    return domain_class(**dirac_params, comm=comm)


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
    """Tests for DataCoupler point-to-point communication across all domains and function spaces."""

    def test_send_receive_scalar(self):
        """Test sending and receiving scalar Data across all domains and function spaces."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # Test each domain configuration
        for domain_module, domain_class, dim, params in DOMAIN_CONFIGS:
            with self.subTest(domain=f"{domain_module}.{domain_class}"):
                domain = create_domain(domain_module, domain_class, params,
                                      domain_array.getDomainComm())

                # Test each function space
                for fs_name, fs_constructor in [
                    ('Solution', Solution),
                    ('ContinuousFunction', ContinuousFunction),
                    ('ReducedSolution', ReducedSolution),
                    ('Function', Function),
                    ('ReducedFunction', ReducedFunction),
                    ('FunctionOnBoundary', FunctionOnBoundary),
                    ('ReducedFunctionOnBoundary', ReducedFunctionOnBoundary),
                    ('DiracDeltaFunctions', DiracDeltaFunctions),
                ]:
                    # Skip unsupported speckley function space combinations
                    if domain_module == 'esys.speckley' and fs_name in [
                        'ReducedSolution', 'FunctionOnBoundary', 'ReducedFunctionOnBoundary'
                    ]:
                        continue

                    with self.subTest(function_space=fs_name):
                        fs = fs_constructor(domain)

                        # All ranks participate
                        if my_domain_idx == 0:
                            # Domain 0 sends
                            data = Scalar(42.0, fs)
                            coupler.send(data, dest_domain_index=1, tag=100)
                        elif my_domain_idx == 1:
                            # Domain 1 receives
                            received = coupler.receive(fs, source_domain_index=0, tag=100)

                            # All ranks in domain 1 verify (sup is collective)
                            received_val = sup(received)
                            self.assertAlmostEqual(received_val, 42.0, places=10,
                                                 msg=f"{domain_module}.{domain_class} {fs_name}")

    def test_send_receive_vector(self):
        """Test sending and receiving vector Data across all domains and function spaces."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # Test each domain configuration
        for domain_module, domain_class, dim, params in DOMAIN_CONFIGS:
            with self.subTest(domain=f"{domain_module}.{domain_class}"):
                domain = create_domain(domain_module, domain_class, params,
                                      domain_array.getDomainComm())

                # Test each function space
                for fs_name, fs_constructor in [
                    ('Solution', Solution),
                    ('ContinuousFunction', ContinuousFunction),
                    ('ReducedSolution', ReducedSolution),
                    ('Function', Function),
                    ('ReducedFunction', ReducedFunction),
                    ('FunctionOnBoundary', FunctionOnBoundary),
                    ('ReducedFunctionOnBoundary', ReducedFunctionOnBoundary),
                    ('DiracDeltaFunctions', DiracDeltaFunctions),
                ]:
                    # Skip unsupported speckley function space combinations
                    if domain_module == 'esys.speckley' and fs_name in [
                        'ReducedSolution', 'FunctionOnBoundary', 'ReducedFunctionOnBoundary'
                    ]:
                        continue

                    with self.subTest(function_space=fs_name):
                        fs = fs_constructor(domain)

                        # All ranks participate
                        if my_domain_idx == 0:
                            # Domain 0 sends vector
                            if dim == 2:
                                data = Vector([1.0, 2.0], fs)
                            else:  # dim == 3
                                data = Vector([1.0, 2.0, 3.0], fs)
                            coupler.send(data, dest_domain_index=1, tag=200)
                        elif my_domain_idx == 1:
                            # Domain 1 receives
                            received = coupler.receive(fs, source_domain_index=0, tag=200)

                            # All ranks verify (sup is collective)
                            if dim == 2:
                                self.assertAlmostEqual(sup(received[0]), 1.0, places=10,
                                                     msg=f"{domain_module}.{domain_class} {fs_name} [0]")
                                self.assertAlmostEqual(sup(received[1]), 2.0, places=10,
                                                     msg=f"{domain_module}.{domain_class} {fs_name} [1]")
                            else:  # dim == 3
                                self.assertAlmostEqual(sup(received[0]), 1.0, places=10,
                                                     msg=f"{domain_module}.{domain_class} {fs_name} [0]")
                                self.assertAlmostEqual(sup(received[1]), 2.0, places=10,
                                                     msg=f"{domain_module}.{domain_class} {fs_name} [1]")
                                self.assertAlmostEqual(sup(received[2]), 3.0, places=10,
                                                     msg=f"{domain_module}.{domain_class} {fs_name} [2]")

    def test_exchange_bidirectional(self):
        """Test bidirectional exchange across all domains and function spaces."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # Test each domain configuration
        for domain_module, domain_class, dim, params in DOMAIN_CONFIGS:
            with self.subTest(domain=f"{domain_module}.{domain_class}"):
                domain = create_domain(domain_module, domain_class, params,
                                      domain_array.getDomainComm())

                # Test each function space
                for fs_name, fs_constructor in [
                    ('Solution', Solution),
                    ('ContinuousFunction', ContinuousFunction),
                    ('ReducedSolution', ReducedSolution),
                    ('Function', Function),
                    ('ReducedFunction', ReducedFunction),
                    ('FunctionOnBoundary', FunctionOnBoundary),
                    ('ReducedFunctionOnBoundary', ReducedFunctionOnBoundary),
                    ('DiracDeltaFunctions', DiracDeltaFunctions),
                ]:
                    # Skip unsupported speckley function space combinations
                    if domain_module == 'esys.speckley' and fs_name in [
                        'ReducedSolution', 'FunctionOnBoundary', 'ReducedFunctionOnBoundary'
                    ]:
                        continue

                    with self.subTest(function_space=fs_name):
                        fs = fs_constructor(domain)

                        # All ranks participate
                        if my_domain_idx == 0:
                            send_data = Scalar(100.0, fs)
                            received = coupler.exchange(send_data, fs, peer_domain_index=1, tag=300)
                            # Verify received data from domain 1
                            received_val = sup(received)
                            self.assertAlmostEqual(received_val, 200.0, places=10,
                                                 msg=f"{domain_module}.{domain_class} {fs_name}")
                        elif my_domain_idx == 1:
                            send_data = Scalar(200.0, fs)
                            received = coupler.exchange(send_data, fs, peer_domain_index=0, tag=300)
                            # Verify received data from domain 0
                            received_val = sup(received)
                            self.assertAlmostEqual(received_val, 100.0, places=10,
                                                 msg=f"{domain_module}.{domain_class} {fs_name}")


@unittest.skipIf(SKIP_MPI_TESTS, "MPI with at least 2 processes required")
class Test_DataCoupler_Collective_Data(unittest.TestCase):
    """Tests for DataCoupler collective operations on Data objects."""

    def test_broadcast_data_scalar(self):
        """Test broadcasting scalar Data across all domains and function spaces."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # Test each domain configuration
        for domain_module, domain_class, dim, params in DOMAIN_CONFIGS:
            with self.subTest(domain=f"{domain_module}.{domain_class}"):
                domain = create_domain(domain_module, domain_class, params,
                                      domain_array.getDomainComm())

                # Test each function space
                for fs_name, fs_constructor in [
                    ('Solution', Solution),
                    ('ContinuousFunction', ContinuousFunction),
                    ('ReducedSolution', ReducedSolution),
                    ('Function', Function),
                    ('ReducedFunction', ReducedFunction),
                    ('FunctionOnBoundary', FunctionOnBoundary),
                    ('ReducedFunctionOnBoundary', ReducedFunctionOnBoundary),
                    ('DiracDeltaFunctions', DiracDeltaFunctions),
                ]:
                    # Skip unsupported speckley function space combinations
                    if domain_module == 'esys.speckley' and fs_name in [
                        'ReducedSolution', 'FunctionOnBoundary', 'ReducedFunctionOnBoundary'
                    ]:
                        continue

                    with self.subTest(function_space=fs_name):
                        fs = fs_constructor(domain)

                        # All ranks participate in broadcast
                        if my_domain_idx == 0:
                            # Domain 0 broadcasts
                            data = Scalar(99.0, fs)
                            result = coupler.broadcast(data, root_domain_index=0)
                        else:
                            # Other domains receive
                            result = coupler.broadcast(function_space=fs, root_domain_index=0)

                        # All ranks verify
                        result_val = sup(result)
                        self.assertAlmostEqual(result_val, 99.0, places=10,
                                             msg=f"{domain_module}.{domain_class} {fs_name}")

    def test_allreduce_data_sum(self):
        """Test all-reduce SUM on Data across all domains and function spaces."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # Test each domain configuration
        for domain_module, domain_class, dim, params in DOMAIN_CONFIGS:
            with self.subTest(domain=f"{domain_module}.{domain_class}"):
                domain = create_domain(domain_module, domain_class, params,
                                      domain_array.getDomainComm())

                # Test each function space
                for fs_name, fs_constructor in [
                    ('Solution', Solution),
                    ('ContinuousFunction', ContinuousFunction),
                    ('ReducedSolution', ReducedSolution),
                    ('Function', Function),
                    ('ReducedFunction', ReducedFunction),
                    ('FunctionOnBoundary', FunctionOnBoundary),
                    ('ReducedFunctionOnBoundary', ReducedFunctionOnBoundary),
                    ('DiracDeltaFunctions', DiracDeltaFunctions),
                ]:
                    # Skip unsupported speckley function space combinations
                    if domain_module == 'esys.speckley' and fs_name in [
                        'ReducedSolution', 'FunctionOnBoundary', 'ReducedFunctionOnBoundary'
                    ]:
                        continue

                    with self.subTest(function_space=fs_name):
                        print(domain_module, fs_name)
                        fs = fs_constructor(domain)

                        # All ranks participate - each domain contributes its index
                        data = Scalar(float(my_domain_idx), fs)
                        result = coupler.allreduce(data, op=MPI.SUM)

                        # All ranks verify
                        expected = sum(range(num_domains))  # 0 + 1 = 1
                        result_val = sup(result)
                        self.assertAlmostEqual(result_val, expected, places=10,
                                             msg=f"{domain_module}.{domain_class} {fs_name}")
                        print("DONE")

@unittest.skipIf(SKIP_MPI_TESTS, "MPI with at least 2 processes required")
class Test_DataCoupler_Collective_Values(unittest.TestCase):
    """Tests for DataCoupler collective operations on scalar values (not Data objects)."""

    def test_broadcast_value_scalar(self):
        """Test broadcasting a scalar value."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # Create a dummy domain for isRootRank check
        domain = create_domain('esys.ripley', 'Rectangle', {'n0': 10, 'n1': 10},
                               domain_array.getDomainComm())

        # All ranks participate in broadcast
        if my_domain_idx == 0 and domain.isRootRank():
            result = coupler.broadcast_value(3.14159, root_domain_index=0)
        else:
            result = coupler.broadcast_value(root_domain_index=0)

        # All ranks verify
        self.assertAlmostEqual(result, 3.14159, places=5)

    def test_allreduce_value_sum(self):
        """Test all-reduce with SUM operation on scalars."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # All ranks participate
        my_value = float(my_domain_idx)
        result = coupler.allreduce_value(my_value, op=MPI.SUM)

        # All ranks verify
        expected = sum(range(num_domains))
        self.assertAlmostEqual(result, expected)

    def test_allreduce_value_max(self):
        """Test all-reduce with MAX operation."""
        num_domains = 2
        domain_array = MPIDomainArray(numDomains=num_domains, comm=MPI.COMM_WORLD)
        my_domain_idx = domain_array.getDomainIndex()
        coupler = DataCoupler(domain_array)

        # All ranks participate
        my_value = float(my_domain_idx)
        result = coupler.allreduce_value(my_value, op=MPI.MAX)

        # All ranks verify
        expected = float(num_domains - 1)
        self.assertAlmostEqual(result, expected)


# Create test suite
def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(Test_MPIDomainArray))
    test_suite.addTest(unittest.makeSuite(Test_DataCoupler_PointToPoint))
    test_suite.addTest(unittest.makeSuite(Test_DataCoupler_Collective_Data))
    test_suite.addTest(unittest.makeSuite(Test_DataCoupler_Collective_Values))
    return test_suite


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
