
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""
Domain Array and Data Coupler for Multi-Domain Simulations
============================================================

This module provides infrastructure for multi-domain parallel simulations
using MPI communicator topology and data exchange between domains.

Classes
-------

- **MPIDomainArray**: MPI communicator topology for domain arrays.
  Manages 2D Cartesian communicator topology with domain and subdomain
  communicators for multi-domain simulations.

- **DataCoupler**: Communication of Data objects between domains.
  Handles sending, receiving, broadcasting, and reduction of esys.escript.Data
  objects between identically distributed domains.

Example Usage
-------------

::

    from mpi4py import MPI
    from esys.escript import Rectangle, Scalar, Solution
    from esys.escript.domaincoupler import MPIDomainArray, DataCoupler

    # Create domain array with 2 domains
    domain_array = MPIDomainArray(numDomains=2, comm=MPI.COMM_WORLD)

    # Create domains (same mesh parameters for identical distribution)
    domain = Rectangle(n0=50, n1=50, comm=domain_array.getDomainComm())

    # Create coupler for data exchange
    coupler = DataCoupler(domain_array)

    # Exchange data between domains
    if domain_array.getDomainIndex() == 0:
        data = Scalar(1.0, Solution(domain))
        coupler.send(data, dest_domain_index=1, tag=100)
    elif domain_array.getDomainIndex() == 1:
        received = coupler.receive(Solution(domain), source_domain_index=0, tag=100)

"""

from mpi4py import MPI
import numpy as np
from esys.escript import Data


__all__ = ['MPIDomainArray', 'DataCoupler']


class MPIDomainArray(object):
    """
    MPI communicator topology for multi-domain simulations.

    Creates a 2D Cartesian communicator topology where:
    - First dimension: different domains (e.g., different physics, ensemble members)
    - Second dimension: subdomains within each domain (for domain decomposition)

    This allows each domain to have its own MPI communicator for parallel operations
    while enabling communication between corresponding ranks across domains.

    Parameters:
    -----------
    numDomains : int
        Number of independent domains
    numSubDomains : int, optional
        Number of MPI ranks per domain. If None, computed as total_ranks // numDomains
    strict : bool, optional
        If True, requires total ranks to exactly match numDomains * numSubDomains.
        Default: True
    comm : MPI.Comm, optional
        Parent MPI communicator. Default: MPI.COMM_WORLD

    Attributes:
    -----------
    numDomains : int
        Number of domains
    numSubDomains : int
        Number of subdomains per domain
    comm : MPI.Comm
        Parent communicator
    domain_comm : MPI.Comm
        Communicator for ranks within same domain (for intra-domain parallelization)
    subdomain_comm : MPI.Comm
        Communicator for corresponding ranks across domains (for inter-domain communication)
    idom : int
        This rank's domain index (0 to numDomains-1)
    isub : int
        This rank's subdomain index (0 to numSubDomains-1)

    Example::

        # 8 total ranks: 2 domains x 4 subdomains each
        domain_array = MPIDomainArray(numDomains=2, comm=MPI.COMM_WORLD)

        # Create domain using the domain communicator
        from esys.ripley import Rectangle
        domain = Rectangle(n0=100, n1=100, comm=domain_array.getDomainComm())

        # Ranks 0-3: domain 0, subdomains 0-3
        # Ranks 4-7: domain 1, subdomains 0-3
        # subdomain_comm connects rank 0 with rank 4, rank 1 with rank 5, etc.
    """

    def __init__(self, numDomains, numSubDomains=None, strict=True, comm=None):
        self.numDomains = numDomains
        if comm is None:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = comm

        if numSubDomains is None:
            numSubDomains = self.comm.Get_size() // self.numDomains
        self.numSubDomains = numSubDomains

        myrank = self.comm.Get_rank()
        size = self.comm.Get_size()

        if strict and not numDomains * numSubDomains == size:
            raise ValueError(f"Requested ranks {numDomains * numSubDomains} does not match available ranks {size}")
        if numDomains * numSubDomains > size:
            raise ValueError(f"Requested ranks {numDomains*numSubDomains} is larger than available ranks {size}")

        cc = comm.Create_cart([numDomains, numSubDomains], periods=[False, True], reorder=False)
        if myrank < cc.Get_size():
           self.group_by_subdomain = {}
           for isub in range(cc.dims[1]):
               dom_group_list = [cc.Get_cart_rank((idom, isub)) for idom in range(cc.dims[0])]
               if myrank in dom_group_list:
                   self.isub = isub
                   self.dom_group_list = dom_group_list
               self.group_by_subdomain[isub] = comm.Get_group().Incl(dom_group_list)
           self.subdomain_comm = self.comm.Create(self.group_by_subdomain[self.isub])

           self.group_by_domain = {}
           for idom in range(cc.dims[0]):
                    sub_group_list = [cc.Get_cart_rank((idom, isub)) for isub in range(cc.dims[1])]
                    if myrank in sub_group_list:
                        self.idom = idom
                        self.sub_group_list = sub_group_list
                    self.group_by_domain[idom] = comm.Get_group().Incl(sub_group_list)
           self.domain_comm = self.comm.Create(self.group_by_domain[self.idom])

    def getSubdomainComm(self):
        """
        Get the subdomain communicator.

        Returns:
        --------
        MPI.Comm
            Communicator connecting corresponding ranks across all domains.
            Size equals numDomains.
        """
        return self.subdomain_comm

    def getDomainComm(self):
        """
        Get the domain communicator.

        Returns:
        --------
        MPI.Comm
            Communicator for all ranks within this domain.
            Size equals numSubDomains. Use this for domain creation.
        """
        return self.domain_comm

    def getSubdomainIndex(self):
        """
        Get this rank's subdomain index within its domain.

        Returns:
        --------
        int
            Subdomain index (0 to numSubDomains-1)
        """
        return self.isub

    def getDomainIndex(self):
        """
        Get this rank's domain index.

        Returns:
        --------
        int
            Domain index (0 to numDomains-1)
        """
        return self.idom


class DataCoupler:
    """
    Handles communication of Data objects between identically distributed domains.

    This class manages point-to-point and collective communication of esys.escript.Data
    objects between domains in an MPIDomainArray. Each rank communicates with its
    corresponding rank in other domains using the subdomain communicator.

    Key Assumption:
        Source and target domains must have the same FunctionSpace structure and
        identical domain decomposition (same number of samples per rank, same distribution).

        This is guaranteed when using MPIDomainArray where all domains are created with
        the same mesh parameters and distributed across communicators of the same size.

    Parameters:
    -----------
    mpi_domain_array : MPIDomainArray
        The domain array managing communicator topology

    Methods:
    --------
    send(data, dest_domain_index, tag=0)
        Send Data to corresponding rank in destination domain

    receive(function_space, source_domain_index, tag=0)
        Receive Data from corresponding rank in source domain

    exchange(send_data, recv_function_space, peer_domain_index, tag=0)
        Bidirectional exchange between two domains

    broadcast(data=None, function_space=None, root_domain_index=0)
        Broadcast Data from root domain to all domains

    allreduce(data, op=MPI.SUM)
        All-reduce Data point-by-point across all domains

    broadcast_value(value=None, root_domain_index=0)
        Broadcast scalar/array from root domain

    allreduce_value(value, op=MPI.SUM)
        All-reduce scalar/array across all domains

    Example::

        from mpi4py import MPI
        from esys.ripley import Rectangle
        from esys.escript import Scalar, Solution
        from esys.escript.domaincoupler import MPIDomainArray, DataCoupler

        # Setup domain array
        domain_array = MPIDomainArray(numDomains=2, comm=MPI.COMM_WORLD)
        domain = Rectangle(n0=50, n1=50, comm=domain_array.getDomainComm())
        coupler = DataCoupler(domain_array)

        # Point-to-point communication
        if domain_array.getDomainIndex() == 0:
            data = Scalar(1.0, Solution(domain))
            coupler.send(data, dest_domain_index=1, tag=100)
        elif domain_array.getDomainIndex() == 1:
            received = coupler.receive(Solution(domain), source_domain_index=0, tag=100)

        # Collective operations
        my_data = Scalar(domain_array.getDomainIndex(), Solution(domain))
        sum_data = coupler.allreduce(my_data, op=MPI.SUM)
    """

    def __init__(self, mpi_domain_array):
        """
        Initialize DataCoupler with an MPIDomainArray.

        Parameters:
        -----------
        mpi_domain_array : MPIDomainArray
            The domain array managing communicator topology
        """
        self.mpi_domain_array = mpi_domain_array
        self.subdomain_comm = mpi_domain_array.getSubdomainComm()
        self.my_domain_index = mpi_domain_array.getDomainIndex()
        self.my_subdomain_index = mpi_domain_array.getSubdomainIndex()

    def _get_peer_rank(self, peer_domain_index):
        """
        Get the subdomain communicator rank for a peer domain.

        Since domains are arranged in a cartesian topology, the peer rank
        in the subdomain communicator corresponds to the peer domain index.
        """
        return peer_domain_index

    def send(self, data, dest_domain_index, tag=0):
        """
        Send a Data object to the corresponding rank in destination domain.

        This performs a rank-to-rank send where each rank sends its local portion
        of the Data to the same subdomain index in the destination domain.

        Parameters:
        -----------
        data : esys.escript.Data
            The Data object to send (will use local portion on this rank)
        dest_domain_index : int
            The index of the destination domain in the MPIDomainArray
        tag : int, optional
            MPI message tag (default: 0)

        Notes:
        ------
        - Data is automatically expanded and resolved before sending
        - Only the local portion on each rank is sent
        - Complex data is supported (sends real and imaginary parts separately)
        """
        if dest_domain_index == self.my_domain_index:
            raise ValueError("Cannot send to same domain (use local copy instead)")

        # Expand and resolve data before sending
        if not data.isExpanded():
            data.expand()
        if data.isLazy():
            data.resolve()

        # Get metadata
        shape = data.getShape()
        is_complex = data.isComplex()
        is_expanded = data.isExpanded()
        num_data_points = data.getNumberOfDataPoints()

        # Allocate numpy array for local data
        # Shape is (num_data_points,) + shape
        if len(shape) == 0:
            # Scalar data
            array_shape = (num_data_points,)
        else:
            # Vector/tensor data
            array_shape = (num_data_points,) + tuple(shape)

        if is_complex:
            local_data = np.empty(array_shape, dtype=np.complex128)
        else:
            local_data = np.empty(array_shape, dtype=np.float64)

        # Copy data point-by-point into numpy array
        for i in range(num_data_points):
            val = data.getTupleForDataPoint(i)
            # For scalar data, getTupleForDataPoint returns (value,) so extract it
            if len(shape) == 0:
                local_data[i] = val[0] if isinstance(val, tuple) else val
            else:
                local_data[i] = val

        # Flatten for sending
        local_data = local_data.flatten()

        # Get destination rank in subdomain communicator
        dest_rank = self._get_peer_rank(dest_domain_index)

        # Send metadata (shape, expanded flag, complex flag)
        metadata = {
            'shape': shape,
            'is_complex': is_complex,
            'is_expanded': is_expanded,
            'data_size': len(local_data)
        }
        self.subdomain_comm.send(metadata, dest=dest_rank, tag=tag)

        # Send data
        if is_complex:
            # Send real and imaginary parts separately
            self.subdomain_comm.Send(local_data.real, dest=dest_rank, tag=tag+1)
            self.subdomain_comm.Send(local_data.imag, dest=dest_rank, tag=tag+2)
        else:
            self.subdomain_comm.Send(local_data, dest=dest_rank, tag=tag+1)

    def receive(self, function_space, source_domain_index, tag=0):
        """
        Receive a Data object from the corresponding rank in source domain.

        This performs a rank-to-rank receive where each rank receives data
        from the same subdomain index in the source domain.

        Parameters:
        -----------
        function_space : esys.escript.FunctionSpace
            The FunctionSpace for the received Data (must be compatible with sender)
        source_domain_index : int
            The index of the source domain in the MPIDomainArray
        tag : int, optional
            MPI message tag (must match sender's tag, default: 0)

        Returns:
        --------
        data : esys.escript.Data
            The received Data object on the specified FunctionSpace

        Notes:
        ------
        - The FunctionSpace must have the same structure as the sender's
        - Each rank receives only its local portion
        - Returns a Data object with the same properties (expanded/constant, complex/real)
        """
        if source_domain_index == self.my_domain_index:
            raise ValueError("Cannot receive from same domain (use local copy instead)")

        # Get source rank in subdomain communicator
        source_rank = self._get_peer_rank(source_domain_index)

        # Receive metadata
        metadata = self.subdomain_comm.recv(source=source_rank, tag=tag)

        shape = metadata['shape']
        is_complex = metadata['is_complex']
        is_expanded = metadata['is_expanded']
        data_size = metadata['data_size']

        # Receive data
        if is_complex:
            real_part = np.empty(data_size, dtype=np.float64)
            imag_part = np.empty(data_size, dtype=np.float64)
            self.subdomain_comm.Recv(real_part, source=source_rank, tag=tag+1)
            self.subdomain_comm.Recv(imag_part, source=source_rank, tag=tag+2)
            local_data = real_part + 1j * imag_part
        else:
            local_data = np.empty(data_size, dtype=np.float64)
            self.subdomain_comm.Recv(local_data, source=source_rank, tag=tag+1)

        # Reshape received data
        if len(shape) == 0:
            # Scalar data
            array_shape = (data_size,)
        else:
            # Vector/tensor data
            array_shape = (data_size // np.prod(shape),) + tuple(shape)

        local_data = local_data.reshape(array_shape)

        # Create empty Data object and populate point-by-point
        received_data = Data(0.0, shape, function_space)
        # Expand to allow point-by-point setting
        received_data.expand()

        num_data_points = received_data.getNumberOfDataPoints()

        for i in range(num_data_points):
            val_tuple = tuple(local_data[i]) if len(shape) > 0 else local_data[i].item()
            received_data.setValueOfDataPoint(i, val_tuple)

        return received_data

    def exchange(self, send_data, recv_function_space, peer_domain_index, tag=0):
        """
        Bidirectional exchange of Data between two domains.

        This is a convenience method for simultaneous send and receive.
        Both domains must use the same tag value.

        Parameters:
        -----------
        send_data : esys.escript.Data
            Data to send
        recv_function_space : esys.escript.FunctionSpace
            FunctionSpace for received data
        peer_domain_index : int
            Index of the peer domain
        tag : int, optional
            Base tag for messages (default: 0)
            Note: Uses tag, tag+1, tag+2 internally

        Returns:
        --------
        received_data : esys.escript.Data
            The received Data

        Notes:
        ------
        - Both domains must call this method with the same tag
        - Operations are ordered by domain index to avoid deadlock
        """
        if peer_domain_index == self.my_domain_index:
            raise ValueError("Cannot exchange with same domain (use local copy instead)")

        # Order operations by domain index to avoid deadlock
        if self.my_domain_index < peer_domain_index:
            # Lower index: send first, then receive
            self.send(send_data, peer_domain_index, tag=tag)
            received_data = self.receive(recv_function_space, peer_domain_index, tag=tag)
        else:
            # Higher index: receive first, then send
            received_data = self.receive(recv_function_space, peer_domain_index, tag=tag)
            self.send(send_data, peer_domain_index, tag=tag)

        return received_data

    def broadcast_value(self, value=None, root_domain_index=0):
        """
        Broadcast a scalar or numpy array from root domain to all other domains.

        The value is replicated across all ranks in all domains, maintaining
        the escript invariant that scalars/arrays are identical across all subdomains.

        Parameters:
        -----------
        value : float, int, numpy.ndarray, or None
            The value to broadcast (only used on rank 0 of root domain)
        root_domain_index : int, optional
            The index of the root domain (default: 0)

        Returns:
        --------
        result : float, int, or numpy.ndarray
            The broadcasted value (same on all ranks of all domains)

        Example::

            # Broadcast a scalar
            if my_domain_idx == 0:
                max_val = coupler.broadcast_value(42.0, root_domain_index=0)
            else:
                max_val = coupler.broadcast_value(root_domain_index=0)

            # Broadcast an array
            if my_domain_idx == 0:
                coeffs = coupler.broadcast_value(np.array([1.0, 2.0, 3.0]), root_domain_index=0)
            else:
                coeffs = coupler.broadcast_value(root_domain_index=0)
        """
        domain_comm = self.mpi_domain_array.getDomainComm()

        # First, broadcast from root rank 0 to all corresponding rank 0s via subdomain comm
        if domain_comm.Get_rank() == 0:
            result = self.subdomain_comm.bcast(value, root=root_domain_index)
        else:
            result = None

        # Then broadcast within each domain to replicate across all subdomains
        result = domain_comm.bcast(result, root=0)

        return result

    def allreduce_value(self, value, op=MPI.SUM):
        """
        All-reduce a scalar or numpy array across all domains.

        Since the value is the same on all ranks within a domain (escript invariant),
        only rank 0 of each domain participates in the reduction, then the result
        is broadcast within each domain.

        Parameters:
        -----------
        value : float, int, or numpy.ndarray
            The value to reduce (should be same on all ranks within domain)
        op : MPI.Op, optional
            MPI reduction operation (default: MPI.SUM)
            Options: MPI.SUM, MPI.MAX, MPI.MIN, MPI.PROD

        Returns:
        --------
        result : float, int, or numpy.ndarray
            Reduced value (same on all ranks of all domains)

        Example::

            # Reduce scalars from all domains
            my_max_temp = Lsup(temperature)  # Same on all ranks within domain
            global_max_temp = coupler.allreduce_value(my_max_temp, op=MPI.MAX)

            # Reduce arrays from all domains
            my_histogram = compute_histogram(data)  # Same on all ranks within domain
            total_histogram = coupler.allreduce_value(my_histogram, op=MPI.SUM)
        """
        domain_comm = self.mpi_domain_array.getDomainComm()

        # Only rank 0 of each domain participates in subdomain reduction
        if domain_comm.Get_rank() == 0:
            if isinstance(value, np.ndarray):
                # Use Allreduce for numpy arrays
                result = np.empty_like(value)
                self.subdomain_comm.Allreduce(value, result, op=op)
            else:
                # Use allreduce for scalars
                result = self.subdomain_comm.allreduce(value, op=op)
        else:
            result = None

        # Broadcast result to all ranks within each domain
        result = domain_comm.bcast(result, root=0)

        return result

    def broadcast(self, data=None, function_space=None, root_domain_index=0):
        """
        Broadcast a Data object from root domain to all other domains.

        This performs a collective broadcast where each rank in the root domain
        broadcasts its local portion to the corresponding rank in all other domains.

        Parameters:
        -----------
        data : esys.escript.Data or None
            The Data object to broadcast (only used on root domain)
        function_space : esys.escript.FunctionSpace or None
            The FunctionSpace for received Data (only used on non-root domains)
        root_domain_index : int, optional
            The index of the root domain (default: 0)

        Returns:
        --------
        result : esys.escript.Data
            The broadcasted Data object (same on all domains)

        Example::

            # Domain 0 broadcasts initial conditions to all ensemble members
            if my_domain_idx == 0:
                initial_temp = Scalar(20.0 + compute_perturbation(), Solution(domain))
                temp = coupler.broadcast(data=initial_temp, root_domain_index=0)
            else:
                temp = coupler.broadcast(function_space=Solution(domain), root_domain_index=0)

        Notes:
        ------
        - This is a collective operation - all domains must call it
        - Root domain must provide 'data', non-root domains must provide 'function_space'
        - All Data/FunctionSpace must have compatible structure
        - Each rank independently broadcasts with corresponding ranks in other domains
        - Works with multi-rank domains via subdomain communicator
        """
        is_root = (self.my_domain_index == root_domain_index)

        if is_root:
            if data is None:
                raise ValueError("Root domain must provide 'data' parameter")

            # Expand and resolve data before broadcasting
            if not data.isExpanded():
                data.expand()
            if data.isLazy():
                data.resolve()

            # Get metadata
            shape = data.getShape()
            is_complex = data.isComplex()
            num_data_points = data.getNumberOfDataPoints()
            function_space = data.getFunctionSpace()

        else:
            if function_space is None:
                raise ValueError("Non-root domains must provide 'function_space' parameter")

            # These will be received via broadcast
            shape = None
            is_complex = None
            num_data_points = None

        # Broadcast metadata from root to all domains
        metadata = {
            'shape': shape if is_root else None,
            'is_complex': is_complex if is_root else None,
            'num_data_points': num_data_points if is_root else None
        }
        metadata = self.subdomain_comm.bcast(metadata, root=root_domain_index)

        shape = metadata['shape']
        is_complex = metadata['is_complex']
        num_data_points = metadata['num_data_points']

        # Check if complex (not yet supported)
        if is_complex:
            raise NotImplementedError("Broadcast not yet implemented for complex Data")

        # Allocate numpy array for data
        if len(shape) == 0:
            array_shape = (num_data_points,)
        else:
            array_shape = (num_data_points,) + tuple(shape)

        local_array = np.empty(array_shape, dtype=np.float64)

        if is_root:
            # Extract local values from Data
            for i in range(num_data_points):
                val = data.getTupleForDataPoint(i)
                if len(shape) == 0:
                    local_array[i] = val[0] if isinstance(val, tuple) else val
                else:
                    local_array[i] = val

        # Broadcast array from root to all domains via subdomain communicator
        self.subdomain_comm.Bcast(local_array, root=root_domain_index)

        # Reconstruct Data object from broadcasted values
        result = Data(0.0, shape, function_space)
        result.expand()

        for i in range(num_data_points):
            val_tuple = tuple(local_array[i]) if len(shape) > 0 else local_array[i].item()
            result.setValueOfDataPoint(i, val_tuple)

        return result

    def allreduce(self, data, op=MPI.SUM):
        """
        All-reduce Data point-by-point across all domains in the subdomain communicator.

        This performs a collective reduction where each data point is reduced independently
        across all domains. Each rank reduces its local portion of Data with the corresponding
        rank in all other domains. The result has the same FunctionSpace and structure as the input.

        Parameters:
        -----------
        data : esys.escript.Data
            The Data object to reduce (local portion on this rank)
        op : MPI.Op, optional
            MPI reduction operation (default: MPI.SUM)
            Options: MPI.SUM, MPI.MAX, MPI.MIN, MPI.PROD

        Returns:
        --------
        result : esys.escript.Data
            Reduced Data object on the same FunctionSpace (available on all domains)

        Example::

            # Sum temperature fields from all ensemble members
            coupler = DataCoupler(domain_array)
            my_temperature = Scalar(value_array, Solution(domain))
            total_temperature = coupler.allreduce(my_temperature, op=MPI.SUM)
            avg_temperature = total_temperature / num_domains

        Notes:
        ------
        - This is a collective operation - all domains must call it
        - All Data objects must have compatible shapes and function spaces
        - Reduction happens point-by-point across domains via subdomain communicator
        - Each rank independently reduces with corresponding ranks in other domains
        - Works with multi-rank domains (each rank reduces its local portion)
        """
        # Expand and resolve data before reduction
        if not data.isExpanded():
            data.expand()
        if data.isLazy():
            data.resolve()

        # Check if complex (not yet supported for allreduce)
        if data.isComplex():
            raise NotImplementedError("All-reduce not yet implemented for complex Data")

        # Get the FunctionSpace for reconstruction
        function_space = data.getFunctionSpace()

        # Each rank independently reduces its local data with corresponding
        # ranks in other domains via the subdomain communicator

        # Get number of data points
        num_data_points = data.getNumberOfDataPoints()
        shape = data.getShape()

        # Allocate numpy array for local data
        if len(shape) == 0:
            # Scalar data
            array_shape = (num_data_points,)
        else:
            # Vector/tensor data
            array_shape = (num_data_points,) + tuple(shape)

        local_array = np.empty(array_shape, dtype=np.float64)

        # Extract all local values point-by-point
        for i in range(num_data_points):
            val = data.getTupleForDataPoint(i)
            # For scalar data, getTupleForDataPoint returns (value,) so extract it
            if len(shape) == 0:
                local_array[i] = val[0] if isinstance(val, tuple) else val
            else:
                local_array[i] = val

        # Allreduce across subdomain communicator
        reduced_array = np.empty_like(local_array)
        self.subdomain_comm.Allreduce(local_array, reduced_array, op=op)

        # Reconstruct Data object from reduced values
        result = Data(0.0, shape, function_space)
        # Expand to allow point-by-point setting
        result.expand()

        for i in range(num_data_points):
            val_tuple = tuple(reduced_array[i]) if len(shape) > 0 else reduced_array[i].item()
            result.setValueOfDataPoint(i, val_tuple)

        return result
