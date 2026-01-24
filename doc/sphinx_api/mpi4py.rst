.. _chap_mpi4py:

==============================
Using mpi4py with escript
==============================

Introduction
============

`mpi4py <https://mpi4py.readthedocs.io/>`_ is a Python package that provides bindings to the
Message Passing Interface (MPI) standard. When escript is built with mpi4py support enabled,
you can use standard MPI programming techniques directly in your Python scripts to control
how escript uses parallel resources.

This chapter explains how to use mpi4py to:

* Pass custom MPI communicators to escript domains
* Run multiple simulations simultaneously on different process groups
* Implement ensemble simulations and parameter sweeps
* Couple multiple physics domains on separate communicators
* Retrieve MPI communicators from escript objects for custom communication

Why use mpi4py with escript?
----------------------------

By default, when escript runs with multiple MPI processes, it automatically uses all available
processes to solve each PDE. This works well for large, high-resolution problems. However,
there are scenarios where you might want more control:

1. **Ensemble simulations**: Run multiple parameter values simultaneously, each on a smaller group of processes
2. **Parameter sweeps**: Explore parameter space by running many simulations in parallel
3. **Multi-physics coupling**: Solve different physics on different process groups
4. **Uncertainty quantification**: Run multiple realizations of stochastic problems
5. **Optimization**: Evaluate multiple trial solutions in parallel

Using mpi4py gives you direct control over how processes are organized and allows you to
implement these patterns using standard MPI techniques.

Requirements
------------

To use mpi4py with escript:

* escript must be built with ``mpi4py=True`` in the build configuration
* mpi4py must be installed and compiled against the same MPI implementation as escript
* Use ``mpi='auto'`` in the build configuration to automatically detect the correct MPI implementation

See the :doc:`installation guide <installation>` for detailed instructions.


MPI Initialization and the pythonMPI Wrapper
============================================

.. important::

   escript requires MPI to be initialized **before** the Python interpreter starts to ensure
   proper initialization of C++ libraries like Trilinos. This section explains why and how to
   do it correctly.

Why MPI Initialization Order Matters
------------------------------------

escript uses the Trilinos library for distributed linear algebra operations. Trilinos contains
C++ code with static objects that capture MPI state when shared libraries are loaded. If MPI
is initialized **after** Python starts, some C++ static objects may capture an uninitialized
MPI state, leading to inconsistent behavior and errors.

**Symptoms of incorrect initialization:**

* Trilinos errors like: ``global index N does not have a corresponding local index in the Directory Map``
* Hangs during interpolation between function spaces
* Crashes or incorrect results when using multiple MPI processes

These errors occur because different parts of the C++ code disagree about how many global indices
exist, leading to off-by-one errors and invalid memory access.

The pythonMPI Wrapper
---------------------

To solve this problem, escript provides a custom Python launcher called ``pythonMPI`` that ensures
MPI is initialized **before** Python starts:

1. The ``pythonMPI`` executable (a C++ program) calls ``MPI_Init_thread()``
2. MPI runtime is fully initialized with ``MPI_COMM_WORLD`` ready
3. ``pythonMPI`` then starts the Python interpreter via ``Py_Main()``
4. Python and all C++ libraries load with MPI already active
5. All C++ static initialization sees a consistent, initialized MPI state

The ``run-escript`` launcher automatically uses ``pythonMPI`` when running with MPI.

Correct Usage: Always Use run-escript
-------------------------------------

**CORRECT:**

.. code-block:: bash

   # For MPI jobs, always use run-escript
   run-escript -n 4 myscript.py

   # For serial jobs, run-escript also works
   run-escript myscript.py

**WRONG:**

.. code-block:: bash

   # DO NOT use mpirun/mpiexec with python3 directly
   mpirun -n 4 python3 myscript.py  # Will cause Trilinos errors!


What Happens with Direct mpirun python3
---------------------------------------

When you run ``mpirun -n 4 python3 myscript.py``, the initialization sequence is incorrect:

1. Python interpreter starts **without** MPI
2. Python imports mpi4py
3. mpi4py calls ``MPI_Init()`` ← MPI initialized **after** Python started
4. Python imports escript modules
5. Some Trilinos C++ static objects were created **before** ``MPI_Init()``
6. Later code uses MPI state from **after** ``MPI_Init()``
7. Inconsistent state causes index calculation errors

**Example error:** For a mesh with 861 nodes (indices 0--860), some code may think there are
862 nodes and try to access index 862, which doesn't exist.

Import Order When Using run-escript
-----------------------------------

When using ``run-escript``, the import order doesn't matter because MPI is already initialized:

.. code-block:: python

   # Both work with run-escript
   from esys.ripley import Rectangle
   from mpi4py import MPI

However, for portability and clarity, it's still good practice to import mpi4py first:

.. code-block:: python

   # Best practice - import mpi4py first
   from mpi4py import MPI
   from esys.escript import *
   from esys.ripley import Rectangle

Summary: Best Practices
-----------------------

* **Always use** ``run-escript`` for MPI programs (required for correct Trilinos initialization)
* **Never use** ``mpirun python3`` directly (will cause Trilinos errors)
* **Always import** mpi4py before escript modules (good practice even with ``run-escript``)
* The ``pythonMPI`` wrapper ensures MPI initializes before Python starts
* This prevents C++ static initialization order issues with Trilinos


Basic Concepts
==============

MPI Communicators
-----------------

An MPI communicator defines a group of processes that can communicate with each other.
By default, all processes belong to ``MPI_COMM_WORLD``, which includes every process in the job.

In mpi4py, you can access this default communicator:

.. code-block:: python

   from mpi4py import MPI

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()  # My process rank (0 to size-1)
   size = comm.Get_size()  # Total number of processes

Splitting Communicators
-----------------------

You can create smaller communicators by splitting ``MPI_COMM_WORLD``:

.. code-block:: python

   # Split into groups based on color
   color = rank % 2  # Processes with same color go in same group
   sub_comm = comm.Split(color, rank)

   # Now sub_comm is a communicator with only part of the processes
   sub_rank = sub_comm.Get_rank()
   sub_size = sub_comm.Get_size()

Processes with the same ``color`` are grouped together into a new communicator.
The second argument to ``Split()`` determines the rank ordering within each new group.


Passing Communicators to escript Domains
========================================

All escript domain constructors accept an optional ``comm`` parameter that allows you to
specify which MPI communicator the domain should use:

.. code-block:: python

   from mpi4py import MPI
   from esys.ripley import Rectangle

   # Create a domain using default MPI_COMM_WORLD
   domain = Rectangle(n0=100, n1=100)

   # Or explicitly pass MPI_COMM_WORLD
   domain = Rectangle(n0=100, n1=100, comm=MPI.COMM_WORLD)

   # Or use a custom communicator
   sub_comm = MPI.COMM_WORLD.Split(color, rank)
   domain = Rectangle(n0=100, n1=100, comm=sub_comm)

The ``comm`` parameter is available in all domain types:

* **Finley**: ``ReadMesh()``, ``ReadGmsh()``, ``Rectangle()``, ``Brick()``
* **Ripley**: ``Rectangle()``, ``Brick()``, ``MultiRectangle()``, ``MultiBrick()``
* **Speckley**: ``Rectangle()``, ``Brick()``
* **Oxley**: ``Rectangle()``, ``Brick()`` (requires Trilinos)


Retrieving Communicators from escript Objects
=============================================

You can retrieve the MPI communicator from any escript domain, function space, or data object:

.. code-block:: python

   # From a domain
   domain_comm = domain.getMPIComm()

   # From a function space
   fs = ContinuousFunction(domain)
   fs_comm = fs.getMPIComm()

   # From a data object
   data = Scalar(1.0, fs)
   data_comm = data.getMPIComm()

   # All return the same communicator
   assert MPI.Comm.Compare(domain_comm, fs_comm) == MPI.IDENT
   assert MPI.Comm.Compare(fs_comm, data_comm) == MPI.IDENT

This allows you to:

* Verify which communicator a domain is using
* Perform custom MPI operations using the domain's communicator
* Pass the communicator to other libraries or functions

.. important::

   When comparing MPI communicators, always use ``MPI.Comm.Compare()``, never use ``==``.
   The comparison returns ``MPI.IDENT`` if the communicators are identical, ``MPI.CONGRUENT``
   if they have the same processes but different contexts, or ``MPI.UNEQUAL`` if they differ.


Simple Example: Ensemble Simulation
===================================

This example demonstrates running multiple parameter values in parallel.
Each process group solves the same problem with a different parameter:

.. code-block:: python

   from mpi4py import MPI
   from esys.escript import *
   from esys.escript.linearPDEs import LinearPDE
   from esys.ripley import Rectangle

   # Get MPI information
   world_comm = MPI.COMM_WORLD
   world_rank = world_comm.Get_rank()
   world_size = world_comm.Get_size()

   # Define parameter values to explore
   parameters = [0.1, 0.5, 1.0, 2.0]
   num_params = len(parameters)

   # Ensure we have enough processes
   if world_size < num_params:
       if world_rank == 0:
           print(f"Need at least {num_params} processes for this example")
       exit(1)

   # Assign each process to a parameter
   param_id = world_rank % num_params
   my_parameter = parameters[param_id]

   # Split communicator: processes with same parameter work together
   processes_per_param = world_size // num_params
   group_comm = world_comm.Split(param_id, world_rank)

   # Create domain with group communicator
   domain = Rectangle(n0=50, n1=50, comm=group_comm)

   # Solve PDE with this parameter value
   pde = LinearPDE(domain)
   x = domain.getX()
   pde.setValue(A=my_parameter*kronecker(2), Y=x[0]*x[1])
   solution = pde.getSolution()

   # Each group can now process its result independently
   if group_comm.Get_rank() == 0:
       result = integrate(solution)
       print(f"Parameter {my_parameter}: Result = {result}")

   # Clean up communicator
   group_comm.Free()

**Running this example:**

.. code-block:: bash

   run-escript -n 8 ensemble_example.py

With 8 processes and 4 parameter values, each parameter gets 2 processes.
All 4 simulations run simultaneously.


Advanced Example: Parameter Sweep
=================================

This example shows a more realistic scenario where different grid resolutions are tested simultaneously:

.. code-block:: python

   from mpi4py import MPI
   from esys.escript import *
   from esys.escript.linearPDEs import LinearPDE
   from esys.ripley import Rectangle
   import time

   world_comm = MPI.COMM_WORLD
   world_rank = world_comm.Get_rank()
   world_size = world_comm.Get_size()

   # Grid resolutions to test
   grid_sizes = [
       (10, 10),   # Coarse
       (20, 20),   # Medium
       (40, 40),   # Fine
       (80, 80),   # Very fine
   ]
   num_grids = len(grid_sizes)

   # Split into groups
   grid_id = world_rank % num_grids
   nx, ny = grid_sizes[grid_id]
   grid_comm = world_comm.Split(grid_id, world_rank)
   grid_rank = grid_comm.Get_rank()

   if world_rank == 0:
       print(f"Running {num_grids} grid resolutions simultaneously")
       print(f"Total processes: {world_size}")

   # Synchronize before starting
   world_comm.Barrier()
   start_time = time.time()

   # Create domain with this grid size
   domain = Rectangle(n0=nx, n1=ny, comm=grid_comm)

   # Solve a test problem
   pde = LinearPDE(domain)
   x = domain.getX()
   pde.setValue(A=kronecker(2), Y=exp(-((x[0]-0.5)**2 + (x[1]-0.5)**2)/0.1))
   solution = pde.getSolution()

   # Compute result
   result = integrate(solution)
   elapsed = time.time() - start_time

   # Leader of each group reports
   if grid_rank == 0:
       procs = grid_comm.Get_size()
       print(f"Grid {nx}x{ny} ({procs} proc): "
             f"time={elapsed:.4f}s, result={result:.6f}")

   grid_comm.Free()


Master-Slave Pattern Example
============================

The master-slave pattern is useful when you have a collection of independent tasks to distribute
among worker processes. One process (the master) manages task distribution and result collection,
while other processes (slaves) compute results.

This example shows how to use a master-slave pattern to solve magnetotelluric (MT) problems
for multiple frequencies:

.. code-block:: python

   from mpi4py import MPI
   from esys.escript import *
   from esys.finley import Rectangle
   from esys.escript.linearPDEs import LinearSinglePDE
   import numpy as np

   comm_world = MPI.COMM_WORLD
   rank = comm_world.Get_rank()
   size = comm_world.Get_size()

   def solve_MT_for_frequency(freq, domain):
       """Solve MT problem using provided domain"""
       # Setup PDE (complex-valued for frequency domain)
       pde = LinearSinglePDE(domain, isComplex=True)
       Mu0 = 4*np.pi*1e-7
       pde.setValue(D=1j*2*np.pi*freq*Mu0)

       # Set resistivity and boundary conditions
       pde.setValue(A=100*kronecker(2))  # 100 Ohm-m
       x = domain.getX()
       pde.setValue(q=whereZero(x[1]-20000), r=1)

       # Solve and return impedance
       Hx = pde.getSolution()
       return integrate(abs(Hx))

   # Create domain once on COMM_SELF (each process has its own)
   if rank != 0:
       domain = Rectangle(n0=200, n1=100, l0=40000, l1=20000,
                         comm=MPI.COMM_SELF)

   if rank == 0:
       # MASTER PROCESS
       frequencies = np.logspace(-1, 1, 10)  # 0.1 to 10 Hz
       results = []
       tasks_sent = 0
       tasks_received = 0

       # Send initial tasks to slaves
       for slave in range(1, min(size, len(frequencies)+1)):
           if tasks_sent < len(frequencies):
               comm_world.send(frequencies[tasks_sent], dest=slave, tag=1)
               tasks_sent += 1

       # Receive results and send more tasks
       while tasks_received < len(frequencies):
           status = MPI.Status()
           result = comm_world.recv(source=MPI.ANY_SOURCE, tag=2,
                                   status=status)
           slave = status.Get_source()
           results.append(result)
           tasks_received += 1

           # Send next task if available
           if tasks_sent < len(frequencies):
               comm_world.send(frequencies[tasks_sent], dest=slave, tag=1)
               tasks_sent += 1

       # Send termination signal
       for slave in range(1, size):
           comm_world.send(None, dest=slave, tag=1)

       print(f"Completed {len(results)} frequencies")

   else:
       # SLAVE PROCESS
       while True:
           freq = comm_world.recv(source=0, tag=1)
           if freq is None:
               break

           result = solve_MT_for_frequency(freq, domain)
           comm_world.send((freq, result), dest=0, tag=2)

   comm_world.Barrier()

**Key features of this pattern:**

* **Dynamic load balancing**: Master assigns tasks as slaves become available
* **Independent domains**: Each slave uses ``MPI.COMM_SELF`` for its own domain
* **Asynchronous communication**: ``MPI.ANY_SOURCE`` allows receiving from any slave
* **Scalability**: Works with any number of slave processes

**Running this example:**

.. code-block:: bash

   run-escript -n 11 mt_master_slave.py

This runs with 1 master and 10 slaves, solving 10 frequencies efficiently.

A complete working example (``MT_master_slave.py``) is provided in the
``doc/examples/usersguide/`` directory.


Multi-Physics Example: Thermo-Mechanical Coupling
=================================================

This example demonstrates time-dependent coupling of thermal and mechanical physics using
DataCoupler for inter-domain communication:

.. code-block:: python

   from mpi4py import MPI
   from esys.escript import *
   from esys.escript.linearPDEs import LinearPDE
   from esys.ripley import Rectangle
   from esys.escript import MPIDomainArray, DataCoupler

   world_comm = MPI.COMM_WORLD
   world_rank = world_comm.Get_rank()
   world_size = world_comm.Get_size()

   # Create domain array: 2 domains (thermal and mechanical)
   num_domains = 2
   domain_array = MPIDomainArray(numDomains=num_domains, comm=world_comm)

   my_domain_idx = domain_array.getDomainIndex()
   domain_comm = domain_array.getDomainComm()
   is_thermal = (my_domain_idx == 0)

   # Create DataCoupler for inter-domain communication
   coupler = DataCoupler(domain_array)

   # Create domain on domain-specific communicator
   domain = Rectangle(n0=80, n1=80, l0=1.0, l1=1.0, comm=domain_comm)
   x = domain.getX()

   # Time stepping parameters
   dt = 0.00025
   n_steps = 200

   if is_thermal:
       # THERMAL SOLVER
       pde = LinearPDE(domain)
       kappa = 1.0
       pde.setValue(A=kappa*kronecker(2), D=1/dt)
       T = 20.0 + 100.0*exp(-50*((x[0]-0.5)**2 + (x[1]-0.5)**2))
       T_old = T

       for step in range(n_steps):
           pde.setValue(Y=T_old/dt)
           T = pde.getSolution()
           coupler.send(interpolate(T, Function(domain)),
                       dest_domain_index=1, tag=step)
           T_old = T
   else:
       # MECHANICAL SOLVER
       pde = LinearPDE(domain, numEquations=2)
       # ... setup mechanical PDE ...

       for step in range(n_steps):
           T_data = coupler.receive(Function(domain),
                                   source_domain_index=0, tag=step)
           # ... solve with thermal loading ...

   world_comm.Barrier()

**Key features:**

* **DataCoupler**: High-level API for inter-domain communication
* **MPIDomainArray**: Manages communicator topology
* **Parallel physics**: Both solvers run simultaneously


Best Practices
==============

Communicator Management
-----------------------

Always free custom communicators when you're done with them:

.. code-block:: python

   sub_comm = MPI.COMM_WORLD.Split(color, rank)
   # ... use sub_comm ...
   sub_comm.Free()  # Important: free when done

Do **not** free ``MPI.COMM_WORLD`` or communicators returned by ``getMPIComm()``,
as these are managed by the system.

Synchronization
---------------

Use barriers when you need to synchronize across groups:

.. code-block:: python

   # Synchronize within a group
   group_comm.Barrier()

   # Synchronize all processes
   MPI.COMM_WORLD.Barrier()

Error Handling
--------------

Exceptions in one process group won't automatically propagate to others.
Consider using MPI collective operations to detect and handle errors:

.. code-block:: python

   try:
       solution = pde.getSolution()
       error_code = 0
   except Exception as e:
       print(f"Rank {rank} error: {e}")
       error_code = 1

   # Check if any process had an error
   global_error = group_comm.allreduce(error_code, op=MPI.MAX)
   if global_error != 0:
       print(f"Group {group_id} had errors")

Load Balancing
--------------

When splitting into groups, ensure each group has sufficient resources:

.. code-block:: python

   # Good: evenly divides processes
   num_groups = 4
   group_id = rank // (world_size // num_groups)

   # Bad: might create groups with different sizes
   group_id = rank % num_groups  # Only works if size divisible by num_groups


Performance Considerations
==========================

Process Group Size
------------------

Each domain needs enough processes to solve efficiently:

* Very small groups (1-2 processes) may not parallelize well
* Very large groups waste resources on small problems
* Experiment to find the optimal group size for your problem

Communication Overhead
----------------------

Minimize communication between process groups:

* Use group-local operations when possible
* Reduce data before sending between groups
* Consider asynchronous communication (``isend/irecv``)


Summary
=======

Using mpi4py with escript provides:

* Direct control over parallel resource allocation
* Standard MPI programming techniques
* Flexibility for complex parallel patterns
* Interoperability with other MPI libraries
* Fine-grained performance tuning

Key capabilities:

* Pass custom communicators to domains via ``comm`` parameter
* Retrieve communicators via ``getMPIComm()``
* Run ensemble simulations and parameter sweeps
* Implement multi-physics coupling
* Use standard MPI collective operations

For more information:

* mpi4py documentation: https://mpi4py.readthedocs.io/
* MPI standard: https://www.mpi-forum.org/
* escript examples directory: ``doc/examples/``
