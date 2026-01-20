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
Simple Ensemble Simulation Example

This example demonstrates running multiple parameter values in parallel.
Each process group solves the same problem with a different parameter.

Usage:
    mpirun -n 8 run-escript ensemble_example.py

With 8 processes and 4 parameter values, each parameter gets 2 processes.
All 4 simulations run simultaneously.
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

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
pde.setValue(A=my_parameter*kronecker(2), Y=x[0]*x[1], q=whereZero(x[0]))
solution = pde.getSolution()

# Each group can now process its result independently
# Note: integrate() is collective, must be called by all processes in group
result = integrate(solution)
if not domain.isRootRank():
    result = None  # Only root rank keeps the result

# Clean up communicator
group_comm.Free()

# Gather all results to rank 0 for display
all_results = world_comm.gather((my_parameter, result), root=0)

if world_rank == 0:
    print("Ensemble simulation results:")
    for param, res in all_results:
        if res is not None:
            print(f"  Parameter {param}: Result = {res}")
    print("Ensemble simulation completed successfully")
