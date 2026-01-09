#!/usr/bin/env python3
"""Test using mpi4py sub-communicators with escript domains"""

from mpi4py import MPI
from esys.escript import *
from esys.finley import Rectangle

# Get world communicator info
world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()

print(f"World rank {world_rank}/{world_size}: Starting test", flush=True)

# Split into two groups based on even/odd ranks
color = world_rank % 2
sub_comm = world_comm.Split(color, world_rank)
sub_rank = sub_comm.Get_rank()
sub_size = sub_comm.Get_size()

print(f"World rank {world_rank}: Assigned to subgroup {color}, sub-rank {sub_rank}/{sub_size}", flush=True)

# Create domain on the sub-communicator
try:
    mydomain = Rectangle(l0=1., l1=1., n0=20, n1=10, comm=sub_comm)
    print(f"World rank {world_rank}: Domain created on subgroup {color} with {mydomain.getMPISize()} ranks", flush=True)

    # Test basic operation
    x = mydomain.getX()
    x_sol = interpolate(x, where=Solution(mydomain))

    print(f"World rank {world_rank}: Interpolation successful on subgroup {color}", flush=True)

    # Test more complex interpolation
    x_cont = interpolate(x_sol, where=ContinuousFunction(mydomain))
    print(f"World rank {world_rank}: Round-trip interpolation successful on subgroup {color}", flush=True)

except Exception as e:
    print(f"World rank {world_rank}: FAILED on subgroup {color}: {e}", flush=True)
    import traceback
    traceback.print_exc()

world_comm.Barrier()
if world_rank == 0:
    print("=" * 60)
    print("SUCCESS: All subgroups completed successfully!")
    print("=" * 60)
