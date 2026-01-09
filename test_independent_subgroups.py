#!/usr/bin/env python3
"""
Test that independent subgroups can run different models simultaneously
This is the use case that SplitWorld was designed for
"""

from mpi4py import MPI
from esys.escript import *
from esys.finley import Rectangle
import sys

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()

# Split into subgroups
color = world_rank % 2
sub_comm = world_comm.Split(color, world_rank)
sub_rank = sub_comm.Get_rank()

sys.stderr.write(f"Rank {world_rank}: Starting in subgroup {color}\n")
sys.stderr.flush()

# Each subgroup runs a DIFFERENT model with different parameters
if color == 0:
    # Subgroup 0: coarse mesh
    domain = Rectangle(l0=1., l1=1., n0=10, n1=10, comm=sub_comm)
    param_value = 1.0
    sys.stderr.write(f"Rank {world_rank}: Subgroup 0 running COARSE model (10x10)\n")
else:
    # Subgroup 1: fine mesh
    domain = Rectangle(l0=1., l1=1., n0=20, n1=20, comm=sub_comm)
    param_value = 2.0
    sys.stderr.write(f"Rank {world_rank}: Subgroup 1 running FINE model (20x20)\n")

sys.stderr.flush()

# Run independent calculations
x = domain.getX()
result = interpolate(x[0] * param_value, where=Solution(domain))

sys.stderr.write(f"Rank {world_rank}: Subgroup {color} computation complete\n")
sys.stderr.flush()

# Synchronize within each subgroup first
sub_comm.Barrier()
sys.stderr.write(f"Rank {world_rank}: Subgroup {color} barrier passed\n")
sys.stderr.flush()

# Then synchronize all ranks
world_comm.Barrier()

if world_rank == 0:
    print("\n" + "=" * 70)
    print("SUCCESS: Multiple independent models run simultaneously on subgroups!")
    print("This replaces the old SplitWorld functionality with mpi4py control.")
    print("=" * 70)