#!/usr/bin/env python3
"""
Demonstrate the MPI initialization order and its impact
"""
import sys
import ctypes

def check_mpi_state(label):
    """Check if MPI is initialized at this point"""
    try:
        mpi_lib = ctypes.CDLL("libmpi.so", mode=ctypes.RTLD_GLOBAL)
        initialized = ctypes.c_int(0)
        mpi_lib.MPI_Initialized(ctypes.byref(initialized))

        status = "INITIALIZED" if initialized.value else "NOT initialized"
        sys.stderr.write(f"[{label}] MPI state: {status}\n")
        sys.stderr.flush()
        return bool(initialized.value)
    except Exception as e:
        sys.stderr.write(f"[{label}] Could not check MPI state: {e}\n")
        sys.stderr.flush()
        return None

# Check 1: Very first thing - before any imports
print("\n" + "="*70)
print("DEMONSTRATION: MPI Initialization Order")
print("="*70 + "\n")

was_init_before_imports = check_mpi_state("STEP 1: Before any imports")

# Check 2: After importing mpi4py
print("\nImporting mpi4py...")
from mpi4py import MPI

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

check_mpi_state(f"STEP 2: After mpi4py import (rank {rank}/{size})")

# Check 3: Before importing escript
sys.stderr.write(f"\n[STEP 3] Rank {rank}: About to import escript/Trilinos...\n")
sys.stderr.flush()

from esys.escript import *
from esys.finley import Rectangle

check_mpi_state(f"STEP 4: After escript import (rank {rank})")

# Now show what this means
print("\n" + "-"*70)
if was_init_before_imports:
    print("RESULT: MPI was initialized BEFORE Python started")
    print("  ✓ This is CORRECT - using pythonMPI wrapper")
    print("  ✓ All C++ libraries initialized with MPI context")
    print("  ✓ Trilinos/Tpetra will work correctly")
else:
    print("RESULT: MPI was initialized AFTER Python started")
    print("  ✗ This is WRONG - using direct mpirun python3")
    print("  ✗ Some C++ static objects created without MPI")
    print("  ✗ Trilinos/Tpetra may fail with index errors")

print("-"*70)

# Try a simple operation to see if it works
sys.stderr.write(f"\n[STEP 5] Rank {rank}: Testing domain creation...\n")
sys.stderr.flush()

try:
    domain = Rectangle(l0=1., l1=1., n0=10, n1=10)
    sys.stderr.write(f"[STEP 5] Rank {rank}: Domain created OK\n")
    sys.stderr.flush()

    # Try interpolation (this is where the bug appears)
    x = domain.getX()
    x_sol = interpolate(x, where=Solution(domain))
    x_cont = interpolate(x_sol, where=ContinuousFunction(domain))

    sys.stderr.write(f"[STEP 5] Rank {rank}: Interpolation OK\n")
    sys.stderr.flush()

    if rank == 0:
        print("\n✓ SUCCESS: Operations completed without errors")

except Exception as e:
    sys.stderr.write(f"[STEP 5] Rank {rank}: FAILED: {e}\n")
    sys.stderr.flush()

    if rank == 0:
        print(f"\n✗ FAILED: {e}")
        print("\nThis failure is likely due to incorrect MPI initialization order")

MPI.COMM_WORLD.Barrier()
