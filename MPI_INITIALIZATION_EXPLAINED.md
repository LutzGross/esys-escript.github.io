# MPI Initialization Order: Detailed Explanation

## The Core Issue

When using escript with MPI, **the order in which MPI is initialized matters** because Trilinos (the linear algebra library) needs MPI to be fully set up before any C++ code in escript/Trilinos loads.

## Two Different Execution Paths

### ✓ CORRECT: `./bin/run-escript -n 2 script.py`

**Execution Timeline:**

```
Time →
─────────────────────────────────────────────────────────────────────

1. mpirun launches pythonMPI executable
   ├─ pythonMPI is a C++ program that embeds Python
   └─ Source: pythonMPI/src/ScriptMPI.cpp

2. pythonMPI calls MPI_Init_thread()
   ├─ MPI runtime fully initialized
   ├─ MPI_COMM_WORLD created and ready
   └─ All ranks know their rank/size

3. pythonMPI calls Py_Main(argc, argv)
   ├─ Python interpreter starts
   ├─ Python starts INSIDE an initialized MPI context
   └─ Environment variables (PYTHONPATH, LD_LIBRARY_PATH) set by run-escript

4. Python imports mpi4py
   ├─ mpi4py checks: is MPI already initialized?
   ├─ Answer: YES! (from step 2)
   ├─ mpi4py just wraps existing MPI_COMM_WORLD
   └─ No second initialization

5. Python imports escript → libescript.so loads
   ├─ C++ shared library loads with dlopen()
   ├─ C++ static constructors run
   ├─ Trilinos libraries load
   └─ ALL C++ code sees MPI is already initialized

6. Create domain: Rectangle(n0=40, n1=20)
   ├─ Calls into finley C++ code
   ├─ Trilinos creates Tpetra::Map for 861 nodes (41×21)
   ├─ Rank 0 gets indices: 0-430 (431 nodes)
   ├─ Rank 1 gets indices: 431-860 (430 nodes)
   └─ Total: 861 nodes, indices 0-860 ✓ CORRECT

7. Interpolate between function spaces
   ├─ Creates coupling matrix between Solution and ContinuousFunction
   ├─ Tpetra builds distributed sparse matrix
   ├─ All global indices 0-860 exist in Maps
   └─ ✓ SUCCESS
```

**Key Point:** MPI is initialized in step 2, **BEFORE** Python and all C++ libraries load.

---

### ✗ WRONG: `mpirun -n 2 python3 script.py`

**Execution Timeline:**

```
Time →
─────────────────────────────────────────────────────────────────────

1. mpirun launches python3 executable
   ├─ python3 is the standard Python interpreter
   └─ NO MPI initialization yet

2. Python interpreter starts
   ├─ Python starts WITHOUT MPI
   ├─ Python's C extensions loading mechanism initialized
   └─ MPI_COMM_WORLD does not exist yet

3. Python imports mpi4py
   ├─ mpi4py checks: is MPI already initialized?
   ├─ Answer: NO
   ├─ mpi4py calls MPI_Init()
   └─ MPI initialized NOW (after Python started)

4. Python imports escript → libescript.so loads
   ├─ C++ shared library loads
   ├─ C++ static constructors run
   ├─ PROBLEM: Some Trilinos static objects may have been
   │   created BEFORE mpi4py called MPI_Init()
   └─ Inconsistent MPI state across different C++ objects

5. Create domain: Rectangle(n0=40, n1=20)
   ├─ Calls into finley C++ code
   ├─ Trilinos creates Tpetra::Map for nodes
   ├─ PROBLEM: Index calculation uses mismatched MPI state
   ├─ Rank 0: calculates some indices
   ├─ Rank 1: calculates other indices
   └─ Total indices don't match! Some off-by-one errors

6. Interpolate between function spaces
   ├─ Creates coupling matrix
   ├─ Code references global index 862
   ├─ But Maps only have indices 0-860!
   └─ ✗ CRASH: "global index 862 does not exist"
```

**Key Point:** MPI is initialized in step 3, **AFTER** Python started. Some C++ static initialization happened without MPI being ready.

## Why Does This Matter?

### C++ Static Initialization Order

C++ allows global objects that are initialized when a shared library loads:

```cpp
// In some Trilinos header/cpp file
static SomeTrilinosObject globalThing;  // Initialized when library loads!
```

**With pythonMPI:**
- Library loads after MPI_Init → globalThing sees initialized MPI ✓

**With direct mpirun python3:**
- Library might load before MPI_Init → globalThing sees uninitialized MPI ✗
- Later code expects initialized MPI
- Mismatch leads to bugs

### The "Index 862" Bug Explained

For a 40×20 Rectangle mesh:
- Nodes: (40+1) × (20+1) = 41 × 21 = **861 nodes**
- Valid indices: **0 through 860**

With 2 MPI ranks, should divide as:
- Rank 0: indices 0-430 (431 nodes)
- Rank 1: indices 431-860 (430 nodes)

But with wrong initialization order:
- Part of code thinks there are 862 nodes (0-861)
- Tries to access index 862
- Index doesn't exist → crash

This happens because:
1. Domain decomposition code read uninitialized MPI state (before MPI_Init)
2. Matrix assembly code read initialized MPI state (after MPI_Init)
3. They disagree on how many indices exist
4. Off-by-one error: thinks max index is 861 instead of 860

## The pythonMPI Wrapper

The `pythonMPI` executable (source: `pythonMPI/src/ScriptMPI.cpp`) is specifically designed to solve this problem:

```cpp
int main(int argc, char **argv)
{
    // STEP 1: Initialize MPI FIRST
    int provided;
    status = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (status != MPI_SUCCESS) {
        std::cerr << "MPI_Init failed" << std::endl;
        return status;
    }

    // STEP 2: NOW start Python
    status = Py_Main(argc, wargv);

    // STEP 3: Clean up MPI when Python exits
    if (status != 0 && status != 2) {
        MPI_Abort(MPI_COMM_WORLD, status);
    } else {
        MPI_Finalize();
    }

    return status;
}
```

This ensures:
- MPI initialized before Python starts
- All Python imports happen with MPI active
- All C++ static initialization sees MPI
- Clean MPI shutdown when Python exits

## How run-escript Calls pythonMPI

The `run-escript` script (source: `bin/run-escript`) does:

```bash
# Set up environment
export PYTHONPATH=$ESCRIPT_ROOT:$PYTHONPATH
export LD_LIBRARY_PATH=$ESCRIPT_ROOT/lib/esys:$LD_LIBRARY_PATH

# Use pythonMPI instead of python3 for MPI jobs
if [ "$MPI_FLAVOUR" != "none" ]; then
    PYTHON_MPI="$ESCRIPT_ROOT/lib/esys/pythonMPI"
    EXEC_CMD="$PYTHON_MPI $@"
else
    EXEC_CMD="python3 $@"
fi

# Launch with mpirun
mpirun -n $NUM_PROCS $EXEC_CMD
```

So `./bin/run-escript -n 2 script.py` actually runs:
```bash
mpirun -n 2 /path/to/pythonMPI script.py
```

NOT:
```bash
mpirun -n 2 python3 script.py
```

## Working with Sub-Communicators

**Good news:** You CAN use mpi4py sub-communicators, as long as you use `run-escript`:

```python
from mpi4py import MPI
from esys.escript import *
from esys.finley import Rectangle

# Split MPI_COMM_WORLD
world = MPI.COMM_WORLD
color = world.Get_rank() % 2
sub_comm = world.Split(color, world.Get_rank())

# Pass sub-communicator to domain
domain = Rectangle(l0=1., l1=1., n0=20, n1=10, comm=sub_comm)

# Each subgroup runs independently!
```

This works because:
- MPI was initialized before Python (by pythonMPI)
- mpi4py can create sub-communicators from MPI_COMM_WORLD
- escript domains accept sub-communicators via `comm=` parameter
- Trilinos works correctly because MPI was initialized first

## Summary

| Method | MPI Init Order | Works? | Reason |
|--------|---------------|--------|---------|
| `mpirun python3 script.py` | After Python starts | ❌ | Trilinos sees inconsistent MPI state |
| `./bin/run-escript script.py` | Before Python starts | ✓ | All C++ code sees initialized MPI |
| `run-escript` + mpi4py subcomms | Before Python starts | ✓ | Sub-comms created after proper init |

**Always use `run-escript` for MPI programs with escript!**

## Analogy

Think of MPI initialization like turning on the electricity in a house:

**pythonMPI (correct):**
1. Turn on electricity (MPI_Init)
2. Move in furniture (Python starts)
3. Plug in appliances (import escript)
✓ All appliances work because electricity was on first

**Direct mpirun python3 (wrong):**
1. Move in furniture (Python starts)
2. Plug in some appliances (import libraries)
3. Turn on electricity (mpi4py MPI_Init)
✗ Some appliances were plugged in before power was on, they're confused!

The appliances (C++ libraries) that were plugged in before electricity (MPI) was on captured the "no power" state and don't work correctly even after power is turned on.
