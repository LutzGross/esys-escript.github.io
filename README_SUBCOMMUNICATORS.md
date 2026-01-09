# Using Sub-Communicators with escript

This document explains how to use mpi4py sub-communicators with escript domains, replacing the old SplitWorld functionality.

## Quick Start

```python
from mpi4py import MPI
from esys.escript import *
from esys.finley import Rectangle

# Split MPI_COMM_WORLD into subgroups
world = MPI.COMM_WORLD
color = world.Get_rank() % 2  # Create 2 groups
sub_comm = world.Split(color, world.Get_rank())

# Create domain on sub-communicator
domain = Rectangle(l0=1., l1=1., n0=20, n1=10, comm=sub_comm)

# Each subgroup runs independently!
```

**IMPORTANT**: Always use `run-escript` to launch your scripts:
```bash
./bin/run-escript -n 4 -t 1 myscript.py
```

**DO NOT** use `mpirun python3` directly - this causes MPI initialization order issues.

## Why MPI Initialization Order Matters

escript uses Trilinos for distributed linear algebra. Trilinos requires MPI to be initialized **before** Python starts to ensure proper C++ static initialization.

- **CORRECT**: `./bin/run-escript -n 4 script.py` → MPI initializes before Python
- **WRONG**: `mpirun -n 4 python3 script.py` → MPI initializes after Python starts → Trilinos errors

See `MPI_INITIALIZATION_EXPLAINED.md` for technical details.

## Example Scripts

### test_subcommunicators.py
Basic test showing sub-communicator usage:
```bash
./bin/run-escript -n 4 -t 1 test_subcommunicators.py
```

### test_independent_subgroups.py
Demonstrates running independent models on different subgroups simultaneously:
```bash
./bin/run-escript -n 4 -t 1 test_independent_subgroups.py
```

### demonstrate_init_order.py
Shows MPI initialization order and its effects:
```bash
./bin/run-escript -n 2 -t 1 demonstrate_init_order.py
```

## Common Patterns

### Ensemble Simulations

Run multiple parameter values in parallel:

```python
from mpi4py import MPI
from esys.escript import *
from esys.ripley import Rectangle

world = MPI.COMM_WORLD
rank = world.Get_rank()

# Define parameters to explore
parameters = [0.1, 0.5, 1.0, 2.0]
param_id = rank % len(parameters)
my_param = parameters[param_id]

# Create subgroup for this parameter
sub_comm = world.Split(param_id, rank)

# Solve with this parameter
domain = Rectangle(n0=50, n1=50, comm=sub_comm)
# ... solve PDE with my_param ...

sub_comm.Free()  # Clean up when done
```

### Master-Slave Pattern

One master distributes tasks to many slaves:

```python
from mpi4py import MPI
from esys.finley import Rectangle

world = MPI.COMM_WORLD
rank = world.Get_rank()

if rank == 0:
    # MASTER: distribute tasks
    for slave in range(1, world.Get_size()):
        task = get_next_task()
        world.send(task, dest=slave, tag=1)

    # Collect results...
else:
    # SLAVE: solve independently
    domain = Rectangle(n0=100, n1=100, comm=MPI.COMM_SELF)
    while True:
        task = world.recv(source=0, tag=1)
        if task is None:
            break
        result = solve_problem(domain, task)
        world.send(result, dest=0, tag=2)
```

### Multi-Physics Coupling

Different physics on different subgroups:

```python
from mpi4py import MPI
from esys.ripley import Rectangle

world = MPI.COMM_WORLD
rank = world.Get_rank()
size = world.Get_size()

# Split into thermal and mechanical groups
is_thermal = (rank < size // 2)
physics_id = 0 if is_thermal else 1
physics_comm = world.Split(physics_id, rank)

# Each physics creates its own domain
domain = Rectangle(n0=80, n1=80, comm=physics_comm)

if is_thermal:
    # Solve heat equation...
    T = solve_thermal(domain)
    # Send temperature to mechanical group
    if rank == 0:
        world.send(T_data, dest=size//2, tag=step)
else:
    # Receive temperature and solve mechanics
    if rank == size//2:
        T_data = world.recv(source=0, tag=step)
    T_data = physics_comm.bcast(T_data, root=0)
    u = solve_mechanical(domain, T_data)

physics_comm.Free()
```

## Best Practices

1. **Always use `run-escript`** for MPI programs
2. **Import mpi4py first** (before escript modules)
3. **Free custom communicators** when done: `sub_comm.Free()`
4. **Don't free `MPI.COMM_WORLD`** or communicators from `getMPIComm()`
5. **Use barriers** to synchronize between groups
6. **Handle errors** with MPI collective operations

## Documentation

- **User Guide**: `doc/user/mpi4py.tex` - Comprehensive examples and patterns
- **Technical Details**: `MPI_INITIALIZATION_EXPLAINED.md` - Why initialization order matters
- **Project Guide**: `CLAUDE.md` - See "Working with MPI" section

## Transitioning from SplitWorld

If you used the old `SplitWorld` API:

**Old (SplitWorld)**:
```python
from esys.escript import SplitWorld
# ... no longer available ...
```

**New (mpi4py)**:
```python
from mpi4py import MPI

# Split communicator
sub_comm = MPI.COMM_WORLD.Split(color, rank)

# Pass to domain
domain = Rectangle(n0=100, n1=100, comm=sub_comm)

# Clean up
sub_comm.Free()
```

## Troubleshooting

### "MPI has not been initialized" error
**Solution**: Import mpi4py before escript, or use `run-escript`

### "global index N does not exist" Trilinos error
**Solution**: Use `run-escript` instead of `mpirun python3`

### Hangs during interpolation with multiple ranks
**Solution**: Use `run-escript` - this is an MPI initialization order issue

### Sub-communicator not working
**Solution**: Ensure you're using `run-escript` launcher, not direct `mpirun python3`

## Further Reading

- mpi4py documentation: https://mpi4py.readthedocs.io/
- MPI standard: https://www.mpi-forum.org/
- escript examples: `doc/examples/usersguide/`
