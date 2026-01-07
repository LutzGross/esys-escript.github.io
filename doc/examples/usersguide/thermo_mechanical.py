##############################################################################
#
# Copyright (c) 2003-2025 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
##############################################################################

"""
Multi-Physics Example: Thermo-Mechanical Coupling

This example demonstrates time-dependent coupling of thermal and mechanical
physics on separate communicators. The thermal solver computes temperature
evolution while the mechanical solver calculates thermal stresses. The
mechanical model runs one time-step behind to use the temperature from the
previous step (staggered coupling).

Usage:
    mpirun -n 8 run-escript thermo_mechanical.py

This uses 4 processes for thermal and 4 for mechanical.
"""

__copyright__ = """Copyright (c) 2003-2025 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

from mpi4py import MPI
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE
from esys.ripley import Rectangle
import numpy as np

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()

# Split processes: first half for thermal, second half for mechanical
is_thermal = (world_rank < world_size // 2)
physics_id = 0 if is_thermal else 1

# Create separate communicators for each physics
physics_comm = world_comm.Split(physics_id, world_rank)

# Create domain on physics-specific communicator
domain = Rectangle(n0=80, n1=80, l0=1.0, l1=1.0, comm=physics_comm)
x = domain.getX()

# Time stepping parameters
dt = 0.01      # Time step
n_steps = 50   # Number of time steps

if is_thermal:
    # THERMAL SOLVER: Transient heat equation
    # du/dt = kappa * nabla^2(u) + Q

    pde = LinearPDE(domain)
    kappa = 1.0  # Thermal diffusivity
    pde.setValue(A=kappa*kronecker(2), D=1.0)

    # Initial temperature: hot spot in center
    T = 20.0 + 100.0*exp(-50*((x[0]-0.5)**2 + (x[1]-0.5)**2))
    T_old = T

    for step in range(n_steps):
        # Heat source (decaying with time)
        Q = 10.0 * exp(-0.1*step*dt) * \
            exp(-50*((x[0]-0.5)**2 + (x[1]-0.5)**2))

        # Backward Euler: (T - T_old)/dt = kappa*nabla^2(T) + Q
        pde.setValue(Y=T_old/dt + Q)
        T = pde.getSolution()

        # Send temperature field to mechanical solver
        if domain.isRootRank():
            # Convert to numpy for efficient MPI transfer
            T_data = interpolate(T, Function(domain))
            T_array = convertToNumpy(T_data)
            world_comm.send(T_array, dest=world_size//2, tag=step)

            if step % 10 == 0:
                T_max = Lsup(T)
                print(f"Thermal step {step}: T_max = {T_max:.2f}")

        T_old = T

else:
    # MECHANICAL SOLVER: Thermal stress (one step behind)
    # Solve: -div(sigma) = 0, sigma = E*epsilon - alpha*E*Delta_T

    pde = LinearPDE(domain, numEquations=2)
    E = 70e9        # Young's modulus (Pa)
    nu = 0.3        # Poisson's ratio
    alpha = 23e-6   # Thermal expansion coefficient (1/K)
    T0 = 20.0       # Reference temperature

    # Elastic stiffness tensor (plane stress)
    C = E/(1-nu**2) * kronecker(2)
    pde.setValue(A=C)

    # Fixed displacement on left edge
    x = domain.getX()
    pde.setValue(q=whereZero(x[0]), r=[0.0, 0.0])

    for step in range(1, n_steps):  # Start at 1 (one step behind)
        # Receive temperature from thermal solver (previous step)
        if domain.isRootRank():
            T_array = world_comm.recv(source=0, tag=step-1)
            # Reconstruct escript Data from numpy array
            T_data = Data(T_array, Function(domain))
        else:
            T_data = None

        # Broadcast temperature within mechanical group
        T_data = physics_comm.bcast(T_data, root=0)

        # Thermal strain causes body force
        Delta_T = T_data - T0
        thermal_stress = alpha * E/(1-nu) * Delta_T

        # Thermal stress acts as body force in equilibrium equation
        pde.setValue(Y=[grad(thermal_stress)[0], grad(thermal_stress)[1]])

        # Solve for displacement
        u = pde.getSolution()

        if domain.isRootRank() and step % 10 == 0:
            u_max = Lsup(length(u))
            print(f"Mechanical step {step}: u_max = {u_max:.2e} m")

physics_comm.Free()
world_comm.Barrier()

if world_rank == 0:
    print("Thermo-mechanical coupling completed successfully")
