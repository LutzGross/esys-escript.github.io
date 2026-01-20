##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
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
physics using DataCoupler for communication between domains. The thermal
solver computes temperature evolution while the mechanical solver calculates
thermal stresses. The mechanical model runs one time-step behind to use the
temperature from the previous step (staggered coupling).

Usage:
    # Standard mode (only rank 0 output visible):
    ./bin/run-escript -n 8 thermo_mechanical.py

    # Verbose mode (see output from all domains):
    ./bin/run-escript -n 8 -a thermo_mechanical.py

This uses 4 processes for thermal domain and 4 for mechanical domain.
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

from mpi4py import MPI
from esys.escript import *
from esys.escript import MPIDomainArray, DataCoupler
from esys.escript.linearPDEs import LinearPDE
from esys.escript.util import identityTensor4
from esys.ripley import Rectangle
from esys.weipa import saveVTK

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()

# Create output directory
mkDir("output")

# Create domain array: 2 domains (thermal and mechanical)
# Domain 0: thermal solver
# Domain 1: mechanical solver
num_domains = 2
domain_array = MPIDomainArray(numDomains=num_domains, comm=world_comm)

my_domain_idx = domain_array.getDomainIndex()
domain_comm = domain_array.getDomainComm()
subdomain_comm = domain_array.getSubdomainComm()
domain_rank = domain_comm.Get_rank()

is_thermal = (my_domain_idx == 0)

# Print initialization from each domain's rank 0
if domain_rank == 0:
    print(f"Domain {my_domain_idx} ({'thermal' if is_thermal else 'mechanical'}) "
          f"initialized with {domain_comm.Get_size()} processes", flush=True)

# Create DataCoupler for inter-domain communication
coupler = DataCoupler(domain_array)

# Create domain on domain-specific communicator
domain = Rectangle(n0=80, n1=80, l0=1.0, l1=1.0, comm=domain_comm)
x = domain.getX()

# Time stepping parameters
dt = 0.00025      # Time step
n_steps = 200  # Number of time steps

if is_thermal:
    # THERMAL SOLVER: Transient heat equation
    # du/dt = kappa * nabla^2(u) + Q

    pde = LinearPDE(domain)
    kappa = 1.0  # Thermal diffusivity
    pde.setValue(A=kappa*kronecker(2), D=1/dt)

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

        # Send temperature field to mechanical solver using DataCoupler
        # Interpolate to Function space for consistent transfer
        T_data = interpolate(T, Function(domain))
        coupler.send(T_data, dest_domain_index=1, tag=step)

        T_max = Lsup(T)
        if domain_rank == 0 and step % 10 == 0:
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
    # C_ijkl = λ δ_ij δ_kl + μ (δ_ik δ_jl + δ_il δ_jk)
    lam = E*nu/(1-nu**2)       # Lamé's first parameter (plane stress)
    mu = E/(2*(1+nu))           # Shear modulus

    # Construct rank-4 elastic stiffness tensor
    d = 2  # 2D
    C = lam * kronecker(d).reshape(d,1,d,1) * kronecker(d).reshape(1,d,1,d) + \
        mu * (identityTensor4(d) + identityTensor4(d).transpose(0,1,3,2))

    pde.setValue(A=C)
    # Fixed displacement on left edge
    x = domain.getX()
    pde.setValue(q=whereZero(x[0])*[1,1], r=[0.0, 0.0])
    for step in range(0, n_steps):
        # Receive temperature from thermal solver (previous step) using DataCoupler
        # DataCoupler handles the receive and distribution across all ranks
        T_data = coupler.receive(Function(domain), source_domain_index=0, tag=step)
        # Thermal stress (prestress due to thermal expansion)
        Delta_T = T_data - T0
        sigma_thermal = (lam + 2./3.*mu) * alpha * Delta_T * kronecker(2)

        # Set thermal stress as prestress (X coefficient)
        pde.setValue(X=sigma_thermal)

        # Solve for displacement
        u = pde.getSolution()
        u_max = Lsup(length(u))
        if domain_rank == 0 and step % 10 == 0:
            print(f"Mechanical step {step}: u_max = {u_max:.2e} m")

        # Save displacement and temperature fields
        if step % 10 == 0:
            saveVTK(f"output/thermal_{step//10:04d}.vtu", u=u, T=T_data)

world_comm.Barrier()

if world_rank == 0:
    print("Thermo-mechanical coupling completed successfully")
