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
Master-Slave MPI Pattern for Magnetotelluric (MT) Problem

This example demonstrates how to use mpi4py with escript to implement a
master-slave workflow. The master process distributes frequency values to
slave processes, which each solve an MT problem independently. Results are
collected by the master and plotted together.

Usage:
    mpirun -n 11 python3 MT_master_slave.py

This runs with 1 master + 10 slaves for 10 different frequencies.
"""

__copyright__ = """Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

from mpi4py import MPI
from esys.escript import *
from esys.finley import Rectangle
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearSinglePDE
import numpy as np
import matplotlib.pyplot as plt

# MPI Setup
comm_world = MPI.COMM_WORLD
rank = comm_world.Get_rank()
size = comm_world.Get_size()

# Constants for MT problem
L0 = 80000   # horizontal extent [m]
L1 = 40000   # depth [m]
NE0 = 400    # number of elements in horizontal direction
NE1 = 200    # number of elements in vertical direction

rho_b = 100  # background resistivity [Ohm m]
rho_a = 0.5  # resistivity in the anomaly [Ohm m]
D = 250.     # depth of the top edge of the anomaly [m]
W = 1000.    # width of the anomaly [m]
H = 2000.    # vertical extent of the anomaly [m]
OFFSET = 1000  # offset of the central vertical axis of the anomaly

Mu0 = 4 * np.pi * 1e-7
h0, h1 = L0/NE0, L1/NE1


def solve_MT_problem(frequency):
    """
    Solve the MT problem for a given frequency on MPI.COMM_SELF.

    This function creates its own escript domain using MPI.COMM_SELF,
    ensuring that each slave process runs independently.

    Parameters
    ----------
    frequency : float
        Frequency in Hz

    Returns
    -------
    tuple
        (x_positions, apparent_resistivity, phase)
    """
    # Create domain using COMM_SELF (single process domain)
    domain = Rectangle(n0=NE0, n1=NE1, l0=L0, l1=L1, comm=MPI.COMM_SELF)

    # Set up PDE (complex-valued)
    pde = LinearSinglePDE(domain, isComplex=True)
    pde.setValue(D=1j*2*np.pi*frequency*Mu0)

    # Define anomaly mask
    X = ReducedFunction(domain).getX()
    m1 = whereNonPositive(X[1] - (L1-D))
    m2 = whereNonNegative(X[1] - (L1-D-H))
    m3 = whereNonNegative(X[0] - (L0/2+OFFSET-W/2))
    m4 = whereNonPositive(X[0] - (L0/2+OFFSET+W/2))
    m = m1 * m2 * m3 * m4

    # Set resistivity distribution
    rho = rho_b * (1-m) + rho_a * m
    rho = interpolate(rho, Function(domain))
    pde.setValue(A=rho*np.eye(2))

    # Set boundary condition: magnetic field = 1 at top
    x = domain.getX()
    mD = whereZero(x[1] - L1)
    pde.setValue(q=mD, r=1)

    # Solve for magnetic field Hx
    Hx = pde.getSolution()

    # Calculate electric field and impedance
    Ey = rho * grad(Hx, ReducedFunction(domain))[1]
    Zyx = Ey / Hx

    # Extract transect along top of domain
    locations_in_transect = [(h0*k + h0/2., L1) for k in range(0, NE0, 2)]
    locator_transect = Locator(where=ReducedFunction(domain),
                               x=locations_in_transect)

    x_ts = np.array([x[0] for x in locator_transect.getX()])
    Zyx_ts = np.array(locator_transect(Zyx))

    # Calculate apparent resistivity and phase
    rho_a_ts = 1./(2*np.pi*frequency*Mu0) * abs(Zyx_ts)**2
    phi_ts = np.angle(Zyx_ts, deg=True)

    return x_ts, rho_a_ts, phi_ts


def master_process():
    """
    Master process: distribute frequencies to slaves and collect results.
    """
    # Define frequencies to compute (logarithmic spacing)
    frequencies = np.logspace(-1, 1, 10)  # 0.1 Hz to 10 Hz
    n_frequencies = len(frequencies)

    # Check we have enough processes
    n_slaves = size - 1
    if n_slaves < 1:
        print("Error: Need at least 2 MPI processes (1 master + 1 slave)")
        print(f"Running with {size} process(es)")
        comm_world.Abort(1)

    print(f"Master: Running MT problem for {n_frequencies} frequencies")
    print(f"Master: Using {n_slaves} slave process(es)")
    print(f"Frequencies: {frequencies}")

    # Storage for results
    results = []

    # Send frequencies to slaves (round-robin distribution)
    tasks_sent = 0
    tasks_received = 0

    # Initially send one task to each slave
    for slave_rank in range(1, min(n_slaves + 1, n_frequencies + 1)):
        if tasks_sent < n_frequencies:
            freq = frequencies[tasks_sent]
            print(f"Master: Sending f={freq:.2f} Hz to slave {slave_rank}")
            comm_world.send(freq, dest=slave_rank, tag=1)
            tasks_sent += 1

    # Receive results and send new tasks
    while tasks_received < n_frequencies:
        # Receive result from any slave
        status = MPI.Status()
        result = comm_world.recv(source=MPI.ANY_SOURCE, tag=2, status=status)
        slave_rank = status.Get_source()
        tasks_received += 1

        freq, x_ts, rho_a_ts, phi_ts = result
        results.append((freq, x_ts, rho_a_ts, phi_ts))
        print(f"Master: Received result for f={freq:.2f} Hz from slave {slave_rank} "
              f"({tasks_received}/{n_frequencies})")

        # Send next task to this slave if available
        if tasks_sent < n_frequencies:
            freq = frequencies[tasks_sent]
            print(f"Master: Sending f={freq:.2f} Hz to slave {slave_rank}")
            comm_world.send(freq, dest=slave_rank, tag=1)
            tasks_sent += 1

    # Send termination signal to all slaves
    for slave_rank in range(1, size):
        comm_world.send(None, dest=slave_rank, tag=1)

    print("Master: All results collected")

    # Sort results by frequency
    results.sort(key=lambda x: x[0])

    # Plot all transects together
    plot_results(results)


def slave_process():
    """
    Slave process: receive frequencies, solve MT problem, send results back.
    """
    while True:
        # Receive frequency from master
        freq = comm_world.recv(source=0, tag=1)

        # Check for termination signal
        if freq is None:
            print(f"Slave {rank}: Received termination signal")
            break

        print(f"Slave {rank}: Solving MT problem for f={freq:.2f} Hz")

        # Solve MT problem
        x_ts, rho_a_ts, phi_ts = solve_MT_problem(freq)

        print(f"Slave {rank}: Finished f={freq:.2f} Hz, sending results")

        # Send results back to master
        result = (freq, x_ts, rho_a_ts, phi_ts)
        comm_world.send(result, dest=0, tag=2)


def plot_results(results):
    """
    Plot apparent resistivity and phase for all frequencies.

    Parameters
    ----------
    results : list of tuples
        Each tuple contains (freq, x_ts, rho_a_ts, phi_ts)
    """
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Color map for different frequencies
    colors = plt.cm.viridis(np.linspace(0, 1, len(results)))

    # Plot apparent resistivity
    for (freq, x_ts, rho_a_ts, phi_ts), color in zip(results, colors):
        ax1.plot(x_ts/1000, rho_a_ts, label=f'{freq:.2f} Hz',
                color=color, linewidth=2)

    ax1.set_xlabel('Distance [km]', fontsize=12)
    ax1.set_ylabel('Apparent Resistivity [Ω·m]', fontsize=12)
    ax1.set_title('MT Apparent Resistivity Transects', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10, ncol=2)
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')

    # Plot phase
    for (freq, x_ts, rho_a_ts, phi_ts), color in zip(results, colors):
        ax2.plot(x_ts/1000, phi_ts, label=f'{freq:.2f} Hz',
                color=color, linewidth=2)

    ax2.set_xlabel('Distance [km]', fontsize=12)
    ax2.set_ylabel('Phase [degrees]', fontsize=12)
    ax2.set_title('MT Phase Transects', fontsize=14, fontweight='bold')
    ax2.legend(loc='best', fontsize=10, ncol=2)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('MT_master_slave_results.png', dpi=150)
    print("Master: Saved plot to MT_master_slave_results.png")

    # Also create a separate plot showing how anomaly signature varies with frequency
    fig2, ax = plt.subplots(figsize=(12, 6))

    # Find the approximate location of the anomaly center
    anomaly_center_km = (L0/2 + OFFSET) / 1000

    for (freq, x_ts, rho_a_ts, phi_ts), color in zip(results, colors):
        # Normalize by background value for better comparison
        rho_normalized = rho_a_ts / rho_b
        ax.plot(x_ts/1000, rho_normalized, label=f'{freq:.2f} Hz',
               color=color, linewidth=2)

    ax.axvline(x=anomaly_center_km, color='red', linestyle='--',
              linewidth=2, alpha=0.5, label='Anomaly center')
    ax.set_xlabel('Distance [km]', fontsize=12)
    ax.set_ylabel('Normalized Resistivity (ρ_a / ρ_background)', fontsize=12)
    ax.set_title('Frequency-Dependent Anomaly Response', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('MT_master_slave_normalized.png', dpi=150)
    print("Master: Saved normalized plot to MT_master_slave_normalized.png")


# Main execution
if __name__ == '__main__':
    if rank == 0:
        # Master process
        master_process()
    else:
        # Slave process
        slave_process()

    # Synchronize all processes before exit
    comm_world.Barrier()

    if rank == 0:
        print("Master: All processes finished successfully")
