##############################################################################
#
# Copyright (c) 2003-2022 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2022 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import *
from esys.finley import Rectangle
from esys.escript.pdetools import Locator
from esys.escript.linearPDEs import LinearSinglePDE
import numpy as np
import matplotlib.pyplot as plt

L0 = 80000   # horizontal extend [m]
L1 = 40000   # depth [m]
NE0 = 400    # number of element in horizontal direction
NE1 = 200    # number of element in vertical direction

rho_b = 100  # background resistivity [Ohm m]
rho_a = 0.5  # resistivity in the anomaly [Ohm m]
D = 250.     # depth of the top edge of the anomaly [m]
W = 1000.    # width of the anomaly [m]
H = 2000.    # vertical extend of the anomaly [m]
OFFSET = 1000 # offset of the central vertical axis of the anomaly from the domain (<0 move to the left, 0> move to the right)
f = 1.       # frequency [Hz]


Mu0=4*np.pi*1e-7
# ... mesh resolution:
h0, h1=L0/NE0, L1/NE1

# create Domain:
domain = Rectangle(n0=NE0, n1=NE1, l0=L0, l1=L1)

# set-up PDE. This one is complex:"
pde = LinearSinglePDE(domain, isComplex=True)
pde.setValue(D=1j*2*np.pi*f*Mu0)

# mask for anomaly boundaries:
X = ReducedFunction(domain).getX()
m1 = whereNonPositive(X[1]-(L1-D))
m2 = whereNonNegative(X[1]-(L1-D-H))
m3 = whereNonNegative(X[0]-(L0/2+OFFSET-W/2))
m4 = whereNonPositive(X[0]-(L0/2+OFFSET+W/2))
#  product of the masks defines mask for the anomaly:
m = m1*m2*m3*m4
# ... create control plot:
m_np = convertToNumpy(m)
x_np = convertToNumpy(m.getX())

plt.figure()
plt.tricontourf(x_np[0], x_np[1], m_np[0], 15)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.colorbar()
plt.savefig('output/mt_anomaly.png')
print("anomaly plot generated.")
# now we can set the resistivity distribution:
rho = rho_b*(1-m)+rho_a*m
rho = interpolate(rho, Function(domain))
pde.setValue(A=rho*np.eye(2))
# set magnetic field to one at the top of the domain:
x = domain.getX()
mD = whereZero(x[1]-L1)
print("SS")
pde.setValue(q=mD, r=1)
# now we can get the solution:
Hx = pde.getSolution()
print("SS")
# electrical field:
Ey = rho*grad(Hx, ReducedFunction(domain))[1]
# impedance
Zyx = Ey/Hx
print("SS")
# =======================================
#
# .. we grab the impedance along a transect
#
locations_in_transect = [( h0*k+h0/2. , L1 ) for k in range(0, NE0, 2)]
locator_transect = Locator(where=ReducedFunction(domain), x=locations_in_transect)
#
# .. horizontal positions and Z_yx along transect
#
x_ts = np.array([ x[0] for x in locator_transect.getX() ])
Zyx_ts=np.array(locator_transect(Zyx))
#
# and save to file:
#
np.savetxt('output/MTData2.csv', (x_ts, Zyx_ts.real, Zyx_ts.imag ), delimiter=',')
#
# ... plot the apparent resistivity and phase:
#
rho_a_ts = 1./(2*np.pi*f*Mu0)*abs(Zyx_ts)**2
phi_ts = np.angle(Zyx_ts, deg=True)
#
# and then we can plot these data:
#
plt.figure()
plt.plot(x_ts, rho_a_ts)
plt.xlabel('x [m]')
plt.ylabel('resistivity [m]')
plt.title(f"Resistivity vs offset  (f={f} Hz)")
plt.savefig('output/mt_resistivity.png')
print("resistivity plot generated.")

plt.figure()
plt.plot(x_ts, phi_ts)
plt.xlabel('x [m]')
plt.ylabel('phase [deg]')
plt.title(f"phase vs offset (f={f} Hz)")
plt.savefig("output/mt_phase.png")
print("phase plot generated.")
