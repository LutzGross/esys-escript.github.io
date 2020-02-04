__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross"

"""
This is a simple example for the TE and TM MT applications (forward PDE only).
It is for layered conductivity (resistivity).
Impedance is calculated over various periods (frequencies) and
plots of apparent resistivity and phase are generated. 

The script will not run under more than one MPI rank.
"""
from esys.escript import *
from esys.ripley import Rectangle
from esys.downunder.apps import MT2DTEModel, MT2DTMModel

from esys.weipa import saveVTK, saveSilo
from esys.escript.pdetools import Locator
from math import ceil,pi
import cmath
import os
import numpy as np
import matplotlib.pylab as plt

#
# ... dimensions of the domain and element size
#
DEPTH=2000.
WIDTH=4000.
NEX=400
NEZ=int(NEX/WIDTH*DEPTH+0.5)
L_AIR=DEPTH*0.2   # thickness of the air layer. The setup does not really need this as there is no horizontal conductivity variation
                  # this is added for illustration purposes.

#
# ... layered conductivity 
#
LAYERS=[L_AIR, 300, 100]   # thickness starting from the top 
SIGMA0=1./100.             # conductivity below the lowest layer if there is still room to reach the full DEPTH
SIGMA=[0., SIGMA0, 1/20.]  # layered setup following Schaa et al. 2016, doi: 10.1088/1742-2132/13/2/S59
# SIGMA=[0., SIGMA0, SIGMA0] #  this is for infinite half space (then app. rho=1/SIGMA, phase=45)

#
#  ... horizontal position of observation stations:
#
OFFSET_STATION=WIDTH*0.3              # from the left and right face of the domain 
NUM_STATION=5                         # number of stations 
DEPTH_STATIONS = L_AIR + DEPTH/NEZ/2  # just below the surface
#
# ... create an array of periods to be processed:
#
PERIODS=np.logspace(-4, 4, num=11, endpoint=True, base=10.0, dtype=float)
# you want to increase the value for `num` to get nicer plots.

#================================================

print("Width [m] = ",WIDTH)
print("Depth [m] = ", DEPTH)
print("Mesh size [m] = ", WIDTH/NEX)
print("Number of cells in x direction = ", NEX)
print("Number of cells in z direction = ", NEZ)
print("Air layer [m] = ", L_AIR)
print("periods [s] = ", PERIODS)

#==============================================

print("generating mesh ...")
domain=Rectangle(NEX,NEZ,l0=WIDTH,l1=DEPTH)

print("generating conductivity ...")
# you can replace this by an interpolation table read from a CSV file.
z=ReducedFunction(domain).getX()[domain.getDim()-1]
sigma=Scalar(SIGMA0, z.getFunctionSpace())
rho=Scalar(1./SIGMA0, z.getFunctionSpace())
z_top=sup(domain.getX()[domain.getDim()-1])
m_top=0.
for l, s in zip(LAYERS, SIGMA ):
    m2=wherePositive(z-z_top+l)
    m=m2-m_top
    sigma=(1-m)*sigma+m*s
    if s > 0:
       rho=(1-m)*rho+m*1./s
    else:
       rho=(1-m)*rho+m*0. # arbitray number as air_layer is backed out in TM mode.
       
    z_top, m_top=z_top-l, m2
    
print("sigma =", sigma)
print("rho =", rho)
#
# ... create Locator to get impedance at position of observations:
#
stationX=np.linspace(OFFSET_STATION, WIDTH-OFFSET_STATION, num=NUM_STATION, endpoint=True)
loc=Locator(ReducedFunction(domain), [ (s, DEPTH-DEPTH_STATIONS) for s in stationX])

print("position of observation stations %s:"%(NUM_STATION//2,), loc.getX()[NUM_STATION//2]) # moved to the next element center!!!
#================================================================================================
FRQ=1./PERIODS
#================================================================================================
print("Start TM mode ...")
model=MT2DTMModel(domain, airLayer=DEPTH-L_AIR)
model.setResistivity(rho, rho_boundary=1/SIGMA0) # rho can be interpolated to the boundary (e.g. when conductivity is given on node) rho_boundary can be omited. 

# collects app. rho and phase for the frequencies:
arho_TM=[]
phase_TM=[]

for frq in FRQ:
    Zyx = model.getImpedance(f=frq)
    arho=model.getApparentResitivity(frq, Zyx)
    phi=model.getPhase(frq, Zyx)
    #saveVTK('sol', rho=rho, Hx=abs(self.Hx), rho=rho, phi=phi) # write to file for visualization.
    arho_TM.append(loc(arho)[NUM_STATION//2])
    phase_TM.append(loc(phi)[NUM_STATION//2])
    print("frequency %s Hz completed. app. rho = %s, phase = %s deg"%(frq, arho_TM[-1], phase_TM[-1]))
    
# ... and let's plot this:
plt.clf()
plt.plot(PERIODS, arho_TM)
plt.xlabel('period [s]')
plt.ylabel('app. resistivity [Ohm m]')
plt.grid(True)
plt.xscale("log")
plt.title("TM Mode: Apparent resistivity vs. period")
plt.savefig("TM_appResistivity.png")
print("app. resistivity image written to file TM_appResistivity.png")
#plt.show()

plt.clf()
plt.plot(PERIODS, phase_TM)
plt.xlabel('period [s]')
plt.ylabel('phase [deg]')
plt.grid(True)
plt.xscale("log")
plt.title("TM Mode: Phase vs. period")
plt.savefig("TM_Phases.png")
print("phase image written to file TM_Phases.png")
#plt.show()

print("End TM mode ...")
#=================================================================================================

#=================================================================================================
print("Start TE mode ...")
model=MT2DTEModel(domain)
model.setConductivity(sigma, sigma_boundary=SIGMA0)

arho_TE=[]
phase_TE=[]
for frq in FRQ:
    Zxy = model.getImpedance(f=frq)
    arho=model.getApparentResitivity(frq, Zxy)
    phi=model.getPhase(frq, Zxy)
    #saveVTK('sol', sigma=sigma, Ex=abs(model.Ex), rho=arho, phi=phi)
    arho_TE.append(loc(arho)[NUM_STATION//2])
    phase_TE.append(loc(phi)[NUM_STATION//2])
    print("frequency %s Hz completed. app. rho = %s, phase = %s deg"%(frq, arho_TE[-1], phase_TE[-1]))
    
    #print(loc.getX()[NUM_STATION//2][0], frq, loc(rho)[NUM_STATION//2], loc(phi)[NUM_STATION//2])

# ... and let's plot this:
plt.clf()
plt.plot(PERIODS, arho_TE)
plt.xlabel('period [s]')
plt.ylabel('app. resistivity [Ohm m]')
plt.grid(True)
plt.xscale("log")
plt.title("TE Mode: Apparent resistivity vs. period")
plt.savefig("TE_appResistivity.png")
print("app. resistivity image written to file TE_appResistivity.png")
#plt.show()

plt.clf()
plt.plot(PERIODS, phase_TE)
plt.xlabel('period [s]')
plt.ylabel('phase [deg]')
plt.grid(True)
plt.xscale("log")
plt.title("TE Mode: Phase vs. period")
plt.savefig("TE_Phases.png")
print("phase image written to file TE_Phases.png")
#plt.show()

print("End TE mode ...")
#=================================================================================================


