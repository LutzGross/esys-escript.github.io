"""
simple seismic solver in frequency domain 
with reflectors and PML.
"""

from esys.escript import *
from esys.finley import Rectangle
from esys.downunder.apps import SonicWaveInFrequencyDomain, PMLCondition
from esys.weipa import saveSilo
from esys.escript.pdetools import Locator

import numpy as np

# horizontal and vertical grid spacing in 
dx=10.
dz=10.

# grid size
NEx=400
NEz=180

# element order (1 or 2):
Order=2

# extend of domain:
Width=NEx*dx
Depth=NEz*dz

# position of geophones:
StationOffset = 60*dx # position of first phone from left boundary
StationSpacing = dx # phone spacing
NStations = int((Width-2*StationOffset)//StationSpacing+0.5)+1 # number of stations

# velocity:
Layers=[     30*dx ,   40*dx,  90*dx ] # thickness of layer from the top down 
Vs        =[ 1500.,    2000.,  3000. ] # corresponding velocities
Vbase=3500                             # base velocity

Frequency = 8.   # frequency
Amplitude= 1.     # and amplitude

# thickness of PML zone:
L_PML=40*dx


print("Domain = %g  x %g m"%(Width, Depth))
print("Grid Cells= %d x %d"%(NEx, NEz))
print("Grid Nodes= %d x %d"%(Order*NEx+1, Order*NEz+1))
print("Grid Spacing = %g  x %g m"%(dx, dz))
print("PML zone = %g m"%(L_PML))
print("Number of Stations = %d"%(NStations,))
print("Station spacing = %g m"%(StationSpacing,))
print("Station offset  = %g m"%(StationOffset,))
print("Reflector Layers [m]: ", Layers)
print("Velocities [m/s]: ", Vs)
print("Base velocity : [m/s]", Vbase)
print("Frequency = %g Hz "%Frequency)
print("Amplitude = %g"%Amplitude)


# generate domain:
stationPositions = [ (StationOffset+s*StationSpacing, Depth ) for s in range(NStations)]
stationTags= [ "P%4.4d"%s for s in range(NStations)]
domain=Rectangle(NEx,NEz,l0=Width,l1=Depth, diracPoints=stationPositions, diracTags=stationTags, order=Order, fullOrder=True)
print("Domain has been generated.")

Dim=domain.getDim()
# setup velocity field:
vp=Scalar(Vbase, Function(domain))
X=vp.getX()
Top=sup(domain.getX()[Dim-1])
for l in range(len(Vs)-1, -1, -1 ):
    mask=wherePositive(X[Dim-1]-(Top-Layers[l]))
    vp=vp*(1-mask)+mask*Vs[l]
print("Velocity field generated:",vp)

# check spacing vs, wavelength:
wavelength=inf(vp/(2*np.pi*Frequency))
print("minimum wavelength   = %g"%wavelength)
assert wavelength > max(dx,dz)/Order


# locator to grap amplitudes at geophones:
geophones=Locator(Solution(domain), stationPositions)


# create PML condition:
pml_condition=PMLCondition(sigma0=vp*100/Frequency, Lleft=[L_PML,L_PML], Lright=[L_PML, None], m=4)
print("PMLConditon generated ...")

# sonic wave model
wave=SonicWaveInFrequencyDomain(domain, vp=vp, pml_condition=pml_condition)
wave.setFrequency(Frequency)
print("model created...")
# set source and get wave response: 
source=Scalar(0j, DiracDeltaFunctions(domain))
source.setTaggedValue(stationTags[NStations//2], Amplitude)
u=wave.getWave(source)
print("model solved...")

# get amplitude and phase
amps=abs(u)
phase=phase(u)/np.pi*180.


# save to file:
saveSilo("solution", velocity=vp, amplitude=amps, phase=phase, pmlzone=pml_condition.getPMLMask(domain))
print("solution saved to file ...")

# get the amplitude and phase at geophones:
import matplotlib.pyplot as plt

fig=plt.gcf()
ax1=plt.gca()

geophoneX=[ x[0] for x in geophones.getX()]
color = 'tab:red'
ax1.set_xlabel('offset')
ax1.set_ylabel('amplitude', color=color)
ax1.plot(geophoneX , geophones(amps), color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('phase [deg]', color=color)  # we already handled the x-label with ax1
ax2.plot(geophoneX , geophones(phase), color=color, marker=".")
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()
