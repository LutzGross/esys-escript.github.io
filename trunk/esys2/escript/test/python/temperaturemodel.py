# $Id$
import numarray
from esys.modelframe import Link,Simulation,ExplicitSimulation
from esys.visualization import WriteVTK
from esys.materials import SimpleMaterialTable
from esys.geometry import RectangularDomain,ScalarConstrainer
from esys.input import InterpolatedTimeProfile,GausseanProfile
from esys.temperature import TemperatureDiffusion

dom=RectangularDomain()
constraints=ScalarConstrainer()
constraints.domain=Link(dom)
constraints.top=1
constraints.bottom=1
   
mt=SimpleMaterialTable()
  
pf=InterpolatedTimeProfile()
pf.t=[0.,0.25,0.5,0.75]
pf.values=[0.,1.,1.,0.]
   
q=GausseanProfile()
q.domain=Link(dom)
q.width=0.05
q.x_c=numarray.array([0.5,0.5,0.5])
q.r=0.01
q.A=Link(pf,"out")
   
tt=TemperatureDiffusion()
tt.domain=Link(dom)
tt.tend=1.
tt.dt=0.1
tt.temperature=0.
tt.density=Link(mt)
tt.c_p=Link(mt)
tt.thermal_permabilty=Link(mt)
tt.reference_temperature=0.
tt.radiation_coefficient=Link(mt)
tt.thermal_source=Link(q,"out")
tt.location_fixed_temperature=Link(constraints,"location_of_constraint")
tt.implicit=True
   
vis=WriteVTK()
vis.scalar=Link(tt,"temperature")
   
s=ExplicitSimulation([dom,constraints,pf,q,Simulation([mt,tt],debug=True),vis],debug=True)
s.writeXML()
s.run()
