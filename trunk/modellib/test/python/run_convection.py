# $Id$

#
#
#    a very simple convection model in model frame:
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

import os
from esys.escript.modelframe import Link,Simulation
from esys.modellib.geometry import RectangularDomain,ScalarConstrainer,VectorConstrainer
from esys.modellib.input import Sequencer,InterpolateOverBox,GaussianProfile,LinearCombination
from esys.modellib.flow import SteadyIncompressibleFlow
from esys.modellib.temperature import TemperatureAdvection
from esys.modellib.materials import SimpleEarthModel,GravityForce
from esys.modellib.visualization import WriteVTK

dom=RectangularDomain()
dom.order=2

temp_constraints=ScalarConstrainer()
temp_constraints.domain=Link(dom)
temp_constraints.top=1
temp_constraints.bottom=1

vel_constraints=VectorConstrainer()
vel_constraints.domain=Link(dom)
vel_constraints.left=[1,0,0]
vel_constraints.right=[1,0,0]
vel_constraints.top=[0,1,0]
vel_constraints.bottom=[0,1,0]
vel_constraints.front=[0,0,1]
vel_constraints.back=[0,0,1]


temp_val=InterpolateOverBox()
temp_val.domain=Link(dom,"domain")
temp_val.right_top_back=Link(dom,"l")
temp_val.value_left_bottom_front=1.
temp_val.value_right_bottom_front=1.
temp_val.value_left_top_front=0.
temp_val.value_right_top_front=0.
temp_val.value_left_bottom_back=1.
temp_val.value_right_bottom_back=1.
temp_val.value_left_top_back=0.
temp_val.value_right_top_back=0.

mat=SimpleEarthModel()
mat.density0=1.
mat.viscocity0=1.
mat.rayleigh_number=10000.
mat.alpha=0.001

temp=TemperatureAdvection(debug=True)
temp.domain=Link(dom)
temp.density=Link(mat,"density0")
temp.heat_capacity=Link(mat,"heat_capacity")
temp.location_fixed_temperature=Link(temp_constraints,"location_of_constraint")
temp.fixed_temperature=Link(temp_val,"out")
temp.safety_factor=0.01
mat.temperature=Link(temp,"temperature")


grav=GravityForce()
grav.domain=Link(dom,"domain")
grav.direction=[0.,-1.,0.]
grav.density=Link(mat,"density")
grav.gravity=Link(mat,"gravity")


vel=SteadyIncompressibleFlow(debug=True)
vel.domain=Link(dom)
vel.internal_force=Link(grav,"gravity_force")
vel.viscosity=Link(mat,"viscosity")
vel.location_prescribed_velocity=Link(vel_constraints,"location_of_constraint")
vel.rel_tol=1.e-6
temp.velocity=Link(vel,"velocity")

sq=Sequencer()
sq.t_end=0.005

vis=WriteVTK()
vis.t=Link(sq)
vis.scalar=Link(temp,"temperature")
vis.vector=Link(vel,"velocity")
vis.stride=5
vis.filename=os.environ['MODELLIB_WORKING_DIR']+"/temp.xml"

per=GaussianProfile()
per.domain=Link(dom)
per.x_c=[0.5,0.5,0.5]
per.A=0.0001
per.width=0.01
per.r=0

lc=LinearCombination()
lc.f0=1.
lc.v0=Link(per,"out")
lc.f1=1.
lc.v1=Link(temp_val,"out")
temp.temperature=Link(lc,"out")

# s=Simulation([dom,sq,temp_constraints,vel_constraints,temp_val,per,lc,mat,grav,Simulation([vel],debug=True),temp,vis],debug=True)
s=Simulation([dom,sq,Simulation([vel],debug=True),temp,vis],debug=True)
s.writeXML()
s.run()
