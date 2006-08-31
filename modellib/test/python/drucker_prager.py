# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

# run by scons build/posix/modellib/test/python/drucker_prager.passed from check out cd.

import os
from esys.modellib.mechanics import DruckerPrager
from esys.escript.modelframe import Link,Simulation
from esys.modellib.geometry import RectangularDomain, VectorConstrainer, UpdateGeometry
from esys.modellib.input import Sequencer, InterpolateOverBox
from esys.modellib.visualization import WriteVTK

try:
   WORKDIR=os.environ['MODELLIB_WORKDIR']
except KeyError:
   WORKDIR='.'


debug=True

dom=RectangularDomain(debug)
dom.l=[5.,1.,1.]
dom.n=[25,5,5]
dom.order=1


sq=Sequencer(debug)
sq.t=0
sq.t_end=3.0
sq.dt_max=10.

iob=InterpolateOverBox(debug)
iob.domain=Link(dom,"domain")
iob.value_left_bottom_front=[1.,0.,0.]
iob.value_right_bottom_front=[0.,0.,0.]
iob.value_left_bottom_back=[1.,0.,0.]
iob.value_right_bottom_back=[0.,0.,0.]
iob.value_left_top_front=[1.,0.,0.]
iob.value_right_top_front=[0.,0.,0.]
iob.value_left_top_back=[1.,0.,0.]
iob.value_right_top_back=[0.,0.,0.]

iob.value_left_bottom_front=[1.,0.]
iob.value_right_bottom_front=[0.,0.]
iob.value_left_bottom_back=[1.,0.]
iob.value_right_bottom_back=[0.,0.]
iob.value_left_top_front=[1.,0.]
iob.value_right_top_front=[0.,0.]
iob.value_left_top_back=[1.,0.]
iob.value_right_top_back=[0.,0.]

m=DruckerPrager(debug)
m.domain=Link(dom,"domain")

cv=VectorConstrainer(debug)
cv.domain=Link(dom,"domain")
cv.value=Link(iob,"out")

cv.left=[True, False, False] 
cv.right= [True, False, False]
cv.bottom= [False, True, False]
cv.top= [False, False, False]
cv.front= [False, False, True]
cv.back= [False, False, False]

m.velocity=Link(iob,"out")
m.prescribed_velocity=Link(cv,"value_of_constraint")
m.location_prescribed_velocity=Link(cv,"location_of_constraint")
m.rel_tol=0.005

m.expansion_coefficient= 0.
m.bulk_modulus=10000.
m.shear_modulus=1.
m.plastic_stress=0.
m.friction_parameter=0.
m.dilatancy_parameter=0.
m.shear_length=m.shear_modulus*2000.


ug=UpdateGeometry(debug)
ug.domain=Link(dom,"domain")
ug.displacement=Link(m,"displacement")

vis=WriteVTK()
vis.t=Link(sq)
vis.vector=Link(m,"velocity")
# vis.scalar=Link(m,"plastic_stress")
vis.tensor=Link(m,"stress")
vis.dt=0.2
vis.filename=WORKDIR+"/temp.xml"

s=Simulation([sq,m,vis],debug=True)
s.run()
