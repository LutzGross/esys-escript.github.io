#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

# run by scons build/posix/modellib/test/python/drucker_prager.passed from check out cd.

import os
from esys.modellib.mechanics import DruckerPrager
from esys.escript.modelframe import Link,Simulation
from esys.modellib.geometry import RectangularDomain, VectorConstrainerOverBox, UpdateGeometry
from esys.modellib.input import Sequencer, InterpolateOverBox
from esys.modellib.visualization import WriteVTK

try:
   WORKDIR=os.environ['MODELLIB_WORKDIR']
except KeyError:
   WORKDIR='.'


dbg=True

dom=RectangularDomain(debug=dbg)
dom.dim=3
dom.l=[0.5,1.,1.]
dom.n=[30,6,6]
dom.order=1
dom.integrationOrder=-1


sq=Sequencer(debug=dbg)
sq.t=0
sq.t_end=0.8
sq.dt_max=100.

iob=InterpolateOverBox(debug=dbg)
iob.domain=Link(dom,"domain")
iob.value_left_bottom_front=[-1.,0.,0.]
iob.value_right_bottom_front=[0.,0.,0.]
iob.value_left_bottom_back=[-1.,0.,0.]
iob.value_right_bottom_back=[0.,0.,0.]
iob.value_left_top_front=[-1.,0.,0.]
iob.value_right_top_front=[0.,0.,0.]
iob.value_left_top_back=[-1.,0.,0.]
iob.value_right_top_back=[0.,0.,0.]

# iob.value_left_bottom_front=[-1.,0.]
# iob.value_right_bottom_front=[0.,0.]
# iob.value_left_bottom_back=[-1.,0.]
# iob.value_right_bottom_back=[0.,0.]
# iob.value_left_top_front=[-1.,0.]
# iob.value_right_top_front=[0.,0.]
# iob.value_left_top_back=[-1.,0.]
# iob.value_right_top_back=[0.,0.]

m=DruckerPrager(debug=dbg)
m.domain=Link(dom,"domain")

cv=VectorConstrainerOverBox(debug=dbg)
cv.domain=Link(dom,"domain")
cv.value=Link(iob,"out")

cv.left=[True, False, False] 
cv.right= [True, True, True]
cv.bottom= [False, False, False]
cv.top= [False, False, False]
cv.front= [False, False, False]
cv.back= [False, False, False]

m.velocity=Link(iob,"out")
m.prescribed_velocity=Link(cv,"value_of_constraint")
m.location_prescribed_velocity=Link(cv,"location_of_constraint")
m.rel_tol=0.01

m.expansion_coefficient= 0.
m.bulk_modulus=1000.
m.shear_modulus=1.
m.plastic_stress=0.
m.friction_parameter=0.
m.dilatancy_parameter=0.
m.shear_length=m.shear_modulus*0.9


ug=UpdateGeometry(debug=dbg)
ug.domain=Link(dom,"domain")
ug.displacement=Link(m,"displacement")

vis=WriteVTK(debug=dbg)
vis.t=Link(sq)
vis.data0=Link(m,"plastic_stress")
vis.data1=Link(m,"velocity")
vis.data2=Link(m,"stress")
vis.dt=0.01
vis.filename=WORKDIR+"/temp.xml"

s=Simulation([sq,cv,m,vis],debug=True)
s.writeXML()
s.run()
