# $Id$

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

# run by scons build/posix/modellib/test/python/drucker_prager.passed from check out cd.

from esys.modellib.mechanics import DruckerPrager
from esys.escript.modelframe import Link,Simulation
from esys.modellib.geometry import RectangularDomain, VectorConstrainer
from esys.modellib.input import Sequencer, InterpolateOverBox

debug=True

dom=RectangularDomain(debug)
dom.order=1


sq=Sequencer(debug)
sq.t=0
sq.t_end=1.
sq.dt_max=0.1

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

m=DruckerPrager(debug)
m.domain=Link(dom,"domain")

cv=VectorConstrainer(debug)
cv.domain=Link(dom,"domain")
cv.value=Link(iob,"out")
cv.left=[True, False, False] 
cv.top= [False, False, False]
cv.bottom= [False, False, False]
cv.front= [False, False, False]
cv.back= [False, False, False]
cv.right= [True, True, True]

m.prescribed_velocity=Link(cv,"value_of_constraint")
m.location_prescribed_velocity=Link(cv,"location_of_constraint")

dom.displacement=Link(m,"displacement")

s=Simulation([sq,dom,cv,m],debug=True)
s.run()
