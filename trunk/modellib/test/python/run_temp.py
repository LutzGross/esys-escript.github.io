
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

#
#   small test problem fro temperture advection:
#
#   T=x0*x1*exp(-t), v=[1,-1]
#

from esys.escript.modelframe import Link,Simulation
from esys.modellib.geometry import RectangularDomain,ScalarConstrainerOverBox
from esys.modellib.input import Sequencer
from esys.modellib.probe import Probe,EvaluateExpression
from esys.modellib.temperature import TemperatureAdvection
import numarray
import os

dom=RectangularDomain()
dom.order=2

sqe=Sequencer()
sqe.t=0
sqe.t_end=0.05

constraints=ScalarConstrainerOverBox()
constraints.domain=Link(dom)
constraints.top=1
constraints.bottom=1
constraints.right=1
constraints.left=1

source=EvaluateExpression()
source.domain=Link(dom)
source.expression="(x[1]-x[0])*exp(-t)-exp(-t)*x[0]*x[1]"
source.t=Link(sqe)
boundaryvalue=EvaluateExpression()
boundaryvalue.domain=Link(dom)
boundaryvalue.expression="x[0]*x[1]*exp(-t)"
boundaryvalue.t=Link(sqe)

tt=TemperatureAdvection()
tt.domain=Link(dom)
tt.temperature=Link(boundaryvalue,"out")
tt.velocity=numarray.array([1,-1])
tt.thermal_source=Link(source,"out")
tt.location_fixed_temperature=Link(constraints,"location_of_constraint")
tt.fixed_temperature=Link(boundaryvalue,"out")
tt.safety_factor=0.1

probe=Probe()
probe.expression="x[0]*x[1]*exp(-t)"
probe.t=Link(sqe)
probe.value=Link(tt,"temperature")


s=Simulation([sqe,constraints,tt,probe],debug=True)
s.writeXML()
s.run()
