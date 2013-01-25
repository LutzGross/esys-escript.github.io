# $Id$
#
#   small test problem fro temperture advection:
#
#   T=x0*x1*exp(-t), v=[1,-1]
#
from escript.modelframe import Link,Simulation
from modellib.geometry import RectangularDomain,ScalarConstrainer
from modellib.input import Sequencer
from modellib.probe import Probe,EvaluateExpression
from modellib.temperature import TemperatureAdvection
import numarray

dom=RectangularDomain()
dom.order=2

sqe=Sequencer()
sqe.t=0
sqe.t_end=0.1

constraints=ScalarConstrainer()
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
tt.velocity=numarray.array([1,-1,0])
tt.thermal_source=Link(source,"out")
tt.location_fixed_temperature=Link(constraints,"location_of_constraint")
tt.fixed_temperature=Link(boundaryvalue,"out")
tt.safety_factor=0.01

probe=Probe()
probe.expression="x[0]*x[1]*exp(-t)"
probe.t=Link(sqe)
probe.value=Link(tt,"temperature")


s=Simulation([sqe,dom,tt,probe],debug=True)
s.writeXML()
s.run()
