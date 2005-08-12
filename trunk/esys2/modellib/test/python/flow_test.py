# $Id$

#
#   small test problem fro temperture advection:
#
#   p=(x0+x1)*t
#
from escript.modelframe import Link,Simulation
from modellib.geometry import RectangularDomain,VectorConstrainer
from modellib.input import Sequencer
from modellib.probe import Probe,EvaluateExpression
from modellib.flow import SteadyIncompressibleFlow

dom=RectangularDomain()
dom.order=2

constraints=VectorConstrainer()
constraints.domain=Link(dom)
constraints.left=[1,0,0]
constraints.right=[1,0,0]
constraints.top=[0,1,0]
constraints.bottom=[0,1,0]
constraints.front=[0,0,1]
constraints.back=[0,0,1]

sqe=Sequencer()
sqe.dt_max=0.3
sqe.t_end=1.

source=EvaluateExpression()
source.domain=Link(dom)
source.t=Link(sqe)
source.expression=["t","t"]

flow=SteadyIncompressibleFlow()
flow.domain=Link(dom)
flow.internal_force=Link(source,"out")
flow.location_prescribed_velocity=Link(constraints,"location_of_constraint")
flow.prescribed_velocity=[0.,0.]

ptest=Probe()
ptest.expression="(x[0]+x[1]-1.)*t"
ptest.t=Link(sqe)
ptest.value=Link(flow,"pressure")

s=Simulation([dom,sqe,constraints,Simulation([source,flow],debug=True),ptest],debug=True)
s.writeXML()
s.run()
