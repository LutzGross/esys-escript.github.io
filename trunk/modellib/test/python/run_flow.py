##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

#
#   small test problem fro temperture advection:
#
#   p=(x0+x1)*t
#

import os
import sys
if sys.version_info >= (3,0):
    from io import StringIO
else:
    #specifically to avoid non-unicode default strings
    #being passed to a python3 style StringIO that expects unicode
    from StringIO import StringIO
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript.modelframe import Link,Simulation
from esys.escript import hasFeature
from esys.modellib.input import Sequencer
from esys.modellib.probe import Probe,EvaluateExpression
from esys.modellib.flow import SteadyIncompressibleFlow

try:
    import esys.finley
    from esys.modellib.geometry import RectangularDomain,VectorConstrainerOverBox
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

# TODO: once Amesos2 can deal with block matrices uncomment
HAVE_DIRECT = hasFeature("PASO_DIRECT") #or hasFeature('trilinos')

#Link() behaves badly inside a TestCase class
def run(dom, stream):
    constraints=VectorConstrainerOverBox()
    constraints.domain=Link(dom)
    constraints.left=[1,0,0]
    constraints.right=[1,0,0]
    constraints.top=[0,1,0]
    constraints.bottom=[0,1,0]
    constraints.front=[0,0,1]
    constraints.back=[0,0,1]

    sqe=Sequencer()
    sqe.dt_max=0.5
    sqe.t_end=1.

    source=EvaluateExpression()
    source.domain=Link(dom)
    source.t=Link(sqe)
    source.expression=["t","t"]

    flow=SteadyIncompressibleFlow()
    flow.domain=Link(dom,"domain")
    flow.internal_force=Link(source,"out")
    flow.location_prescribed_velocity=Link(constraints,"location_of_constraint")
    flow.prescribed_velocity=[0.,0.]

    ptest=Probe()
    ptest.expression="(x[0]+x[1]-1.)*t"
    ptest.t=Link(sqe)
    ptest.value=Link(flow,"pressure")

    s=Simulation([sqe,constraints,Simulation([flow],debug=True),ptest],debug=True)
    s.writeXML(stream)
    s.run()

class Test_RunFlow(unittest.TestCase):
    def setUp(self):
        import sys
        self.old = sys.stdout
        sys.stdout = StringIO()
    
    def tearDown(self):
        import sys
        sys.stdout = self.old

    @unittest.skipIf(not HAVE_FINLEY, "Finley module not available")
    @unittest.skipIf(not HAVE_DIRECT, "Direct solver not available")
    def test_order2(self):
        dom=RectangularDomain()
        dom.order=2
        run(dom, sys.stdout)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
