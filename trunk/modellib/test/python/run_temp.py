
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

#
#   small test problem fro temperture advection:
#
#   T=x0*x1*exp(-t), v=[1,-1]
#

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript.modelframe import Link,Simulation
from esys.modellib.geometry import RectangularDomain,ScalarConstrainerOverBox
from esys.modellib.input import Sequencer
from esys.modellib.probe import Probe,EvaluateExpression
from esys.modellib.temperature import TemperatureAdvection
import numpy
import os
import sys
if sys.version_info >= (3,0):
    from io import StringIO
else:
    #specifically to avoid non-unicode default strings
    #being passed to a python3 style StringIO that expects unicode
    from StringIO import StringIO

#Link() behaves badly inside a TestCase class
def run(dom, stream):
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
    tt.velocity=numpy.array([1,-1])
    tt.thermal_source=Link(source,"out")
    tt.location_fixed_temperature=Link(constraints,"location_of_constraint")
    tt.fixed_temperature=Link(boundaryvalue,"out")
    tt.safety_factor=0.1

    probe=Probe()
    probe.expression="x[0]*x[1]*exp(-t)"
    probe.t=Link(sqe)
    probe.value=Link(tt,"temperature")


    s=Simulation([sqe,constraints,tt,probe],debug=True)
    s.writeXML(stream)
    s.run()

class Test_RunTemp(unittest.TestCase):
    def setUp(self):
        import sys
        self.old = sys.stdout
        sys.stdout = StringIO()
    
    def tearDown(self):
        import sys
        sys.stdout = self.old
    
    def test_order2(self):
        dom=RectangularDomain()
        dom.order=2
        run(dom, sys.stdout)

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

