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

from __future__ import division, print_function

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

#
#   Testing reading domains
#
#   Tests for finley and dudley reading .fly and .gmsh
#

import os
from subprocess import PIPE, Popen
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
try:
    import esys.dudley
    HAVE_DUDLEY = True
except ImportError:
    HAVE_DUDLEY = False

try:
    import esys.finley
    from esys.modellib.geometry import *
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

from esys.escript import getMPISizeWorld, hasFeature
from esys.escript.modelframe import DataSource
from esys.pycad.gmsh import *
from esys.pycad import *

# Putting it in modellib dir
try:
     MODELLIB_WORKDIR=os.environ['MODELLIB_WORKDIR']
except KeyError:
     MODELLIB_WORKDIR='.'

GMSH = hasFeature('gmsh')
mpiSize = getMPISizeWorld()

@unittest.skipIf(not GMSH, "gmsh not available")
@unittest.skipIf(mpiSize > 1, "not tested with more than 1 MPI rank")
class Test_domainReaders(unittest.TestCase):
    def domain_family(self, dommodule, f):
        dom=RectangularDomain(dommodule, parameters=["fish","dummy"], debug=True)
        # need to write to both .fly and .gmsh
        dom.domain().write(f+"dr.fly")
        r1=DomainReader(dommodule)
        r1.source=DataSource(uri=f+"dr.fly", fileformat="fly")
        r1.domain()

        des=Design()
        b=Volume(Brick(Point(0,0,0), Point(1,1,1)))
        des.addItems(b)
        des.setMeshFileName(os.path.join(MODELLIB_WORKDIR,"TESTgmsh_test.msh"))
        esys.finley.MakeDomain(des)
        
        r2=DomainReader(dommodule)
        r2.source=DataSource(uri=os.path.join(MODELLIB_WORKDIR,"TESTgmsh_test.msh"), fileformat="gmsh")
        r2.domain()

    @unittest.skipIf(not HAVE_FINLEY, "Finley module not available")
    def test_domain_families(self):
        self.domain_family(None,os.path.join(MODELLIB_WORKDIR,'TESTnone'))
        if HAVE_DUDLEY:
            self.domain_family(esys.dudley,os.path.join(MODELLIB_WORKDIR,'TESTDud'))
        self.domain_family(esys.finley,os.path.join(MODELLIB_WORKDIR,'TESTfin'))

    @unittest.skipIf(not HAVE_FINLEY, "Finley module not available")
    def test_finleyReader(self):
        self.domain_family(esys.finley,os.path.join(MODELLIB_WORKDIR,'TESTfin'))
        rf=FinleyReader(parameters=["fish","dummy"], debug=True)
        rf.source=DataSource(uri=os.path.join(MODELLIB_WORKDIR,"TESTfindr.fly"), fileformat="fly")
        rf.domain()


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

