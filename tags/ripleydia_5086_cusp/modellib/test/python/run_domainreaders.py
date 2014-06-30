
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

from __future__ import print_function

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
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
import esys.finley
import esys.dudley

from esys.escript import getMPISizeWorld
from esys.escript.modelframe import DataSource
from esys.modellib.geometry import *
from esys.pycad.gmsh import *
from esys.pycad import *

# Putting it in modellib dir
try:
     MODELLIB_WORKDIR=os.environ['MODELLIB_WORKDIR']
except KeyError:
     MODELLIB_WORKDIR='.'

GMSH = None
try:
    p=Popen(['gmsh', '-info'], stderr=PIPE)
    _,e=p.communicate()
    if e.split().count("MPI"):
        GMSH = 'm'
    else:
        GMSH = 's'
except OSError:
    pass
    
@unittest.skipIf(GMSH is None, "gmsh not available")
@unittest.skipIf(esys.escript.getEscriptParamInt("MPIBUILD",0)>0,
        "not tested with MPI builds")
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

    def test_domain_families(self):
        self.domain_family(None,os.path.join(MODELLIB_WORKDIR,'TESTnone'))
        self.domain_family(esys.dudley,os.path.join(MODELLIB_WORKDIR,'TESTDud'))
        self.domain_family(esys.finley,os.path.join(MODELLIB_WORKDIR,'TESTfin'))
   
    def test_finleyReader(self):
        self.domain_family(esys.finley,os.path.join(MODELLIB_WORKDIR,'TESTfin'))
        rf=FinleyReader(parameters=["fish","dummy"], debug=True)
        rf.source=DataSource(uri=os.path.join(MODELLIB_WORKDIR,"TESTfindr.fly"), fileformat="fly")
        rf.domain()


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

