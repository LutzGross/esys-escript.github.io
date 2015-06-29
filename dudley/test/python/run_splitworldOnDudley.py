
##############################################################################
#
# Copyright (c) 2015 by The University of Queensland
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2015 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.dudley import Rectangle, Brick, ReadMesh, ReadGmsh
from test_splitworld import Test_SplitWorld, sw_testing


mpisize=getMPISizeWorld()
NE=4 # number elements, must be even

class Test_SplitOnDudley(Test_SplitWorld):
  def setUp(self):
    self.domainpars=[Rectangle, NE, NE]
    
  def tearDown(self):
    del self.domainpars
    
class Test_dudley_sw_2D(sw_testing):
    def setUp(self):
        from esys.dudley import Rectangle
        self.domain_ctr=Rectangle
        self.domain_vec=(6,6)
        self.domain_dict={}

    def tearDown(self):
        del self.domain_ctr
        del self.domain_vec


class Test_dudley_sw_3D(sw_testing):
    def setUp(self):
        from esys.dudley import Brick
        self.domain_ctr=Brick
        self.domain_vec=(6,6,6)
        self.domain_dict={}
        
    def tearDown(self):
        del self.domain_ctr
        del self.domain_vec
    



if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
