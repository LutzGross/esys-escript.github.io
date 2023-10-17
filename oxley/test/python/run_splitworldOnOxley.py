
##############################################################################
#
# Copyright (c) 2015-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2015-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.oxley import Rectangle, Brick
from test_splitworld import Test_SplitWorld, sw_testing, sw_testmany


mpisize=getMPISizeWorld()
NE=4 # number elements, must be even

#TODO
# class Test_SplitOnOxley(Test_SplitWorld):
#   def setUp(self):
#     self.domainpars=[Rectangle, NE, NE]
    
#   def tearDown(self):
#     del self.domainpars
    
# class Test_Oxley_sw_2D(sw_testing):
#     def setUp(self):
#         from esys.oxley import Rectangle
#         self.domain_ctr=Rectangle
#         self.domain_vec=(6,6)
#         self.domain_dict={}

#     def tearDown(self):
#         del self.domain_ctr
#         del self.domain_vec


# class Test_Oxley_sw_3D(sw_testing):
#     def setUp(self):
#         from esys.oxley import Brick
#         self.domain_ctr=Brick
#         self.domain_vec=(6,6,6)
#         self.domain_dict={}
        
#     def tearDown(self):
#         del self.domain_ctr
#         del self.domain_vec
    
# class Test_SplitOnMultiOxley(Test_SplitWorld):
#   def setUp(self):
#     self.domainpars=[MultiRectangle, NE, NE]
    
#   def tearDown(self):
#     del self.domainpars

#TODO  
# class Test_multioxley_sw_2D(sw_testing):
#     def setUp(self):
#         from esys.oxley import MultiRectangle
#         self.domain_ctr=MultiRectangle
#         self.domain_vec=(6,6)
#         self.domain_dict={}

#     def tearDown(self):
#         del self.domain_ctr
#         del self.domain_vec


# class Test_multioxley_sw_3D(sw_testmany):
#     def setUp(self):
#         from esys.oxley import MultiBrick
#         self.domain_ctr=MultiBrick
#         self.domain_vec=(6,6,6)
#         self.domain_dict={}
        
#     def tearDown(self):
#         del self.domain_ctr
#         del self.domain_vec


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
