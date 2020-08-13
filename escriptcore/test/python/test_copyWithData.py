
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

import esys.escriptcore.utestselect as unittest
from esys.escript import *
from esys.escriptcore.testing import *

# class Test_copyWithMask(unittest.TestCase):
    #  TODO


    # def test_copyWithMask_sss(self):
    #     d=Scalar(10.0,Function(self.domain))
    #     c=Scalar(0.0,Function(self.domain))
    #     m=Scalar(1.0,Function(self.domain))
    #     d.copyWithMask(c,m)
    #     tup=d.toListOfTuples()
    #     self.assertTrue(max(tup)==0,"copyWithMask")

    # def test_copyWithMask_vss(self):
    #     d=Vector([10.0,10.0],Function(self.domain))
    #     c=Scalar(0.0,Function(self.domain))
    #     m=Scalar(1.0,Function(self.domain))
    #     d.copyWithMask(c,m)
    #     tup=d.toListOfTuples()
    #     self.assertTrue(max(max(tup))==0,"copyWithMask")

    # def test_copyWithMask_vvs(self):
    #     d=Vector([10.0,10.0],Function(self.domain))
    #     c=Vector([0.0,0.0],Function(self.domain))
    #     m=Scalar(1.0,Function(self.domain))
    #     d.copyWithMask(c,m)
    #     tup=d.toListOfTuples()
    #     self.assertTrue(max(max(tup))==0,"copyWithMask")

    # def test_copyWithMask_vvv(self):
    #     d=Vector([10.0,10.0],Function(self.domain))
    #     c=Vector([0.0,0.0],Function(self.domain))
    #     m=Vector([1.0,1.0],Function(self.domain))
    #     d.copyWithMask(c,m)
    #     tup=d.toListOfTuples()
    #     self.assertTrue(max(max(tup))==0,"copyWithMask")

    # def test_copyWithMask_vvv(self):
    #     d=Vector([10.0,10.0],Function(self.domain))
    #     c=Scalar(0.0,Function(self.domain))
    #     m=Vector([1.0,1.0],Function(self.domain))
    #     d.copyWithMask(c,m)
    #     tup=d.toListOfTuples()
    #     self.assertTrue(max(max(tup))==0,"copyWithMask")


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
