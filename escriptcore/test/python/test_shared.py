
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


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Test to ensure that modification of shared Data does not occur
"""

import esys.escriptcore.utestselect as unittest
from esys.escript import *

class Test_Shared(unittest.TestCase):
  def test_setToZero(self):
        d=Data(42)
        e=d.delay()
        d.setToZero()
        self.assertTrue(Lsup(e-42)<=self.tol)
        
  def test_copyConstr(self):
        d=Data(42)
        e=Data(d)
        d+=17
        self.assertTrue(Lsup(e-42)<=self.tol)
        
  # This should not fail (even in the old code) but it doesn't hurt to check
  def test_Copy(self):
        d=Data(42)
        e=d.copy()
        d+=17
        self.assertTrue(Lsup(e-42)<=self.tol)
        
  def  test_eqops(self):
        d=Data(42)
        e=d.delay()
        d+=17
        self.assertTrue(Lsup(e-42)<=self.tol)
        d=Data(42)
        e=d.delay()
        d-=1
        self.assertTrue(Lsup(e-42)<=self.tol)
        d*=3
        d=Data(42)
        e=d.delay()     
        d/=2
        self.assertTrue(Lsup(e-42)<=self.tol)
        
  def test_setItem(self):
        d=Data(42)
        e=d.delay()
        d[tuple()]=17
        self.assertTrue(Lsup(e-42)<=self.tol)


  def test_setTaggedValue(self):
        d=Data(42,self.domain.getX().getFunctionSpace())        # doesn't really matter which non-NULL FS we use
        d.tag()
        self.domain.setTagMap("TestTag",2)
        e=d.delay()
        d.setTaggedValue("TestTag",17)
        e.resolve()
        self.assertTrue(str(e)!=str(d)) 
        e=d.delay()
        d.setTaggedValue(1,12)
        e.resolve()
        self.assertTrue(str(e)!=str(d))


        
        
