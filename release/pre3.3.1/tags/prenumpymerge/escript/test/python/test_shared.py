
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
Test to ensure that modification of shared Data does not occur
"""

import unittest
from esys.escript import *

class Test_Shared(unittest.TestCase):
  def test_setToZero(self):
	d=Data(42)
	e=d.delay()
	d.setToZero()
	self.failUnless(Lsup(e-42)<=self.tol)
	
  def test_copyConstr(self):
	d=Data(42)
	e=Data(d)
	d+=17
	self.failUnless(Lsup(e-42)<=self.tol)
	
  # This should not fail (even in the old code) but it doesn't hurt to check
  def test_Copy(self):
  	d=Data(42)
	e=d.copy()
	d+=17
	self.failUnless(Lsup(e-42)<=self.tol)
	
  def  test_eqops(self):
	d=Data(42)
	e=d.delay()
	d+=17
	self.failUnless(Lsup(e-42)<=self.tol)
	d=Data(42)
	e=d.delay()
	d-=1
	self.failUnless(Lsup(e-42)<=self.tol)
	d*=3
	d=Data(42)
	e=d.delay()	
	d/=2
	self.failUnless(Lsup(e-42)<=self.tol)
	
  def test_setItem(self):
	d=Data(42)
	e=d.delay()
	d[tuple()]=17
	self.failUnless(Lsup(e-42)<=self.tol)


  def test_setTaggedValue(self):
	d=Data(42,self.domain.getX().getFunctionSpace())	# doesn't really matter which non-NULL FS we use
	d.tag()
	self.domain.setTagMap("TestTag",2)
	e=d.delay()
	d.setTaggedValue("TestTag",17)
	e.resolve()
	self.failUnless(str(e)!=str(d))	
	e=d.delay()
	d.setTaggedValue(1,12)
	e.resolve()
	self.failUnless(str(e)!=str(d))


	
	