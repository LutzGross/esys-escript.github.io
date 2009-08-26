
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import unittest
import sys
from esys.escript import *
from esys.escript.util import EPSILON

# This test assumes that samples are in order in the object, ie the last point of sample x is 
# immediately before the first point of sample x+1
class DataAccessTestCase(unittest.TestCase):
    #This is a very basic test - it only contains one value.
    def testtoListOfTuplesScalarOnNullDomain(self):
	inp=42.0
	d=Data(inp)
	t=d.toListOfTuples(scalarastuple=True)[0]
	self.failUnless(type(t)==type((1.0,)), "Did not return tuple for scalar data")
	self.failUnless(abs(inp-t[0])<EPSILON, "Did not return correct value")
	t=d.toListOfTuples(scalarastuple=False)[0]
	self.failUnless(type(t)==float, "Did not return non-tuple when asked")
	self.failUnless(abs(inp-t)<EPSILON, "Did not return correct (non-tuple)")

    #Test for one point per sample
    def testtoListOfTuples_SinglePPS(self):
	fs=getTestDomainFunctionSpace(1,5)
	inp=42.0
	d=Data(inp,fs,True)
	for x in range(5):
	    d.setValueOfDataPoint(x,43+x)
	t=d.toListOfTuples(scalarastuple=True)
	self.failUnless(len(t)==5,"Returned list has the wrong length")
	total=0
	for x in range(5):
	    total+=t[x][0]-43-x
	self.failUnless(abs(total)<EPSILON,"Returned list has wrong elements")
	inp=[[1,2],[3,4]]
	d=Data(inp,fs,True)
	for x in range(5):
	    d.setValueOfDataPoint(x,((x,x+1),(x+2,x+3)))
	t=d.toListOfTuples(scalarastuple=True)
	ok=True
	for x in range(5):
	    if t[x]!=((x,x+1),(x+2,x+3)): ok=False
	self.failUnless(ok, "Returned matrix does not match")
	
    def testtoListOfTuples_MultiPPS(self):
	fs=getTestDomainFunctionSpace(3,5)
	inp=0
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,x)
	ok=True
	t=d.toListOfTuples(scalarastuple=True)
	for x in range(15):
	    if t[x]!=(x,): ok=False
	self.failUnless(ok,"Returned scalar does not match")
	inp=(0,0)
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,(2*(x/2),(2*(x/2)+1)))
	t=d.toListOfTuples(scalarastuple=True)
	ok=True
	for x in range(15):
	    if t[x]!=(2*(x/2),2*(x/2)+1): ok=False
	self.failUnless(ok,"Returned vector does not match")
	# Now we try Matricies
	inp=((0,0),(0,0))
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,((x,x+1),(x+2,x+3)))
	t=d.toListOfTuples(scalarastuple=True)
	ok=True
	for x in range(15):
	    if t[x]!=((x,x+1),(x+2,x+3)): ok=False
	self.failUnless(ok,"Returned matrix does not match")
	#Now 3-Tensors
	inp=(((0,0),(0,0)),((0,0),(0,0)))
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,(((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))))
	t=d.toListOfTuples(scalarastuple=True)
	ok=True
	for x in range(15):
	    if t[x]!=(((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))): ok=False
	self.failUnless(ok,"Returned 3-Tensor does not match")
	#Now 4-Tensors
	inp=((((0,0),(0,0)),((0,0),(0,0))),(((0,0),(0,0)),((0,0),(0,0))))
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,((((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))),(((9+x,9+x+1),(9+x+2,9+x+3)),((9+x+4,9+x+5),(9+x+6,9+x+7)))))
	t=d.toListOfTuples()
	ok=True
	for x in range(15):
	    if t[x]!=((((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))),(((9+x,9+x+1),(9+x+2,9+x+3)),((9+x+4,9+x+5),(9+x+6,9+x+7)))): ok=False
	self.failUnless(ok,"Returned 4-Tensor does not match")	
	
    # This test sets values then gets them.
    # Strictly speaking, it is not a complete test because it cannot tell if functions are broken
    # in a symmetric manner	
    def testToFromTupleTogether(self):
	fs=getTestDomainFunctionSpace(3,5)
	inp=0
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,x)
	ok=True
	for x in range(15):
	    if d.getTupleForDataPoint(x)!=(x,): ok=False
	self.failUnless(ok,"Returned scalar does not match")
	inp=(0,0)
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,(2*(x/2),(2*(x/2)+1)))
	ok=True
	for x in range(15):
	    if d.getTupleForDataPoint(x)!=(2*(x/2),2*(x/2)+1): ok=False
	self.failUnless(ok,"Returned vector does not match")
	# Now we try Matricies
	inp=((0,0),(0,0))
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,((x,x+1),(x+2,x+3)))
	ok=True
	for x in range(15):
	    if d.getTupleForDataPoint(x)!=((x,x+1),(x+2,x+3)): ok=False
	self.failUnless(ok,"Returned matrix does not match")
	#Now 3-Tensors
	inp=(((0,0),(0,0)),((0,0),(0,0)))
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,(((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))))
	ok=True
	for x in range(15):
	    if d.getTupleForDataPoint(x)!=(((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))): ok=False
	self.failUnless(ok,"Returned 3-Tensor does not match")
	#Now 4-Tensors
	inp=((((0,0),(0,0)),((0,0),(0,0))),(((0,0),(0,0)),((0,0),(0,0))))
	d=Data(inp,fs,True)
	for x in range(15):
	    d.setValueOfDataPoint(x,((((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))),(((9+x,9+x+1),(9+x+2,9+x+3)),((9+x+4,9+x+5),(9+x+6,9+x+7)))))
	ok=True
	for x in range(15):
	    if d.getTupleForDataPoint(x)!=((((x,x+1),(x+2,x+3)),((x+4,x+5),(x+6,x+7))),(((9+x,9+x+1),(9+x+2,9+x+3)),((9+x+4,9+x+5),(9+x+6,9+x+7)))): ok=False
	self.failUnless(ok,"Returned 4-Tensor does not match")		
	
	
if __name__ == "__main__":
    if getMPISizeWorld() == 1: 
       suite = unittest.TestSuite()
       suite.addTest(unittest.makeSuite(DataAccessTestCase))
       s=unittest.TextTestRunner(verbosity=2).run(suite)
       if not s.wasSuccessful(): sys.exit(1)
    else:
        print "run_data_access is not executed as more than one processor is used."
