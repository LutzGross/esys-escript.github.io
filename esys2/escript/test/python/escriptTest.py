import sys
import unittest
import os

esys_root=os.getenv('ESYS_ROOT')
sys.path.append(esys_root+'/finley/lib')
sys.path.append(esys_root+'/escript/lib')
sys.path.append(esys_root+'/escript/py_src')


import escript
import finley
import numarray

class escriptTestCase(unittest.TestCase):
  def testPow(self):
    """Test the pow function."""
    myFuncSpac=escript.FunctionSpace()
    myVector=numarray.array([[1,2],[3,4]])
    myData=escript.Data(myVector,myFuncSpac,True)
    print myData**3
  def testFunctionSpace(self):
      """Test the creation of FunctionSpace objects."""
      print
      mesh=finley.Brick(1,1,1,1,1.,1.,1.,1,1,1,1,1)
      cFunc=escript.ContinuousFunction(mesh)
  def testDataOperations(self):
    """Test the operations that may be performed on Data."""
    myFuncSpac=escript.FunctionSpace()
    myVector=numarray.array([[1,2],[3,4]])
    myData=escript.Data(myVector,myFuncSpac,True)
    myData2=3+myData
    print myData2
    myList=[[1,2],[3,4]]
    myData3=myList+myData
    print myData3
    myData3=myData+myData2+myData3
    print myData3
    myData4=myList-myData
    print myData4
    myData5=0-myData
    print myData5
    myData6=1/myData
    print myData6
    myData7=5*myData
    print myData7
    myData8=(3*myData)/3+2*myData-2*myData
    print myData8
    myData9=myData.sin()
    print myData9
    myData9=myData.cos()
    print myData9
    print myData9.wherePositive()
    print myData9.whereNegative()
#
#   Lsup test for DataExpanded
    assert (myData.Lsup()==4)
    assert (myData.inf()==1)
    assert (myData.Lsup()==myData.sup())
    assert ((myData*-1).sup()!=(myData*-1).Lsup())
#
#   Lsup test for DataConstant
    myData10=escript.Data(myVector,myFuncSpac,False)
    myData11=-1.0*myData10
    assert(myData10.Lsup()==myData11.Lsup())
    
suite=unittest.TestSuite()
suite.addTest(unittest.makeSuite(escriptTestCase))
unittest.TextTestRunner(verbosity=2).run(suite)




