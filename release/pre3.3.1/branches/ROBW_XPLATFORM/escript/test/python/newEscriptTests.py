
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""

import sys
import unittest

import numarray

from esys.escript import *

class DataTestCase(unittest.TestCase):

   def powTest(self):
      a=escript.Data([1,2],3,3)
      print a
      a=a**3
      print a

   def dataOperationTest(self):
      #
      # Create expanded data with a vector at each node
      myFuncSpac=escript.FunctionSpace()
      myVector=numarray.array([[1,2],[3,4]])
      myData=escript.Data(myVector,myFuncSpac,True)
      myData2=myData
      #
      # Test operator+
      myData3=myData+myData2
      print myData3

   def binaryOperationTest(self):
      a=escript.Data([1,2])
      print a
      a-=1
      print a
      
   def unaryOperationTest(self):
      constData=escript.Data([1,2])
      result=sin(constData)
      print result

   def dataExpansionTest(self):
      constData=escript.Data([1,2])
      print constData
      constData.expand()
      print constData
#
#     check expansion occurs automaticaly when needed
      constData=escript.Data(1)
      myFuncSpac=escript.FunctionSpace()
      expandedData=escript.Data(2,myFuncSpac,True)
#
#     this operation should cause automatic expansion
      try:
         resultData=constData+expandedData
         assert (False,'Failed shape mismatch exception test.')
      except Exception, e:
         print e
         print 'Passed shape mismatch exception test.'
#
#     do the expansion test
      expandedData=escript.Data(2,myFuncSpac,True)
      resultData=constData+expandedData
      print resultData

   def runTest(self):
      self.binaryOperationTest()
      self.unaryOperationTest()
      self.dataExpansionTest()
      self.dataOperationTest()
      #
      # Create constant data
      myData=escript.Data(1)
      myFuncSpac=escript.FunctionSpace()
      #
      # Create expanded data with a vector at each node
      a=numarray.array([1,2])
      myData2=escript.Data(a,myFuncSpac,True)
      try:
         myData2+='fred'
         assert(False, 'Failed illegal argument to operator+= test.')
      except Exception, e:
         print e
         print 'Passed illegal argument to operator+= test.'
#
#     Test adding a scalar to a 1 dimensional Data
      myData3=escript.Data(a,myFuncSpac,True)
      myData3+=2.0
      print myData2
      myData2+=2.0
      print myData2
      myData+=1
      try:
         myData+=myData2
         assert(False, 'Failed shape mismatch on += exception test.')
      except Exception, e:
         print e
         print 'Passed shape mismatch on += exception test.'
      try:
         myData3=myData+myData2
         assert(False, 'Failed currently illegal operation test.')
      except Exception, e:
         print e
         print 'Passed currently illegal operation test.'
      myData3=escript.Data(2.1,myFuncSpac,True)
      print myData3
      myData3+=1.0
      print myData3
      try:
         myData3+=a
         assert (False, 'Failed shape mismatch exception test.')
      except Exception, e:
         print e
         print 'Passed shape mismatch exception test.'
      newData=myData3+1.0
      print newData
      print myData3
      try:
          myData4=escript.Data(3.2,myFuncSpac,True)
          assert (False, 'Failed 0 size exception test.')
      except Exception, e:
          print e
          print 'Passed 0 size exception Test'
          
      myData5=escript.Data(a,myFuncSpac,True)
      myData5=myData5+3.1
      print myData5
 
dataTestCase=DataTestCase()
dataTestCase.runTest()

