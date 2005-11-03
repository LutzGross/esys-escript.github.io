# $Id$

"""
Test suite for the util.py module

The tests must be linked with a function space class object in the setUp method:

   from esys.finley import Rectangle
   class Test_UtilOnFinley(Test_Util):
       def setUp(self):
           self.functionspace = FunctionOnBoundary(Rectangle(10,10,2))
   suite = unittest.TestSuite()
   suite.addTest(unittest.makeSuite(Test_UtilOnFinley))
   unittest.TextTestRunner(verbosity=2).run(suite)

The test assumes that that the functionspace has samples with tags equal to an different from 1.
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision"
__date__="$Date$"

import unittest
import numarray
from esys.escript import *

class Test_Util(unittest.TestCase):
   tol=1.e-7
#=========================================================
#  constants
#=========================================================
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_1(self):
      val=kronecker(d=1)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0],1.0,"wrong value for (0,0)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_2(self):
      val=kronecker(d=2)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.failUnlessEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.failUnlessEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.failUnlessEqual(val[1,1],1.0,"wrong value for (1,1)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_3(self):
      val=kronecker(d=3)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.failUnlessEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.failUnlessEqual(val[0,2],0.0,"wrong value for (0,2)")
      self.failUnlessEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.failUnlessEqual(val[1,1],1.0,"wrong value for (1,1)")
      self.failUnlessEqual(val[1,2],0.0,"wrong value for (1,2)")
      self.failUnlessEqual(val[2,0],0.0,"wrong value for (2,0)")
      self.failUnlessEqual(val[2,1],0.0,"wrong value for (2,1)")
      self.failUnlessEqual(val[2,2],1.0,"wrong value for (2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_kronecker_domain(self):
      val=kronecker(d=self.functionspace)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val.shape,(self.functionspace.getDim(),self.functionspace.getDim()),"wrong shape.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_1(self):
      val=identityTensor(d=1)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0],1.0,"wrong value for (0,0)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_2(self):
      val=identityTensor(d=2)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.failUnlessEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.failUnlessEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.failUnlessEqual(val[1,1],1.0,"wrong value for (1,1)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_3(self):
      val=identityTensor(d=3)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0],1.0,"wrong value for (0,0)")
      self.failUnlessEqual(val[0,1],0.0,"wrong value for (0,1)")
      self.failUnlessEqual(val[0,2],0.0,"wrong value for (0,2)")
      self.failUnlessEqual(val[1,0],0.0,"wrong value for (1,0)")
      self.failUnlessEqual(val[1,1],1.0,"wrong value for (1,1)")
      self.failUnlessEqual(val[1,2],0.0,"wrong value for (1,2)")
      self.failUnlessEqual(val[2,0],0.0,"wrong value for (2,0)")
      self.failUnlessEqual(val[2,1],0.0,"wrong value for (2,1)")
      self.failUnlessEqual(val[2,2],1.0,"wrong value for (2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor_domain(self):
      val=identityTensor(d=self.functionspace)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val.shape,(self.functionspace.getDim(),self.functionspace.getDim()),"wrong shape.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_1(self):
      val=identityTensor4(d=1)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_2(self):
      val=identityTensor4(d=2)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
      self.failUnlessEqual(val[0,0,0,1],0.0,"wrong value for (0,0,0,1)")
      self.failUnlessEqual(val[0,0,1,0],0.0,"wrong value for (0,0,1,0)")
      self.failUnlessEqual(val[0,0,1,1],0.0,"wrong value for (0,0,1,1)")
      self.failUnlessEqual(val[0,1,0,0],0.0,"wrong value for (0,1,0,0)")
      self.failUnlessEqual(val[0,1,0,1],1.0,"wrong value for (0,1,0,1)")
      self.failUnlessEqual(val[0,1,1,0],0.0,"wrong value for (0,1,1,0)")
      self.failUnlessEqual(val[0,1,1,1],0.0,"wrong value for (0,1,1,1)")
      self.failUnlessEqual(val[1,0,0,0],0.0,"wrong value for (1,0,0,0)")
      self.failUnlessEqual(val[1,0,0,1],0.0,"wrong value for (1,0,0,1)")
      self.failUnlessEqual(val[1,0,1,0],1.0,"wrong value for (1,0,1,0)")
      self.failUnlessEqual(val[1,0,1,1],0.0,"wrong value for (1,0,1,1)")
      self.failUnlessEqual(val[1,1,0,0],0.0,"wrong value for (1,1,0,0)")
      self.failUnlessEqual(val[1,1,0,1],0.0,"wrong value for (1,1,0,1)")
      self.failUnlessEqual(val[1,1,1,0],0.0,"wrong value for (1,1,1,0)")
      self.failUnlessEqual(val[1,1,1,1],1.0,"wrong value for (1,1,1,1)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_3(self):
      val=identityTensor4(d=3)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val[0,0,0,0],1.0,"wrong value for (0,0,0,0)")
      self.failUnlessEqual(val[0,0,0,1],0.0,"wrong value for (0,0,0,1)")
      self.failUnlessEqual(val[0,0,0,2],0.0,"wrong value for (0,0,0,2)")
      self.failUnlessEqual(val[0,0,1,0],0.0,"wrong value for (0,0,1,0)")
      self.failUnlessEqual(val[0,0,1,1],0.0,"wrong value for (0,0,1,1)")
      self.failUnlessEqual(val[0,0,1,2],0.0,"wrong value for (0,0,1,2)")
      self.failUnlessEqual(val[0,0,2,0],0.0,"wrong value for (0,0,2,0)")
      self.failUnlessEqual(val[0,0,2,1],0.0,"wrong value for (0,0,2,1)")
      self.failUnlessEqual(val[0,0,2,2],0.0,"wrong value for (0,0,2,2)")
      self.failUnlessEqual(val[0,1,0,0],0.0,"wrong value for (0,1,0,0)")
      self.failUnlessEqual(val[0,1,0,1],1.0,"wrong value for (0,1,0,1)")
      self.failUnlessEqual(val[0,1,0,2],0.0,"wrong value for (0,1,0,2)")
      self.failUnlessEqual(val[0,1,1,0],0.0,"wrong value for (0,1,1,0)")
      self.failUnlessEqual(val[0,1,1,1],0.0,"wrong value for (0,1,1,1)")
      self.failUnlessEqual(val[0,1,1,2],0.0,"wrong value for (0,1,1,2)")
      self.failUnlessEqual(val[0,1,2,0],0.0,"wrong value for (0,1,2,0)")
      self.failUnlessEqual(val[0,1,2,1],0.0,"wrong value for (0,1,2,1)")
      self.failUnlessEqual(val[0,1,2,2],0.0,"wrong value for (0,1,2,2)")
      self.failUnlessEqual(val[0,2,0,0],0.0,"wrong value for (0,2,0,0)")
      self.failUnlessEqual(val[0,2,0,1],0.0,"wrong value for (0,2,0,1)")
      self.failUnlessEqual(val[0,2,0,2],1.0,"wrong value for (0,2,0,2)")
      self.failUnlessEqual(val[0,2,1,0],0.0,"wrong value for (0,2,1,0)")
      self.failUnlessEqual(val[0,2,1,1],0.0,"wrong value for (0,2,1,1)")
      self.failUnlessEqual(val[0,2,1,2],0.0,"wrong value for (0,2,1,2)")
      self.failUnlessEqual(val[0,2,2,0],0.0,"wrong value for (0,2,2,0)")
      self.failUnlessEqual(val[0,2,2,1],0.0,"wrong value for (0,2,2,1)")
      self.failUnlessEqual(val[0,2,2,2],0.0,"wrong value for (0,2,2,2)")
      self.failUnlessEqual(val[1,0,0,0],0.0,"wrong value for (1,0,0,0)")
      self.failUnlessEqual(val[1,0,0,1],0.0,"wrong value for (1,0,0,1)")
      self.failUnlessEqual(val[1,0,0,2],0.0,"wrong value for (1,0,0,2)")
      self.failUnlessEqual(val[1,0,1,0],1.0,"wrong value for (1,0,1,0)")
      self.failUnlessEqual(val[1,0,1,1],0.0,"wrong value for (1,0,1,1)")
      self.failUnlessEqual(val[1,0,1,2],0.0,"wrong value for (1,0,1,2)")
      self.failUnlessEqual(val[1,0,2,0],0.0,"wrong value for (1,0,2,0)")
      self.failUnlessEqual(val[1,0,2,1],0.0,"wrong value for (1,0,2,1)")
      self.failUnlessEqual(val[1,0,2,2],0.0,"wrong value for (1,0,2,2)")
      self.failUnlessEqual(val[1,1,0,0],0.0,"wrong value for (1,1,0,0)")
      self.failUnlessEqual(val[1,1,0,1],0.0,"wrong value for (1,1,0,1)")
      self.failUnlessEqual(val[1,1,0,2],0.0,"wrong value for (1,1,0,2)")
      self.failUnlessEqual(val[1,1,1,0],0.0,"wrong value for (1,1,1,0)")
      self.failUnlessEqual(val[1,1,1,1],1.0,"wrong value for (1,1,1,1)")
      self.failUnlessEqual(val[1,1,1,2],0.0,"wrong value for (1,1,1,2)")
      self.failUnlessEqual(val[1,1,2,0],0.0,"wrong value for (1,1,2,0)")
      self.failUnlessEqual(val[1,1,2,1],0.0,"wrong value for (1,1,2,1)")
      self.failUnlessEqual(val[1,1,2,2],0.0,"wrong value for (1,1,2,2)")
      self.failUnlessEqual(val[1,2,0,0],0.0,"wrong value for (1,2,0,0)")
      self.failUnlessEqual(val[1,2,0,1],0.0,"wrong value for (1,2,0,1)")
      self.failUnlessEqual(val[1,2,0,2],0.0,"wrong value for (1,2,0,2)")
      self.failUnlessEqual(val[1,2,1,0],0.0,"wrong value for (1,2,1,0)")
      self.failUnlessEqual(val[1,2,1,1],0.0,"wrong value for (1,2,1,1)")
      self.failUnlessEqual(val[1,2,1,2],1.0,"wrong value for (1,2,1,2)")
      self.failUnlessEqual(val[1,2,2,0],0.0,"wrong value for (1,2,2,0)")
      self.failUnlessEqual(val[1,2,2,1],0.0,"wrong value for (1,2,2,1)")
      self.failUnlessEqual(val[1,2,2,2],0.0,"wrong value for (1,2,2,2)")
      self.failUnlessEqual(val[2,0,0,0],0.0,"wrong value for (2,0,0,0)")
      self.failUnlessEqual(val[2,0,0,1],0.0,"wrong value for (2,0,0,1)")
      self.failUnlessEqual(val[2,0,0,2],0.0,"wrong value for (2,0,0,2)")
      self.failUnlessEqual(val[2,0,1,0],0.0,"wrong value for (2,0,1,0)")
      self.failUnlessEqual(val[2,0,1,1],0.0,"wrong value for (2,0,1,1)")
      self.failUnlessEqual(val[2,0,1,2],0.0,"wrong value for (2,0,1,2)")
      self.failUnlessEqual(val[2,0,2,0],1.0,"wrong value for (2,0,2,0)")
      self.failUnlessEqual(val[2,0,2,1],0.0,"wrong value for (2,0,2,1)")
      self.failUnlessEqual(val[2,0,2,2],0.0,"wrong value for (2,0,2,2)")
      self.failUnlessEqual(val[2,1,0,0],0.0,"wrong value for (2,1,0,0)")
      self.failUnlessEqual(val[2,1,0,1],0.0,"wrong value for (2,1,0,1)")
      self.failUnlessEqual(val[2,1,0,2],0.0,"wrong value for (2,1,0,2)")
      self.failUnlessEqual(val[2,1,1,0],0.0,"wrong value for (2,1,1,0)")
      self.failUnlessEqual(val[2,1,1,1],0.0,"wrong value for (2,1,1,1)")
      self.failUnlessEqual(val[2,1,1,2],0.0,"wrong value for (2,1,1,2)")
      self.failUnlessEqual(val[2,1,2,0],0.0,"wrong value for (2,1,2,0)")
      self.failUnlessEqual(val[2,1,2,1],1.0,"wrong value for (2,1,2,1)")
      self.failUnlessEqual(val[2,1,2,2],0.0,"wrong value for (2,1,2,2)")
      self.failUnlessEqual(val[2,2,0,0],0.0,"wrong value for (2,2,0,0)")
      self.failUnlessEqual(val[2,2,0,1],0.0,"wrong value for (2,2,0,1)")
      self.failUnlessEqual(val[2,2,0,2],0.0,"wrong value for (2,2,0,2)")
      self.failUnlessEqual(val[2,2,1,0],0.0,"wrong value for (2,2,1,0)")
      self.failUnlessEqual(val[2,2,1,1],0.0,"wrong value for (2,2,1,1)")
      self.failUnlessEqual(val[2,2,1,2],0.0,"wrong value for (2,2,1,2)")
      self.failUnlessEqual(val[2,2,2,0],0.0,"wrong value for (2,2,2,0)")
      self.failUnlessEqual(val[2,2,2,1],0.0,"wrong value for (2,2,2,1)")
      self.failUnlessEqual(val[2,2,2,2],1.0,"wrong value for (2,2,2,2)")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_identityTensor4_domain(self):
      val=identityTensor4(d=self.functionspace)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val.shape,(self.functionspace.getDim(),self.functionspace.getDim(),self.functionspace.getDim(),self.functionspace.getDim()),"wrong shape.")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_1(self):
      val=unitVector(i=0,d=1)
      self.failUnlessEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_2(self):
      val=unitVector(i=0,d=2)
      self.failUnlessEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
      self.failUnlessEqual(val[1],0.0,"wrong value for 1 in the 0 vector")
      val=unitVector(i=1,d=2)
      self.failUnlessEqual(val[0],0.0,"wrong value for 0 in the 1 vector")
      self.failUnlessEqual(val[1],1.0,"wrong value for 1 in the 1 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_3(self):
      val=unitVector(i=0,d=3)
      self.failUnlessEqual(val[0],1.0,"wrong value for 0 in the 0 vector")
      self.failUnlessEqual(val[1],0.0,"wrong value for 1 in the 0 vector")
      self.failUnlessEqual(val[2],0.0,"wrong value for 2 in the 0 vector")
      val=unitVector(i=1,d=3)
      self.failUnlessEqual(val[0],0.0,"wrong value for 0 in the 1 vector")
      self.failUnlessEqual(val[1],1.0,"wrong value for 1 in the 1 vector")
      self.failUnlessEqual(val[2],0.0,"wrong value for 2 in the 1 vector")
      val=unitVector(i=2,d=3)
      self.failUnlessEqual(val[0],0.0,"wrong value for 0 in the 2 vector")
      self.failUnlessEqual(val[1],0.0,"wrong value for 1 in the 2 vector")
      self.failUnlessEqual(val[2],1.0,"wrong value for 2 in the 2 vector")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_unitVector_domain(self):
      val=unitVector(i=0,d=self.functionspace)
      self.failUnless(isinstance(val,numarray.ArrayType),"wrong type of result.")
      self.failUnlessEqual(val.shape,(self.functionspace.getDim(),),"wrong shape.")
#=========================================================================
#   global reduction operations (these functions have no symbolic version)
#=========================================================================
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_float_rank0(self):
      arg=0.479077251703
      ref=0.479077251703
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_array_rank0(self):
      arg=0.352800421569
      ref=0.352800421569
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_array_rank1(self):
      arg=numarray.array([0.58364106865247445, 0.19224319360367659])
      ref=0.583641068652
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_array_rank2(self):
      arg=numarray.array([[0.67017748174064184, 0.27426357568545234, 0.99809374007262508, 0.0068333566694658288, -0.27942939334057559], [-0.41062296082648619, -0.036816602223561423, -0.50580074937952246, 0.93227848108675948, -0.061517050082725788], [0.36561750746233845, 0.41114839130078873, 0.52258027672142848, -0.16534372330544111, 0.20772668552253304], [0.821900382760401, -0.84255628577421948, -0.69396587198625026, -0.57918798921236458, -0.72171447032975466]])
      ref=0.998093740073
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_array_rank3(self):
      arg=numarray.array([[[0.058641541671277109, -0.90451294682583527], [-0.24359709498927606, -0.57748558070477163]], [[0.035804882991596898, 0.62841290637910441], [-0.28533080116748288, -0.97886508166774955]], [[0.31393622401598642, -0.43905852202615403], [-0.86251727012547685, 0.028980168735740941]], [[0.64796855283921229, -0.48583038861071492], [0.18113352051559328, -0.41145930584343637]], [[0.039393878628251944, 0.8768398562091233], [-0.17607723439655953, -0.88597401556177768]], [[-0.015710131346685419, -0.1460065558640582], [0.97739538148461858, -0.96991499683153215]]])
      ref=0.978865081668
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_array_rank4(self):
      arg=numarray.array([[[[-0.78982105495266985, -0.63220666806337955, 0.17512704833203396, 0.87905338838606228], [0.23198845220692799, 0.039972475068823288, 0.012169097173197807, 0.44188750570302382], [0.76550090905039703, 0.31693138180972436, 0.1053031163896907, -0.35110494947362092]], [[0.98438531435704557, 0.64830270240860366, 0.17255823042313656, -0.89376135594562944], [-0.24129437029981871, -0.142955980423916, 0.16557185988864043, 0.97182386578689162], [0.68369373227893937, -0.86170550438838256, 0.30120477894454822, 0.38702330237685523]]], [[[0.77139284396922037, 0.20032741426304668, 0.57845916425558697, -0.29867163908832151], [-0.068269410287581156, 0.5940891737261742, 0.076472990825278808, -0.099092183170674364], [-0.052727700907511776, 0.86303703635283835, -0.87561628108225542, 0.98706354430335175]], [[0.59243014649382819, 0.1550040875984271, -0.2755507051420949, -0.0013143184448647371], [0.49341486033505921, 0.47331310491746503, -0.79931467469262252, -0.90673470029976722], [-0.032268150780954796, 0.296035852616644, 0.51579882806939303, 0.46437108203184607]]], [[[-0.54940019219066349, 0.063961557315018069, 0.58950734587654585, -0.98334853918198539], [-0.3624096661573355, 0.41744569348555416, 0.30209950686844023, 0.51268273249278518], [0.18884359916930848, -0.71707023426140903, -0.30560603652072227, 0.50521867139895282]], [[0.48925658559264695, -0.22791551552340583, -0.0018172920946346593, -0.35038144063572618], [-0.92608233760416425, -0.58447575161042908, 0.6419293813902982, -0.9165521427783867], [0.32116313637555338, 0.64441081354246466, 0.57516697859586241, -0.30456483792192746]]]])
      ref=0.987063544303
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_constData_rank0(self):
      arg=Data(0.196366308048,self.functionspace)
      ref=0.196366308048
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_constData_rank1(self):
      arg=Data(numarray.array([-0.013183241788205846, 0.30081447346639489]),self.functionspace)
      ref=0.300814473466
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_constData_rank2(self):
      arg=Data(numarray.array([[0.5711180583492661, 0.70931947195628298, -0.66895311699347904, 0.96005746113679025, 0.73085528334644767], [-0.63904611175106618, 0.2843691804450883, 0.44023994297671054, 0.74230048057601272, 0.32582591826440876], [0.058605148358656045, 0.17856553839104938, 0.92397360311332144, -0.96449976222010503, -0.2936728605307215], [-0.54599501106213921, 0.76941479487476183, 0.071247548913826231, 0.19101147233175308, -0.1697403800152153]]),self.functionspace)
      ref=0.96449976222
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_constData_rank3(self):
      arg=Data(numarray.array([[[0.72501574945437675, 0.097442689963523454], [0.81136110912526904, 0.30019286779005516]], [[-0.49590270028453376, 0.89873757442194169], [-0.77574675514072333, 0.00090692035026496143]], [[0.30313499990678294, -0.22304437168798286], [0.26434595235717628, 0.56043553186944139]], [[-0.82536121216538372, 0.017266274277504934], [0.15087851023611853, 0.85422443819044291]], [[-0.85528228633213454, 0.21599153787828373], [-0.8320606477196939, 0.8359530516934528]], [[-0.32478507656272382, 0.11549647741760993], [-0.87438785398253049, 0.58454806081387956]]]),self.functionspace)
      ref=0.898737574422
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_constData_rank4(self):
      arg=Data(numarray.array([[[[-0.1675633776917147, 0.33827478137880718, -0.93890402023643449, 0.65020437341791437], [0.29507018847480526, 0.98622587753446878, 0.12652012726316597, -0.31134038684685694], [-0.046095794370747178, 0.52924578464459526, -0.6479404156998898, -0.50108997075395512]], [[-0.89461015899273821, -0.079360196866752331, 0.24950542226018069, 0.6689436082056277], [0.92392213781413735, 0.3873078097702356, 0.19593123983162242, -0.24092882483013001], [-0.64621424798001881, 0.9822743623774457, 0.89791841241748926, 0.61910184653693512]]], [[[-0.93993640130694156, 0.86452728798536005, 0.094496916350070848, 0.59825417382728907], [0.55042390382543216, 0.83625046124041091, -0.59865905280251042, 0.60081510989738351], [0.96300656863917622, 0.45676715577013183, 0.96765574240961594, 0.35217159943804499]], [[-0.44344990079738578, -0.62540931368504271, 0.47046830875624712, 0.56727920796684694], [0.68754074058706793, -0.20419202844112316, -0.0095491803785341389, 0.013299778291189002], [0.17824394120278897, -0.27714200037108694, -0.2616405339148673, -0.32155257707876661]]], [[[0.47113927793594357, -0.99306136743656892, 0.30468996581271934, -0.55558797016447881], [0.83216176170936151, 0.016003159554198287, 0.50607924358488665, -0.44441953149310631], [0.81919419287951278, -0.65849894919643503, 0.91897977494732008, -0.52338741357416407]], [[0.71408966944475138, -0.49347702658095161, 0.35676281330171133, 0.87268025092466872], [0.38401738326898771, -0.66323897612738114, 0.57309433517459518, 0.72101582669934583], [-0.0289954568811297, 0.55204032281174009, 0.51120867863750807, -0.67373936301915327]]]]),self.functionspace)
      ref=0.993061367437
      res=Lsup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_taggedData_rank0(self):
      arg=Data(0.860813503322,self.functionspace)
      arg.setTaggedValue(1,0.860813503322)
      res=Lsup(arg)
      ref=0.860813503322
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_taggedData_rank1(self):
      arg=Data(numarray.array([-0.54932912559284452, 0.29396676960376178]),self.functionspace)
      arg.setTaggedValue(1,[0.98025125990414441, -0.070257235982443378])
      res=Lsup(arg)
      ref=0.980251259904
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_taggedData_rank2(self):
      arg=Data(numarray.array([[-0.86895475746708128, 0.45103596542916824, 0.89202718384469271, 0.66165880808530297, 0.71929063210904221], [-0.054750345740449236, -0.26270085023397649, -0.44869339310367407, 0.84127602579890803, 0.4084040169910117], [-0.80258081555352101, 0.71946204694435134, -0.97606916814646971, -0.88087380297928397, 0.91540441306863141], [0.53133024472568935, -0.60623654813712635, 0.82280414663810242, 0.64010933901991374, 0.62566314353300356]]),self.functionspace)
      arg.setTaggedValue(1,[[-0.84852153765445437, 0.13244202632711666, -0.64133508534494599, -0.73706953458433633, 0.55834403408867184], [0.27998214461793847, 0.31446145164831063, -0.63410404784852048, 0.2813747329563423, 0.41221195047082393], [-0.79513090436643696, 0.92563876120768263, 0.80602538500705001, 0.21092919617246042, -0.21449414451693616], [0.50885151366468984, -0.53247783698745965, -0.98502684901017235, 0.36104863911630503, -0.68481313205160554]])
      res=Lsup(arg)
      ref=0.98502684901
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_taggedData_rank3(self):
      arg=Data(numarray.array([[[-0.0289179522354861, -0.53873332225554216], [-0.87313007813556509, -0.47147149825784584]], [[0.046403579054177468, -0.66499184318911042], [-0.14945300648197457, 0.33023752485562841]], [[0.73609028529612153, 0.62400582710031194], [0.18047782954118574, 0.98299132707347403]], [[0.97452943106570422, -0.80052218344822124], [0.90989474269184356, 0.74467116925414456]], [[-0.40975095375636039, 0.35721815590834538], [-0.023117827122894896, 0.38726163442133732]], [[0.35214474483480052, 0.79626235681759927], [0.072063982160859297, -0.13255981975702369]]]),self.functionspace)
      arg.setTaggedValue(1,[[[0.35540494212277429, -0.40452986468200347], [0.92646378475498059, -0.60976701230157482]], [[-0.17488076275939557, 0.1383489038535719], [0.87222102776136068, 0.05521388649844039]], [[-0.45846974731683399, 0.84645585780786292], [-0.36620222778926448, 0.8758026265447818]], [[0.55804586547774848, -0.19954715807059986], [0.51849302021923482, 0.29871500421422281]], [[0.98995968883285035, -0.78797081527577162], [0.3108075746688399, 0.5474101080348186]], [[-0.74670637823709085, 0.16925394395842575], [-0.76281911656771095, 0.79574461985041189]]])
      res=Lsup(arg)
      ref=0.989959688833
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_taggedData_rank4(self):
      arg=Data(numarray.array([[[[-0.17143955149003642, -0.0089504254973518105, 0.53089483578382191, 0.32381408252152899], [-0.77296292864801019, -0.089110721047862773, 0.46285355679515883, 0.62221785170960708], [-0.67474474335134915, 0.34682047265462157, 0.9384448074012548, 0.5682850087388498]], [[0.52652198841191855, 0.52394055634009762, 0.41923325950497858, -0.48507989905455706], [-0.90073796582592069, -0.40217266425438258, -0.60530063652424215, 0.68062718938448441], [-0.59931923693732347, -0.79549982795384744, -0.70714734772804722, -0.46042778371080284]]], [[[0.58538756755140686, 0.98385384505005846, 0.7777811719634411, -0.64306377174574281], [0.72961354558815694, 0.10696472171933968, -0.11372282342784068, 0.87929133681948879], [-0.67126196529672244, -0.64730190047646907, 0.64547629928395711, 0.50361974274373145]], [[0.96265942240931546, -0.20746026072477042, 0.47323657518133921, -0.78443796621025053], [0.61977887181220659, -0.0192018581010025, -0.0016015804221325425, 0.25446656696052594], [0.19964691572203019, -0.44122579240360293, 0.89836148991232134, -0.97914559157737457]]], [[[-0.32474221003830039, 0.50501185871734799, 0.081832847990893409, 0.49226411509256796], [-0.58561709012191865, -0.97753141368900409, 0.50702769958783778, 0.46965610524959978], [0.19394052487354463, 0.32118138849740641, 0.48348630138165749, -0.61570132061285632]], [[0.78997938799317668, 0.48729593728848108, 0.86690961213187001, -0.55317005853484491], [-0.38985400756166189, -0.79197087340853445, 0.150444446088422, 0.30366473850354492], [0.16673919050825758, -0.28616432413953641, -0.49042930009947883, -0.80964116966434485]]]]),self.functionspace)
      arg.setTaggedValue(1,[[[[-0.72676590268291097, 0.98420782971554899, -0.58004995296952444, 0.37649505647780179], [-0.36963117451708949, -0.38478644500667469, 0.1606599749645139, 0.26146427896482427], [-0.99755391430668583, 0.96243322443760793, -0.34748898506056713, 0.28223401802658166]], [[0.41892282572460227, -0.068327589700850844, -0.92249969532644394, -0.2927104302765704], [0.63237889769391709, -0.61446924102341649, -0.9271255632289408, 0.72693928120951368], [-0.099138333893530106, -0.93278471458000989, 0.16805036953472618, 0.13406769552186848]]], [[[0.1322308020971239, -0.15094779056740282, 0.48419178200868274, -0.90259173990902308], [0.088806733010250438, -0.44134645109664827, 0.50169033175317468, -0.16413576472992863], [0.10447947060273766, 0.59946651445651744, -0.28648625172498821, -0.26114646276357711]], [[-0.17647875332717788, -0.95243401465773969, 0.066994364736289391, 0.76072295812282875], [-0.29974152935779652, -0.87018574916912828, -0.40027227651920905, -0.27566894336852044], [-0.87505794257603342, 0.53786153286888583, -0.23579951775243324, 0.29461110217796826]]], [[[-0.031292782596848978, 0.19001451946176351, 0.51137483078731094, 0.35855090738394124], [-0.62796181019314523, -0.017622867812650655, -0.20994152673731148, 0.21972116995451207], [-0.53419638828850147, 0.61964526276926013, 0.83633801914948402, -0.22627427949817003]], [[-0.25275677187826617, 0.92174213140825789, 0.29387486254521544, 0.2851840648022741], [0.99521823294639589, 0.30976825827796484, 0.39585066725930163, -0.037512976967312373], [-0.0098417329760405181, -0.72834591016301697, -0.2368701950529164, -0.075161686057492183]]]])
      res=Lsup(arg)
      ref=0.997553914307
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_expandedData_rank0(self):
      arg=Data(0.907507663119,self.functionspace,True)
      arg.setTaggedValue(1,0.907507663119)
      res=Lsup(arg)
      ref=0.907507663119
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_expandedData_rank1(self):
      arg=Data(numarray.array([0.64842023599463228, 0.7950449048151853]),self.functionspace,True)
      arg.setTaggedValue(1,[-0.53356920308809608, -0.30392740367264159])
      res=Lsup(arg)
      ref=0.795044904815
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_expandedData_rank2(self):
      arg=Data(numarray.array([[0.72772455370935374, -0.47120637916343311, -0.82075870555018726, -0.13186541688079845, 0.14046405940683471], [0.10791450121551649, 0.043686402172190775, 0.4587164889316806, -0.16569163575870283, 0.64477762072726041], [0.62628246294878309, -0.46538827310659792, 0.58089235621217727, -0.745300901798017, -0.1072674226756638], [0.93074707226494824, 0.17195108914746116, 0.77645205158899833, -0.55814650711894975, -0.68929261213084403]]),self.functionspace,True)
      arg.setTaggedValue(1,[[-0.25856012658041871, -0.64132706411711782, -0.90934384682634128, 0.13076992241789931, 0.23456684296051011], [0.54052140856785269, 0.78868044275368643, 0.20419986484049479, 0.64782383948156319, 0.12883249735345115], [-0.44575026820636654, -0.86972644697707824, 0.74937006939719653, 0.64390867433141019, 0.57227950445890885], [-0.59430616226417832, -0.77932115906125854, 0.60864641730739155, 0.909083762068297, -0.5444504332265574]])
      res=Lsup(arg)
      ref=0.930747072265
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_expandedData_rank3(self):
      arg=Data(numarray.array([[[0.31609964814729086, -0.59322042077218184], [0.62671563648606643, -0.99924882354010092]], [[-0.81285219908299133, -0.91261945476639861], [-0.66394864058744174, -0.070011911370653657]], [[-0.4798784324091383, -0.017929635369934749], [0.87935995021589952, 0.73748462709583618]], [[-0.89673095768516986, 0.44255562883781563], [-0.33009427566326166, 0.89415170271508537]], [[-0.070411620428932897, 0.34854312339042304], [-0.54088506672687542, 0.57996368816677069]], [[0.98447862226498417, 0.31010343079927294], [0.18832525314748882, 0.46594491838516161]]]),self.functionspace,True)
      arg.setTaggedValue(1,[[[0.036725320623094637, -0.95874317021810596], [-0.66138152872168576, -0.52870418789420959]], [[0.59182952544611012, -0.31117513196914603], [0.17957160072325573, -0.93130315483187909]], [[0.33548100103066814, -0.6503677585469938], [-0.15995741665912955, -0.79138987042982367]], [[0.12353483100690976, -0.72197260504479188], [-0.35933752275788389, -0.46752695895022667]], [[0.7611306903818762, -0.88807332904594882], [0.91131651925077373, 0.81255438802194258]], [[-0.56892641899978824, 0.0010213385566093525], [0.40194539652124472, -0.006585080723547021]]])
      res=Lsup(arg)
      ref=0.99924882354
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_Lsup_expandedData_rank4(self):
      arg=Data(numarray.array([[[[0.4996033665120716, -0.074132654028136846, 0.2908837888281397, 0.79804117547820219], [0.25334186456470942, -0.10892196087517791, 0.90411872189757747, -0.21372334914190327], [0.65549677694560438, 0.97979738684954532, -0.29737844854710138, 0.93659843329597914]], [[0.22668444631484608, -0.92481913749465305, -0.9980815386745403, -0.022502768980955601], [-0.21769755240946398, 0.77216348666766876, -0.19843685166706204, 0.54270333879579558], [-0.11274856721131221, -0.29600869223659299, 0.1458222910080329, -0.83739782177046851]]], [[[0.86722932032155531, 0.39888432468517876, -0.8991679849590255, -0.088621935923834272], [-0.58464076321585412, -0.09564380294323116, 0.18232602464536307, 0.7910046931530843], [0.15923450234841652, -0.39331159996226872, 0.18298541750645669, 0.99889484861795563]], [[0.36793558813747418, -0.64593764280975363, 0.048503028175158613, -0.8304805399530264], [0.072019074767407432, -0.066883567381289311, -0.55849542620276127, -0.32521841292447484], [0.83256632210175896, -0.52124955617445723, -0.0047287521832242163, 0.84184001532121422]]], [[[-0.81375499823702158, 0.12901434959756353, -0.51500727423485215, 0.52626362435118912], [-0.47602620905811044, 0.81525173294213982, 0.023145745277130203, 0.5818744103097242], [-0.26074066195347489, 0.62737248392572886, 0.24246935026650718, 0.86155298917514145]], [[0.40180649524587109, -0.13468267908829512, -0.66576279256576543, -0.97664336021962694], [-0.81183732113700424, -0.10477655696019839, -0.90212494842448732, 0.50784279020015499], [0.29352929816605089, 0.10640245030222295, -0.16640870997460122, 0.91371366707232826]]]]),self.functionspace,True)
      arg.setTaggedValue(1,[[[[0.6465819200939884, 0.21755363919340609, 0.73973535907059662, -0.36114669719855241], [0.16087051863228452, -0.96295075570141986, 0.93283344105282717, 0.8498346294196879], [-0.62871563035312805, 0.028501977366871101, -0.76183996578150004, -0.42396762024023338]], [[0.45139869884828876, 0.9085099092003921, 0.90516675545818392, -0.2797423591331305], [-0.012176715080714828, 0.40935600076669765, -0.010156279663344314, -0.45527372677880185], [-0.56697253731600106, -0.88076957837377901, -0.43708368779823448, -0.98115855535334329]]], [[[-0.6812131434679467, -0.75007359636996074, 0.52195871968240559, 0.74207747673309732], [0.53576769134014213, -0.19432873205999046, -0.87970756195132904, -0.36970944422105911], [0.18377651984390431, -0.30093639418244189, 0.30640551056952825, -0.95779743159891284]], [[0.3069655071293016, 0.42532244942656305, 0.27182877898608804, 0.89926151593228765], [-0.94227360921249192, 0.17309985459832045, -0.067341615594060267, -0.24017528169767255], [0.72377020653147883, -0.60287365881872312, 0.17612550848074338, -0.89499327987049915]]], [[[-0.8985020338092089, -0.27805317471704494, -0.096352572725887375, 0.26107376060313836], [-0.98264038134460852, -0.40101944215897967, 0.80787105663414827, -0.91046803072373206], [-0.11056024039811629, -0.35146855329949345, 0.62639019941990948, 0.029258586953160748]], [[0.5190564210634494, 0.25178673168519605, -0.095466912631134937, -0.66223610619416462], [0.86572944270431917, -0.0070686656495086986, -0.56404011740509774, 0.5156978630437381], [0.15106963134402651, 0.12900511640038159, 0.6022471822104567, 0.48643914768022012]]]])
      res=Lsup(arg)
      ref=0.998894848618
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_float_rank0(self):
      arg=0.870743835413
      ref=0.870743835413
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_array_rank0(self):
      arg=0.469212543992
      ref=0.469212543992
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_array_rank1(self):
      arg=numarray.array([0.8163530200305178, 0.7844191729334391])
      ref=0.816353020031
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_array_rank2(self):
      arg=numarray.array([[-0.52338785100595131, 0.72523140545134046, -0.23883623210726745, 0.29329903553233394, 0.77300897701720128], [0.56646202903477705, -0.67833617682539948, 0.71280801753916911, 0.108973189514324, -0.86675353843929437], [0.37080584846118247, 0.61817009100957776, -0.20780655998890807, 0.085315295987765438, -0.73527023026482174], [-0.97586476277122935, 0.14501540684207481, 0.57375473938632338, 0.08516777342367865, -0.22644451314946301]])
      ref=0.773008977017
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_array_rank3(self):
      arg=numarray.array([[[0.16933183602716984, -0.42964457496769226], [-0.63714228263554573, -0.28513481547494179]], [[0.72479530443714335, -0.097501515360919111], [-0.28611653510816737, -0.58370472731498535]], [[-0.18432738416022554, 0.79010596522300558], [-0.65367387441910196, 0.90861652898349976]], [[0.56004415223670123, 0.20178156913861489], [0.90730594499457595, 0.91196305151516754]], [[-0.46179421349599847, -0.54555869532019163], [0.36014998847058499, -0.70585188726413306]], [[0.49988705904335418, -0.52181171665742077], [0.14475259007357621, -0.94336078709637383]]])
      ref=0.911963051515
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_array_rank4(self):
      arg=numarray.array([[[[0.085213773984220698, -0.9837351364458633, -0.23428780807935823, -0.54350985593703971], [-0.021654619752513593, -0.58267295960777465, -0.95135334805945893, -0.82172163939108089], [0.92970460708060565, 0.12721615089598348, 0.36243089042630872, 0.50434042290503855]], [[0.20617957875725845, 0.13657289785878368, 0.7942986246389756, -0.92641374730818571], [0.30431241116181762, -0.13563881945622858, 0.37921852742514695, -0.39275408991812211], [-0.016546769015328033, 0.50932041928102878, -0.27343457607173893, -0.0076289641375255624]]], [[[0.97189015970083137, -0.71286035174080009, 0.10143028964095313, -0.41856292303563181], [-0.97563458262665792, 0.79724078659701769, -0.70932749973904685, 0.018497784992804522], [-0.86374122662275021, 0.047715471334789816, -0.95453593058418518, 0.54562170290354794]], [[0.40249406070198157, -0.54609432046574469, -0.22682900899629166, 0.98394939138178539], [0.11049172557176901, 0.42172241721325388, 0.71050000578192951, 0.35353993854713206], [0.35412886303451896, -0.98169699399727617, 0.04758881049644037, 0.96971205948237493]]], [[[0.44672925866052249, -0.51476498049696828, 0.56442823935318875, -0.39769248164928039], [-0.40340965812893304, -0.87947712857546945, 0.55927022788356706, -0.89941016574892707], [-0.43878304559423431, 0.20915357555548764, -0.76548553334601799, -0.67202522557876954]], [[-0.56749721271516318, -0.10244683680777245, 0.17727779142251943, -0.57219284260940473], [-0.17044718853145446, 0.91117482665936023, -0.30039592703806584, -0.73813808369358713], [0.63771084365736663, -0.61427668096170129, 0.34365587989446378, -0.11877233448104674]]]])
      ref=0.983949391382
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_constData_rank0(self):
      arg=Data(0.165371505685,self.functionspace)
      ref=0.165371505685
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_constData_rank1(self):
      arg=Data(numarray.array([-0.89603386749185288, -0.68712608295212729]),self.functionspace)
      ref=-0.687126082952
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_constData_rank2(self):
      arg=Data(numarray.array([[0.60272743728014655, 0.26067131689446787, -0.48488892003697837, -0.54328378217335027, -0.96627165443113894], [0.38861396631681999, -0.14210447298121753, -0.84480805358330624, -0.25397651427390566, 0.25670041011662192], [-0.062982523786134337, -0.149708363807598, -0.63332360725934489, -0.49175302564011525, -0.97647588301352473], [0.52022334705669038, -0.69717039787412727, -0.28284586409251511, 0.99642563937215467, -0.67058148736338885]]),self.functionspace)
      ref=0.996425639372
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_constData_rank3(self):
      arg=Data(numarray.array([[[-0.98552936119306023, -0.58995212270861552], [0.51743177430155907, 0.68576837981065508]], [[-0.61618414432919089, -0.12325580790677049], [0.32387395300714172, -0.95456083598524333]], [[0.89779642579517049, 0.98676270760314266], [0.71959629907181966, -0.9949078223284622]], [[-0.81547040114414271, 0.10033634427970006], [-0.21591232734408217, -0.68608679705274822]], [[0.30423138886571999, 0.34122142527426802], [-0.4344532377066066, -0.31076903154305779]], [[-0.46111628105416602, -0.18371998359850483], [0.63606993796228117, -0.10703087143670587]]]),self.functionspace)
      ref=0.986762707603
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_constData_rank4(self):
      arg=Data(numarray.array([[[[0.89488524023952776, 0.8669396968091807, -0.45769331537553315, -0.89738195395075349], [0.39077366764066168, -0.71075932241646922, -0.51096889323130612, 0.87130290122807663], [0.17079807940685177, -0.6198922248677714, -0.41261700237404653, -0.30627765803812368]], [[-0.0659839868001868, -0.54367942190111385, 0.79254440140135607, -0.28895269724305006], [0.2554732744127961, -0.0076696085190677277, 0.6846601760873452, 0.38696598742090393], [-0.77125424651939789, 0.63909999575689591, -0.87840142433528379, 0.41711809974302594]]], [[[-0.99322035791310692, 0.27097830567855352, -0.4253855401144222, 0.15768186455727529], [-0.49181115516922302, -0.36126829134959304, 0.52357599944776667, 0.91209852597809005], [0.069076441159411361, -0.18292686281510551, -0.6497679800515983, 0.022610374934600719]], [[0.28755759348156507, -0.08281906224050295, 0.76036900801429907, 0.54802231074240826], [-0.033682724326368874, -0.7414032372095134, -0.86699767465780231, 0.40592904057808044], [0.51593363738292841, -0.72087130860034332, 0.35672334112134374, -0.090721746979026463]]], [[[-0.54866684145444511, -0.96738751715836968, -0.21201752332220436, -0.099425492405464277], [-0.76528700517307313, -0.85955622688708644, -0.10492266626403257, 0.69311319310724762], [-0.33886276086664902, -0.6846128323156393, 0.05873264876508566, 0.88498228323799433]], [[-0.28404277561384639, -0.63570388064518468, -0.67775264818658387, 0.20825454125346576], [-0.84788984114351473, 0.037932422136330635, 0.021981819447397299, -0.2709264612684219], [-0.64072476278735468, 0.46126191894728197, -0.37456096950035489, 0.85599593427453957]]]]),self.functionspace)
      ref=0.912098525978
      res=sup(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_taggedData_rank0(self):
      arg=Data(0.649634736435,self.functionspace)
      arg.setTaggedValue(1,0.649634736435)
      res=sup(arg)
      ref=0.649634736435
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_taggedData_rank1(self):
      arg=Data(numarray.array([-0.91775874675364899, -0.44518660348335226]),self.functionspace)
      arg.setTaggedValue(1,[-0.91030878048996744, -0.36380755471992954])
      res=sup(arg)
      ref=-0.36380755472
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_taggedData_rank2(self):
      arg=Data(numarray.array([[-0.086523811957357255, -0.42774060865319008, 0.35170422393945833, -0.24772487756509576, 0.85056214851195744], [-0.31877386278938702, 0.28773006836605641, 0.44932528141129979, -0.56724416115855569, 0.80633095352264816], [-0.53652434245562564, -0.26697043576159984, -0.88767305488188519, -0.029691195696610828, 0.67899103041623876], [0.92484508322389836, 0.18625473102022339, -0.27285903116359256, 0.63921542460538938, -0.9221199231145456]]),self.functionspace)
      arg.setTaggedValue(1,[[0.65738567253805136, -0.6778218395330815, 0.40806699669092361, 0.34540048412849589, -0.11704616494950493], [0.38512651510421825, 0.74221788961938562, 0.95314896284964967, -0.040871082359481337, 0.73045537711619035], [0.10490367249326416, -0.24457205097868751, 0.23569203929084925, -0.4833470179537227, 0.13727107062761412], [-0.34956075762939753, -0.40510846111177878, -0.60113099618774268, -0.8694743269747125, 0.8300938895447072]])
      res=sup(arg)
      ref=0.95314896285
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_taggedData_rank3(self):
      arg=Data(numarray.array([[[0.50181406381755478, 0.12696291412348604], [0.45580252908683372, -0.69381529096260919]], [[-0.44711719663105876, -0.49939480373343814], [0.49163505815095565, 0.70755602744123269]], [[-0.32952519817455106, -0.56163922933132882], [0.020194551395573912, 0.74013090992747421]], [[-0.13639620097611793, -0.5402749306510144], [-0.71348777995694368, -0.07149424731352183]], [[-0.81298057712035066, -0.12510197890662789], [-0.30874509775629533, -0.58120893128076712]], [[0.86654409796596377, -0.50673775089683915], [-0.12239574780538409, 0.81691821472857717]]]),self.functionspace)
      arg.setTaggedValue(1,[[[-0.57151921079320256, 0.95258628100636189], [-0.77681734612402287, 0.95978727375688866]], [[0.76493561669286803, 0.51125225923442486], [-0.24392383855124655, 0.014944647833669666]], [[-0.93836576145690231, -0.044479851632975853], [-0.30511938835638897, 0.0091738439461943599]], [[0.71921984284603702, -0.22105010374862499], [-0.78589399511594116, -0.8895142672649694]], [[0.6735135460868733, -0.56646772685337399], [0.73605625715117484, -0.68735959525940049]], [[0.38440898374441201, 0.87186026279634277], [-0.59320035048324327, -0.87430848491656854]]])
      res=sup(arg)
      ref=0.959787273757
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_taggedData_rank4(self):
      arg=Data(numarray.array([[[[0.65078702842768221, -0.73900875462993776, -0.65984653996722176, -0.97539191158251093], [0.99720731720088596, -0.19242629799409516, -0.16320360691406655, -0.043083875914333492], [0.59108757025274139, 0.88522336612612706, 0.74088460755604268, 0.80676177752104627]], [[0.82498945289642323, 0.23110039892583623, -0.13089698664065774, -0.400241064232296], [-0.67523171276117, -0.72699269118844501, 0.53208807038398409, -0.49611668096456474], [-0.20482996861953673, 0.83614314521551858, -0.82232394050773272, -0.71904026227092044]]], [[[0.34003140388599551, 0.31905118515154363, 0.51883291903272832, -0.76898011700894209], [-0.00054352334236540401, -0.42459032663080087, 0.72331023817772211, -0.50647033418441123], [0.76195696335445318, -0.22666337528546276, -0.40740135145605061, -0.81651333897992373]], [[-0.9819499586934739, -0.95616911480063527, -0.2900474363658565, -0.16636609538100222], [0.66435123810546681, -0.30523563374854135, -0.21355210817886849, -0.74243246288718034], [0.82105586828161425, -0.93524621329362057, 0.36308161224720026, 0.76840492117538539]]], [[[-0.22872302699935965, 0.047300841535470761, -0.93216772922157576, -0.40541971639813301], [0.22921538079079262, 0.1958937804621288, 0.27198494374967619, 0.55888236453433282], [0.97179646000620856, 0.53691052552944973, -0.13695340427959324, -0.10499588580374764]], [[-0.4821448283439782, -0.96593591454069139, -0.95284453814535297, -0.35311046977499583], [0.4857495870271038, -0.66584818036412852, -0.04957796396188785, 0.28223147859593767], [-0.2171789936962405, -0.016212404435609562, 0.16258357268143042, -0.6781166044152207]]]]),self.functionspace)
      arg.setTaggedValue(1,[[[[0.80768589260978896, -0.22702941364663487, -0.72896606223385407, 0.11413326357409237], [0.40571664072592251, -0.9311405487001907, 0.80450361552618688, 0.56480640933432991], [-0.33782609968052979, 0.39512515837757123, 0.73591694462004398, 0.24768959674737778]], [[0.56569618270183164, -0.93779341659685, 0.64969642196738708, -0.77336556098096976], [0.41311175212178153, 0.056953746906872826, 0.25300968955971204, -0.35019321925911262], [-0.8863122403302417, -0.89705763853428677, 0.3060047535556778, -0.92592167036041095]]], [[[0.45649019646929379, -0.29843125838513096, -0.20714508244855367, 0.246705639144563], [-0.32477747411703084, 0.30488751585973306, 0.53390827733820756, -0.84339943975046583], [-0.12373671376305295, -0.1640054521913612, -0.87414144472044897, -0.0021211693404443732]], [[0.38273461960506805, 0.70999974995969617, 0.22361687978370237, -0.098549468178965371], [0.81724211804899904, -0.88965787513620009, -0.39375673885748119, 0.69490920308416371], [-0.65400552410197754, -0.82376412930222931, 0.046545365304690778, -0.21012949605343434]]], [[[-0.27643919361415081, -0.44880727610691973, -0.57364607821939151, -0.14013355075911993], [-0.99302452440223621, 0.70400083788517742, 0.29183091261608896, 0.57457780218190213], [0.20084128112884403, 0.98904695078892235, 0.87503585272015294, -0.26131933340055569]], [[0.94633204265993198, -0.73295197510079446, 0.56975658926329098, -0.83390352955122538], [0.2617682886960544, 0.14649180808291562, -0.29972426982703726, -0.015848496464521356], [-0.96680270201151153, -0.79982829582732196, -0.29552300849179347, 0.66620264190912515]]]])
      res=sup(arg)
      ref=0.997207317201
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_expandedData_rank0(self):
      arg=Data(0.842459260157,self.functionspace,True)
      arg.setTaggedValue(1,0.842459260157)
      res=sup(arg)
      ref=0.985691469761
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_expandedData_rank1(self):
      arg=Data(numarray.array([0.47024430020620023, -0.40410868427962021]),self.functionspace,True)
      arg.setTaggedValue(1,[0.34568516056640197, -0.43342673126146103])
      res=sup(arg)
      ref=0.470244300206
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_expandedData_rank2(self):
      arg=Data(numarray.array([[-0.88351146515427814, 0.87526615684362929, 0.80242059303502922, -0.6048006248682789, 0.42108402102305953], [0.11385226207304888, -0.66232044805160473, 0.69874162972474063, 0.22449470583280284, -0.19986159203927545], [0.17102700533232951, -0.9790459454968905, 0.4374092200111821, -0.4343902463997138, 0.72993175063985705], [-0.42636794385375065, -0.88420160799308434, -0.073357336228458081, -0.96071213693656698, -0.28591564459422769]]),self.functionspace,True)
      arg.setTaggedValue(1,[[-0.66507554751323883, 0.90956185314332916, 0.39373039994051529, 0.20185922272970869, 0.56272045770025914], [-0.81603670050480082, -0.98308666677421708, 0.081972418155825233, 0.98933334281872432, -0.67280110519862579], [-0.67384516784043069, 0.19168471924963759, -0.59938113619504896, 0.22266455997006829, -0.93324865777959265], [0.93134589933544842, -0.10311385886401503, -0.2331509870020978, 0.37315747180467018, 0.73551730577813057]])
      res=sup(arg)
      ref=0.989333342819
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_expandedData_rank3(self):
      arg=Data(numarray.array([[[-0.69867440642271283, -0.37489565319806606], [-0.66485463581725601, 0.18238313710104448]], [[-0.78388494323459157, 0.24027472943325923], [-0.4377158808146262, -0.34262284842949864]], [[-0.8539094490313166, 0.62109328106735684], [-0.20715549218162899, -0.33320905765908693]], [[-0.25702886172082651, 0.94651536183641283], [-0.2229696935051404, 0.84771132395539794]], [[0.52244706935442076, -0.89344888833386804], [0.064322796922618108, 0.36966382541443532]], [[0.55710175300577691, -0.22780561634551688], [0.89548548132750594, -0.77561238906399077]]]),self.functionspace,True)
      arg.setTaggedValue(1,[[[-0.86933240372621312, 0.32078825553642276], [0.36873095289245739, 0.20578150283627239]], [[-0.13143164944725427, 0.89920065296927199], [0.79295913398164197, 0.78838875458878954]], [[0.068510791358111334, 0.87754319283313054], [-0.96880448620091553, -0.012058603734139028]], [[0.44680342907934811, -0.52293412412648466], [0.51117158355443149, -0.6794979840301234]], [[0.55644022070667498, 0.082838767459920026], [-0.64861425420762142, -0.20781169747814943]], [[0.52302454039265078, -0.71078210239352546], [0.67348959224612859, 0.18668184392186915]]])
      res=sup(arg)
      ref=0.946515361836
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_sup_expandedData_rank4(self):
      arg=Data(numarray.array([[[[-0.97992610933200397, -0.25896245628971326, -0.70483438344031946, -0.99381143197943023], [-0.94915117581748709, -0.41900816047650036, -0.73354223723879386, 0.75059333884969881], [-0.76848095939834238, -0.75674819943701244, 0.065477641668864495, -0.42345850130975426]], [[-0.86758200423510234, -0.91179084363212071, -0.95127824910458103, 0.44686792515820439], [0.24843166983163689, -0.78829756520109728, -0.29421842880871774, 0.081312723661290498], [-0.73676127398786795, 0.91442315487631975, -0.64075473615377643, -0.68228346762450576]]], [[[-0.52913624371899504, 0.18611492696208209, 0.87647965087074042, 0.81733477764270401], [-0.45994626037000352, -0.72068350777641998, -0.67722337882139305, 0.53688427407036721], [-0.62155086776441348, -0.16407031198819877, 0.27545162813408264, 0.82180316351563087]], [[0.86221555518143278, -0.038942303505028031, 0.99073092777336114, -0.59207542384567113], [-0.42324210154238218, 0.76649013040311842, 0.86242470461616949, 0.49956796372115009], [-0.3983839684487791, 0.50100130531149967, -0.61285700166494195, -0.77581662400471862]]], [[[0.393847951204221, -0.72453857578029801, 0.97194727576353768, -0.95039782864739597], [-0.33546269250215266, 0.4831646699512071, -0.36084136685304102, -0.21835205500655297], [-0.82932637046105562, 0.70351683765540218, -0.058619290451321637, -0.2054802590132303]], [[0.94447147902776929, -0.42312399627562503, -0.11683182283655591, -0.40970199303413191], [-0.057628087834810326, -0.92418858637343759, 0.17632737701884516, -0.45930384359491305], [-0.47612419752878798, -0.96256404924809247, 0.10415265869578216, -0.63202348546194909]]]]),self.functionspace,True)
      arg.setTaggedValue(1,[[[[-0.97434762405813102, 0.7014239623245393, 0.62448079624150288, -0.37769806924712546], [-0.85193729675303742, -0.77375392434251022, 0.36142120851055592, 0.26580596564769365], [-0.87276378893603579, 0.38973819311591718, 0.98354255063199636, -0.66795317131402343]], [[0.90231744676982428, -0.27585106370360957, -0.32108129691288001, 0.60440834831736634], [-0.045905880902579188, -0.29199624499256793, 0.025945636614535861, -0.098321374439206854], [-0.56344512008582659, -0.3766256629210043, -0.16641839175601203, -0.47355537095296341]]], [[[-0.98963277387073156, -0.84954014836916003, 0.58413699452656576, -0.90574771672091781], [-0.076165236495137068, -0.33185579541568311, 0.71508582816996036, 0.092982625872773284], [-0.59803076383019627, 0.96534061564337792, -0.86101074809006883, 0.4514321077629706]], [[0.68169942993010979, -0.46839019201575827, 0.86748862862295062, 0.077239140299438347], [-0.41039932807903745, 0.77351932859252215, -0.91858786520587499, -0.83803209999478323], [-0.66149540952220676, 0.21690925014134543, 0.49329666114093262, 0.22750525795569843]]], [[[0.49476796496001074, 0.54374663478072427, -0.64963120626875592, 0.57704965092442939], [0.74520801893525035, -0.18704889796101698, -0.73119506283113278, 0.30233417986347821], [-0.40842036405150561, -0.58450470869605797, 0.084020496164487923, -0.58538663622202569]], [[0.75611417989143481, -0.52196566347687212, 0.47144797786301607, -0.99505702762998172], [-0.0043333507755138889, -0.8740422871734248, 0.69683145050592943, -0.015830354510735223], [-0.047302931197539255, -0.033238606395723957, -0.95126827803180114, 0.88954001433619978]]]])
      res=sup(arg)
      ref=0.990730927773
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_float_rank0(self):
      arg=0.857535693433
      ref=0.857535693433
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_array_rank0(self):
      arg=0.170725403135
      ref=0.170725403135
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_array_rank1(self):
      arg=numarray.array([-0.20582799927627726, 0.0065475369467946631])
      ref=-0.205827999276
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_array_rank2(self):
      arg=numarray.array([[-0.00088705991533410966, -0.78884753663192009, 0.51871980588812661, -0.58204851803513313, 0.14241101940826861], [0.79094574969805964, 0.79672216617995617, -0.5690477768894624, 0.51157272654417052, 0.18066938665191556], [0.32364745994353683, 0.4748425103497671, 0.66679519455306924, -0.69644515487912217, -0.37494171775165297], [-0.18679695095262239, -0.78312977298360509, 0.044885312146701661, -0.44016241609550066, -0.49756845096624081]])
      ref=-0.788847536632
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_array_rank3(self):
      arg=numarray.array([[[-0.8081125483464171, -0.39512178653285135], [-0.88509761616456495, 0.27821695074525521]], [[-0.12065704909614361, 0.68332883926843135], [0.3403814721074454, -0.32879966956330042]], [[-0.7045229028656752, -0.1577338131736914], [-0.5966034105188045, 0.73456332700634985]], [[0.87843639147943997, 0.94490362642776882], [-0.45552277927474183, -0.9768135246661469]], [[-0.65451540143458864, -0.2796668209543185], [-0.085396852552783953, 0.83466766003117332]], [[0.43465138886082078, 0.61441480296663342], [0.92555078046558337, -0.24612584648713121]]])
      ref=-0.976813524666
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_array_rank4(self):
      arg=numarray.array([[[[-0.05536021838962446, -0.0098764348632232046, -0.52953548051740618, 0.58523070076080819], [0.10613609818504877, 0.67324697282212087, 0.54615663777618906, -0.67572666188479613], [-0.14390410848091539, 0.54011444546397591, -0.85609073171970373, 0.077099187121419277]], [[-0.20493469025659716, 0.90647730634569368, -0.71749107354722064, -0.12332697315517271], [0.63551267588618598, -0.60802528409862266, 0.052319255834022638, -0.95394697709081688], [-0.88612629669117959, 0.32248454322667519, 0.0072741938614420132, -0.69013368898879235]]], [[[-0.10462047858216028, -0.30097012474135698, -0.050016775782701028, 0.54822125876578376], [0.84395749034726886, -0.53249513893193168, -0.88004100855031275, -0.80779542570577179], [-0.79677629667791683, 0.95096027764472169, 0.63643207567783144, 0.97757008271555401]], [[-0.65415697736192047, -0.97050764835238645, -0.84814693021942777, -0.43855897064286542], [-0.37135308255264543, 0.041120751125186095, 0.036995657114785807, -0.35706630152349828], [-0.0030591331649565401, 0.48192500000712779, 0.18102011879743984, -0.78573775232009435]]], [[[-0.31965876602783605, 0.10351748464331689, 0.067424791069963907, -0.049524027182576535], [-0.5213817364489115, 0.027521683153738818, -0.24734661576641237, 0.24321699964875232], [-0.83947613904690699, 0.77162806253216987, -0.90740945316368071, -0.3420545897410685]], [[0.91845344502663262, -0.70878381509801414, 0.90861837177726379, -0.4013061463136427], [-0.18540388033546473, 0.9254510240675875, 0.30634230347058677, -0.97817133509804033], [-0.43975591131244984, 0.30020642565479139, 0.36841633323637479, 0.3066739733421715]]]])
      ref=-0.978171335098
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_constData_rank0(self):
      arg=Data(0.0114629834279,self.functionspace)
      ref=0.0114629834279
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_constData_rank1(self):
      arg=Data(numarray.array([-0.13734485813185704, -0.54812466656634307]),self.functionspace)
      ref=-0.548124666566
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_constData_rank2(self):
      arg=Data(numarray.array([[-0.89509582806410015, 0.10780316690435621, -0.93073500763754335, -0.38534759506991545, -0.6935160645644014], [-0.056672310128515813, 0.6285075027787359, 0.73632355512072167, -0.60238897825476267, 0.77403597203864094], [-0.5930215600641755, 0.72623233579382429, -0.32117191475695361, -0.081104170523293773, 0.62137628665436373], [0.2669734570396014, -0.65030905665614136, -0.53589374176691495, -0.48334830355881309, -0.89125004784938211]]),self.functionspace)
      ref=-0.930735007638
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_constData_rank3(self):
      arg=Data(numarray.array([[[0.94493948950092999, 0.13737629155757691], [-0.81599535906086107, -0.35832278646984816]], [[-0.53692839435234041, 0.53679590218669571], [0.038856705021854232, -0.18543838436402926]], [[0.19718168292863836, -0.55405958298613656], [0.16477861536800242, 0.17787953041277582]], [[0.51547288009005165, 0.35889372726595203], [-0.033476505587150873, -0.42142418570614026]], [[0.80507204877418204, -0.79581688922832838], [-0.85909254497735588, 0.66095083521227149]], [[0.46206420953978222, 0.53654696439305005], [0.57618105395776831, -0.22241758047110038]]]),self.functionspace)
      ref=-0.859092544977
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_constData_rank4(self):
      arg=Data(numarray.array([[[[0.50544713768476202, 0.96922321849050874, -0.81524480218696649, -0.36499730379849193], [-0.48131882706974372, 0.026812357207576465, 0.090903267401989618, -0.24742363369877829], [-0.51631372893805438, 0.30410275437953183, -0.75149566289642533, -0.19930300338453599]], [[0.82034878499482788, -0.70904661587698792, -0.27637223434426073, -0.34818734117560401], [0.11686048779802416, -0.76746266142163178, -0.75578186306174833, 0.14509316330390232], [0.1590050723141736, 0.69684384552537937, -0.58747105640080832, -0.28640840371441523]]], [[[0.14956532194045669, 0.081514192262221119, 0.32061383569406399, -0.2444346881437609], [0.79564139071785278, -0.5456680167461434, 0.24722978802719742, 0.28286130725068315], [0.10385207763921711, -0.064749181840278336, 0.21325254547672734, -0.71875644540473838]], [[0.58552496009870802, 0.35472373485671338, -0.18411162994671826, 0.71609038134967773], [-0.20966804574945064, -0.49286619989346314, 0.85116051808632553, -0.94417114370961075], [-0.40434528979823714, 0.62250343758157611, 0.64860074098639742, 0.0043146814280992096]]], [[[-0.14242849200713259, 0.42551908502898095, 0.7691157770973962, -0.37595641162856674], [0.026655444032149589, -0.82186407521644167, 0.40285091480648783, -0.53328831035315982], [-0.12887729257054481, 0.75610663428133451, 0.022049613835531723, 0.59949338706293043]], [[-0.34506254315071772, 0.019719877473602043, 0.10216765908478709, 0.022681548062032153], [0.2228614880408597, 0.26944547311401901, -0.10122095357202965, -0.51019076850180589], [-0.081439546799124463, 0.18829632566943544, 0.12366885442775377, 0.73651436499107814]]]]),self.functionspace)
      ref=-0.94417114371
      res=inf(arg)
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_taggedData_rank0(self):
      arg=Data(0.162086575852,self.functionspace)
      arg.setTaggedValue(1,0.162086575852)
      res=inf(arg)
      ref=0.162086575852
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_taggedData_rank1(self):
      arg=Data(numarray.array([0.97502548554439095, 0.14468929449768342]),self.functionspace)
      arg.setTaggedValue(1,[-0.80902425002058509, -0.89805781018804365])
      res=inf(arg)
      ref=-0.898057810188
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_taggedData_rank2(self):
      arg=Data(numarray.array([[-0.5694816711808135, 0.80659676528259427, 0.30198257784452243, -0.056794149786052461, 0.28747138649972137], [-0.55351234371619573, -0.7889070920065997, 0.88098536867480282, 0.18532519445126572, 0.4953885498039805], [0.034981671345569287, 0.32837659010563613, -0.12027778324017335, 0.6186028529626495, 0.94436199061580317], [-0.32534943311134756, 0.54392552535517247, -0.07184818553543626, -0.32050443715783694, -0.85778834873938736]]),self.functionspace)
      arg.setTaggedValue(1,[[-0.39066692582632312, 0.56443807842368665, 0.50983381732205602, 0.8533352755557535, 0.26857966396123456], [-0.43675811505896145, 0.061780174738994775, -0.70572251236028349, 0.58277425693964768, -0.77149002637252218], [0.21410554898576928, 0.63314655619690563, 0.83857171132417307, -0.64841751958506944, 0.49689361559212686], [-0.77395862547728433, -0.49190916492680103, -0.96727611195835639, 0.75749173365014144, 0.31952116010595022]])
      res=inf(arg)
      ref=-0.967276111958
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_taggedData_rank3(self):
      arg=Data(numarray.array([[[-0.42552857556236279, -0.35412978699074116], [0.32347364437155268, -0.06660424254888464]], [[0.69470606979783533, -0.27025063802848237], [0.60488222118615087, 0.75136007820911765]], [[-0.18334541260777715, 0.46826845351595714], [0.50971438744602593, -0.75173109388533477]], [[-0.83925390214063711, 0.10222838155029201], [0.97240290055902889, 0.61070698842426729]], [[-0.096661305892300042, 0.3060232400934193], [0.44355296215710527, -0.42328090263660423]], [[-0.37229736098865907, -0.61446651581066591], [-0.59861049707188863, -0.083231793539160881]]]),self.functionspace)
      arg.setTaggedValue(1,[[[0.098193175227128116, 0.70584535076979682], [0.64798160196689181, -0.4003221305355702]], [[-0.04270365807543719, -0.81717683614201619], [-0.40362171906510702, -0.64589178916948042]], [[-0.34373454409447479, -0.77185875470824516], [-0.0065539211815313081, 0.83912029516917275]], [[0.50979364877951405, 0.43268906334074786], [0.97762994631388556, 0.097902836960438]], [[0.25841131938252593, -0.075670175812033058], [-0.57350008458063262, -0.54227522655426319]], [[0.023586725726449931, 0.51624902917426008], [-0.88765839495000032, 0.61903659822536805]]])
      res=inf(arg)
      ref=-0.88765839495
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_taggedData_rank4(self):
      arg=Data(numarray.array([[[[-0.28454409743308506, -0.5668179780431859, 0.94184808892377703, 0.92323038327044005], [-0.62255144042695854, 0.10076051065495695, -0.75730451262334242, 0.92747176366592243], [-0.27419623999345366, -0.2073834606691658, -0.80951479323876696, 0.7135757150367612]], [[-0.97314628579045248, 0.68637741632500915, 0.53335906255833199, 0.96086938773033026], [-0.94945339539011386, -0.81392220879720334, -0.85288194032365472, -0.16477859877633083], [-0.21874528606061738, 0.02516499720401133, -0.049189220512075194, 0.037621105958556278]]], [[[-0.047849568456701519, 0.17313004158353174, -0.73652898602638794, 0.54596512107217787], [0.53827755889057149, 0.82310276090363499, -0.28442357567638377, -0.27179900237071242], [-0.19399624326093612, -0.12349387443529247, -0.56735519289734504, -0.68464632829023375]], [[0.33637651162122606, 0.22878413048564927, -0.70679532991608829, -0.25607454879877101], [-0.29261869752827607, -0.210849838585774, -0.15002781188439474, -0.40403456108229485], [0.98194953410530106, -0.10120885200799257, -0.73197515875287644, -0.79609802011910569]]], [[[-0.25327135865022421, -0.57325273532977317, 0.81059261361542845, -0.046095046141370499], [0.96531384124283903, -0.61153728454098477, -0.97457818352675707, -0.62195519514120923], [0.69039154876488285, 0.92574092049825851, -0.16814749709359544, -0.10824422727782279]], [[-0.54851863859404215, -0.6147932835008143, -0.12084391085229162, 0.57476183090130717], [0.23951711064034398, -0.54814262685597792, 0.94406155895646182, 0.5501815266421175], [0.21626599181010975, -0.9661525769220054, 0.85862231711120418, 0.97115736784479956]]]]),self.functionspace)
      arg.setTaggedValue(1,[[[[-0.44268265295180487, -0.94661739524160704, -0.58229956703493468, -0.83249621866194934], [-0.68913939538381319, -0.75804956419646286, 0.64749957799330526, 0.32616508244527531], [-0.14277223200775602, 0.74824266545401041, -0.75701815640062908, 0.27359871030382887]], [[-0.48033424095007193, -0.067929127807738965, 0.45865972622825102, -0.21732884173055611], [-0.84531291000880504, 0.83240578985878577, -0.25004340292634253, 0.21302689907942329], [0.46335645072330633, -0.10408082112144545, 0.96768118368872202, -0.14088569164618203]]], [[[0.28258575728428736, -0.45763983811002018, -0.29020353090870232, -0.35290974430439914], [0.8698073989711228, 0.99501431931359741, 0.36030715281910175, -0.27073779707934853], [-0.86095060819411717, -0.85818974539956372, 0.88449786126964258, -0.4780560878654454]], [[-0.12884556589949114, 0.46614981239408104, -0.82753211385973646, 0.71952728123859067], [-0.43573431868291923, 0.41690370693945611, -0.024901756894041061, 0.14934111352857715], [0.78047901800597619, 0.26373552686712709, -0.95380580583786689, 0.63288455897048768]]], [[[-0.94242811017558603, -0.2101877288188847, 0.91797691709730911, 0.45617540058153594], [-0.93413158750269942, -0.70238302965184563, 0.36501953522261488, -0.43956770565509551], [-0.056597755835613439, 0.41357120112496037, -0.60521324615373961, -0.2181851605833629]], [[-0.21130451442872822, 0.53289932629760295, -0.72866050065430343, -0.90245795340674695], [0.1171412044571345, -0.91170721460310178, 0.58826690518723557, -0.39254304440341792], [-0.60338068872893746, -0.40393966609543219, -0.69793709205007581, -0.50271106301761193]]]])
      res=inf(arg)
      ref=-0.974578183527
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_expandedData_rank0(self):
      arg=Data(0.97331285569,self.functionspace,True)
      arg.setTaggedValue(1,0.97331285569)
      res=inf(arg)
      ref=0.911344578019
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_expandedData_rank1(self):
      arg=Data(numarray.array([-0.025782649803792301, 0.15017595667012174]),self.functionspace,True)
      arg.setTaggedValue(1,[-0.69996944983865683, 0.23286168827412723])
      res=inf(arg)
      ref=-0.699969449839
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_expandedData_rank2(self):
      arg=Data(numarray.array([[-0.11924617230337731, -0.69281409701569552, 0.8819805848672797, 0.67606672557817737, -0.15719846909292801], [-0.91958498979699277, 0.35156143886023572, -0.66768679685578647, -0.48673444352737993, -0.33274286155542288], [-0.011930486621410052, -0.035661622752792166, 0.8109198510155371, 0.2680383999767102, 0.38679270088753159], [-0.64107995512896321, 0.58836944145799008, -0.55856121453754959, 0.63346354980998365, -0.74990799644898765]]),self.functionspace,True)
      arg.setTaggedValue(1,[[0.075603512293304975, 0.65374459813382813, -0.68190855977115938, 0.023184369691352069, -0.15500904033575869], [-0.28865053906142646, 0.94679590620919218, -0.53701383780009748, -0.5413294941411595, -0.1693286769782183], [-0.34870739946353968, -0.92230175409787241, -0.95867600429924127, -0.0012605712136730673, 0.33180267589091184], [-0.60152160703395885, 0.3397221706065463, 0.58860878406913142, 0.59589494790106934, -0.96157487116725071]])
      res=inf(arg)
      ref=-0.961574871167
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_expandedData_rank3(self):
      arg=Data(numarray.array([[[0.76729109013175933, 0.93446692332588199], [-0.85883010569426133, 0.94426468618035608]], [[0.93376991286356659, -0.31027020688815887], [0.68776369494826284, -0.71422083083074273]], [[0.88895851954297633, -0.0041028262794304826], [-0.56213269421968981, -0.53305350433881493]], [[0.91080353792869873, 0.87460732654624618], [-0.43411934331432445, -0.66454826732696648]], [[-0.7458989571542316, 0.0054889915153857327], [-0.43786572648758848, 0.45559010284736901]], [[0.97274627384999923, 0.28719064869132427], [-0.068944574283213989, 0.92462015394940389]]]),self.functionspace,True)
      arg.setTaggedValue(1,[[[0.45653716082021401, -0.23702471933521396], [-0.45177153258555913, -0.41351543065221219]], [[0.24913453932538854, -0.29954039408048372], [0.47106013214852771, -0.86663080347446608]], [[-0.75781905181266151, 0.18918623668675538], [-0.66636479853389385, 0.74776214919244066]], [[-0.1199599974824086, 0.90306995612278795], [0.77248429149482201, 0.094536885699270545]], [[-0.5505821492977705, -0.42454185319972981], [-0.035064708067998662, -0.75303714911308206]], [[-0.63966405281134198, -0.5545794126646626], [-0.010523374154368881, -0.6779302443539712]]])
      res=inf(arg)
      ref=-0.866630803474
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")
   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   def test_inf_expandedData_rank4(self):
      arg=Data(numarray.array([[[[0.59301561766072242, 0.6053289847447505, 0.1720524601196014, 0.4833679091408507], [-0.09424534157102582, 0.91976800441974427, -0.11950945037648264, -0.10547294628619919], [0.045302052907340462, -0.26699801365763864, -0.59100789693502542, 0.8220898248447317]], [[-0.71764117646989489, -0.31097048109028447, -0.36880919201656037, -0.12380648590370025], [0.71171672214402615, -0.87716399369053821, 0.53673002158727501, -0.29732656272582969], [-0.71737153055394787, -0.18861795164190798, -0.55325559014387649, 0.50828072891961451]]], [[[0.7139879661021804, -0.54655602578694507, -0.62738810494945163, -0.71453639471562314], [0.25287437655577372, 0.47571671295036144, -0.86354401602062025, -0.45354287599804666], [0.44863440392806031, -0.64112646326775291, 0.91877609169151908, 0.98277738739744636]], [[-0.23627452345538003, 0.26824331838587279, -0.11255839194272887, 0.62327509528972769], [0.87773514238000572, 0.62819546928589953, 0.53463187163744053, -0.11780048176113112], [0.17743656841726674, -0.47418242507357755, 0.37012543459814751, 0.84217962384585232]]], [[[-0.51615332690558602, -0.32180746309362296, -0.86871747985681558, -0.89347967000489548], [0.79236651642797762, -0.82769899986382356, -0.27503369025221813, 0.70491420567015894], [-0.46319485284379014, -0.71191635366918904, -0.4935498645964489, 0.42096986386791024]], [[-0.95241673224815648, -0.99507534300684708, -0.76759978752396862, -0.1621109099401965], [0.011495974208982407, 0.90547796535927083, -0.79651097416963368, 0.44382172569512224], [-0.11210429012307821, 0.28303519643554576, -0.59783370185364437, 0.75478113456483809]]]]),self.functionspace,True)
      arg.setTaggedValue(1,[[[[0.70392145790103089, -0.7532260844479588, -0.53422583924595046, 0.29798767877629717], [0.87516135156666075, -0.40276494310086908, -0.070320230959511232, 0.13669067815714553], [0.26690559342737163, -0.0067218741792103298, 0.79272133073277207, -0.90080452301155178]], [[-0.90606808543363915, 0.94927254364127123, -0.15130300474521019, -0.34632615506118936], [0.29908319964670294, -0.71717585217034374, -0.11325472101126977, 0.60911487311528911], [0.055579283773870003, -0.42895625110409941, 0.22788310989168781, 0.50000680777627848]]], [[[-0.56683598589959683, 0.05720680440798942, 0.27949696950971137, 0.028472427888231566], [-0.51619860801944828, 0.73067372155966415, 0.81910302859539108, 0.97736576074434089], [-0.68851989268805447, 0.93303896549540433, -0.21083269397032622, -0.93122611728597304]], [[-0.13573136829275501, -0.30858456571228721, 0.86534434039603658, 0.79114631697533278], [-0.32870753097283112, 0.68472493656951161, -0.071089197922887593, -0.032936485145198757], [0.40875380316397925, -0.27280027645573401, -0.27155299841136848, -0.81500786568603067]]], [[[0.080473457255322733, -0.54931257643900566, 0.19082193004877501, -0.52744293877995374], [0.44447571282767639, -0.0060974808128133606, -0.87810788913565485, -0.92670043304972327], [0.45760340069120065, -0.7504981866011311, -0.95431384972577082, 0.4644557658783679]], [[0.18449401878090566, -0.73315009595487135, -0.59404794479416467, 0.076556365298426687], [0.65441299213898763, -0.91649003079073932, -0.082437357750317686, 0.88150531806538579], [-0.45829076922310796, 0.91236619416957021, -0.91128142664652279, 0.9524897339484899]]]])
      res=inf(arg)
      ref=-0.995075343007
      self.failUnless(isinstance(res,float),"wrong type of result.")
      self.failUnless(abs(res-ref)<=self.tol*abs(ref),"wrong result")

