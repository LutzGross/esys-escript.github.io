#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
Test suite for the linearPDE  and pdetools test on finley

@remark:

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Lutz Gross, l.gross@uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"


import unittest
from test_linearPDEs import Test_Poisson,Test_LinearPDE, Test_LinearPDE_noLumping
from test_assemblage import Test_assemblage_2Do1, Test_assemblage_2Do2, Test_assemblage_3Do1, Test_assemblage_3Do2, \
                            Test_assemblage_2Do1_Contact,Test_assemblage_2Do2_Contact, Test_assemblage_3Do1_Contact, Test_assemblage_3Do2_Contact
from test_pdetools import Test_pdetools, Test_pdetools_noLumping
from esys.escript import *
from esys.finley import Rectangle,Brick,JoinFaces, ReadMesh
import sys


try:
     FINLEY_TEST_DATA=os.environ['FINLEY_TEST_DATA']
except KeyError:
     FINLEY_TEST_DATA='.'

FINLEY_TEST_MESH_PATH=FINLEY_TEST_DATA+"/data_meshes/"

NE=6 # number of element in each spatial direction (must be even)

class Test_LinearPDEOnFinleyHex2DOrder1(Test_LinearPDE):
# class Test_LinearPDEOnFinleyHex2DOrder1(Test_LinearPDE,Test_pdetools,Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,1)
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyHex2DOrder2(Test_LinearPDE,Test_pdetools,Test_assemblage_2Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2)
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyHex3DOrder1(Test_LinearPDE,Test_pdetools,Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,1)

class Test_LinearPDEOnFinleyHex3DOrder2(Test_LinearPDE,Test_pdetools,Test_assemblage_3Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Brick(NE,NE,NE,2)
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet2DOrder1(Test_LinearPDE,Test_pdetools,Test_assemblage_2Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_2D_order1.fly",optimize=False)
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet2DOrder2(Test_LinearPDE_noLumping,Test_pdetools_noLumping,Test_assemblage_2Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_2D_order2.fly",optimize=False)
   def tearDown(self):
        del self.domain

class Test_LinearPDEOnFinleyTet3DOrder1(Test_LinearPDE,Test_pdetools,Test_assemblage_3Do1):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_3D_order1.fly",optimize=False)

class Test_LinearPDEOnFinleyTet3DOrder2(Test_LinearPDE,Test_pdetools,Test_assemblage_3Do2):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = ReadMesh(FINLEY_TEST_MESH_PATH+"tet_3D_order2.fly",optimize=False)
   def tearDown(self):
        del self.domain

class Test_PoissonOnFinley(Test_Poisson):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2)
   def tearDown(self):
        del self.domain


class Test_AssemblePDEwithFinley_2Do1_Contact(Test_assemblage_2Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do2_Contact(Test_assemblage_2Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do1_Contact(Test_assemblage_3Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do2_Contact(Test_assemblage_3Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain


class Test_AssemblePDEwithFinley_2Do1_Contact_withElementsOnFace(Test_assemblage_2Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1,useElementsOnFace=True)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1,useElementsOnFace=True)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do2_Contact_withElementsOnFace(Test_assemblage_2Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2,useElementsOnFace=True)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2,useElementsOnFace=True)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do1_Contact_withElementsOnFace(Test_assemblage_3Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do2_Contact_withElementsOnFace(Test_assemblage_3Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
       x1 = ContinuousFunction(d1).getX()
       ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
       ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       d2.setX(d2.getX()+[0.5,0.,0.])
       self.domain = JoinFaces([d1,d2],optimize=False)
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True:
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyHex2DOrder1))
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyHex2DOrder2))
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyHex3DOrder1))
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyHex3DOrder2))

      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyTet2DOrder1))
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyTet2DOrder2))
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyTet3DOrder1))
      suite.addTest(unittest.makeSuite(Test_LinearPDEOnFinleyTet3DOrder2))
  
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_2Do1_Contact))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_2Do2_Contact))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_3Do1_Contact))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_3Do2_Contact))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_2Do1_Contact_withElementsOnFace))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_2Do2_Contact_withElementsOnFace))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_3Do1_Contact_withElementsOnFace))
      suite.addTest(unittest.makeSuite(Test_AssemblePDEwithFinley_3Do2_Contact_withElementsOnFace))
   else:
      pass

   s=unittest.TextTestRunner(verbosity=2).run(suite)
   if not s.wasSuccessful(): sys.exit(1)

