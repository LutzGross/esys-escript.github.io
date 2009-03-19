
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
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
Test suite for the linearPDE  and pdetools test on finley

@remark:

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import os

import unittest
from test_linearPDEs import Test_Poisson,Test_LinearPDE, Test_LinearPDE_noLumping, Test_TransportPDE, Test_Helmholtz, Test_LameEquation
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

FINLEY_TEST_MESH_PATH=os.path.join(FINLEY_TEST_DATA,"data_meshes")

NE=8 # number of element in each spatial direction (must be even)

class Test_LameOnFinley(Test_LameEquation):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2,useFullElementOrder=True)
   def tearDown(self):
        del self.domain
class Test_PoissonOnFinley(Test_Poisson):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2,useFullElementOrder=True)
   def tearDown(self):
        del self.domain

class Test_HelmholtzOnFinley(Test_Helmholtz):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
        self.domain = Rectangle(NE,NE,2,useFullElementOrder=True)
   def tearDown(self):
        del self.domain


class Test_AssemblePDEwithFinley_2Do1_Contact(Test_assemblage_2Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do1_Contact.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do2_Contact(Test_assemblage_2Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do2_Contact.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do1_Contact(Test_assemblage_3Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do1_Contact.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do2_Contact(Test_assemblage_3Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do2_Contact.fly"))
   def tearDown(self):
        del self.domain


class Test_AssemblePDEwithFinley_2Do1_Contact_withElementsOnFace(Test_assemblage_2Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1,useElementsOnFace=True)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=1,useElementsOnFace=True)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do1_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_2Do2_Contact_withElementsOnFace(Test_assemblage_2Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2,useElementsOnFace=True)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Rectangle(n0=int(NE/2),n1=NE,l0=0.5,order=2,useElementsOnFace=True)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_2Do2_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do1_Contact_withElementsOnFace(Test_assemblage_3Do1_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=1,useElementsOnFace=True)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do1_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

class Test_AssemblePDEwithFinley_3Do2_Contact_withElementsOnFace(Test_assemblage_3Do2_Contact):
   RES_TOL=1.e-7
   ABS_TOL=1.e-8
   def setUp(self):
       # d1 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
       # x1 = ContinuousFunction(d1).getX()
       # ContinuousFunction(d1).setTags(1,Scalar(1,ContinuousFunction(d1)))
       # d2 = Brick(n0=int(NE/2),n1=NE,n2=NE,l0=0.5,order=2,useElementsOnFace=True)
       # ContinuousFunction(d2).setTags(2,Scalar(1,ContinuousFunction(d2)))
       # d2.setX(d2.getX()+[0.5,0.,0.])
       # self.domain = JoinFaces([d1,d2],optimize=False)
       self.domain=ReadMesh(os.path.join(FINLEY_TEST_MESH_PATH,"mesh_3Do2_Contact_withElementsOnFace.fly"))
   def tearDown(self):
        del self.domain

if __name__ == '__main__':
   suite = unittest.TestSuite()
   if True :
      suite.addTest(unittest.makeSuite(Test_PoissonOnFinley))
      suite.addTest(unittest.makeSuite(Test_HelmholtzOnFinley))
      suite.addTest(unittest.makeSuite(Test_LameOnFinley))
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

