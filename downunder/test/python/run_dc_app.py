__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"


import unittest
import sys
import numpy as np
from esys.escript import *
from esys.escript.pdetools import Locator
from esys.ripley import Brick
from esys.downunder.apps import DCResistivityModel, DCResistivityModelNoPrimary
from esys.downunder import *
from esys.escriptcore.testing import *
import esys.escriptcore.utestselect as unittest
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions

class TestDCResistivityApp(unittest.TestCase):

    dx = 10         # grid line spacing in [m]
    NEx = 50      # number of nodes in the x direction
    NEy = 50      # number of nodes in the y direction
    NEz = 50      # number of nodes in the z direction
    H0 = 600       # height [m] of transect above bottom of domain (will be locked to grid)
    sig_p = 1.     # primary conductivity
    sig_2 = 100.   # conductivity in ball    xdim = 1000.
    TESTTOL = 1e-3
    POSx = 100
    POSy = 70
    NEGx = 100
    NEGy = 110
    POS_node = ( POSx*dx, POSy*dx, NEz*dx)
    NEG_node = ( NEGx*dx, NEGy*dx, NEz*dx)
    Z0=100           # vertical position of circle below transect [m]
    Lx = dx*NEx
    Ly = dx*NEy
    Lz = dx*NEz
    c=[Lx/2.,Ly/2., H0-Z0] # circle center
    R=50.            # radius




    def setUp(self):
        self.domain=Brick(n0=self.NEx, n1=self.NEy, n2=self.NEz, l0=self.Lx, l1=self.Ly, l2=self.Lz,
              diracPoints= [self.POS_node, self.NEG_node], diracTags = ['e0','e1'])
 
    def tearDown(self):
        del self.domain  

    def test_pde_Primary(self):

        model1=DCResistivityModel(self.domain,sigma0=self.sig_p,useFastSolver = True)
        model1.setPrimaryPotentialForHalfSpace(sources= [self.POS_node, self.NEG_node], 
                                      charges=[1.0, -1.0] )

        analyticprimary = model1.setPrimaryPotentialForHalfSpace(sources= [self.POS_node, self.NEG_node], 
                                      charges=[1.0, -1.0] ) 
        primaryans1 = model1.getPrimaryPotential()
        model2=DCResistivityModel(self.domain,sigma0=self.sig_p,useFastSolver = True)
        src=Scalar(0., DiracDeltaFunctions(self.domain))
        src.setTaggedValue('e0', 1.0)
        src.setTaggedValue('e1', -1.0)
        model2.setPrimaryPotential(source=src) 
        primaryans2 = model2.getPrimaryPotential()
        abserror=abs(primaryans1-primaryans2)
        self.assertLessEqual(sup(abserror),self.TESTTOL*10.,"primaries")
        sigma = Scalar(self.sig_p,ContinuousFunction(self.domain))
        x=self.domain.getX()
        d=length(x-self.c)
        sphereCond=sigma+(self.sig_2-self.sig_p)*whereNegative(d-self.R)    # 0 for d>R and 1 for d<R
        model1.setConductivity(sphereCond)
        model2.setConductivity(sphereCond)
        us1=model1.getSecondaryPotential()
        us2=model1.getSecondaryPotential()
        abserror=abs(us1-us2)
        self.assertLessEqual(sup(abserror),self.TESTTOL*1.,"secondaries")
        ut1=model1.getPotential()
        ut2=model2.getPotential()
        model3=DCResistivityModelNoPrimary(self.domain,source=src,sigma=sphereCond,useFastSolver = True)
        ut3=model3.getPotential()
        abserror=abs(ut1-ut2)
        self.assertLessEqual(sup(abserror),self.TESTTOL*10.,"total with primaries")    
        abserror=abs(ut1-ut3)
        self.assertLessEqual(sup(abserror),self.TESTTOL*10.,"total with analytic primary")    
        abserror=abs(ut2-ut3)
        self.assertLessEqual(sup(abserror),self.TESTTOL*10.,"total with FE primary")    
                        
        
        
        
        
        

################################
if __name__ == '__main__':

    mySuite=unittest.TestSuite()
    mySuite.addTest(unittest.makeSuite(TestDCResistivityApp)) 
    runner=unittest.TextTestRunner(verbosity=2)
    runner.run(mySuite)

    
