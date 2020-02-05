__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import sys
import numpy as np
import esys.escript.unitsSI as U
from esys.escript import *
from esys.escript.pdetools import Locator
from esys.ripley import Brick
from esys.downunder.apps import GravityModel

class TestGravityApp(unittest.TestCase):


    xdim = 1000.
    ydim = 1000.
    zdim = 1000.
    NEX = 100
    NEZ = 100
    TESTTOL = 1e-3
    

    def setUp(self):
        self.domain = Brick(self.NEX,self.NEX,self.NEZ,l0=self.xdim,l1=self.ydim,l2=self.zdim)

    def tearDown(self):
        del self.domain  

    def test_pde_answer(self):
        model = GravityModel(self.domain,fixBase=True)
        x=self.domain.getX()
        xmin=inf(x[0])
        xmax=sup(x[0])
        ymin=inf(x[1])
        ymax=sup(x[1])
        zmin=inf(x[2])
        zmax=sup(x[2])
        xp = 2.*np.pi/(xmax-xmin)
        yp = 2.*np.pi/(ymax-ymin)
        zp = 2.*np.pi/(zmax-zmin)
        Dens = (xp**2+yp**2+zp**2)*cos(xp*x[0])*cos(yp*x[1])*sin(zp*x[2])/(4.0*np.pi*U.Gravitational_Constant) 
        model.setDensity(Dens)
        outdens= model.getDensity()
        densdiff=abs(Dens-outdens)
        self.assertLessEqual(sup(densdiff),self.TESTTOL)
        actualanswer = -cos(xp*x[0])*cos(yp*x[1])*sin(zp*x[2])
        abserror = abs(actualanswer-model.getGravityPotential())
        biggesterror = sup(abserror)
        self.assertLessEqual(biggesterror,self.TESTTOL)
        gz=model.getzGravity()
        actualgz= -zp*cos(xp*x[0])*cos(yp*x[1])*cos(zp*x[2])
        errorgz=abs(gz-actualgz)
        self.assertLessEqual(sup(errorgz),self.TESTTOL*100.)


if __name__ == '__main__':

    mySuite=unittest.TestSuite()
    mySuite.addTest(unittest.makeSuite(TestGravityApp)) 
    runner=unittest.TextTestRunner(verbosity=2)
    runner.run(mySuite)
