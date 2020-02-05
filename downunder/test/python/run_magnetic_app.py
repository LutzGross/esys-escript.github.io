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
from esys.ripley import Brick,Rectangle
from esys.downunder.apps import MagneticModel2D,MagneticModel3D

class TestMagneticApp2D(unittest.TestCase):


    xdim = 1000.
    zdim = 1000.
    NEX = 100
    NEZ = 100
    TESTTOL = 1e-3
    

    def setUp(self):
        self.domain = Rectangle(self.NEX,self.NEZ,l0=self.xdim,l1=self.zdim)

    def tearDown(self):
        del self.domain  

    def test_pde_answer(self):
        model = MagneticModel2D(self.domain,fixVert=True)
        x=self.domain.getX()
        xmin=inf(x[0])
        xmax=sup(x[0])
        ymin=inf(x[1])
        ymax=sup(x[1])
        xp = 2.*np.pi/(xmax-xmin)
        yp = 2.*np.pi/(ymax-ymin)
        Bh = [1., 0.] 
        k = ((xp**2+yp**2)/xp)*cos(xp*x[0])*cos(yp*x[1]) 
        model.setSusceptibility(k)
        outk= model.getSusceptibility()
        kdiff=abs(k-outk)
        self.assertLessEqual(sup(kdiff),self.TESTTOL)
        model.setBackgroundMagneticField(Bh)
        actualanswer = sin(xp*x[0])*cos(yp*x[1])
        abserror = abs(actualanswer-model.getAnomalyPotential())
        biggesterror = sup(abserror)
        self.assertLessEqual(biggesterror,self.TESTTOL)

class TestMagneticApp3D(unittest.TestCase):


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
        model = MagneticModel3D(self.domain,fixVert=True)
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
        Bh = [1., 0., 0.] 
        k = ((xp**2+yp**2+zp**2)/xp)*cos(xp*x[0])*sin(yp*x[1])*cos(zp*x[2]) 
        model.setSusceptibility(k)
        outk= model.getSusceptibility()
        kdiff=abs(k-outk)
        self.assertLessEqual(sup(kdiff),self.TESTTOL)
        model.setBackgroundMagneticField(Bh)
        actualanswer = sin(xp*x[0])*sin(yp*x[1])*cos(zp*x[2])
        abserror = abs(actualanswer-model.getAnomalyPotential())
        biggesterror = sup(abserror)
        self.assertLessEqual(biggesterror,self.TESTTOL)


if __name__ == '__main__':

    mySuite=unittest.TestSuite()
    mySuite.addTest(unittest.makeSuite(TestMagneticApp2D)) 
    mySuite.addTest(unittest.makeSuite(TestMagneticApp3D)) 
    runner=unittest.TextTestRunner(verbosity=2)
    runner.run(mySuite)
