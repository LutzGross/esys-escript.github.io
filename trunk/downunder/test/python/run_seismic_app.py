__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross"

from esys.escriptcore.testing import *
import esys.escriptcore.utestselect as unittest
import sys, scipy.special
from esys.escript import *
from esys.finley import Rectangle
from esys.downunder.apps import SonicWaveInFrequencyDomain, PMLCondition
from esys.weipa import saveSilo
from esys.escript.pdetools import Locator

import numpy as np

NO_TRILINOS = not hasFeature("trilinos")

class TestSeismic2D(unittest.TestCase):
    # grid spacing in 
    DX=10.
    NEx=200
    Order=2
    VP=2000
    NE_PML=40  # PML thickness in term of number of elements
    NE_STATION0=70 # first station offset in term of number of elements
    Amplitude=1
    S0=1000
    TESTTOL=1e-2
    
    def setUp(self):
        Width=self.DX*self.NEx
        Center=self.NEx//2*self.DX
        self.domain=Rectangle(self.NEx,self.NEx,l0=Width,l1=Width, diracPoints=[(Center,Center)], diracTags=[ "src" ], order=self.Order, fullOrder=True)
        numStation=(self.NEx//2-3-self.NE_STATION0-1)
        stations=[ ( (self.NE_STATION0 + i) * self.DX, Center) for i in range(numStation)] 
        self.loc=Locator(Solution(self.domain), stations )
        self.source=Scalar(0j, DiracDeltaFunctions(self.domain))
        self.source.setTaggedValue("src", self.Amplitude)
        self.sourceX=(Center, Center)
    def tearDown(self):
        del self.domain
        del self.loc
        del self.source
        del self.sourceX
                 
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")           
    def test_RadialWaveLowF(self):
        Frequency=0.01
        k=Frequency/self.VP*2*np.pi
        assert k*self.DX <1
        L_PML=self.DX*self.NE_PML
        pml_condition=PMLCondition(sigma0=self.S0, Lleft=[L_PML,L_PML], Lright=[L_PML, L_PML], m=4)
        # sonic wave model
        wave=SonicWaveInFrequencyDomain(self.domain, vp=self.VP, pml_condition=pml_condition)
        wave.setFrequency(Frequency)
        # set source and get wave response: 
        u=wave.getWave(self.source)
        
        r=np.array([sqrt((x[0]-self.sourceX[0])**2+(x[1]-self.sourceX[1])**2) for x in self.loc.getX() ])
        uu2=-scipy.special.hankel2(0, k*r)*1j/4
        uu=np.array(self.loc(u))
        errorA=max(abs(uu2-uu))
        A=max(abs(uu2))
        self.assertLessEqual(errorA, A*self.TESTTOL)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")        
    def test_RadialWaveMediumF(self):
        Frequency=1.
        k=Frequency/self.VP*2*np.pi
        assert k*self.DX <1
        L_PML=self.DX*self.NE_PML
        pml_condition=PMLCondition(sigma0=self.S0, Lleft=[L_PML,L_PML], Lright=[L_PML, L_PML], m=4)
        # sonic wave model
        wave=SonicWaveInFrequencyDomain(self.domain, vp=self.VP, pml_condition=pml_condition)
        wave.setFrequency(Frequency)
        # set source and get wave response: 
        u=wave.getWave(self.source)
        
        r=np.array([sqrt((x[0]-self.sourceX[0])**2+(x[1]-self.sourceX[1])**2) for x in self.loc.getX() ])
        uu2=-scipy.special.hankel2(0, k*r)*1j/4
        uu=np.array(self.loc(u))
        errorA=max(abs(uu2-uu))
        A=max(abs(uu2))
        self.assertLessEqual(errorA, A*self.TESTTOL)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")        
    def test_RadialWaveHighF(self):
        Frequency=10.
        k=Frequency/self.VP*2*np.pi
        assert k*self.DX <1
        L_PML=self.DX*self.NE_PML
        pml_condition=PMLCondition(sigma0=self.S0, Lleft=[L_PML,L_PML], Lright=[L_PML, L_PML], m=4)
        # sonic wave model
        wave=SonicWaveInFrequencyDomain(self.domain, vp=self.VP, pml_condition=pml_condition)
        wave.setFrequency(Frequency)
        # set source and get wave response: 
        u=wave.getWave(self.source)
        
        r=np.array([sqrt((x[0]-self.sourceX[0])**2+(x[1]-self.sourceX[1])**2) for x in self.loc.getX() ])
        uu2=-scipy.special.hankel2(0, k*r)*1j/4
        uu=np.array(self.loc(u))
        errorA=max(abs(uu2-uu))
        A=max(abs(uu2))
        self.assertLessEqual(errorA, A*self.TESTTOL)
        

        
if __name__ == '__main__':

    mySuite=unittest.TestSuite()
    mySuite.addTest(unittest.makeSuite(TestSeismic2D)) 
    runner=unittest.TextTestRunner(verbosity=2)
    runner.run(mySuite)
