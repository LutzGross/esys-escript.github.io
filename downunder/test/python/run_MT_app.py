__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross"

# import unittest
import sys
import numpy as np
from esys.escript import *
from esys.escript.pdetools import Locator
from esys.ripley import Rectangle
from esys.downunder.apps import MT2DTEModel, MT2DTMModel
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *

NO_TRILINOS = not hasFeature("trilinos")

class TestMT2DTE(unittest.TestCase):
    
    DEPTH=3000.
    WIDTH=1000.
    SIGMA0=1./100.
    NEX=100
    NEZ=300
    STATION_OFFSET=500.
    
    TRUE_PHASE=45.
    TRUE_RHO=100.
    USEFASTSOLVER=False # True works not yet
    
    def setUp(self):
        self.domain=Rectangle(self.NEX,self.NEZ,l0=self.WIDTH,l1=self.DEPTH)
        self.loc=Locator(ReducedFunction(self.domain), [ (self.STATION_OFFSET,self.DEPTH-self.DEPTH/self.NEZ/2)] )
    
    def tearDown(self):
        del self.domain
        del self.loc
    
    def runModel(self, model, PERIODS, TOL):
        for frq in PERIODS:
            Zxy = model.getImpedance(f=frq)
            self.assertIsInstance(Zxy, Data)
            self.assertEqual(Zxy.getShape(), ())
            self.assertEqual(Zxy.getFunctionSpace(), ReducedFunction(self.domain))
            
            
            rho=model.getApparentResitivity(frq, Zxy)
            self.assertIsInstance(rho, Data)
            self.assertEqual(rho.getShape(), ())
            self.assertEqual(rho.getFunctionSpace(), ReducedFunction(self.domain))
            
            
            phi=model.getPhase(frq, Zxy)
            self.assertIsInstance(phi, Data)
            self.assertEqual(phi.getShape(), ())
            self.assertEqual(phi.getFunctionSpace(), ReducedFunction(self.domain))

            #print(frq, self.loc(rho)[0], self.loc(phi)[0])
            self.assertAlmostEqual(self.loc(rho)[0], self.TRUE_RHO, delta=TOL*self.TRUE_RHO)
            self.assertAlmostEqual(self.loc(phi)[0], self.TRUE_PHASE, delta=TOL*self.TRUE_PHASE)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_SigmaFloatSigmaBFloat(self):
        """
        conductivity set as float, radiation condition
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        model.setConductivity(self.SIGMA0, sigma_boundary=self.SIGMA0)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_SigmaFloatNoSigmaB(self):
        """
        conductivity set as float, no value on boundary, radiation condition 
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        model.setConductivity(self.SIGMA0)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_SigmaDataNoSigmaB(self):
        """
        conductivity set as Scalar, no value on boundary, radiation condition 
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        sigma=Scalar(self.SIGMA0, ContinuousFunction(self.domain))
        model.setConductivity(sigma)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")    
    def test_SigmaDataNoSigmaBFails(self):
        """
        conductivity set as Scalar, no value on boundary but interpolation is not possible
        """        
        model=MT2DTEModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        sigma=Scalar(self.SIGMA0, Function(self.domain))
        self.assertRaises(RuntimeError, model.setConductivity, *(sigma, ))
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_SigmaFloatSigmaBFloatFixedButtom(self):
        """
        conductivity set as float, fixed bottom -> TRUE values are not matched exactly
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        model.setConductivity(self.SIGMA0, sigma_boundary=self.SIGMA0)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_SigmaFloatNoSigmaBFixedButtom(self):
        """
        conductivity set as float, no value on boundary, fixed bottom -> TRUE values are not matched exactly
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        model.setConductivity(self.SIGMA0)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_SigmaDataNoSigmaBFixedButtom(self):
        """
        conductivity set as Scalar, no value on boundary, fixed bottom -> TRUE values are not matched exactly
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        sigma=Scalar(self.SIGMA0, ContinuousFunction(self.domain))
        model.setConductivity(sigma)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")    
    def test_SigmaDataNoSigmaBFixedButtom2(self):
        """
        conductivity set as Scalar, no value on boundary, interpolation is not possible but fixed bottom -> TRUE values are not matched exactly
        """        
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTEModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        sigma=Scalar(self.SIGMA0, Function(self.domain))
        model.setConductivity(sigma)
        self.runModel(model, PERIODS, TOL=1.e-1)        

class TestMT2DTMNoAirLayer(unittest.TestCase):
    
    DEPTH=3000.
    WIDTH=1000.
    RHO0=100.
    NEX=100
    NEZ=300
    STATION_OFFSET=500.
    
    TRUE_PHASE=45.
    TRUE_RHO=100.
    USEFASTSOLVER=False # True works not yet
    
    def setUp(self):
        self.domain=Rectangle(self.NEX,self.NEZ,l0=self.WIDTH,l1=self.DEPTH)
        self.loc=Locator(ReducedFunction(self.domain), [ (self.STATION_OFFSET,self.DEPTH-self.DEPTH/self.NEZ/2)] )
    
    def tearDown(self):
        del self.domain
        del self.loc
    
    def runModel(self, model, PERIODS, TOL):
        for frq in PERIODS:
            Zyx = model.getImpedance(f=frq)
            self.assertIsInstance(Zyx, Data)
            self.assertEqual(Zyx.getShape(), ())
            self.assertEqual(Zyx.getFunctionSpace(), ReducedFunction(self.domain))
            
            
            rho=model.getApparentResitivity(frq, Zyx)
            self.assertIsInstance(rho, Data)
            self.assertEqual(rho.getShape(), ())
            self.assertEqual(rho.getFunctionSpace(), ReducedFunction(self.domain))
            
            
            phi=model.getPhase(frq, Zyx)
            self.assertIsInstance(phi, Data)
            self.assertEqual(phi.getShape(), ())
            self.assertEqual(phi.getFunctionSpace(), ReducedFunction(self.domain))

            #print(frq, self.loc(rho)[0], self.loc(phi)[0])
            self.assertAlmostEqual(self.loc(rho)[0], self.TRUE_RHO, delta=TOL*self.TRUE_RHO)
            self.assertAlmostEqual(self.loc(phi)[0], self.TRUE_PHASE, delta=TOL*self.TRUE_PHASE)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoFloatRhoBFloat(self):
        """
        resistivity set as float, radiation condition
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        model.setResistivity(self.RHO0, rho_boundary=self.RHO0)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoFloatNoRhoB(self):
        """
        resistivity set as float, no value on boundary, radiation condition 
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        model.setResistivity(self.RHO0)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoDataNoRhoB(self):
        """
        resistivity set as Scalar, no value on boundary, radiation condition 
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        rho=Scalar(self.RHO0, ContinuousFunction(self.domain))
        model.setResistivity(rho)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")    
    def test_RhoDataNoRhoBFails(self):
        """
        resistivity set as Scalar, no value on boundary but interpolation is not possible
        """        
        model=MT2DTMModel(self.domain, fixBottom=False, useFastSolver=self.USEFASTSOLVER)
        rho=Scalar(self.RHO0, Function(self.domain))
        self.assertRaises(RuntimeError, model.setResistivity, *(rho, ))
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoFloatRhoBFloatFixedButtom(self):
        """
        resistivity set as float, fixed bottom -> TRUE values are not matched exactly
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        model.setResistivity(self.RHO0, rho_boundary=self.RHO0)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoFloatNoRhoBFixedButtom(self):
        """
        resistivity set as float, no value on boundary, fixed bottom -> TRUE values are not matched exactly
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        model.setResistivity(self.RHO0)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoDataNoRhoBFixedButtom(self):
        """
        resistivity set as Scalar, no value on boundary, fixed bottom -> TRUE values are not matched exactly
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        rho=Scalar(self.RHO0, ContinuousFunction(self.domain))
        model.setResistivity(rho)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")    
    def test_RhoDataNoRhoBFixedButtom2(self):
        """
        resistivity set as Scalar, no value on boundary, interpolation is not possible but fixed bottom -> TRUE values are not matched exactly
        """        
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        model=MT2DTMModel(self.domain, fixBottom=True, useFastSolver=self.USEFASTSOLVER)
        rho=Scalar(self.RHO0, Function(self.domain))
        model.setResistivity(rho)
        self.runModel(model, PERIODS, TOL=1.e-1)  
        
class TestMT2DTMWithAirLayer(unittest.TestCase):
    
    DEPTH=3000.
    WIDTH=1000.
    RHO0=100.
    NEX=100
    NEZ=300
    STATION_OFFSET=1000.
    AIR_LAYER=300.
    
    TRUE_PHASE=45.
    TRUE_RHO=100.
    USEFASTSOLVER=False # True works not yet
    
    def setUp(self):
        self.domain=Rectangle(self.NEX,self.NEZ,l0=self.WIDTH,l1=self.DEPTH)
        self.loc=Locator(ReducedFunction(self.domain), [ (self.STATION_OFFSET,self.DEPTH-self.AIR_LAYER-self.DEPTH/self.NEZ/2)] )
        self.airLayerMask=whereNonNegative(self.domain.getX()[1]-self.DEPTH+self.AIR_LAYER)
        self.airLayerMaskCenter=whereNonNegative(ReducedFunction(self.domain).getX()[1]-self.DEPTH+self.AIR_LAYER)
        
    
    def tearDown(self):
        del self.domain
        del self.loc
    
    def runModel(self, model, PERIODS, TOL):
        for frq in PERIODS:
            Zyx = model.getImpedance(f=frq)
            self.assertIsInstance(Zyx, Data)
            self.assertEqual(Zyx.getShape(), ())
            self.assertEqual(Zyx.getFunctionSpace(), ReducedFunction(self.domain))
            
            
            rho=model.getApparentResitivity(frq, Zyx)
            self.assertIsInstance(rho, Data)
            self.assertEqual(rho.getShape(), ())
            self.assertEqual(rho.getFunctionSpace(), ReducedFunction(self.domain))
            
            
            phi=model.getPhase(frq, Zyx)
            self.assertIsInstance(phi, Data)
            self.assertEqual(phi.getShape(), ())
            self.assertEqual(phi.getFunctionSpace(), ReducedFunction(self.domain))

            #print(frq, self.loc(rho)[0], self.loc(phi)[0])
            self.assertAlmostEqual(self.loc(rho)[0], self.TRUE_RHO, delta=TOL*self.TRUE_RHO)
            self.assertAlmostEqual(self.loc(phi)[0], self.TRUE_PHASE, delta=TOL*self.TRUE_PHASE)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")            
    def test_RhoBFloat(self):
        """
        resistivity set on elements, radiation condition
        """
        PERIODS=np.logspace(-2, 2, num=5, endpoint=True, base=10.0, dtype=float)
        rho=self.airLayerMaskCenter*9999999+(1-self.airLayerMaskCenter)*self.RHO0
        model=MT2DTMModel(self.domain, fixBottom=False, airLayer=self.DEPTH-self.AIR_LAYER, useFastSolver=self.USEFASTSOLVER)
        self.assertEqual(0, Lsup(model.airLayer-self.airLayerMask))
        model.setResistivity(rho, rho_boundary=self.RHO0)
        self.runModel(model, PERIODS, TOL=1.e-3)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoBDataFixedButtom(self):
        """
        resistivity set on elements, no value on boundary but interpolation is not possible fixed bottom 
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        rho=self.airLayerMaskCenter*9999999+(1-self.airLayerMaskCenter)*self.RHO0
        model=MT2DTMModel(self.domain, fixBottom=True, airLayer=self.DEPTH-self.AIR_LAYER, useFastSolver=self.USEFASTSOLVER)
        self.assertEqual(0, Lsup(model.airLayer-self.airLayerMask))
        model.setResistivity(rho)
        self.runModel(model, PERIODS, TOL=1.e-1)
            
    @unittest.skipIf(NO_TRILINOS, "requires Trilinos")
    def test_RhoBDataFailed(self):
        """
        resistivity set on elements, no value on boundary but interpolation is not possible fixed bottom 
        """
        PERIODS=np.logspace(2, 3, num=3, endpoint=True, base=10.0, dtype=float)
        rho=self.airLayerMaskCenter*9999999+(1-self.airLayerMaskCenter)*self.RHO0
        model=MT2DTMModel(self.domain, fixBottom=False, airLayer=self.DEPTH-self.AIR_LAYER, useFastSolver=self.USEFASTSOLVER)
        self.assertEqual(0, Lsup(model.airLayer-self.airLayerMask))
        self.assertRaises(RuntimeError, model.setResistivity, *(rho, ))
        
            
if __name__ == '__main__':

    mySuite=unittest.TestSuite()
    mySuite.addTest(unittest.makeSuite(TestMT2DTE)) 
    mySuite.addTest(unittest.makeSuite(TestMT2DTMNoAirLayer))
    mySuite.addTest(unittest.makeSuite(TestMT2DTMWithAirLayer))
    runner=unittest.TextTestRunner(verbosity=2)
    runner.run(mySuite)
