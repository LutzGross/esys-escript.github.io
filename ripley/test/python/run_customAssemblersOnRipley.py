
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.ripley import Rectangle, Brick
from esys.escript.linearPDEs import LameEquation, LinearPDESystem, WavePDE
#TODO from esys.downunder import HTIWave, VTIWave, Ricker

EXPANDED, SCALAR, CONSTANT = range(3)

class RipleyWaveAssemblerTestBase(unittest.TestCase):
    def generate_fast_HTI_PDE_solution(self):
        pde = WavePDE(self.domain, [("c11", self.c11),
                    ("c23", self.c23), ("c13", self.c13), ("c33", self.c33),
                    ("c44", self.c44), ("c66", self.c66)])
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()
        dim = pde.getDim()
        X = self.domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(self.domain))
        du = grad(X*X)
#        D = Scalar(2500., Function(self.domain))*kronecker(dim)
        D = 2500.*kronecker(dim)

        pde.setValue(D=D, y_dirac=y, du=du)
        return pde.getSolution()

    def generate_slow_HTI_PDE_solution(self):
        pde = LinearPDESystem(self.domain)
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()

        dim = pde.getDim()
        X = self.domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(self.domain))
        du = grad(X*X)
        D = 2500.*kronecker(dim)
        pde.setValue(X=pde.createCoefficient('X'))
        sigma = pde.getCoefficient('X')
        if dim == 3:
            e11=du[0,0]
            e22=du[1,1]
            e33=du[2,2]

            sigma[0,0]=self.c11*e11+self.c13*(e22+e33)
            sigma[1,1]=self.c13*e11+self.c33*e22+self.c23*e33
            sigma[2,2]=self.c13*e11+self.c23*e22+self.c33*e33

            s=self.c44*(du[2,1]+du[1,2])
            sigma[1,2]=s
            sigma[2,1]=s

            s=self.c66*(du[2,0]+du[0,2])
            sigma[0,2]=s
            sigma[2,0]=s

            s=self.c66*(du[0,1]+du[1,0])
            sigma[0,1]=s
            sigma[1,0]=s

        else:
            e11=du[0,0]
            e22=du[1,1]
            sigma[0,0]=self.c11*e11+self.c13*e22
            sigma[1,1]=self.c13*e11+self.c33*e22

            s=self.c66*(du[1,0]+du[0,1])
            sigma[0,1]=s
            sigma[1,0]=s

        pde.setValue(D=D, X=-sigma, y_dirac=y)
        return pde.getSolution()

    def generate_fast_VTI_PDE_solution(self):
        pde = WavePDE(self.domain, [("c11", self.c11),
                    ("c12", self.c12), ("c13", self.c13), ("c33", self.c33),
                    ("c44", self.c44), ("c66", self.c66)])
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()
        dim = pde.getDim()

        X = self.domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(self.domain))
        du = grad(X*X)
        D = 2500.*kronecker(dim)

        pde.setValue(D=D, y_dirac=y, du=du)
        return pde.getSolution()

    def generate_slow_VTI_PDE_solution(self):
        pde = LinearPDESystem(self.domain)
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()
        dim = pde.getDim()

        dim = pde.getDim()
        X = self.domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(self.domain))
        du = grad(X*X)
        D = 2500.*kronecker(dim)
        pde.setValue(X=pde.createCoefficient('X'))
        sigma = pde.getCoefficient('X')
        
        if dim == 3:
            e11=du[0,0]
            e22=du[1,1]
            e33=du[2,2]
            sigma[0,0]=self.c11*e11+self.c12*e22+self.c13*e33
            sigma[1,1]=self.c12*e11+self.c11*e22+self.c13*e33
            sigma[2,2]=self.c13*(e11+e22)+self.c33*e33

            s=self.c44*(du[2,1]+du[1,2])
            sigma[1,2]=s
            sigma[2,1]=s             

            s=self.c44*(du[2,0]+du[0,2])
            sigma[0,2]=s
            sigma[2,0]=s

            s=self.c66*(du[0,1]+du[1,0])
            sigma[0,1]=s
            sigma[1,0]=s
        else:
            e11=du[0,0]
            e22=du[1,1]
            sigma[0,0]=self.c11*e11+self.c13*e22
            sigma[1,1]=self.c13*e11+self.c33*e22
            s=self.c44*(du[1,0]+du[0,1])
            sigma[0,1]=s
            sigma[1,0]=s

        pde.setValue(D=D, X=-sigma, y_dirac=y)
        return pde.getSolution()

    def run_HTI_assembly(self):
        model = HTIWave(self.domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True,
                disable_fast_assemblers=True)
        self.assertFalse(model.fastAssembler) #ensure the arg is obeyed

        model = HTIWave(self.domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True)
        self.assertTrue(model.fastAssembler) #ensure fast is actually used

        slow = self.generate_slow_HTI_PDE_solution()
        fast = self.generate_fast_HTI_PDE_solution()

        self.assertLess(Lsup(fast - slow), 1e-12) #comparison between them

    def run_VTI_assembly(self):
        model = VTIWave(self.domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True,
                disable_fast_assemblers=True)
        self.assertFalse(model.fastAssembler) #ensure the arg is obeyed

        model = VTIWave(self.domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True)
        self.assertTrue(model.fastAssembler) #ensure fast is actually used
        
        slow = self.generate_slow_VTI_PDE_solution()
        fast = self.generate_fast_VTI_PDE_solution()

        self.assertLess(Lsup(fast - slow), 1e-12) #comparison between them

    def test_Function_params(self):
        self.V_p = Data(2500., (), Function(self.domain))
        self.V_s = Data(1250., (), Function(self.domain))
        self.c11 = Data(11., (), Function(self.domain))
        self.c12 = Data(12., (), Function(self.domain))
        self.c13 = Data(13., (), Function(self.domain))
        self.c23 = Data(23., (), Function(self.domain))
        self.c33 = Data(33., (), Function(self.domain))
        self.c44 = Data(44., (), Function(self.domain))
        self.c66 = Data(66., (), Function(self.domain))
        for i in [self.V_p, self.V_s, self.c11, self.c12, self.c13, self.c23,
                self.c33, self.c44, self.c66]:
            i.expand()
        self.run_HTI_assembly()
        self.run_VTI_assembly()
    
    def test_ReducedFunction_params(self):
        self.V_p = Data(2500., (), ReducedFunction(self.domain))
        self.V_s = Data(1250., (), ReducedFunction(self.domain))
        self.c11 = Data(11., (), ReducedFunction(self.domain))
        self.c12 = Data(12., (), ReducedFunction(self.domain))
        self.c13 = Data(13., (), ReducedFunction(self.domain))
        self.c23 = Data(23., (), ReducedFunction(self.domain))
        self.c33 = Data(33., (), ReducedFunction(self.domain))
        self.c44 = Data(44., (), ReducedFunction(self.domain))
        self.c66 = Data(66., (), ReducedFunction(self.domain))
        for i in [self.V_p, self.V_s, self.c11, self.c12, self.c13, self.c23,
                self.c33, self.c44, self.c66]:
            i.expand()
        with self.assertRaises(ValueError) as e:
            self.run_HTI_assembly()
        self.assertTrue("mismatching function spaces" in str(e.exception))
        with self.assertRaises(ValueError) as e:
            self.run_VTI_assembly()
        self.assertTrue("mismatching function spaces" in str(e.exception))

    def test_Constant_params(self):
        self.V_p = Scalar(2500., Function(self.domain))
        self.V_s = Scalar(1250., Function(self.domain))
        self.c11 = Scalar(11., ReducedFunction(self.domain))
        self.c12 = Scalar(12., ReducedFunction(self.domain))
        self.c13 = Scalar(13., ReducedFunction(self.domain))
        self.c23 = Scalar(23., ReducedFunction(self.domain))
        self.c33 = Scalar(33., ReducedFunction(self.domain))
        self.c44 = Scalar(44., ReducedFunction(self.domain))
        self.c66 = Scalar(66., ReducedFunction(self.domain))
        with self.assertRaises(ValueError) as e:
            self.run_HTI_assembly()
        self.assertTrue("mismatching function spaces" in str(e.exception))
        with self.assertRaises(ValueError) as e:
            self.run_VTI_assembly()
        self.assertTrue("mismatching function spaces" in str(e.exception))
@unittest.skip("Ripley wave solver 2D skipping")
class Test_RipleyWaveAssembler2D(RipleyWaveAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(20,20,l0=100.,l1=100., diracTags=["source"],
                diracPoints=[(0,0)])
        self.wavelet = Ricker(100.)
        
    def tearDown(self):
        del self.domain

@unittest.skip("Ripley wave solver #D skipping")
class Test_RipleyWaveAssembler3D(RipleyWaveAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(10,10,10,l0=100.,l1=100., l2=100.,
                diracTags=["source"], diracPoints=[(0,0,0)])
        self.wavelet = Ricker(100.)

    def tearDown(self):
        del self.domain


class RipleyLameAssemblerTestBase(unittest.TestCase): #requires subclassing
    def run_lame(self, fast, test_type, mu=3, lamb=50):
        d=self.domain.getDim()
        mypde=LameEquation(self.domain, useFast=fast)
        cf=ContinuousFunction(self.domain)
        x=cf.getX()
        u_ex=x
        msk=Vector(0.,cf)
        for i in range(d):
            msk[i]=whereZero(x[i])
        if test_type != CONSTANT:
            mu = Scalar(mu, cf)
            lamb = Scalar(lamb, cf)
            if test_type == EXPANDED:
                mu.expand()
                lamb.expand()
        mypde.setValue(q=msk,r=u_ex,lame_mu=mu,lame_lambda=lamb,f=(2*3+50*d)*FunctionOnBoundary(self.domain).getNormal())
        return mypde.getSolution()

    def test_lameExpanded(self):
        #check default and lame assemblers agree
        default = self.run_lame(False, EXPANDED)
        lame = self.run_lame(True, EXPANDED)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "expanded data"%self.domain.getDim())
        #reverse order, ensure default assembler still operational
        lame = self.run_lame(True, EXPANDED, mu=7, lamb=40)
        default = self.run_lame(False, EXPANDED, mu=7, lamb=40)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "expanded data"%self.domain.getDim())

    def test_lameScalar(self):
        #check default and lame assemblers agree
        default = self.run_lame(False, SCALAR)
        lame = self.run_lame(True, SCALAR)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "scalar data"%self.domain.getDim())
        #reverse order, ensure default assembler still operational
        lame = self.run_lame(True, SCALAR, mu=7, lamb=40)
        default = self.run_lame(False, SCALAR, mu=7, lamb=40)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "scalar data"%self.domain.getDim())

    def test_lameConstant(self):
        #check default and lame assemblers agree
        default = self.run_lame(False, CONSTANT)
        lame = self.run_lame(True, CONSTANT)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "constant data"%self.domain.getDim())
        #reverse order, ensure default assembler still operational
        lame = self.run_lame(True, CONSTANT, mu=7, lamb=40)
        default = self.run_lame(False, CONSTANT, mu=7, lamb=40)
        self.assertLess(Lsup(default - lame), 1e-8,
                "Default and Lame %dDassembler solutions differ for "\
                "constant data"%self.domain.getDim())

class Test_RipleyLameAssemblers2D(RipleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = Rectangle(20,20)

    def tearDown(self):
        del self.domain

class Test_RipleyLameAssemblers3D(RipleyLameAssemblerTestBase):
    def setUp(self):
        self.domain = Brick(10,10,10)

    def tearDown(self):
        del self.domain


if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

