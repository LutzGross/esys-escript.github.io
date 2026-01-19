
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
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
from esys.speckley import Rectangle, Brick
from esys.escript.linearPDEs import LameEquation, LinearPDESystem, WavePDE, LinearSinglePDE
from esys.escript.pdetools import Ricker
#TODO from esys.downunder import HTIWave, VTIWave

EXPANDED, SCALAR, CONSTANT = range(3)

class SpeckleyWaveAssemblerTestBase(unittest.TestCase):
    def generate_fast_HTI_PDE_solution(self, domain):
        pde = WavePDE(domain, [("c11", self.c11),
                    ("c23", self.c23), ("c13", self.c13), ("c33", self.c33),
                    ("c44", self.c44), ("c66", self.c66)])
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()
        dim = pde.getDim()
        X = domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(domain))
        du = grad(X*X)
        D = 2500.*kronecker(dim)

        pde.setValue(D=D, y_dirac=y, du=du)
        return pde.getSolution()

    def generate_slow_HTI_PDE_solution(self, domain):
        pde = LinearPDESystem(domain)
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()

        dim = pde.getDim()
        X = domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(domain))
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

    def generate_fast_VTI_PDE_solution(self, domain):
        pde = WavePDE(domain, [("c11", self.c11),
                    ("c12", self.c12), ("c13", self.c13), ("c33", self.c33),
                    ("c44", self.c44), ("c66", self.c66)])
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()
        dim = pde.getDim()

        X = domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(domain))
        du = grad(X*X)
        D = 2500.*kronecker(dim)

        pde.setValue(D=D, y_dirac=y, du=du)
        return pde.getSolution()

    def generate_slow_VTI_PDE_solution(self, domain):
        pde = LinearPDESystem(domain)
        pde.getSolverOptions().setSolverMethod(SolverOptions.HRZ_LUMPING)
        pde.setSymmetryOn()
        dim = pde.getDim()

        dim = pde.getDim()
        X = domain.getX()
        y = Vector([2.,3.,4.][:dim], DiracDeltaFunctions(domain))
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
    @unittest.skip("is failing")
    def run_HTI_assembly(self, domain):
        model = HTIWave(domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True,
                disable_fast_assemblers=True)
        self.assertFalse(model.fastAssembler) #ensure the arg is obeyed

        model = HTIWave(domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True)
        self.assertTrue(model.fastAssembler) #ensure fast is actually used

        slow = self.generate_slow_HTI_PDE_solution(domain)
        fast = self.generate_fast_HTI_PDE_solution(domain)

        self.assertLess(Lsup(fast - slow), 1e-12*Lsup(slow)) #comparison between them
    @unittest.skip("is failing")
    def run_VTI_assembly(self, domain):
        model = VTIWave(domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True,
                disable_fast_assemblers=True)
        self.assertFalse(model.fastAssembler) #ensure the arg is obeyed

        model = VTIWave(domain, self.V_p, self.V_s, self.wavelet, "source",
                source_vector=[0,0,1], eps=0., gamma=0., delta=0.,
                rho=2000., absorption_zone=None, lumping=True)
        self.assertTrue(model.fastAssembler) #ensure fast is actually used

        slow = self.generate_slow_VTI_PDE_solution(domain)
        fast = self.generate_fast_VTI_PDE_solution(domain)

        self.assertLess(Lsup(fast - slow), 1e-12*Lsup(slow)) #comparison between them

    def test_Function_params(self):
        for domain in self.domains:
            self.V_p = Data(2500., (), Function(domain))
            self.V_s = Data(1250., (), Function(domain))
            self.c11 = Data(11., (), Function(domain))
            self.c12 = Data(12., (), Function(domain))
            self.c13 = Data(13., (), Function(domain))
            self.c23 = Data(23., (), Function(domain))
            self.c33 = Data(33., (), Function(domain))
            self.c44 = Data(44., (), Function(domain))
            self.c66 = Data(66., (), Function(domain))
            for i in [self.V_p, self.V_s, self.c11, self.c12, self.c13, self.c23,
                    self.c33, self.c44, self.c66]:
                i.expand()
            with self.assertRaises(RuntimeError) as e:
                self.run_HTI_assembly(domain)
            self.assertTrue("C tensor elements must be reduced" in str(e.exception))
            with self.assertRaises(RuntimeError) as e:
                self.run_VTI_assembly(domain)
            self.assertTrue("C tensor elements must be reduced" in str(e.exception))

    def test_ReducedFunction_params(self):
        for domain in self.domains:
            self.V_p = Data(2500., (), ReducedFunction(domain))
            self.V_s = Data(1250., (), ReducedFunction(domain))
            self.c11 = Data(11., (), ReducedFunction(domain))
            self.c12 = Data(12., (), ReducedFunction(domain))
            self.c13 = Data(13., (), ReducedFunction(domain))
            self.c23 = Data(23., (), ReducedFunction(domain))
            self.c33 = Data(33., (), ReducedFunction(domain))
            self.c44 = Data(44., (), ReducedFunction(domain))
            self.c66 = Data(66., (), ReducedFunction(domain))
            for i in [self.V_p, self.V_s, self.c11, self.c12, self.c13, self.c23,
                    self.c33, self.c44, self.c66]:
                i.expand()
            self.run_HTI_assembly(domain)
            self.run_VTI_assembly(domain)

    def test_Constant_params(self):
        for domain in self.domains:
            self.V_p = Scalar(2500., ReducedFunction(domain))
            self.V_s = Scalar(1250., ReducedFunction(domain))
            self.c11 = Scalar(11., ReducedFunction(domain))
            self.c12 = Scalar(12., ReducedFunction(domain))
            self.c13 = Scalar(13., ReducedFunction(domain))
            self.c23 = Scalar(23., ReducedFunction(domain))
            self.c33 = Scalar(33., ReducedFunction(domain))
            self.c44 = Scalar(44., ReducedFunction(domain))
            self.c66 = Scalar(66., ReducedFunction(domain))
            self.run_HTI_assembly(domain)
            self.run_VTI_assembly(domain)

class Test_SpeckleyWaveAssembler2D(SpeckleyWaveAssemblerTestBase):
    def setUp(self):
        self.domains = []
        for order in range(2,11):
            self.domains.append(Rectangle(order,10,10,l0=100,l1=100,diracTags=["source"],
                    diracPoints=[(0,0)]))
        self.wavelet = Ricker(100.)

    def tearDown(self):
        del self.domains

class Test_SpeckleyWaveAssembler3D(SpeckleyWaveAssemblerTestBase):
    def setUp(self):
        self.domains = []
        for order in range(2,11):
            self.domains.append(Brick(order, 4,4,4,
                    diracTags=["source"], diracPoints=[(0,0,0)]))
        self.wavelet = Ricker(100.)

    def tearDown(self):
        del self.domains

class Test_Complex_Assembler(unittest.TestCase):
    def test_complex_params_Rectangle(self):
        domain=Rectangle(order=2,n0=10,n1=10)
        pde = LinearSinglePDE(domain, isComplex=True)
        pde.setValue(D=1j)
        pde.setValue(Y=1.0)
        self.assertTrue(Lsup(pde.getSolution())==1.0, "Failed test_complex_params_Rectangle")
        del domain

    def test_complex_params_Brick(self):
        domain=Brick(order=2,n0=10,n1=10,n2=10)
        pde = LinearSinglePDE(domain, isComplex=True)
        pde.setValue(D=1j)
        pde.setValue(Y=1.0)
        self.assertTrue(Lsup(pde.getSolution())==1.0, "Failed test_complex_params_Brick")
        del domain

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
