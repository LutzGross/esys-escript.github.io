
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

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
solvers for the stokes problem

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.escript import *
from esys.escript.pdetools import SaddlePointProblem
from esys.escript.linearPDEs import LinearPDE
from esys.finley import Rectangle
from esys.weipa import saveVTK

class SimpleStokesProblem(SaddlePointProblem):
      """
      simple example of saddle point problem
      """
      def __init__(self,domain):
         super(SimpleStokesProblem, self).__init__(self)

         self.__pde_u=LinearPDE(domain)
         self.__pde_u.setSymmetryOn()
         self.__pde_u.setValue(A=identityTensor4(dom))

         self.__pde_p=LinearPDE(domain)
         self.__pde_p.setReducedOrderOn()
         self.__pde_p.setSymmetryOn()
         self.__pde_p.setValue(D=1.)

      def initialize(self,f=Data(),fixed_u_mask=Data()):
         self.__pde_u.setValue(q=fixed_u_mask,Y=f)
      def inner(self,p0,p1):
         return integrate(p0*p1,Function(self.__pde_p.getDomain()))

      def solve_f(self,u,p,tol=1.e-8):
         self.__pde_u.setTolerance(tol)
         self.__pde_u.setValue(X=grad(u)+p*kronecker(self.__pde_u.getDomain()))
         return  self.__pde_u.getSolution()
      def solve_g(self,u,tol=1.e-8):
         self.__pde_p.setTolerance(tol)
         self.__pde_p.setValue(X=-u) 
         dp=self.__pde_p.getSolution()
         return  dp

class StokesProblem(SaddlePointProblem):
      """
      simple example of saddle point problem
      """
      def __init__(self,domain):
         super(StokesProblem, self).__init__(self)
         self.domain=domain
         self.__pde_u=LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_u.setSymmetryOn()

         self.__pde_p=LinearPDE(domain)
         self.__pde_p.setReducedOrderOn()
         self.__pde_p.setSymmetryOn()

      def initialize(self,f=Data(),fixed_u_mask=Data(),eta=1):
         self.eta=eta
         A =self.__pde_u.createCoefficientOfGeneralPDE("A")
         for i in range(self.domain.getDim()):
           for j in range(self.domain.getDim()):
             A[i,j,j,i] += self.eta
             A[i,j,i,j] += self.eta
         self.__pde_p.setValue(D=1./self.eta)
         self.__pde_u.setValue(A=A,q=fixed_u_mask,Y=f)

      def inner(self,p0,p1):
         return integrate(p0*p1,Function(self.__pde_p.getDomain()))

      def solve_f(self,u,p,tol=1.e-8):
         self.__pde_u.setTolerance(tol)
         g=grad(u)
         self.__pde_u.setValue(X=self.eta*symmetric(g)+p*kronecker(self.__pde_u.getDomain()))
         return  self.__pde_u.getSolution()

      def solve_g(self,u,tol=1.e-8):
         self.__pde_p.setTolerance(tol)
         self.__pde_p.setValue(X=-u) 
         dp=self.__pde_p.getSolution()
         return  dp

NE=50
dom=Rectangle(NE,NE,order=2)
# prop=SimpleStokesProblem(dom)
prop=StokesProblem(dom)
x=dom.getX()
mask=(whereZero(x[0])+whereZero(x[0]-1.)+whereZero(x[1]-1.))*unitVector(0,dom)+(whereZero(x[1]-1.)+whereZero(x[1]))*unitVector(1,dom)
u0=Vector(0.,Solution(dom))
u0[0]=x[1]*whereZero(x[1]-1.)
p0=Scalar(0,ReducedSolution(dom))
# prop.initialize(fixed_u_mask=mask)
prop.initialize(fixed_u_mask=mask,eta=10.)
u,p=prop.solve(u0,p0,tolerance=0.01)
# saveVTK("stokes.vtu",u=u,p=p,m=mask,u0=u0)

eta=whereNegative(x[1]-0.5)*1.e6+whereNonNegative(x[1]-0.5)
prop.initialize(fixed_u_mask=mask,eta=eta)
u,p=prop.solve(u0,p0,tolerance=0.01,tolerance_u=0.1,accepted_reduction=0.8)
saveVTK("stokes.vtu",u=u,p=p,m=mask,u0=u0)
          
# vim: expandtab shiftwidth=4:
