# $Id$

"""
calculation of the stress distribution around a fault from the slip on the fault

e.g. use slip_stress_mesh.py to generate mesh

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, Louise Kettle"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"

from esys.escript import *
from esys.escript.pdetools import SaddlePointProblem
from esys.escript.linearPDEs import LinearPDE
from esys.finley import ReadMesh


rho=0.
lam_lmbd=1.
lam_mu=1.
g=9.81

class SlippingFault(SaddlePointProblem):
      """
      simple example of saddle point problem
      """
      def __init__(self,domain):
         super(SlippingFault, self).__init__(self)
         self.domain=domain
         self.__pde_u=LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_u.setSymmetryOn()

      def initialize(self,density=1.,lmbd=1., mu=1., traction=Data(),fixed_u_mask=Data(), slip=0.):
         d=self.domain.getDim()
         self.slip=slip
         A =self.__pde_u.createCoefficientOfGeneralPDE("A")
         for i in range(self.domain.getDim()):
           for j in range(self.domain.getDim()):
             A[i,j,j,i] += mu
             A[i,j,i,j] += mu
             A[i,i,j,j] += lmbd
         self.__pde_u.setValue(A=A,q=fixed_u_mask,Y=-kronecker(Function(self.domain))[d-1]*g*density,y=traction)

      def inner(self,p0,p1):
         return integrate(p0*p1,FunctionOnContactZero(self.domain))

      def solve_f(self,u,p,tol=1.e-8):
         self.__pde_u.setTolerance(tol)
         self.__pde_u.setValue(y_contact=-p)
         return  self.__pde_u.getSolution()

      def solve_g(self,u,tol=1.e-8):
         dp=-(self.slip-jump(u))
         return  dp


s=numarray.array([0.,1.,1.])
dom=ReadMesh("meshfault3D.fly")
prop=SlippingFault(dom)
d=dom.getDim()
x=dom.getX()[d-1]
mask=whereZero(x-inf(x))*numarray.ones((d,))
u0=Vector(0.,Solution(dom))
p0=Vector(1.,FunctionOnContactZero(dom))
prop.initialize(fixed_u_mask=mask,slip=Data(s,FunctionOnContactZero(dom)), density=rho,lmbd=lam_lmbd, mu=lam_mu)
u,p=prop.solve(u0,p0,iter_max=50,tolerance=0.01)
saveVTK("dis.xml",u=u)
saveVTK("fault.xml",sigma=p,s=jump(u))
