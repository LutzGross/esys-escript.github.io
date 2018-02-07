
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
calculation of the stress distribution around a fault from the slip on the fault

e.g. use slip_stress_mesh.py to generate mesh

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, Louise Kettle"

from esys.escript import *
from esys.escript.pdetools import SaddlePointProblem
from esys.escript.linearPDEs import LinearPDE
from esys.finley import ReadMesh
from esys.weipa import saveVTK


rho=0.
lam_lmbd=1.7e11
lam_mu=1.7e11
g=9.81
fstart =  [50000.0, 40000.0, 10909.09090909091]
fend =  [50000.0, 60000.0,   19090.909090909092]




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
         return integrate(inner(p0,p1),FunctionOnContactZero(self.domain))

      def solve_f(self,u,p,tol=1.e-8):
         self.__pde_u.setTolerance(tol)
         self.__pde_u.setValue(y_contact=-p)
         # print "p:",inf(p),sup(p)
         # print "u:",inf(u),sup(u)
         self.__pde_u.setValue(y_contact=-p)
         return  self.__pde_u.getSolution()

      def solve_g(self,u,tol=1.e-8):
         dp=Vector(0.,FunctionOnContactZero(self.domain))
         h=FunctionOnContactZero(self.domain).getSize()
         # print jump(u)-self.slip
         dp[0]=(self.slip[0]-jump(u[0]))*lam_mu/h
         dp[1]=(self.slip[1]-jump(u[1]))*lam_mu/h
         dp[2]=(self.slip[2]-jump(u[2]))*lam_mu/h
         return  dp


dom=ReadMesh("meshfault3D.fly",integrationOrder=-1)
prop=SlippingFault(dom)
d=dom.getDim()
x=dom.getX()[0]
# x=dom.getX()[d-1]
mask=whereZero(x-inf(x))*numpy.ones((d,))
x=FunctionOnContactZero(dom).getX()
s=numpy.array([-100000.,1.,1.])
for i in range(3):
     d=fend[i]-fstart[i]
     if d>0:
         q=(x[i]-fstart[i])/d
         s=q*(1-q)*4*s
     elif d<0:
         q=(x[i]-fend[i])/d
         s=q*(1-q)*4*s
u0=Vector(0.,Solution(dom))
p0=Vector(1.,FunctionOnContactZero(dom))
prop.initialize(fixed_u_mask=mask,slip=Data(s,FunctionOnContactZero(dom)), density=rho,lmbd=lam_lmbd, mu=lam_mu)
u,p=prop.solve(u0,p0,iter_max=50,tolerance=0.13,accepted_reduction=1.)
saveVTK("dis.vtu",u=u)
saveVTK("fault.vtu",sigma=p,s=jump(u))
