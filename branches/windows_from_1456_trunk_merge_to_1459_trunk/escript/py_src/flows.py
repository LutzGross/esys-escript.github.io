# $Id:$
#
#######################################################
#
#       Copyright 2008 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
Some models for flow

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2008 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision:$"
__date__="$Date:$"

from escript import *
import util
from linearPDEs import LinearPDE
from pdetools import HomogeneousSaddlePointProblem

class StokesProblemCartesian(HomogeneousSaddlePointProblem):
      """
      solves 

          -(eta*(u_{i,j}+u_{j,i}))_j - p_i = f_i
                u_{i,i}=0

          u=0 where  fixed_u_mask>0
          eta*(u_{i,j}+u_{j,i})*n_j=surface_stress 

      if surface_stress is not give 0 is assumed. 

      typical usage:

            sp=StokesProblemCartesian(domain)
            sp.setTolerance()
            sp.initialize(...)
            v,p=sp.solve(v0,p0)
      """
      def __init__(self,domain,**kwargs):
         HomogeneousSaddlePointProblem.__init__(self,**kwargs)
         self.domain=domain
         self.vol=util.integrate(1.,Function(self.domain))
         self.__pde_u=LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_u.setSymmetryOn()
         self.__pde_u.setSolverMethod(preconditioner=LinearPDE.ILU0)
            
         self.__pde_prec=LinearPDE(domain)
         self.__pde_prec.setReducedOrderOn()
         self.__pde_prec.setSymmetryOn()

         self.__pde_proj=LinearPDE(domain)
         self.__pde_proj.setReducedOrderOn()
         self.__pde_proj.setSymmetryOn()
         self.__pde_proj.setValue(D=1.)

      def initialize(self,f=Data(),fixed_u_mask=Data(),eta=1,surface_stress=Data()):
        self.eta=eta
        A =self.__pde_u.createCoefficientOfGeneralPDE("A")
	self.__pde_u.setValue(A=Data())
        for i in range(self.domain.getDim()):
		for j in range(self.domain.getDim()):
			A[i,j,j,i] += 1. 
			A[i,j,i,j] += 1.
	self.__pde_prec.setValue(D=1./self.eta)
        self.__pde_u.setValue(A=A*self.eta,q=fixed_u_mask,Y=f,y=surface_stress)

      def B(self,arg):
         d=util.div(arg)
         self.__pde_proj.setValue(Y=d)
         self.__pde_proj.setTolerance(self.getSubProblemTolerance())
         return self.__pde_proj.getSolution(verbose=self.show_details)

      def inner(self,p0,p1):
         s0=util.interpolate(p0,Function(self.domain))
         s1=util.interpolate(p1,Function(self.domain))
         return util.integrate(s0*s1)

      def getStress(self,u):
         mg=util.grad(u)
         return 2.*self.eta*util.symmetric(mg)

      def solve_A(self,u,p):
         """
         solves Av=f-Au-B^*p (v=0 on fixed_u_mask)
         """
         self.__pde_u.setTolerance(self.getSubProblemTolerance())
         self.__pde_u.setValue(X=-self.getStress(u)+p*util.kronecker(self.domain))
         return  self.__pde_u.getSolution(verbose=self.show_details)

      def solve_prec(self,p):
         self.__pde_prec.setTolerance(self.getSubProblemTolerance())
         self.__pde_prec.setValue(Y=p)
         q=self.__pde_prec.getSolution(verbose=self.show_details)
         return q
      def stoppingcriterium(self,Bv,v,p):
          n_r=util.sqrt(self.inner(Bv,Bv))
          n_v=util.Lsup(v)
          if self.verbose: print "PCG step %s: L2(div(v)) = %s, Lsup(v)=%s"%(self.iter,n_r,n_v) 
          self.iter+=1
          if n_r <= self.vol**(1./2.-1./self.domain.getDim())*n_v*self.getTolerance():
              if self.verbose: print "PCG terminated after %s steps."%self.iter
              return True
          else:
              return False
