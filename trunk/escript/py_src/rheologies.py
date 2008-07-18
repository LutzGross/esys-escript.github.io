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
from pdetools import HomogeneousSaddlePointProblem,Projector

class PlateMantelModel(HomogeneousSaddlePointProblem):
      """
      This class implements the reology of a self-consistent plate-mantel model proposed in
      U{Hans-Bernd Muhlhaus<emailto:h.muhlhaus@uq.edu.au>}  and U{Klaus Regenauer-Lieb<mailto:klaus.regenauer-lieb@csiro.au>}: 
      I{Towards a self-consistent plate mantle model that includes elasticity: simple benchmarks and application to basic modes of convection}, 
      see U{doi: 10.1111/j.1365-246X.2005.02742.x<http://www3.interscience.wiley.com/journal/118661486/abstract?CRETRY=1&SRETRY=0>}.

      typical usage:

            sp=PlateMantelModel(domain,stress=0,v=0)
            sp.setTolerance()
            sp.initialize(...)
            v,p=sp.solve(v0,p0)
      """
      def __init__(self, domain, stress=0, v=0, p=0, t=0, useJaumannStress=True, **kwargs):
         """
         initializes PlateMantelModel 

         @param domain: problem domain
         @type domain: L{domain}
         @param stress: initial stress
         @param v: initial velocity field
         @param t: initial time
         @param useJaumannStress: C{True} if Jaumann stress is used
         """
         HomogeneousSaddlePointProblem.__init__(self,**kwargs)
         self.__domain=domain
         self.__t=t
         self.setSmall()
         self.__vol=util.integrate(1.,Function(self.__domain))
         self.__useJaumannStress=useJaumannStress
         #=======================
         # state variables:
         #
         if isinstance(stress,Data):
            self.__stress=util.interpolate(stress,Function(domain))
         else:
            self.__stress=Data(stress,(domain.getDim(),domain.getDim()),Function(domain))
         self.__stress-=util.trace(self.__stress)*(util.kronecker(domain)/domain.getDim())
         if isinstance(v,Data):
            self.__v=util.interpolate(v,Solution(domain))
         else:
            self.__v=Data(v,(domain.getDim(),),Solution(domain)) 
         if isinstance(p,Data):
            self.__p=util.interpolate(p,ReducedSolution(domain))
         else:
            self.__p=Data(p,(),ReducedSolution(domain))
         self.__tau=util.sqrt(0.5)*util.length(self.__stress)
         self.__D=util.symmetric(util.grad(self.__v))
         #=======================
         #  parameters
         #
         self.__mu=None
         self.__eta_N=None
         self.__eta_0=None
         self.__tau_0=None
         self.__n=None
         self.__eta_Y=None
         self.__tau_Y=None
         self.__n_Y=None
         #=======================
         # solme value to track:
         #
	 self.__mu_eff=None
	 self.__h_eff=None
	 self.__eta_eff=None
         #=======================
         # PDE related stuff 
         self.__pde_u=LinearPDE(domain,numEquations=self.getDomain().getDim(),numSolutions=self.getDomain().getDim())
         self.__pde_u.setSymmetryOn()
         # self.__pde_u.setSolverMethod(preconditioner=LinearPDE.ILU0)
            
         self.__pde_prec=LinearPDE(domain)
         self.__pde_prec.setReducedOrderOn()
         self.__pde_prec.setSymmetryOn()

         self.__pde_proj=LinearPDE(domain)
         self.__pde_proj.setReducedOrderOn()
         self.__pde_proj.setSymmetryOn()
         self.__pde_proj.setValue(D=1.)

      def useJaumannStress(self):
          """ 
          return C{True} if Jaumann stress is included.
          """
          return self.__useJaumannStress
      def setSmall(self,small=util.sqrt(util.EPSILON)):
          """
          sets small value
      
          @param small: positive small value
          """
          self.__small=small
      def getSmall(self):
          """
          returns small value
          @rtype: positive float
          """
          return self.__small
      def getDomain(self):
          """
          returns the domain
          """
          return self.__domain
      def getStress(self):
          """
          returns current stress
          """
          return self.__stress
      def getStretching(self):
          """
          return stretching
          """
          return self.__D
      def getTau(self):
          """
          returns current second stress deviatoric invariant
          """
          return self.__tau
      def getPressure(self):
          """
          returns current pressure
          """
          return self.__p
      def getVelocity(self):
          """
          returns current velocity
          """
          return self.__v
      def getTime(self):
          """
          returns current time
          """
          return self.__t
      def getVolume(self):
          """
          returns domain volume
          """
          return self.__vol

      def getEtaEffective(self):
          """
          returns effective shear viscocity
          """
          return self.__eta_eff
      def getMuEffective(self):
          """
          returns effective shear modulus
          """
          return self.__mu_eff
      def getHEffective(self):
          """
          returns effective h
          """
          return self.__h_eff
      def getMechanicalPower(self):
          """
          returns the locl mechanical power M{D_ij Stress_ij}
          """
          return util.inner(self.getStress(),self.getStretching())
          
      def initialize(self, mu=None, eta_N=None, eta_0=None, tau_0=None, n=None, eta_Y=None, tau_Y=None, n_Y=None, F=None, q=None):
          """ 
          sets the model parameters.

          @param mu: shear modulus. If C{mu==None} ...
          @param eta_N: Newtonian viscosity. 
          @param eta_0: power law viscovity for tau=tau0
          @param tau_0: reference tau for power low. If C{tau_0==None} constant viscosity is assumed.
          @param n: power of power law. if C{n=None}, constant viscosity is assumed.
          @param eta_Y: =None
          @param tau_Y: =None
          @param n_Y: =None
          @param F: external force.
          @param q: location of constraints
          """
          if mu != None: self.__mu=mu
          if eta_N != None: self.__eta_N=eta_N
          if eta_0 != None: self.__eta_0=eta_0
          if tau_0 != None: self.__tau_0=tau_0
          if n != None: self.__n=n
          if eta_Y != None: self.__eta_Y=eta_Y
          if tau_Y != None: self.__tau_Y=tau_Y
          if n_Y != None: self.__n_Y=n_Y
          if F != None: self.__pde_u.setValue(Y=F)
          if q != None: self.__pde_u.setValue(q=q)

      def updateEffectiveCoefficients(self, dt, tau):
          """
          updates the effective coefficints depending on tau and time step size dt.
          """
          h_vis=None
          h_yie=None
          # calculate eta_eff = ( 1/eta_N + 1/eta_vis + 1/eta_yie)^{-1}
          #    with   eta_vis = eta_0 * (tau/tau_0) ^ (1-n)
          #           eta_yie = eta_Y * (tau/tau_Y) ^ (1-n_Y)
          s=Scalar(0.,Function(self.getDomain()))
          if self.__eta_N!=None:
              s+=1/self.__eta_N
          if self.__eta_0!=None and self.__tau_0 != None and self.__n!=None:
              if self.__tau_0 !=None and  self.__n!=None:
                 eta_vis=self.__eta_0*(tau/self.__tau_0)**(1-self.__n)
                 s+=1./(eta_vis+self.__eta_0*self.getSmall())
                 h_vis=eta_vis/(self.__n-1)
          if self.__tau_Y!=None and self.__eta_Y!=None and self.__n_Y!=None:
               eta_yie=self.__eta_Y*(tau/self.__tau_Y)**(1-self.__n_Y)
               s+=1/(eta_yie+self.getSmall()*self.__eta_Y)
               h_yie=self.__eta_Y/(self.__n_Y-1)
          self.__eta_eff=1/s
          # calculate eta_eff = ( 1/h_vis + 1/h_yie)^{-1}
          #    with   h_vis = eta_vis/(n-1)
          #           h_yie = eta_yie/(n_Y-1)
          if h_vis == None: 
             if h_yie==None:
                self__h_eff=None
             else:
                self__h_eff=h_yie
          else:
             if h_yie==None:
                self__h_eff=h_vis
             else:
                self__h_eff=1/((1./h_vis)+(1./h_yie))
          # calculate mu_eff = ( 1/mu + dt/eta_eff)^{-1} = mu*eta_eff/(mu*dt+eta_eff)
          if self.__mu == None:
             self.__mu_eff=self.__eta_eff/dt
          else:
             self.__mu_eff=1./(1./self.__mu+dt/self.__eta_eff) 

      def update(self,dt,max_inner_iter=20, verbose=False, show_details=False, tol=10., solver="PCG"):
          """
          updates stress, velocity and pressure for time increment dt

          @param dt: time increment
          @param max_inner_iter: maximum number of iteration steps in the incompressible solver
          @param verbose: prints some infos in the incompressible solve
          @param show_details: prints some infos while solving PDEs
          @param tol: tolerance for the time step
          """
          stress_last=self.getStress()
          # we should use something like FGMRES to merge the two iterations!
          e=10.*tol
          # get values from last time step and use them as initial guess:
          v=self.__v
          stress=self.__stress
          p=self.__p
          tau=self.__tau
          while e > tol: # revise stress calculation if this is used.
              #
              #  update the effective coefficients:
              #
              self.updateEffectiveCoefficients(dt,tau)
              eta_eff=self.getEtaEffective()
              mu_eff=self.getMuEffective()
              h_eff=self.getHEffective()
              #
              #   create some temporary variables:
              #
	      self.__pde_u.setValue(A=Data()) # save memory!
              k3=util.kronecker(Function(self.getDomain()))
              k3Xk3=util.outer(k3,k3)
              self.__f3=mu_eff*dt
              A=self.__f3*(util.swap_axes(k3Xk3,0,3)+util.swap_axes(k3Xk3,1,3))
              print "mueff=",util.inf(mu_eff),util.sup(mu_eff)
              if h_eff == None:
                 s=0
                 self.__f0=0
                 self.__f1=mu_eff/eta_eff*dt
              else:
                 s=mu_eff*dt/(h_eff+mu_eff*dt)
                 Lsup_tau=Lsup(tau)
                 if Lsup>0:
                    self.__f0=mu_eff*s*dt/(tau+Lsup_tau*self.getSmall())**2
                    A+=util.outer((-self.__f0)*stress,stress)
                 else:
                    self.__f0=0.
                 self.__f1=mu_eff/eta_eff*(1-s)*dt

              if self.useJaumannStress():
                 self.__f2=mu_eff/eta_eff*dt**2
                 sXk3=util.outer(stress,self.__f2*(-k3/2))

                 A+=util.swap_axes(sXk3,0,3)
                 A+=util.swap_axes(sXk3,1,3)
                 A-=util.swap_axes(sXk3,0,2)
                 A-=util.swap_axes(sXk3,1,2)
              else:
                 self.__f2=0
              self.__pde_u.setValue(A=A)
	      self.__pde_prec.setValue(D=1/mu_eff) 

              print "X f0:",util.inf(self.__f0), util.sup(self.__f0)
              print "X f1:",util.inf(self.__f1), util.sup(self.__f1)
              print "X f2:",util.inf(self.__f2), util.sup(self.__f2)
              print "X f3:",util.inf(self.__f3), util.sup(self.__f3)

              v_old=v
              v,p=self.solve(v,p,max_iter=max_inner_iter, verbose=verbose, show_details=show_details, solver=solver)
              # update stress
              stress=stress_last+dt*self.getStressChange(v)
              stress-=util.trace(stress)*(util.kronecker(self.getDomain())/self.getDomain().getDim())
              tau=util.sqrt(0.5)*util.length(stress)
              # calculate error:
              e=util.Lsup(v_old-v)/util.Lsup(v)
          # state variables:
          self.__t+=dt
          self.__v=v
          self.__p=p
          self.__stress=stress
          self.__tau=tau
          self.__D=util.symmetric(util.grad(v))

      def getStressChange(self,v):
          """
          returns the stress change due to a given velocity field v
          """
          stress=self.getStress()
          g=util.grad(v)
          D=util.symmetric(g)
          W=util.nonsymmetric(g)
          U=util.nonsymmetric(util.tensor_mult(W,stress))
          dstress=2*self.__f3*D-self.__f0*util.inner(stress,D)*stress-self.__f1*stress-2*self.__f2*U
          return dstress

      def B(self,arg):
         """
         div operator
         """
         d=util.div(arg)
         self.__pde_proj.setValue(Y=d)
         self.__pde_proj.setTolerance(self.getSubProblemTolerance())
         return self.__pde_proj.getSolution(verbose=self.show_details)

      def solve_prec(self,p):
         """
         preconditioner
         """
	 #proj=Projector(domain=self.getDomain(), reduce = True, fast=False)
         self.__pde_prec.setTolerance(self.getSubProblemTolerance())
         self.__pde_prec.setValue(Y=p)
         q=self.__pde_prec.getSolution(verbose=self.show_details)
         return q

      def inner(self,p0,p1):
         """
         inner product for pressure
         """
         s0=util.interpolate(p0,Function(self.getDomain()))
         s1=util.interpolate(p1,Function(self.getDomain()))
         return util.integrate(s0*s1)

      def solve_A(self,u,p):
         """
         solves Av=f-Au-B^*p (v=0 on fixed_u_mask)
         """
         self.__pde_u.setTolerance(self.getSubProblemTolerance())
         self.__pde_u.setValue(X=-self.getStressChange(u)-p*util.kronecker(self.getDomain()))
         return  self.__pde_u.getSolution(verbose=self.show_details)


      def stoppingcriterium(self,Bv,v,p):
          n_r=util.sqrt(self.inner(Bv,Bv))
          n_v=util.sqrt(util.integrate(util.length(util.grad(v))**2))
          if self.verbose: print "PCG step %s: L2(div(v)) = %s, L2(grad(v))=%s"%(self.iter,n_r,n_v) 
          if self.iter == 0: self.__n_v=n_v;
          self.__n_v, n_v_old =n_v, self.__n_v
          self.iter+=1
          if self.iter>1 and n_r <= n_v*self.getTolerance() and abs(n_v_old-self.__n_v) <= n_v * self.getTolerance():
              if self.verbose: print "PCG terminated after %s steps."%self.iter
              return True
          else:
              return False
