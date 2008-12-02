
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

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

from escript import *
import util
from linearPDEs import LinearPDE
from pdetools import Defect, NewtonGMRES, ArithmeticTuple

class IncompressibleIsotropicKelvinFlow(Defect):
      """
      This class implements the reology of an isotropic kelvin material. 

      @note: this model has been used in the self-consistent plate-mantel model proposed in
      U{Hans-Bernd Muhlhaus<emailto:h.muhlhaus@uq.edu.au>}  and U{Klaus Regenauer-Lieb<mailto:klaus.regenauer-lieb@csiro.au>}: 
      I{Towards a self-consistent plate mantle model that includes elasticity: simple benchmarks and application to basic modes of convection}, 
      see U{doi: 10.1111/j.1365-246X.2005.02742.x<http://www3.interscience.wiley.com/journal/118661486/abstract?CRETRY=1&SRETRY=0>}.

      typical usage:

            sp=PlateMantelModel(domain,stress=0,v=0)
            sp.setTolerance()
            sp.initialize(...)
            v,p=sp.solve(v0,p0)
      """
      def __init__(self, domain, stress=0, v=0, p=0, t=0, numMaterials=1, useJaumannStress=True, **kwargs):
         """
         initializes PlateMantelModel 

         @param domain: problem domain
         @type domain: L{domain}
         @param stress: initial deviatoric stress
         @param v: initial velocity field
         @param p: initial pressure
         @param t: initial time
         @param useJaumannStress: C{True} if Jaumann stress is used (not supported yet)
         """
         if numMaterials<1:
            raise ValueError,"at least one material must be defined."
         super(IncompressibleIsotropicKelvinFlow, self).__init__(**kwargs)
         self.__domain=domain
         self.__t=t
         self.__vol=util.integrate(1.,Function(self.__domain))
         self.__useJaumannStress=useJaumannStress
         self.__numMaterials=numMaterials
         self.__eta_N=[None for i in xrange(self.__numMaterials)]
         self.__tau_t=[None for i in xrange(self.__numMaterials)]
         self.__power=[1 for i in xrange(self.__numMaterials)]
         self.__tau_Y=None
         self.__friction=None
         self.__mu=None
         self.__v_boundary=Vector(0,Solution(self.__domain))
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
         #=======================
         # PDE related stuff 
         self.__pde_v=LinearPDE(domain,numEquations=self.getDomain().getDim(),numSolutions=self.getDomain().getDim())
         self.__pde_v.setSymmetryOn()
         self.__pde_v.setSolverMethod(preconditioner=LinearPDE.RILU)
            
         self.__pde_p=LinearPDE(domain)
         self.__pde_p.setReducedOrderOn()
         self.__pde_p.setSymmetryOn()

         self.setTolerance()
         self.setSmall()
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

      def setTolerance(self,tol=1.e-4):
          """
          sets the tolerance
          """
          self.__pde_v.setTolerance(tol**2)
          self.__pde_p.setTolerance(tol**2)
          self.__tol=tol

      def getTolerance(self):
          """
          returns the set tolerance
          @rtype: positive float
          """
          return self.__tol

      def getDomain(self):
          """
          returns the domain
          """
          return self.__domain
      def getTime(self):
          """
          returns current time
          """
          return self.__t

      def setPowerLaw(self,id,eta_N, tau_t=None, power=1):
          """
          sets the power-law parameters for material q
          """
          if id<0 or id>=self.__numMaterials:
              raise ValueError,"Illegal material id = %s."%id
          self.__eta_N[id]=eta_N
          self.__power[id]=power
          self.__tau_t[id]=tau_t

      def setPowerLaws(self,eta_N, tau_t, power):
          """
          set the parameters of the powerlaw for all materials
          """
          if len(eta_N)!=self.__numMaterials or len(tau_t)!=self.__numMaterials or len(power)!=self.__numMaterials:
              raise ValueError,"%s materials are expected."%self.__numMaterials
          for i in xrange(self.__numMaterials):
               self.setPowerLaw(i,eta_N[i],tau_t[i],power[i])
      def setDruckerPragerLaw(self,tau_Y=None,friction=0):
          """
          set the parameters for the Drucker-Prager model
          """
          self.__tau_Y=tau_Y
          self.__friction=friction
          
      def setElasticShearModulus(self,mu=None):
          """ 
          sets the elastic shere modulus
          """
          self.__mu=mu
      def setExternals(self, F=None, f=None, q=None, v_boundary=None):
          """ 
          sets 

          @param F: external force.
          @param f: surface force
          @param q: location of constraints
          """
          if F != None: self.__pde_v.setValue(Y=F)
          if f != None: self.__pde_v.setValue(y=f)
          if q != None: self.__pde_v.setValue(q=q)
          if v_boundary != None: self.__v_boundary=v_boundary

      def bilinearform(self,arg0,arg1):
        s0=util.deviatoric(util.symmetric(util.grad(arg0[0])))
        s1=util.deviatoric(util.symmetric(util.grad(arg1[0])))
        # s0=util.interpolate(arg0[0],Function(self.getDomain()))
        # s1=util.interpolate(arg1[0],Function(self.getDomain()))
        p0=util.interpolate(arg0[1],Function(self.getDomain()))
        p1=util.interpolate(arg1[1],Function(self.getDomain()))
        a=util.integrate(self.__p_weight**2*util.inner(s0,s1))+util.integrate(p0*p1)
        return a

      def getEtaEff(self,strain, pressure):
          if self.__mu==None: 
                eps=util.length(strain)*util.sqrt(2)
          else:
                eps=util.length(strain+self.__stress/((2*self.__dt)*self.__mu))*util.sqrt(2)
          p=util.interpolate(pressure,eps.getFunctionSpace())
          if self.__tau_Y!= None:
             tmp=self.__tau_Y+self.__friction*p
             m=util.wherePositive(eps)*util.wherePositive(tmp)
             eta_max=m*tmp/(eps+(1-m)*util.EPSILON)+(1-m)*util.DBLE_MAX
          else:
             eta_max=util.DBLE_MAX
          # initial guess:
          tau=util.length(self.__stress)/util.sqrt(2)
          # start the iteration:
          cc=0
          TOL=1e-7
          dtau=util.DBLE_MAX
          print "tau = ", tau, "eps =",eps
          while cc<10 and dtau>TOL*util.Lsup(tau):
             eta_eff2,eta_eff_dash=self.evalEtaEff(tau,return_dash=True)
             eta_eff=util.clip(eta_eff2-eta_eff_dash*tau/(1-eta_eff_dash*eps),maxval=eta_max)
             tau, tau2=eta_eff*eps, tau
             dtau=util.Lsup(tau2-tau)
             print "step ",cc,dtau, util.Lsup(tau)
             cc+=1
          return eta_eff

      def getEtaCharacteristic(self):
          a=0
          for i in xrange(self.__numMaterials):
            a=a+1./self.__eta_N[i]
          return 1/a
             
      def evalEtaEff(self, tau, return_dash=False):
         a=Scalar(0,tau.getFunctionSpace())  # =1/eta
         if return_dash: a_dash=Scalar(0,tau.getFunctionSpace()) # =(1/eta)'
         s=util.Lsup(tau)
         if s>0:
            m=util.wherePositive(tau)
            tau2=s*util.EPSILON*(1.-m)+m*tau
            for i in xrange(self.__numMaterials):
                 eta_N=self.__eta_N[i]
                 tau_t=self.__tau_t[i]
                 if tau_t==None:
                    a+=1./eta_N
                 else:
                    power=1.-1./self.__power[i]
                    c=1./(tau_t**power*eta_N)
                    a+=c*tau2**power
                    if return_dash: a_dash+=power*c*tau2**(power-1.)
         else:
            for i in xrange(self.__numMaterials):
                 eta_N=self.__eta_N[i]
                 power=1.-1./self.__power[i]
                 a+=util.whereZero(power)/eta_N
         if self.__mu!=None: a+=1./(self.__dt*self.__mu)
         out=1/a
         if return_dash:
             return out,-out**2*a_dash
         else:
             return out
             
      def eval(self,arg):
         v=arg[0]
         p=arg[1]
         D=self.getDeviatoricStrain(v)
         eta_eff=self.getEtaEff(D,p)
         print "eta_eff=",eta_eff
         # solve for dv
         self.__pde_v.setValue(A=Data()) # save memory!
         k3=util.kronecker(Function(self.getDomain()))
         k3Xk3=util.outer(k3,k3)
         self.__pde_v.setValue(A=eta_eff*(util.swap_axes(k3Xk3,0,3)+util.swap_axes(k3Xk3,1,3)),X=-eta_eff*D+p*util.kronecker(self.getDomain()))
         dv=self.__pde_v.getSolution(verbose=self.__verbose)
         print "resistep dv =",dv
         # solve for dp
         v2=v+dv
         self.__pde_p.setValue(D=1/eta_eff,Y=util.div(v2))
         dp=self.__pde_p.getSolution(verbose=self.__verbose)
         print "resistep dp =",dp
         return ArithmeticTuple(dv,dp)

      def update(self,dt, iter_max=100, inner_iter_max=20, verbose=False):
          """
          updates stress, velocity and pressure for time increment dt

          @param dt: time increment
          @param max_inner_iter: maximum number of iteration steps in the incompressible solver
          @param verbose: prints some infos in the incompressible solve
          @param show_details: prints some infos while solving PDEs
          @param tol: tolerance for the time step
          """
          self.__verbose=verbose
          self.__dt=dt
          tol=self.getTolerance()
          # set the initial velocity:
          m=util.wherePositive(self.__pde_v.getCoefficient("q"))
          v_new=self.__v*(1-m)+self.__v_boundary*m
          # and off we go:
          x=ArithmeticTuple(v_new, self.__p)
          # self.__p_weight=util.interpolate(1./self.getEtaCharacteristic(),Function(self.__domain))**2
          self.__p_weight=self.getEtaCharacteristic()
          # self.__p_weight=util.interpolate(1./self.getEtaCharacteristic()**2,self.__p.getFunctionSpace())
          atol=self.norm(x)*self.__tol
          x_new=NewtonGMRES(self, x, iter_max=iter_max,sub_iter_max=inner_iter_max, atol=atol,rtol=0., verbose=verbose)
          self.__v=x_new[0]
          self.__p=x_new[1]
          1/0
          # self.__stress=self.getUpdatedStress(...)
          self.__t+=dt
          return self.__v, self.__p

      #=========================================================================================

      def getNewDeviatoricStress(self,D,eta_eff=None):
         if eta_eff==None: eta_eff=self.evalEtaEff(self.__stress,D,self.__p)
         s=(2*eta_eff)*D
         if self.__mu!=None: s+=eta_eff/(self.__dt*self.__mu)*self.__last_stress
         return s

      def getDeviatoricStress(self):
          """
          returns current stress
          """
          return self.__stress
      def getDeviatoricStrain(self,velocity=None):
          """
          return strain
          """
          if velocity==None: velocity=self.getVelocity()
          return util.deviatoric(util.symmetric(util.grad(velocity)))

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

      def getTau(self,stress=None):
          """
          returns current second stress deviatoric invariant
          """
          if stress==None: stress=self.getDeviatoricStress()
          return util.sqrt(0.5)*util.length(stress)

      def getGammaDot(self,strain=None):
          """
          returns current second stress deviatoric invariant
          """
          if strain==None: strain=self.getDeviatoricStrain()
          return util.sqrt(2)*util.length(strain)

