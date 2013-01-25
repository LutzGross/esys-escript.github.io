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
__url__="https://launchpad.net/escript-finley"

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

class PowerLaw(object):
    """
    this implements the power law for a composition of a set of materials where the viscosity eta of each material is given by a 
    power law relationship of the form

    M{eta=eta_N*(tau/tau_t)**(1.-1./power)}

    where tau is equivalent stress and eta_N, tau_t and power are given constant. Moreover an elastic component can be considered. 
    Moreover tau meets the Drucker-Prager type yield condition
 
    M{tau <= tau_Y + friction * pressure}

    where gamma_dot is the equivalent. 
    """
    def __init__(self, numMaterials=1,verbose=False):
         """
         initializes a power law
   
         @param numMaterials: number of materials
         @type numMaterials: C{int}
         @param verbose: if C{True} some informations are printed.
         @type verbose: C{bool}
         """
         if numMaterials<1:
            raise ValueError,"at least one material must be defined."
         self.__numMaterials=numMaterials
         self.__eta_N=[None for i in xrange(self.__numMaterials)]
         self.__tau_t=[1. for i in xrange(self.__numMaterials)]
         self.__power=[1. for i in xrange(self.__numMaterials)]
         self.__tau_Y=None
         self.__friction=None
         self.__mu=None
         self.__verbose=verbose
         self.setEtaTolerance()
    #===========================================================================
    def getNumMaterials(self):
         """
         returns the numebr of materials
         
         @return: number of materials
         @rtype: C{int}
         """
         return self.__numMaterials
    def validMaterialId(self,id=0):
         """
         checks if a given material id is valid
   
         @param id: a material id
         @type id: C{int}
         @return: C{True} is the id is valid
         @rtype: C{bool}
         """
         return 0<=id and id<self.getNumMaterials()
    def setEtaTolerance(self,rtol=util.sqrt(util.EPSILON)):
         """
         sets the relative tolerance for the effectice viscosity.
 
         @param rtol: relative tolerance
         @type rtol: positive C{float}
         """
         if rtol<=0:
             raise ValueError,"rtol needs to positive."
         self.__rtol=rtol
    def getEtaTolerance(self):
         """
         returns the relative tolerance for the effectice viscosity.
 
         @return: relative tolerance
         @rtype rtol: positive C{float}
         """
         return self.__rtol
    #===========================================================================
    def setDruckerPragerLaw(self,tau_Y=None,friction=None):
          """
          Sets the parameters for the Drucker-Prager model.

          @param tau_Y: yield stress
          @param friction: friction coefficient
          """
          self.__tau_Y=tau_Y
          self.__friction=friction
    def getFriction(self):
         """
         returns the friction coefficient

         @return: friction coefficient
         """
         return self.__friction
    def getTauY(self):
         """
         returns the yield stress

         @return: the yield stress
         """
         return self.__tau_Y
    #===========================================================================
    def getElasticShearModulus(self):
        """
        returns the elastic shear modulus.

        @return: elastic shear modulus
        """
        return self.__mu
    def setElasticShearModulus(self,mu=None):
        """
        Sets the elastic shear modulus.

        @param mu: elastic shear modulus
        """
        self.__mu=mu
    #===========================================================================
    def getPower(self, id=None):
         """
         returns the power in the power law

         @param id:  if present, the power for material C{id} is returned.
         @type id: C{int}
         @return: the list of the powers for all matrials is returned. If C{id} is present only the power for material C{id} is returned.
         """
         if id == None:
            return self.__power
         else:
            if self.validMaterialId(id):
              return self.__power[id]
            else:
              raise ValueError,"Illegal material id %s."%id
    def getEtaN(self, id=None):
         """
         returns the viscosity

         @param id:  if present, the viscosity for material C{id} is returned.
         @type id: C{int}
         @return: the list of the viscosities for all matrials is returned. If C{id} is present only the viscosity for material C{id} is returned.
         """
         if id == None:
            return self.__eta_N
         else:
            if self.validMaterialId(id):
              return self.__eta_N[id]
            else:
             raise ValueError,"Illegal material id %s."%id
    def getTauT(self, id=None):
         """
         returns the transition stress

         @param id:  if present, the transition stress for material C{id} is returned.
         @type id: C{int}
         @return: the list of the transition stresses for all matrials is returned. If C{id} is present only the transition stress for material C{id} is returned.
         """
         if id == None:
            return self.__tau_t
         else:
            if self.validMaterialId(id):
              return self.__tau_t[id]
            else:
              raise ValueError,"Illegal material id %s."%id
    
    def setPowerLaw(self,eta_N, id=0, tau_t=1, power=1):
          """
          Sets the power-law parameters for material id
          
          @param id: material id
          @type id: C{int}
          @param eta_N: viscosity for tau=tau_t
          @param tau_t: transition stress
          @param power: power law coefficient
          """
          if self.validMaterialId(id):
             self.__eta_N[id]=eta_N
             self.__power[id]=power
             self.__tau_t[id]=tau_t
          else:
              raise ValueError,"Illegal material id %s."%id

    def setPowerLaws(self,eta_N, tau_t, power):
          """
          Sets the parameters of the power-law for all materials.

          @param eta_N: list of viscosities for tau=tau_t
          @param tau_t: list of transition stresses
          @param power: list of power law coefficient
          """
          if len(eta_N)!=self.__numMaterials or len(tau_t)!=self.__numMaterials or len(power)!=self.__numMaterials:
              raise ValueError,"%s materials are expected."%self.__numMaterials
          for i in xrange(self.__numMaterials):
               self.setPowerLaw(id=i, eta_N=eta_N[i],tau_t=tau_t[i],power=power[i])

    #===========================================================================
    def getEtaEff(self,gamma_dot, eta0=None, pressure=None,dt=None, iter_max=10):
         """
         returns the effective viscosity eta_eff such that 

         M{tau=eta_eff * gamma_dot}

         by solving a non-linear problem for tau.

         @param gamma_dot: equivalent strain gamma_dot
         @param eta0: initial guess for the effective viscosity (e.g from a previous time step). If not present, an initial guess is calculated.
         @param pressure: pressure used to calculate yield condition
         @param dt: time step size. only needed if elastic component is considered.
         @type dt: positive C{float} if present
         @param iter_max: maximum number of iteration steps.
         @type iter_max: C{int}
         @return: effective viscosity. 
         """
         SMALL=1./(util.DBLE_MAX/100.)
         numMaterial=self.getNumMaterials()
         s=[1.-1./p for p in self.getPower() ]
         eta_N=self.getEtaN()
         tau_t=self.getTauT()
         mu=self.getElasticShearModulus()
         fric=self.getFriction()
         tau_Y=self.getTauY()
         if eta0==None:
             theta=0.
             for i in xrange(numMaterial): 
                  inv_eta_i=0**s[i]/eta_N[i]
                  theta=theta+inv_eta_i
             if util.inf(theta)<=0: 
                 raise ValueError,"unable to set positive initial guess for eta_eff. Most likely no power law with power 1 set."
             eta_eff=1./theta
         else:
             if util.inf(eta0)<=0:
                 raise ValueError,"initial guess for eta_eff is not positive."
             eta_eff=eta0

         if mu !=None and dt == None:
             raise ValueError,"Time stepsize dt must be given."
         if dt !=None:
             if dt<=0: raise ValueError,"time step size must be positive."
         if tau_Y==None and fric==None:
             eta_max=None
         else:
            if fric == None:
                eta_max=tau_Y/(gamma_dot+SMALL*util.whereZero(gamma_dot))
            else:
                if tau_Y==None: tau_Y==0
                if util.inf(fric)<=0: 
                    raise ValueError,"if friction present it needs to be positive."
                eta_max=fric*util.clip(tau_Y/fric+pressure,minval=0)/(gamma_dot+SMALL*util.whereZero(gamma_dot))
         rtol=self.getEtaTolerance()
         iter =0
         converged=False
         tau=eta_eff*gamma_dot
         if self.__verbose: print "Start calculation of eta_eff (tolerance = %s)\ninitial max eta_eff = %s, tau = %s."%(rtol,util.Lsup(eta_eff),util.Lsup(tau))
         while not converged:
             if iter>max(iter_max,1):
                raise RuntimeError,"tolerance not reached after %s steps."%max(iter_max,1)
             #===========================================
             theta=0. # =1/eta
             omega=0. # = tau*theta'= eta'*tau/eta**2
             if mu !=None: theta=1./(dt*mu)
             for i in xrange(numMaterial):
                  inv_eta_i=(tau/tau_t[i])**s[i]/eta_N[i]
                  theta=theta+inv_eta_i
                  omega=omega+s[i]*inv_eta_i
             #===========================================
             eta_eff, eta_eff_old=util.clip(eta_eff*(theta+omega)/(eta_eff*theta**2+omega),maxval=eta_max), eta_eff
             tau=eta_eff*gamma_dot
             d=util.Lsup(eta_eff-eta_eff_old)
             l=util.Lsup(eta_eff)
             iter+=1
             if self.__verbose: print "step %s: correction = %s, max eta_eff = %s, max tau= %s"%(iter, d, l,util.Lsup(tau))
             converged= d<= rtol* l
         return eta_eff

#====================================================================================================================================

class IncompressibleIsotropicKelvinFlow(Defect):
      """
      This class implements the rheology of an isotropic Kelvin material.

      Typical usage::

          sp = IncompressibleIsotropicKelvinFlow(domain, stress=0, v=0)
          sp.setTolerance()
          sp.initialize(...)
          v,p = sp.solve(v0, p0)

      @note: This model has been used in the self-consistent plate-mantle model
             proposed in U{Hans-Bernd Muhlhaus<emailto:h.muhlhaus@uq.edu.au>}
             and U{Klaus Regenauer-Lieb<mailto:klaus.regenauer-lieb@csiro.au>}:
             I{Towards a self-consistent plate mantle model that includes elasticity: simple benchmarks and application to basic modes of convection},
             see U{doi: 10.1111/j.1365-246X.2005.02742.x<http://www3.interscience.wiley.com/journal/118661486/abstract?CRETRY=1&SRETRY=0>}.

      """
      def __init__(self, domain, stress=0, v=0, p=0, t=0, numMaterials=1, useJaumannStress=True, **kwargs):
         """
         Initializes the model.

         @param domain: problem domain
         @type domain: L{domain}
         @param stress: initial deviatoric stress
         @param v: initial velocity field
         @param p: initial pressure
         @param t: initial time
         @param useJaumannStress: C{True} if Jaumann stress is used
                                  (not supported yet)
         """
         super(IncompressibleIsotropicKelvinFlow, self).__init__(**kwargs)
         self.__domain=domain
         self.__t=t
         self.__vol=util.integrate(1.,Function(self.__domain))
         self.__useJaumannStress=useJaumannStress
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
          Returns C{True} if Jaumann stress is included.
          """
          return self.__useJaumannStress

      def setSmall(self,small=util.sqrt(util.EPSILON)):
          """
          Sets a small value to be used.

          @param small: positive small value
          """
          self.__small=small

      def getSmall(self):
          """
          Returns small value.
          @rtype: positive float
          """
          return self.__small

      def setTolerance(self,tol=1.e-4):
          """
          Sets the tolerance.
          """
          self.__pde_v.setTolerance(tol**2)
          self.__pde_p.setTolerance(tol**2)
          self.__tol=tol

      def getTolerance(self):
          """
          Returns the set tolerance.
          @rtype: positive float
          """
          return self.__tol

      def getDomain(self):
          """
          Returns the domain.
          """
          return self.__domain

      def getTime(self):
          """
          Returns current time.
          """
          return self.__t

      def setExternals(self, F=None, f=None, q=None, v_boundary=None):
          """
          Sets externals.

          @param F: external force
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

      def update(self, dt, iter_max=100, inner_iter_max=20, verbose=False):
          """
          Updates stress, velocity and pressure for time increment dt.

          @param dt: time increment
          @param inner_iter_max: maximum number of iteration steps in the
                                 incompressible solver
          @param verbose: prints some infos in the incompressible solver
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

      #========================================================================

      def getNewDeviatoricStress(self,D,eta_eff=None):
         if eta_eff==None: eta_eff=self.evalEtaEff(self.__stress,D,self.__p)
         s=(2*eta_eff)*D
         if self.__mu!=None: s+=eta_eff/(self.__dt*self.__mu)*self.__last_stress
         return s

      def getDeviatoricStress(self):
          """
          Returns current stress.
          """
          return self.__stress

      def getDeviatoricStrain(self,velocity=None):
          """
          Returns strain.
          """
          if velocity==None: velocity=self.getVelocity()
          return util.deviatoric(util.symmetric(util.grad(velocity)))

      def getPressure(self):
          """
          Returns current pressure.
          """
          return self.__p

      def getVelocity(self):
          """
          Returns current velocity.
          """
          return self.__v

      def getTau(self,stress=None):
          """
          Returns current second stress deviatoric invariant.
          """
          if stress==None: stress=self.getDeviatoricStress()
          return util.sqrt(0.5)*util.length(stress)

      def getGammaDot(self,strain=None):
          """
          Returns current second stress deviatoric invariant.
          """
          if strain==None: strain=self.getDeviatoricStrain()
          return util.sqrt(2)*util.length(strain)

