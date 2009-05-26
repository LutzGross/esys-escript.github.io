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
from flows import StokesProblemCartesian
from pdetools import MaxIterReached

class PowerLaw(object):
    """
    this implements the power law for a composition of a set of materials where the viscosity eta of each material is given by a 
    power law relationship of the form

    M{eta=eta_N*(tau/tau_t)**(1./power-1.)}

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
    def getEtaEff(self,gamma_dot, eta0=None, pressure=None,dt=None, iter_max=30):
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
         s=[p-1. for p in self.getPower() ]
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
         if self.__verbose: print "PowerLaw: Start calculation of eta_eff (tolerance = %s)\nPowerLaw: initial max eta_eff = %s, tau = %s."%(rtol,util.Lsup(eta_eff),util.Lsup(tau))
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
             if self.__verbose: print "PowerLaw: step %s: correction = %s, max eta_eff = %s, max tau= %s"%(iter, d, l,util.Lsup(tau))
             converged= d<= rtol* l
         if self.__verbose: print "PowerLaw: Start calculation of eta_eff finalized after %s steps."%iter
         return eta_eff

#====================================================================================================================================
class Rheology(object):
      """
      General framework to implement a rheology
      """
      def __init__(self, domain, stress=None, v=None, p=None, t=0, verbose=True):
         """
         Initializes the rheology

         @param domain: problem domain
         @type domain: L{Domain}
         @param stress: initial (deviatoric) stress
         @type stress: a tensor value/field of order 2
         @param v: initial velocity field
         @type stress: a vector value/field
         @param p: initial pressure
         @type p: a scalar value/field
         @param t: initial time
         @type t: C{float}
         """
         self.__domain=domain
         self.__t=t
         self.__verbose=verbose
         #=======================
         #
         # state variables:
         #
         if stress == None: stress=Tensor(0.,Function(self.__domain))
         if v == None: v=Vector(0.,Solution(self.__domain))
         if p == None: p=Vector(0.,ReducedSolution(self.__domain))
         self.setDeviatoricStress(stress)
         self.setVelocity(v)
         self.setPressure(p)
         self.setDeviatoricStrain()
         self.setTime(t)
         #=============================================================
         self.setExternals(F=Data(), f=Data(), fixed_v_mask=Data(), v_boundary=Data())
         
      def getDomain(self):
          """
          returns the domain.

          @return: the domain
          @rtype: L{Domain}
          """
          return self.__domain

      def getTime(self):
          """
          Returns current time.

          @return: current time
          @rtype: C{float}
          """
          return self.__t   

      def setExternals(self, F=None, f=None, fixed_v_mask=None, v_boundary=None):
          """
          sets external forces and velocity constraints

          @param F: external force
          @type F: vector value/field 
          @param f: surface force
          @type f: vector value/field on boundary
          @param fixed_v_mask: location of constraints maked by positive values
          @type fixed_v_mask: vector value/field 
          @param v_boundary: value of velocity at location of constraints
          @type v_boundary: vector value/field 
          @note: Only changing parameters need to be specified.
          """
          if F != None: self.__F=F
          if f != None: self.__f=f
          if fixed_v_mask != None: self.__fixed_v_mask=fixed_v_mask
          if v_boundary != None: self.__v_boundary=v_boundary 
          
      def getForce(self):
          """
          Returns the external force

          @return:  external force
          @rtype: L{Data}
          """
          return self.__F

      def getSurfaceForce(self):
          """
          Returns the surface force

          @return:  surface force
          @rtype: L{Data}
          """
          return self.__f

      def getVelocityConstraint(self):
          """
          Returns the constraint for the velocity as a pair of the 
          mask of the location of the constraint and the values.

          @return: the locations of fixed velocity and value of velocities at these locations
          @rtype: C{tuple} of L{Data}s
          """
          return self.__fixed_v_mask, self.__v_boundary       

      def checkVerbose(self):
          """
          Returns True if verbose is switched on

          @return: value of verbosity flag
          @rtype: C{bool}
          """
          return self.__verbose

      def setTime(self,t=0.):
          """
          Updates current time.

          @param t: new time mark
          @type t: C{float}
          """
          self.__t=t
      #=======================================================================================
      def getStress(self):
          """
          Returns current stress. 

          @return: current stress
          @rtype: L{Data} of rank 2
          """
          s=self.getDeviatoricStress()
          p=self.getPressure()
          k=util.kronecker(self.getDomain())
          return s-p*(k/trace(k))
            
      def getDeviatoricStress(self):
          """
          Returns current deviatoric stress.

          @return: current deviatoric stress
          @rtype: L{Data} of rank 2
          """
          return self.__stress

      def setDeviatoricStress(self, stress):
          """
          Sets the current deviatoric stress

          @param stress: new deviatoric stress
          @type stress: L{Data} of rank 2
          """
          dom=self.getDomain()
          s=util.interpolate(stress,Function(dom))
          self.__stress=util.deviatoric(s)

      def getPressure(self):
          """
          Returns current pressure.

          @return: current stress
          @rtype: scalar L{Data} 
          """
          return self.__p

      def setPressure(self, p):
          """
          Sets current pressure.
          @param p: new deviatoric stress
          @type p: scalar L{Data}
          """
          self.__p=util.interpolate(p,ReducedSolution(self.getDomain()))

      def getVelocity(self):
          """
          Returns current velocity.

          @return: current velocity
          @rtype: vector L{Data} 
          """
          return self.__v

      def setVelocity(self, v):
          """
          Sets current velocity.

          @param v: new current velocity
          @type v: vector L{Data} 
          """
          self.__v=util.interpolate(v,Solution(self.getDomain()))

      def setDeviatoricStrain(self, D=None):
          """
          set deviatoric strain 

          @param D: new deviatoric strain. If D is not present the current velocity is used.
          @type D: L{Data} of rank 2
          """
          if D==None: D=util.deviatoric(util.symmetric(util.grad(2.*self.getVelocity())))
          self.__D=util.deviatoric(util.interpolate(D,Function(self.getDomain())))

      def getDeviatoricStrain(self):
          """
          Returns deviatoric strain of current velocity. 

          @return: deviatoric strain
          @rtype: L{Data}  of rank 2
          """
          return self.__D

      def getTau(self):
          """
          Returns current second invariant of deviatoric stress

          @return: second invariant of deviatoric stress
          @rtype: scalar L{Data}
          """
          s=self.getDeviatoricStress()
          return util.sqrt(0.5)*util.length(s)

      def getGammaDot(self):
          """
          Returns current second invariant of deviatoric strain

          @return: second invariant of deviatoric strain
          @rtype: scalar L{Data}
          """
          s=self.getDeviatoricStrain()
          return util.sqrt(2)*util.length(s)

      def setTolerance(self,tol=1.e-4):
          """
          Sets the tolerance used to terminate the iteration on a time step.
          See the implementation of the rheology for details.

          @param tol: relative tolerance to terminate iteration on time step.
          @type tol: positive C{float}
          """
          if tol<=0.:
              raise ValueError,"tolerance must be non-negative."
          self.__tol=tol

      def getTolerance(self):
          """
          Returns the set tolerance for terminate the iteration on a time step.

          @rtype: positive C{float}
          """
          return self.__tol

      #=======================================================================
      def setFlowTolerance(self, tol=1.e-4):
          """
          Sets the relative tolerance for the flow solver

          @param tol: desired relative tolerance for the flow solver
          @type tol: positive C{float}
          @note: Typically this method is overwritten by a subclass.
          """
          pass
      def getFlowTolerance(self):
          """
          Returns the relative tolerance for the flow solver

          @return: tolerance of the flow solver
          @rtype: C{float}
          @note: Typically this method is overwritten by a subclass.
          """
          pass
      def setFlowSubTolerance(self, tol=1.e-8):
          """
          Sets the relative tolerance for the subsolver of the flow solver

          @param tol: desired relative tolerance for the subsolver
          @type tol: positive C{float}
          @note: Typically this method is overwritten by a subclass.
          """
          pass
      def getFlowSubTolerance(self):
          """
          Returns the relative tolerance for the subsolver of the flow solver

          @return: tolerance of the flow subsolver
          @rtype: C{float}
          @note: Typically this method is overwritten by a subclass.
          """
          pass


#====================================================================================================================================

class IncompressibleIsotropicFlowCartesian(PowerLaw,Rheology):
      """
      This class implements the rheology of an isotropic Kelvin material.

      Typical usage::

          sp = IncompressibleIsotropicFlow(domain, stress=0, v=0)
          sp.setTolerance()
          sp.initialize(...)
          v,p = sp.solve()

      @note: This model has been used in the self-consistent plate-mantle model
             proposed in U{Hans-Bernd Muhlhaus<emailto:h.muhlhaus@uq.edu.au>}
             and U{Klaus Regenauer-Lieb<mailto:klaus.regenauer-lieb@csiro.au>}:
             I{Towards a self-consistent plate mantle model that includes elasticity: simple benchmarks and application to basic modes of convection},
             see U{doi: 10.1111/j.1365-246X.2005.02742.x<http://www3.interscience.wiley.com/journal/118661486/abstract>}
      """
      def __init__(self, domain, stress=0, v=0, p=0, t=0, numMaterials=1, verbose=True):
         """
         Initializes the model.

         @param domain: problem domain
         @type domain: L{Domain}
         @param stress: initial (deviatoric) stress
         @type stress: a tensor value/field of order 2
         @param v: initial velocity field
         @type stress: a vector value/field
         @param p: initial pressure
         @type p: a scalar value/field
         @param t: initial time
         @type t: C{float}
         @param numMaterials: number of materials
         @type numMaterials: C{int}
         @param verbose: if C{True} some informations are printed.
         @type verbose: C{bool}         
         """
         PowerLaw. __init__(self, numMaterials,verbose)
         Rheology. __init__(self, domain, stress, v, p, t, verbose)
         self.__solver=StokesProblemCartesian(self.getDomain(),verbose=verbose)
         self.__eta_eff=None
         self.setTolerance()
         self.setFlowTolerance()
         self.setFlowSubTolerance()

      def update(self, dt, iter_max=100, inner_iter_max=20, verbose=False, usePCG=True):
          """
          Updates stress, velocity and pressure for time increment dt.

          @param dt: time increment
          @param inner_iter_max: maximum number of iteration steps in the
                                 incompressible solver
          @param verbose: prints some infos in the incompressible solver
          """
          if self.checkVerbose(): print "IncompressibleIsotropicFlowCartesian: start iteration for t = %s."%(self.getTime()+dt,)
          v_last=self.getVelocity()
          s_last=self.getDeviatoricStress()
          F=self.getForce()
          f=self.getSurfaceForce()
          mask_v,v_b=self.getVelocityConstraint()
          mu=self.getElasticShearModulus()
          #=========================================================================
          #
          #   we use velocity and pressure from the last time step as initial guess:
          #
          v=v_last
          p=self.getPressure()
          #
          #  calculate eta_eff  if we don't have one or elasticity is present.
          #
          if self.__eta_eff == None or  mu!=None: 
             D=self.__getDeviatoricStrain(v)
             if mu==None:
                 gamma=util.sqrt(2.)*util.length(D)
             else:
                 gamma=util.sqrt(2.)*util.length(D+s_last/(2*dt*mu))
             if self.__eta_eff == None:
                 eta0=None
             else:
                  eta0=self.__eta_eff
             eta_eff=self.getEtaEff(gamma, pressure=p,dt=dt, eta0=eta0, iter_max=iter_max)
             if self.checkVerbose(): print "IncompressibleIsotropicFlowCartesian: eta_eff has been initialied."
          else:
             eta_eff = self.__eta_eff
          iter=0 
          converged=False
          while not converged:
             #
             #   intialize the solver 
             #
             if mu==None:          
                stress0=Data()
             else:
                stress0=-(eta_eff/(dt*mu))*s_last
             self.__solver.initialize(f=F,fixed_u_mask=mask_v,eta=eta_eff,surface_stress=f,stress=stress0)
             # 
             # get a new velcocity and pressure:
             #
             if mask_v.isEmpty() or v_b.isEmpty():
                v0=v
             else:
                v0=v_b*mask_v+v*(1.-mask_v)
             v,p=self.__solver.solve(v0,p,show_details=False, 
                                          verbose=self.checkVerbose(),max_iter=inner_iter_max,usePCG=usePCG)
             # 
             #   update eta_eff:
             #
             D=self.__getDeviatoricStrain(v)
             if mu==None:
                 gamma=util.sqrt(2.)*util.length(D)
             else:
                 gamma=util.sqrt(2.)*util.length(D+s_last/(2*dt*mu))
             eta_eff_old ,eta_eff=eta_eff, self.getEtaEff(gamma, pressure=p,dt=dt, eta0=eta_eff, iter_max=iter_max)
             if self.checkVerbose(): print "IncompressibleIsotropicFlowCartesian: eta_eff has been updated."
             #
             # check the change on eta_eff:
             #
             diff=util.Lsup(eta_eff_old-eta_eff)
             n=util.Lsup(eta_eff)
             if self.checkVerbose(): print "IncompressibleIsotropicFlowCartesian: step %s: max. change in eta_eff is %s."%(iter,diff)
             converged = diff <= self.getTolerance()* n
             iter+=1
             if iter >= iter_max:
                 raise MaxIterReached,"maximum number of iteration steps on time step %e reached."%(self.getTime()+dt)
          #
          #   finally we can update the return values:
          #
          self.setPressure(p)
          self.setVelocity(v)
          self.setDeviatoricStrain(D)
          if mu==None:          
              stress=(2*eta_eff)*D
          else:
              stress=(2.*eta_eff)*(D+s_last/(2*dt*mu))
          self.setDeviatoricStress(stress)
          self.__eta_eff = eta_eff 
          self.setTime(self.getTime()+dt)
          if self.checkVerbose(): print "IncompressibleIsotropicFlowCartesian: iteration on time step %s completed after %s steps."%(self.getTime(),iter)
          return self.getVelocity(), self.getPressure()

      def __getDeviatoricStrain(self, v):
          """
          Returns deviatoric strain of velocity v:
          """
          return util.deviatoric(util.symmetric(util.grad(v)))

      def setFlowTolerance(self, tol=1.e-6):
          """
          Sets the relative tolerance for the flow solver. See L{StokesProblemCartesian.setTolerance} for details.

          @param tol: desired relative tolerance for the flow solver
          @type tol: positive C{float}
          """
          self.__solver.setTolerance(tol)
      def getFlowTolerance(self):
          """
          Returns the relative tolerance for the flow solver

          @return: tolerance of the flow solver
          @rtype: C{float}
          """
          return self.__solver.getTolerance()
      def setFlowSubTolerance(self, tol=1.e-12):
          """
          Sets the relative tolerance for the subsolver of the flow solver. See L{StokesProblemCartesian.setSubProblemTolerance} for details

          @param tol: desired relative tolerance for the subsolver
          @type tol: positive C{float}
          """
          self.__solver.setSubProblemTolerance(tol)
      def getFlowSubTolerance(self):
          """
          Returns the relative tolerance for the subsolver of the flow solver

          @return: tolerance of the flow subsolver
          @rtype: C{float}
          """
          return self.__solver.getSubProblemTolerance()

