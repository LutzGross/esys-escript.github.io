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
Some models for flow

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from . import escriptcpp as escore
from . import util
from .flows import StokesProblemCartesian
from .pdetools import MaxIterReached

class PowerLaw(object):
    """
    this implements the power law for a composition of a set of materials where the viscosity eta of each material is given by a 
    power law relationship of the form

    *eta=eta_N*(tau/tau_t)**(1./power-1.)*

    where tau is equivalent stress and eta_N, tau_t and power are given constant. Moreover an elastic component can be considered. 
    Moreover tau meets the Drucker-Prager type yield condition
 
    *tau <= tau_Y + friction * pressure*

    where gamma_dot is the equivalent. 
    """
    def __init__(self, numMaterials=1,verbose=False):
         """
         initializes a power law
   
         :param numMaterials: number of materials
         :type numMaterials: ``int``
         :param verbose: if ``True`` some information is printed.
         :type verbose: ``bool``
         """
         if numMaterials<1:
            raise ValueError("at least one material must be defined.")
         self.__numMaterials=numMaterials
         self.__eta_N=[None for i in range(self.__numMaterials)]
         self.__tau_t=[1. for i in range(self.__numMaterials)]
         self.__power=[1. for i in range(self.__numMaterials)]
         self.__tau_Y=None
         self.__friction=None
         self.__mu=None
         self.__verbose=verbose
         self.setEtaTolerance()
    #===========================================================================
    def getNumMaterials(self):
         """
         returns the numebr of materials
         
         :return: number of materials
         :rtype: ``int``
         """
         return self.__numMaterials
    def validMaterialId(self,id=0):
         """
         checks if a given material id is valid
   
         :param id: a material id
         :type id: ``int``
         :return: ``True`` is the id is valid
         :rtype: ``bool``
         """
         return 0<=id and id<self.getNumMaterials()
    def setEtaTolerance(self,rtol=1.e-4):
         """
         sets the relative tolerance for the effectice viscosity.
 
         :param rtol: relative tolerance
         :type rtol: positive ``float``
         """
         if rtol<=0:
             raise ValueError("rtol needs to positive.")
         self.__rtol=rtol
    def getEtaTolerance(self):
         """
         returns the relative tolerance for the effectice viscosity.
 
         :return: relative tolerance
         :rtype: positive ``float``
         """
         return self.__rtol
    #===========================================================================
    def setDruckerPragerLaw(self,tau_Y=None,friction=None):
          """
          Sets the parameters for the Drucker-Prager model.

          :param tau_Y: yield stress
          :param friction: friction coefficient
          """
          self.__tau_Y=tau_Y
          self.__friction=friction
    def getFriction(self):
         """
         returns the friction coefficient

         :return: friction coefficient
         """
         return self.__friction
    def getTauY(self):
         """
         returns the yield stress

         :return: the yield stress
         """
         return self.__tau_Y
    #===========================================================================
    def getElasticShearModulus(self):
        """
        returns the elastic shear modulus.

        :return: elastic shear modulus
        """
        return self.__mu
    def setElasticShearModulus(self,mu=None):
        """
        Sets the elastic shear modulus.

        :param mu: elastic shear modulus
        """
        self.__mu=mu
    #===========================================================================
    def getPower(self, id=None):
         """
         returns the power in the power law

         :param id:  if present, the power for material ``id`` is returned.
         :type id: ``int``
         :return: the list of the powers for all matrials is returned. If ``id`` is present only the power for material ``id`` is returned.
         """
         if id is None:
            return self.__power
         else:
            if self.validMaterialId(id):
              return self.__power[id]
            else:
              raise ValueError("Illegal material id %s."%id)
    def getEtaN(self, id=None):
         """
         returns the viscosity

         :param id:  if present, the viscosity for material ``id`` is returned.
         :type id: ``int``
         :return: the list of the viscosities for all matrials is returned. If ``id`` is present only the viscosity for material ``id`` is returned.
         """
         if id is None:
            return self.__eta_N
         else:
            if self.validMaterialId(id):
              return self.__eta_N[id]
            else:
             raise ValueError("Illegal material id %s."%id)
    def getTauT(self, id=None):
         """
         returns the transition stress

         :param id:  if present, the transition stress for material ``id`` is returned.
         :type id: ``int``
         :return: the list of the transition stresses for all matrials is returned. If ``id`` is present only the transition stress for material ``id`` is returned.
         """
         if id is None:
            return self.__tau_t
         else:
            if self.validMaterialId(id):
              return self.__tau_t[id]
            else:
              raise ValueError("Illegal material id %s."%id)
    
    def setPowerLaw(self,eta_N, id=0, tau_t=1, power=1):
          """
          Sets the power-law parameters for material id
          
          :param id: material id
          :type id: ``int``
          :param eta_N: viscosity for tau=tau_t
          :param tau_t: transition stress
          :param power: power law coefficient
          """
          if self.validMaterialId(id):
             self.__eta_N[id]=eta_N
             self.__power[id]=power
             self.__tau_t[id]=tau_t
          else:
              raise ValueError("Illegal material id %s."%id)

    def setPowerLaws(self,eta_N, tau_t, power):
          """
          Sets the parameters of the power-law for all materials.

          :param eta_N: list of viscosities for tau=tau_t
          :param tau_t: list of transition stresses
          :param power: list of power law coefficient
          """
          if len(eta_N)!=self.__numMaterials or len(tau_t)!=self.__numMaterials or len(power)!=self.__numMaterials:
              raise ValueError("%s materials are expected."%self.__numMaterials)
          for i in range(self.__numMaterials):
               self.setPowerLaw(id=i, eta_N=eta_N[i],tau_t=tau_t[i],power=power[i])

    #===========================================================================
    def getEtaEff(self,gamma_dot, eta0=None, pressure=None,dt=None, iter_max=30):
         """
         returns the effective viscosity eta_eff such that 

         *tau=eta_eff * gamma_dot*

         by solving a non-linear problem for tau.

         :param gamma_dot: equivalent strain gamma_dot
         :param eta0: initial guess for the effective viscosity (e.g from a previous time step). If not present, an initial guess is calculated.
         :param pressure: pressure used to calculate yield condition
         :param dt: time step size. only needed if elastic component is considered.
         :type dt: positive ``float`` if present
         :param iter_max: maximum number of iteration steps.
         :type iter_max: ``int``
         :return: effective viscosity. 
         """
         if pressure is None:
            p2 = None
         else:
            p2=(abs(pressure)+pressure)/2.
         SMALL=1./(util.DBLE_MAX/100.)
         numMaterial=self.getNumMaterials()
         s=[p-1. for p in self.getPower() ]
         eta_N=self.getEtaN()
         tau_t=self.getTauT()
         mu=self.getElasticShearModulus()
         fric=self.getFriction()
         tau_Y=self.getTauY()
         if eta0 is None:
             theta=0.
             for i in range(numMaterial): 
                  inv_eta_i=0**s[i]/eta_N[i]
                  theta=theta+inv_eta_i
             if util.inf(theta)<=0: 
                 raise ValueError("unable to set positive initial guess for eta_eff. Most likely no power law with power 1 set.")
             eta_eff=1./theta
         else:
             if util.inf(eta0)<=0:
                 raise ValueError("initial guess for eta_eff is not positive.")
             eta_eff=eta0

         if mu !=None:
             if dt is None: raise ValueError("Time stepsize dt must be given.")
             if dt<=0: raise ValueError("Time step size must be positive.")
         if tau_Y is None and fric is None:
             eta_max=None
         else:
            if fric is None or p2 is None:
                eta_max=tau_Y/(gamma_dot+SMALL*util.whereZero(gamma_dot))
            else:
                if tau_Y is None: tau_Y==0
                if util.inf(fric)<=0: 
                    raise ValueError("if friction present it needs to be positive.")
                eta_max=fric*util.clip(tau_Y/fric+p2,minval=0)/(gamma_dot+SMALL*util.whereZero(gamma_dot))
         rtol=self.getEtaTolerance()
         iter =0
         converged=False
         tau=eta_eff*gamma_dot
         if self.__verbose: print(("PowerLaw: Start calculation of eta_eff (tolerance = %s)\nPowerLaw: initial max eta_eff = %s, tau = %s."%(rtol,util.Lsup(eta_eff),util.Lsup(tau))))
         while not converged:
             if iter>max(iter_max,1):
                raise RuntimeError("tolerance not reached after %s steps."%max(iter_max,1))
             #===========================================
             theta=0. # =1/eta
             omega=0. # = tau*theta'= eta'*tau/eta**2
             if mu !=None: theta=1./(dt*mu)
             for i in range(numMaterial):
                  inv_eta_i=(tau/tau_t[i])**s[i]/eta_N[i]
                  theta=theta+inv_eta_i
                  omega=omega+s[i]*inv_eta_i
             #===========================================
             eta_eff, eta_eff_old=util.clip(eta_eff*(theta+omega)/(eta_eff*theta**2+omega),maxval=eta_max), eta_eff
             tau=eta_eff*gamma_dot
             d=util.Lsup(eta_eff-eta_eff_old)
             l=util.Lsup(eta_eff)
             iter+=1
             if self.__verbose: print(("PowerLaw: step %s: correction = %s, max eta_eff = %s, max tau= %s"%(iter, d, l,util.Lsup(tau))))
             converged= d<= rtol* l
         if self.__verbose: print(("PowerLaw: Start calculation of eta_eff finalized after %s steps."%iter))
         return eta_eff

#====================================================================================================================================
class Rheology(object):
      """
      General framework to implement a rheology
      """
      def __init__(self, domain, stress=None, v=None, p=None, t=0, verbose=True):
         """
         Initializes the rheology

         :param domain: problem domain
         :type domain: `Domain`
         :param stress: initial (deviatoric) stress
         :type stress: a tensor value/field of order 2
         :param v: initial velocity field
         :type v: a vector value/field
         :param p: initial pressure
         :type p: a scalar value/field
         :param t: initial time
         :type t: ``float``
         """
         self.__domain=domain
         self.__t=t
         self.__verbose=verbose
         #=======================
         #
         # state variables:
         #
         if stress is None: stress=Tensor(0.,escore.Function(self.__domain))
         if v is None: v=Vector(0.,escore.Solution(self.__domain))
         if p is None: p=Vector(0.,escore.ReducedSolution(self.__domain))
         self.setStatus(t, v, p, stress)
         self.setExternals(F=escore.Data(), f=escore.Data(), fixed_v_mask=escore.Data(), v_boundary=escore.Data(), restoration_factor=0)
         
      def getDomain(self):
          """
          returns the domain.

          :return: the domain
          :rtype: `Domain`
          """
          return self.__domain

      def getTime(self):
          """
          Returns current time.

          :return: current time
          :rtype: ``float``
          """
          return self.__t   

      def setExternals(self, F=None, f=None, fixed_v_mask=None, v_boundary=None, restoration_factor=None):
          """
          sets external forces and velocity constraints

          :param F: external force
          :type F: vector value/field 
          :param f: surface force
          :type f: vector value/field on boundary
          :param fixed_v_mask: location of constraints maked by positive values
          :type fixed_v_mask: vector value/field 
          :param v_boundary: value of velocity at location of constraints
          :type v_boundary: vector value/field 
          :param restoration_factor: factor for normal restoration force
          :type restoration_factor: scalar values/field
          :note: Only changing parameters need to be specified.
          """
          if F is not None: self.__F=F
          if f is not None: self.__f=f
          if fixed_v_mask is not None: self.__fixed_v_mask=fixed_v_mask
          if v_boundary is not None: self.__v_boundary=v_boundary 
          if restoration_factor is not None: self.__restoration_factor=restoration_factor
          
      def getForce(self):
          """
          Returns the external force

          :return:  external force
          :rtype: `Data`
          """
          return self.__F

      def getSurfaceForce(self):
          """
          Returns the surface force

          :return:  surface force
          :rtype: `Data`
          """
          return self.__f

      def getVelocityConstraint(self):
          """
          Returns the constraint for the velocity as a pair of the 
          mask of the location of the constraint and the values.

          :return: the locations of fixed velocity and value of velocities at these locations
          :rtype: ``tuple`` of `Data` s
          """
          return self.__fixed_v_mask, self.__v_boundary       

      def getRestorationFactor(self):
          """
          Returns the restoring force factor

          :return:  restoring force factor
          :rtype: `float` or `Data`
          """
          return self.__restoration_factor
          

      def checkVerbose(self):
          """
          Returns True if verbose is switched on

          :return: value of verbosity flag
          :rtype: ``bool``
          """
          return self.__verbose

      def setTime(self,t=0.):
          """
          Updates current time.

          :param t: new time mark
          :type t: ``float``
          """
          self.__t=t
      #=======================================================================================
      def getStress(self):
          """
          Returns current stress. 

          :return: current stress
          :rtype: `Data` of rank 2
          """
          s=self.getDeviatoricStress()
          p=self.getPressure()
          k=util.kronecker(self.getDomain())
          return s-p*(k/trace(k))
            
      def getDeviatoricStress(self):
          """
          Returns current deviatoric stress.

          :return: current deviatoric stress
          :rtype: `Data` of rank 2
          """
          return self.__stress

      def setDeviatoricStress(self, stress):
          """
          Sets the current deviatoric stress

          :param stress: new deviatoric stress
          :type stress: `Data` of rank 2
          """
          dom=self.getDomain()
          s=util.interpolate(stress,escore.Function(dom))
          self.__stress=util.deviatoric(s)

      def getPressure(self):
          """
          Returns current pressure.

          :return: current stress
          :rtype: scalar `Data` 
          """
          return self.__p

      def setPressure(self, p):
          """
          Sets current pressure.
          :param p: new deviatoric stress
          :type p: scalar `Data`
          """
          self.__p=util.interpolate(p,escore.ReducedSolution(self.getDomain()))

      def getVelocity(self):
          """
          Returns current velocity.

          :return: current velocity
          :rtype: vector `Data` 
          """
          return self.__v

      def setVelocity(self, v):
          """
          Sets current velocity.

          :param v: new current velocity
          :type v: vector `Data` 
          """
          self.__v=util.interpolate(v,escore.Solution(self.getDomain()))
      def setStatus(self,t, v, p, stress):
          """
          Resets the current status given by pressure p and velocity v.
    
          :param t: new time mark
          :type t: `float`
          :param v: new current velocity
          :type v: vector `Data`
          :param p: new deviatoric stress
          :type p: scalar `Data`
          :param stress: new deviatoric stress
          :type stress: `Data` of rank 2
          """
          self.setDeviatoricStress(stress)
          self.setVelocity(v)
          self.setPressure(p)
          self.setDeviatoricStrain()
          self.setGammaDot()
          self.setTime(t)

      def setDeviatoricStrain(self, D=None):
          """
          set deviatoric strain 

          :param D: new deviatoric strain. If ``D`` is not present the current velocity is used.
          :type D: `Data` of rank 2
          """
          if D is None: 
              self.__D=self.getDeviatoricStrain(self.getVelocity())
          else:
              self.__D=util.deviatoric(util.interpolate(D,escore.Function(self.getDomain())))

      def getDeviatoricStrain(self, v=None):
          """
          Returns deviatoric strain of current velocity or if ``v`` is present the 
          deviatoric strain of velocity ``v``:

          :param v: a velocity field
          :type v: `Data` of rank 1
          :return: deviatoric strain of the current velocity field or if ``v`` is present the deviatoric strain of velocity ``v``
          :rtype: `Data`  of rank 2
          """
          if v is None:
             return self.__D
          else:
             return util.deviatoric(util.symmetric(util.grad(v)))

      def getTau(self):
          """
          Returns current second invariant of deviatoric stress

          :return: second invariant of deviatoric stress
          :rtype: scalar `Data`
          """
          s=self.getDeviatoricStress()
          return util.sqrt(0.5)*util.length(s)

      def setGammaDot(self, gammadot=None):
          """
          set the second invariant of deviatoric strain rate. If ``gammadot`` is not present zero is used.

          :param gammadot: second invariant of deviatoric strain rate. 
          :type gammadot: `Data` of rank 1
          """
          if gammadot is None:
               self.__gammadot = escore.Scalar(0.,escore.Function(self.getDomain()))
          else:
               self.__gammadot=gammadot
          
      def getGammaDot(self, D=None):
          """
          Returns current second invariant of deviatoric strain rate or if ``D`` is present the second invariant of ``D``.

          :param D: deviatoric strain rate tensor
          :type D: `Data`  of rank 0
          :return: second invariant of deviatoric strain
          :rtype: scalar `Data`
          """
          if D is None: 
              return self.__gammadot
          else:
              return util.sqrt(2.)*util.length(D)
           

#====================================================================================================================================

class IncompressibleIsotropicFlowCartesian(PowerLaw,Rheology, StokesProblemCartesian):
     """
     This class implements the rheology of an isotropic Kelvin material.

     Typical usage::

          sp = IncompressibleIsotropicFlowCartesian(domain, stress=0, v=0)
          sp.initialize(...)
          v,p = sp.solve()

     :note: This model has been used in the self-consistent plate-mantle model
             proposed in `Hans-Bernd Muhlhaus <mailto:h.muhlhaus@uq.edu.au>`_
             and `Klaus Regenauer-Lieb <mailto:klaus.regenauer-lieb@csiro.au>`_:
             "Towards a self-consistent plate mantle model that includes elasticity: simple benchmarks and application to basic modes of convection",
             see `doi: 10.1111/j.1365-246X.2005.02742.x <http://www3.interscience.wiley.com/journal/118661486/abstract>`_
     """
     def __init__(self, domain, stress=0, v=0, p=0, t=0, numMaterials=1, verbose=True):
         """
         Initializes the model.

         :param domain: problem domain
         :type domain: `Domain`
         :param stress: initial (deviatoric) stress
         :type stress: a tensor value/field of order 2
         :param v: initial velocity field
         :type v: a vector value/field
         :param p: initial pressure
         :type p: a scalar value/field
         :param t: initial time
         :type t: ``float``
         :param numMaterials: number of materials
         :type numMaterials: ``int``
         :param verbose: if ``True`` some information is printed.
         :type verbose: ``bool``         
         """
         PowerLaw. __init__(self, numMaterials,verbose=verbose)
         Rheology. __init__(self, domain, stress, v, p, t,verbose=verbose)
         StokesProblemCartesian.__init__(self,domain,verbose=verbose)
         self.__eta_eff=None

     def getCurrentEtaEff(self):
          """
          returns the effective viscosity used in the last iteration step of the last time step.
          """
          return self.__eta_eff


     def updateStokesEquation(self, v, p):
         """
         updates the underlying Stokes equation to consider dependencies from ``v`` and ``p``
         """
         dt=self.__dt
         mu=self.getElasticShearModulus()
         F=self.getForce()
         f=self.getSurfaceForce()
         mask_v,v_b=self.getVelocityConstraint()
         s_last=self.getDeviatoricStress()
         #
         #  calculate eta_eff if we don't have one or elasticity is present.
         #
         if mu is None:
             gamma=self.getGammaDot(self.getDeviatoricStrain(v))
         else:
             gamma=self.getGammaDot(self.getDeviatoricStrain(v)+s_last/(2*dt*mu))

         self.__eta_eff_save=self.getEtaEff(gamma, pressure=p,dt=dt, eta0=self.__eta_eff_save, iter_max=self.__eta_iter_max)

         if self.checkVerbose(): print("IncompressibleIsotropicFlowCartesian: eta_eff has been updated.")

         if mu is None:          
             stress0=escore.Data()
         else:
             stress0=-(self.__eta_eff_save/(dt*mu))*s_last
 
         self.setStokesEquation(eta=self.__eta_eff_save,stress=stress0)


     def initialize(self, F=None, f=None, fixed_v_mask=None, v_boundary=None, restoration_factor=None):
          """
          sets external forces and velocity constraints

          :param F: external force
          :type F: vector value/field 
          :param f: surface force
          :type f: vector value/field on boundary
          :param fixed_v_mask: location of constraints maked by positive values
          :type fixed_v_mask: vector value/field 
          :param v_boundary: value of velocity at location of constraints
          :type v_boundary: vector value/field 
          :param restoration_factor: factor for normal restoration force
          :type restoration_factor: scalar values/field
          :note: Only changing parameters need to be specified.
          """
          self.setExternals(F, f, fixed_v_mask, v_boundary, restoration_factor)

     def update(self, dt, iter_max=10, eta_iter_max=20, verbose=False, usePCG=True, max_correction_steps=50):
          """
          Updates stress, velocity and pressure for time increment dt.

          :param dt: time increment
          :param iter_max: maximum number of iteration steps in the incompressible solver
          :param eta_iter_max: maximum number of iteration steps in the incompressible solver
          :param verbose: prints some infos in the incompressible solver
          """
          mu=self.getElasticShearModulus()
          if mu is not None:
             if not dt > 0.:
                 raise ValueError("dt must be positive.")
          else:
             dt=max(0,dt)
          self.__dt=dt
          self.__eta_iter_max=max(eta_iter_max,1)
          v_last=self.getVelocity() 
          s_last=self.getDeviatoricStress()
          mask_v,v_b=self.getVelocityConstraint()
          p_last=self.getPressure()
          self.__eta_eff_save=self.getCurrentEtaEff()

          self.setStokesEquation(f=self.getForce(),fixed_u_mask=mask_v,surface_stress=self.getSurfaceForce(), restoration_factor=self.getRestorationFactor())

          if self.checkVerbose(): print(("IncompressibleIsotropicFlowCartesian: start iteration for t = %s."%(self.getTime()+dt,)))
          # 
          # get a new velcocity and pressure:
          #
          if mask_v.isEmpty():
               v0=v_last
          else:
              if v_b.isEmpty():
                 v0=v_last*(1.-mask_v)
              else:
                 v0=v_b*mask_v+v_last*(1.-mask_v)

          v,p=self._solve(v0,p_last,verbose=self.checkVerbose(),max_iter=iter_max,usePCG=usePCG, max_correction_steps=max_correction_steps)
          #
          #   finally we can update the return values:
          #
          self.setPressure(p)
          self.setVelocity(v)
          self.setDeviatoricStrain(self.getDeviatoricStrain(v))
          if mu is None:
             D=self.getDeviatoricStrain(v)
          else:
             D=self.getDeviatoricStrain(v)+s_last/(2*dt*mu)
          gamma=self.getGammaDot(D)
          self.setGammaDot(gamma)
          self.__eta_eff = self.getEtaEff(self.getGammaDot(), pressure=p,dt=dt, eta0=self.__eta_eff_save, iter_max=self.__eta_iter_max)
          self.setDeviatoricStress(2.*self.__eta_eff*D)
          self.setTime(self.getTime()+dt)
          if self.checkVerbose(): print(("IncompressibleIsotropicFlowCartesian: iteration on time step %s completed."%(self.getTime(),)))
          return self.getVelocity(), self.getPressure()


