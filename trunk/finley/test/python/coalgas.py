######################################################
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
"""
Gas in Coal Seam (fully coupled version)
"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript.linearPDEs import LinearPDE
from esys.escript import unitsSI as U
from esys.escript.pdetools import Locator
from esys.escript import *
from esys.weipa import saveVTK
import math

USE_NODAL_WELL = False or True

class MaterialProperty(object):
   """
   generic class for material properties depending on one (or several) status 
   variable(s)
   """
   def __init__(self, *args, **kwargs):
       """
       set up
       """
       pass
   def __call__(self,  *args, **kwargs):
      """
      return value of material property for given status
      """
      return self.getValue( *args, **kwargs)
      
      
   def getValue(self, *args, **kwargs): 
      """
      return value of material property for given status
      
      :remark: needs to be overwritten
      """
      raise NotImplementedError

class MaterialPropertyWithDifferential(MaterialProperty):
   """
   generic class for material properties depending on one (or several) status 
   variable(s) where the  derivative of the property with respect to the 
   status variables is available
   """
   def getValueDifferential(self,  *args, **kwargs):
      """
      returns the value of the derivative of material property for given status variable
      
      :remark: needs to be overwritten
      """
      raise NotImplementedError

class Porosity(MaterialPropertyWithDifferential):
    """
    defines porosity phi as function of pressure 
    
         phi = phi_p_ref /(1 + X + X**2/2 ) with X= C * (p-p_ref) 
         
    phi_p_ref is claculted from the initial porosity phi_0 at  pressure p_0
    """
    def __init__(self, phi_0, p_0, p_ref=1.*U.atm , C=4.0e-5/U.bar):
      """
      """
      self.phi_p_ref=1 # will be overwritten later
      self.p_ref=p_ref
      self.C=C
      # reset  phi_p_ref to get phi(p_0)=phi_0    
      self.phi_p_ref=phi_0/self.getValue(p_0)
      
    def getValue(self, p):
      """
      returns the porosity for given pressure p
      """
      X= self.C * ( p - self.p_ref )      
      return  self.phi_p_ref * (1. + X * (1. + X/2 ) )
      
    def getValueDifferential(self, p):
      """
      returns the porosity for given pressure p
      """
      X= self.C * ( p - self.p_ref )      
      return  self.phi_p_ref * self.C * (1. + X )
          
      
class WaterDensity(MaterialPropertyWithDifferential):
    """
    set water density as 
       
          rho = rho_surf (1 + X + X*X/2) with X= C * ( p - p_ref )
          
    with rho_surf =rho_s/B_ref * gravity
    """
    def __init__(self, B_ref=1., p_ref = 1.*U.atm, gravity=1.,  C = 0./U.bar, rho_s= 999.014 * U.kg/U.m**3):
      self.rho_surf = rho_s * gravity
      self.__rho_0 = self.rho_surf/B_ref
      self.C=C
      self.p_ref=p_ref
    
    def getValue(self, p):
      """
      returns the density for given pressure p
      """
      X= self.C * ( p - self.p_ref )
      return self.__rho_0 * (1+ X * (1+ X/2) )  
      
    def getValueDifferential(self,  p):
      """
      """
      X= self.C * ( p - self.p_ref )
      return self.__rho_0 * self.C * (1+ X) 
      
    def getFormationFactor(self, p):
      return self.rho_surf/self.getValue(p)



class WaterViscosity(MaterialProperty):
  """
  defines viscosity of water
  
  mu=mu_ref /(1 + X + X*X/2)
  
  with X=-C*(p-p_ref)
  """

  def __init__(self, mu_ref = 0.3*U.cPoise, p_ref = 1.*U.atm , C = 0./U.bar):
    """
    set up
    """
    self.mu_ref = mu_ref
    self.p_ref = p_ref
    self.C=C
  def getValue(self, p):
      """
      returns the viscosity for given pressure p
      """
      X= -self.C * ( p - self.p_ref )
      return self.mu_ref/(1+ X * (1+ X/2) )  
      
class GasDensity(MaterialPropertyWithDifferential):
    """
    set water density as 
       
          rho = gravity * rho_air_s /B(p) 
          
    where B is given by an interpolation table
    """
    def __init__(self, p, B, gravity=1., rho_air_s=1.2041*U.kg/U.m**3):
      self.rho_surf =rho_air_s * gravity
      self.tab=InterpolationTable(x=p, y=B)
    
    def getValue(self, p):
      """
      returns the density for given pressure p
      """
      return self.rho_surf/self.getFormationFactor(p)
      
    def getValueDifferential(self,  p):
      """
      """
      B    = self.getFormationFactor(p)
      dBdp = self.getFormationFactorDifferential(p)
      return  -self.rho_surf * dBdp/(B * B)
      
    def getFormationFactor(self, p):
      return self.tab.getValue(p)
      
    def getFormationFactorDifferential(self, p):
      return self.tab.getValueDifferential(p)


class InterpolationTable(MaterialPropertyWithDifferential):
   """
   a simple 1D interpolation table for escript Data object with non-equidistant nodes
   """
   def __init__(self,x=[], y=[], obeyBounds=True):
      """
      set up interpolation table. obeyBounds is set an exception is thrown if
      the interpolation argument is below min(x) or larger than max(x). Otherwise
      the value for x is set to y[0] for 
      """
      MaterialPropertyWithDifferential.__init__(self)
      if len(x) != len(y):
         raise ValueError("length of interpolation nodes and value lists need to be identical.")
      if len(x) < 1 :
         raise ValueError("length of interpolation nodes a list needs to at least one.")
      
      x_ref=x[0]
      for i in range(1,len(x)):
         if x_ref >= x[i]:
            raise ValueError("interpolation nodes need to be given in increasing order.")
         x_ref=x[i]
      self.__x = x 
      self.__y = y
      self.__obeyBounds = obeyBounds
      
   def getValue(self, x):
      """
      returns the interpolated values of x
      """
      X=self.__x
      Y=self.__y
      
      x0=X[0]
      m0=whereNegative(x-x0)
      if self.__obeyBounds:
         if sup(m0) > 0:
            raise ValueError("interpolation argument out of range [%e, %e]"%(X[0],X[-1]))
      out=self.__x[0]
      for i in range(1,len(X)):
          z=(Y[i]-Y[i-1])/(X[i]-X[i-1]) * (x-X[i-1]) + Y[i-1]
          out = out * m0 + z * (1-m0)
          m0=whereNonPositive(x-X[i])
          
      if self.__obeyBounds:
            if inf(m0) < 1 :
               raise ValueError("interpolation argument out of range [%e, %e]"%(X[0],X[-1]))
      else:
            out = out * m0 + X[-1] * (1-m0)
      return out

   def getValueDifferential(self, x):
      X=self.__x
      Y=self.__y

      x0=X[0]
      m0=whereNegative(x-x0)
      if self.__obeyBounds:
         if sup(m0) > 0:
            raise ValueError("interpolation argument out of range [%e, %e]"%(X[0],X[-1]))
      out=0.
      for i in range(1,len(X)):
          z=(Y[i]-Y[i-1])/(X[i]-X[i-1])
          out = out * m0 + z * (1-m0)
          m0=whereNonPositive(x-X[i])
     
      if self.__obeyBounds:
            if inf(m0) < 1:
               raise ValueError("interpolation argument out of range [%e, %e]"%(X[0],X[-1]))
      else:
            out = out * m0
      return out   
      

class Well(object):
   """
   generic well
   
   :var WATER: phase identifier for water well
   :var GAS: phase identifier for gas well
   """
   WATER="Water"
   GAS="Gas"
   def __init__(self, name, domain, Q=[0.], schedule=[0.], phase="Water", BHP_limit=1.*U.atm, X0=[0.,0.,0.], *args, **kwargs):
       """
       set up well 
       """
       if not len(schedule) == len(Q):
           raise ValueError("length of schedule and Q must match.")
       self.__schedule=schedule
       self.__Q = Q
       self.__phase=phase
       self.__BHP_limit=BHP_limit
       self.name=name
       self.domain=domain
       self.locator=Locator(DiracDeltaFunctions(self.domain),X0[:self.domain.getDim()])
       self.X0=self.locator.getX()
       
   def getLocation(self):
       return self.X0      

   def getProductivityIndex(self):
       """
       returns the productivity index of the well.
       typically a Gaussian profile around the location of the well.
       
       :note: needs to be overwritten
       """
       raise NotImplementedError
   
   def getFlowRate(self,t):
      """
      returns the flow rate
      """
      out=0.
      for i in range(len(self.__schedule)):
         if t <= self.__schedule[i]:
             out=self.__Q[i]
             break
      return out
         
   def getBHPLimit(self):
      """
      return bottom-hole pressure limit
      
      :note: needs to be overwritten
      """
      return self.__BHP_limit
      
   def getPhase(self):
      """
      returns the pahse the well works on
      """
      return self.__phase
      
class VerticalPeacemanWell(Well):
   """
   defines a well using Peaceman's formula
   """
   def __init__(self,name, domain, schedule = [0.], BHP_limit=1.*U.atm, Q=0, r=10*U.cm, X0=[0.,0.,0.], D=[1.*U.km,1.*U.km, 1.*U.m], 
                     perm=[1.*U.cPoise,1.*U.cPoise, 1.*U.cPoise], phase=Well.WATER, s=0, decay_factor = 5):
                     # reset_factor=0.1):
       """
       set up well 
       
       :param BHP: ottom-hole pressure
       :param Q: flow rate
       :param r: well radius
       :param X: location 
       :param D: dimension of well block
       :param perm: permeabilities
       :param s: skin factor
       """
       Well.__init__(self, name, domain, Q=Q, schedule=schedule, phase=phase,BHP_limit=BHP_limit, X0=X0)
       r_el=0.28 * sqrt( sqrt(perm[1]/perm[0]) * D[0]**2 +  sqrt(perm[0]/perm[1]) * D[1]**2 )\
                         / ( (perm[1]/perm[0])**0.25 + (perm[1]/perm[0])**0.25 )

       r_el=0.11271331774384821*D[0]
       self.__PI = 2 * math.pi * D[2] * sqrt(perm[1]*perm[0]) / (log(r_el/r) + s)
       self.__r = r 

       self.__D=D
       self.r_el=r_el
 
       if self.domain.getDim() == 3:
           self.geo_fac=1
       else:
           self.geo_fac=1./D[2]


       
   def getProductivityIndex(self):
       """
       returns the productivity index of the well.
       """
       return self.__PI
     

class DualPorosity(object):
   """
   generic dual porosity model 
   """
   def __init__(self, domain, phi_f=None, phi_m=None, phi=None,
                      S_fg=None, S_mg=None, 
                      perm_f_0=None, perm_f_1=None, perm_f_2=None,
                      perm_m_0=None, perm_m_1=None, perm_m_2=None,
                      k_w =None, k_g=None, mu_w =None, mu_g =None,
                      rho_w =None, rho_g=None, 
                      wells=[], g=9.81*U.m/U.sec**2):
      """
      set up
      
      :param domain: domain
      :type domain: `esys.escript.Domain`
      :param phi_f: porosity of the fractured rock as function of the gas pressure
      :type phi_f: `MaterialPropertyWithDifferential`
      :param phi_m: porosity of the coal matrix as function of the gas pressure, if None phi_m = phi-phi_f
      :type phi_m: `MaterialPropertyWithDifferential` or None
      :param phi: total porosity if None phi=1.
      :type phi: scalar or None
      :param S_fg: gas saturation in the fractured rock as function of the capillary pressure p_fc=S_fg- p_f
      :type S_fg: `MaterialPropertyWithDifferential`
      :param S_mg: gas saturation in the coal matrix as function of the capillary pressure p_mc=p_mg-p_mw
      :type S_mg: `MaterialPropertyWithDifferential`
      :param perm_f_0: permeability in the x_0 direction in the fractured media
      :type perm_f_0: scalar 
      :param perm_f_1: permeability in the x_1 direction in the fractured media
      :type perm_f_1: scalar
      :param perm_f_2: permeability in the x_2 direction in the fractured media
      :type perm_f_2: scalar
      :param perm_m_0: permeability in the x_0 direction in the coal matrix
      :type perm_m_0: scalar 
      :param perm_m_1: permeability in the x_1 direction in the coal matrix
      :type perm_m_1: scalar
      :param perm_m_2: permeability in the x_2 direction in the coal matrix
      :type perm_m_2: scalar
      :param k_w: relative permeability of water as function of water saturation
      :type k_w: `MaterialProperty`
      :param k_g: relative permeability of gas as function of gas saturation
      :type k_g: `MaterialProperty`
      :param mu_w: viscosity of water as function of water pressure
      :type mu_w: `MaterialProperty`
      :param mu_g: viscosity of gas as function of gas pressure
      :type mu_g: `MaterialProperty`
      :param rho_w: density of water as function of water pressure
      :type rho_w: `MaterialPropertyWithDifferential`
      :param rho_g: density of gas as function of gas pressure
      :type rho_g: `MaterialPropertyWithDifferential`
      :param wells : list of wells
      :type wells: list of `Well`
      """

      self.domain = domain
      self.phi_f = phi_f
      self.phi_m = phi_m
      self.phi = phi
      self.S_fg = S_fg
      self.S_mg = S_mg
      self.perm_f = [ perm_f_0, perm_f_1, perm_f_2 ]
      self.perm_m = [ perm_m_0, perm_m_1, perm_m_2 ]
      self.k_w = k_w
      self.k_g = k_g
      self.mu_w = mu_w
      self.mu_g = mu_g
      self.rho_w = rho_w
      self.rho_g = rho_g    
      self.wells=wells
      self.t =0
      self.g=g
      
      self.__iter_max=1
      self.__rtol=1.e-4
      self.verbose=False
      self.XI=0.5
   def setIterationControl(self, iter_max=None, rtol=None, verbose=None, xi=None):
     """
     sets parameters to control iteration process
     """
     if iter_max !=None: self.__iter_max=iter_max
     if rtol !=None: self.__rtol = rtol
     if verbose !=None: self.verbose=verbose
     if xi !=None: self.XI=xi
   
   def update(self, dt): 
         self.u, self.u_old =tuple([ v.copy() for v in self.u ] ), self.u
         n=0
         converged=False
         while n < self.__iter_max and not converged:
            u=self.solvePDE(dt)
            if self.verbose: print("iteration step %d:"%n)
            converged=True
            for i in range(len(u)):
               if isinstance(u[i], Data):
                 norm_u=Lsup(u[i])
                 norm_e=Lsup(u[i]-self.u[i])
               else:
                 norm_e=0.
                 norm_u=1.
            
               if norm_u > 0:
                   rerr=norm_e/norm_u
               else:
                   rerr=norm_e
               if norm_e>self.__rtol * norm_u + 1.e-10: converged=False
               if self.verbose: print("   comp %i: change = %e (value = %e)"%(i, norm_e,norm_u))
            n+=1
            self.u=u
         print("iteration completed.")
         self.t+=dt


class PorosityOneHalfModel(DualPorosity):
      """
      Model for gas and water in double prosity model tracking water and gas 
      pressure in the fractured  rock and gas concentration in the coal matrix.
      This corresponds to the coal bed methan model in the ECLIPSE code.
      """ 

      def __init__(self, domain, phi_f=None, phi=None, L_g=None, 
                         perm_f_0=None, perm_f_1=None, perm_f_2=None,
                         k_w =None, k_g=None, mu_w =None, mu_g =None,
                         rho_w =None, rho_g=None, sigma=0, A_mg=0, f_rg=1.,   
                           wells=[], g=9.81*U.m/U.sec**2):
         """
         set up
         
         :param domain: domain
         :type domain: `esys.escript.Domain`
         :param phi_f: porosity of the fractured rock as function of the gas pressure
         :type phi_f: `MaterialPropertyWithDifferential`
         :param phi: total porosity if None phi=1.
         :type phi: scalar or None
         :param L_g: gas adsorption as function of gas pressure
         :type L_g: `MaterialPropertyWithDifferential`
         :param S_fg: gas saturation in the fractured rock as function of the capillary pressure p_fc=S_fg- p_f
         :type S_fg: `MaterialPropertyWithDifferential`
         :param perm_f_0: permeability in the x_0 direction in the fractured media
         :type perm_f_0: scalar 
         :param perm_f_1: permeability in the x_1 direction in the fractured media
         :type perm_f_1: scalar
         :param perm_f_2: permeability in the x_2 direction in the fractured media
         :type perm_f_2: scalar
         :param k_w: relative permeability of water as function of water saturation
         :type k_w: `MaterialProperty`
         :param k_g: relative permeability of gas as function of gas saturation
         :type k_g: `MaterialProperty`
         :param mu_w: viscosity of water as function of water pressure
         :type mu_w: `MaterialProperty`
         :param mu_g: viscosity of gas as function of gas pressure
         :type mu_g: `MaterialProperty`
         :param rho_w: density of water as function of water pressure
         :type rho_w: `MaterialPropertyWithDifferential`
         :param rho_g: density of gas as function of gas pressure
         :type rho_g: `MaterialPropertyWithDifferential`
         :param wells : list of wells
         :type wells: list of `Well`
         :param sigma: shape factor for gas matrix diffusion 
         :param A_mg: diffusion constant for gas adsorption
         :param f_rg: gas re-adsorption factor
         """
   
         DualPorosity.__init__(self, domain,
                              phi_f=phi_f, phi_m=None, phi=phi,
                              S_fg=None, S_mg=None, 
                              perm_f_0=perm_f_0, perm_f_1=perm_f_1, perm_f_2=perm_f_2,
                              perm_m_0=None , perm_m_1=None, perm_m_2=None, 
                              k_w =k_w, k_g=k_g, mu_w =mu_w, mu_g =mu_g,
                              rho_w =rho_w, rho_g=rho_g, 
                              wells=wells, g=g)
         self.L_g=L_g
         self.sigma = sigma 
         self.A_mg = A_mg
         self.f_rg  = f_rg
         self.__pde=LinearPDE(self.domain, numEquations=3, numSolutions =3)
         
    
      def getPDEOptions(self):
         """
         returns the `SolverOptions` of the underlying PDE
         """
         return self.__pde.getSolverOptions()
         
      def getState(self): 
         return self.u

      def getOldState(self): 
         return self.u_old

      def setInitialState(self, p_top=1.*U.atm, p_bottom= 1.*U.atm, S_fg=0,  c_mg=None):
            """
            sets the initial state
            
            :param p: pressure
            :param S_fg: gas saturation in fractured rock 
            :param c_mg: gas concentration in coal matrix at standart conditions. if not given it is calculated 
                        using the  gas adsorption curve.
            """    
            if self.domain.getDim() == 2:
               p_init=Scalar((p_top + p_bottom) /2, Solution(self.domain))
            else:
               z=self.u.getDomain().getX()[0]
               z_top=sup(z)
               z_bottom=inf(z)
               p_init=(p_bottom-p_top)/(z_bottom-z_top)*(z-z_top) + p_top

            S_fg_init=Scalar(S_fg, Solution(self.domain))

            if c_mg == None:
              c_mg_init=self.L_g(p_init)
            else:
              c_mg_init=Scalar(c_mg, Solution(self.domain))

            q_gas=Scalar(0., DiracDeltaFunctions(self.domain))
            q_water=Scalar(0., DiracDeltaFunctions(self.domain))
            BHP=interpolate(p_init, DiracDeltaFunctions(self.domain))

            self.u=(p_init,S_fg_init, c_mg_init, BHP, q_gas, q_water)

      def solvePDE(self, dt):
         
         p_f, S_fg, c_mg, BHP_old, q_gas_old, q_water_old =self.getState() 
         p_f_old, S_fg_old, c_mg_old, BHP, q_gas, q_water =self.getOldState()

         S_fw=1-S_fg

         if self.verbose: 
              print("p_f range = ",inf(p_f),sup(p_f)) 
              print("S_fg range = ",inf(S_fg),sup(S_fg))
              print("S_fw range = ",inf(S_fw),sup(S_fw))
              print("c_mg range = ",inf(c_mg),sup(c_mg))
              print("BHP =",BHP)
              print("q_gas =",q_gas)
              print("q_water =",q_water)

         k_fw=self.k_w(S_fw)
         if self.verbose: print("k_fw range = ",inf(k_fw),sup(k_fw)) 


         k_fg=self.k_g(S_fg)
         if self.verbose: print("k_fg range = ",inf(k_fg),sup(k_fg)) 

         mu_fw=self.mu_w(p_f)
         if self.verbose: print("mu_fw range = ",inf(mu_fw),sup(mu_fw)) 

         mu_fg=self.mu_g(p_f)
         if self.verbose: print("mu_fg range = ",inf(mu_fg),sup(mu_fg)) 
         

         phi_f   =self.phi_f.getValue(p_f)
         dphi_fdp=self.phi_f.getValueDifferential(p_f)
         if self.verbose: print("phi_f range = ",inf(phi_f),sup(phi_f)," (slope %e,%e)"%(inf(dphi_fdp), sup(dphi_fdp))) 
         
         rho_fw         = self.rho_w.getValue(p_f)
         drho_fwdp        = self.rho_w.getValueDifferential(p_f)
         if self.verbose: print("rho_fw range = ",inf(rho_fw),sup(rho_fw)," (slope %e,%e)"%(inf(drho_fwdp), sup(drho_fwdp))) 

         rho_fg = self.rho_g.getValue(p_f)
         rho_g_surf = self.rho_g.rho_surf
         drho_fgdp = self.rho_g.getValueDifferential(p_f)
         if self.verbose: 
              print("rho_fg range = ",inf(rho_fg),sup(rho_fg)," (slope %e,%e)"%(inf(drho_fgdp), sup(drho_fgdp))) 
              print("rho_fg surf = ",rho_g_surf)
              
         L_g = self.L_g(p_f)
         dL_gdp = self.L_g.getValueDifferential(p_f)
         if self.verbose: print("L_g range = ",inf(L_g),sup(L_g)," (slope %e,%e)"%(inf(dL_gdp), sup(dL_gdp))) 
                  
         A_fw = rho_fw * k_fw/mu_fw 
         A_fg = rho_fg * k_fg/mu_fg
         
         m = whereNegative(L_g - c_mg) 
         H = self.sigma * self.A_mg * (m + (1-m) * self.f_rg * S_fg )
         
         E_fpp= S_fw * (  rho_fw * dphi_fdp + phi_f  * drho_fwdp )
         E_fps=  -  phi_f * rho_fw 
         E_fsp= S_fg * ( rho_fg * dphi_fdp + phi_f * drho_fgdp )
         E_fss= phi_f * rho_fg 
        
        
         
         F_fw=0.
         F_fg=0.
         D_fw=0.
         D_fg=0.

         prod_index=Scalar(0., DiracDeltaFunctions(self.domain))
         geo_fac=Scalar(0., DiracDeltaFunctions(self.domain))
         for I in self.wells:
             prod_index.setTaggedValue(I.name, I.getProductivityIndex() )
             geo_fac.setTaggedValue(I.name, I.geo_fac)

         F_fp_Y = A_fw * prod_index * BHP  * geo_fac
         F_fs_Y = A_fg * prod_index * BHP * geo_fac
         D_fpp =  A_fw * prod_index * geo_fac
         D_fsp =  A_fg * prod_index * geo_fac

         
         if self.domain.getDim() == 3:
            F_fp_X = ( A_fw * self.g * rho_fw * self.perm_f[2] ) * kronecker(self.domain)[2]
            F_fs_X = ( A_fg * self.g * rho_fg * self.perm_f[2] ) * kronecker(self.domain)[2]
         else:
            F_fp_X = 0 * kronecker(self.domain)[1]
            F_fs_X = 0 * kronecker(self.domain)[1]
            
         F_mg_Y = H * L_g


         D=self.__pde.createCoefficient("D")
         A=self.__pde.createCoefficient("A")
         Y=self.__pde.createCoefficient("Y")
         X=self.__pde.createCoefficient("X") 
         
         d_dirac=self.__pde.createCoefficient("d_dirac")
         y_dirac=self.__pde.createCoefficient("y_dirac")


        
         D[0,0]=E_fpp
         D[0,1]=E_fps
         d_dirac[0,0]=dt * D_fpp
         
         D[1,0]=E_fsp
         D[1,1]=E_fss 
         D[1,2]= rho_g_surf
         d_dirac[1,0]=dt * D_fsp
         
         D[0,2]= -dt * H * dL_gdp * 0
         D[2,2]= 1 + dt * H
         
         
         for i in range(self.domain.getDim()):
            A[0,i,0,i] = dt * A_fw * self.perm_f[i]
            A[1,i,1,i] = dt * A_fg * self.perm_f[i]

         X[0,:]=  dt * F_fp_X
         X[1,:]=  dt * F_fs_X

         Y[0] = E_fpp *  p_f_old + E_fps * S_fg_old 
         Y[1] = E_fsp *  p_f_old + E_fss * S_fg_old + rho_g_surf * c_mg_old 
         Y[2] = c_mg_old                                                    + dt * ( F_mg_Y -  H * dL_gdp * p_f *0 )
         
         y_dirac[0] =dt * F_fp_Y
         y_dirac[1] =dt * F_fs_Y 
         
         print("HHH D[0,0] = ",D[0,0])
         print("HHH D[0,1] = ",D[0,1])
         print("HHH D[0,2] = ",D[0,2])
         print("HHH D[1,0] = ",D[1,0])
         print("HHH D[1,1] = ",D[1,1])
         print("HHH D[1,2] = ",D[1,2])
         print("HHH D[2,0] = ",D[2,0])
         print("HHH D[2,1] = ",D[2,1])
         print("HHH D[2,2] = ",D[2,2])
         print("HHH A_fw = ",A_fw)
         print("HHH A_fg = ",A_fg)
         print("HHH A[0,0,0,0] = ",A[0,0,0,0])
         print("HHH A[0,1,0,1] = ",A[0,1,0,1])
         print("HHH A[1,0,1,0] = ",A[1,0,1,0])
         print("HHH A[1,1,1,1] = ",A[1,1,1,1])
         print("HHH Y[0] ",Y[0])
         print("HHH Y[1] ",Y[1])
         print("HHH Y[2] ",Y[2])
         print("HHH F_fp_Y ",F_fp_Y)
         print("HHH F_fs_Y ",F_fs_Y)
         print("HHH F_mg_Y ",F_mg_Y)
         print("HHH H = ",H)

         self.__pde.setValue(A=A, D=D, X=X, Y=Y, d_dirac=d_dirac , y_dirac=y_dirac)
         
         u2 = self.__pde.getSolution()
         # this is deals with round-off errors to maintain physical meaningful values
         # we should do this in a better way to detect values that are totally wrong
         p_f=u2[0]
         S_fg=clip(u2[1],minval=0, maxval=1)
         c_mg=clip(u2[2],minval=0)
         


         q=Scalar(0., DiracDeltaFunctions(self.domain))
         BHP_limit=Scalar(0., DiracDeltaFunctions(self.domain))
         water_well_mask=Scalar(0., DiracDeltaFunctions(self.domain))
         for I in self.wells:
             q.setTaggedValue(I.name, I.getFlowRate(self.t+dt))
             BHP_limit.setTaggedValue(I.name, I.getBHPLimit())
             if I.getPhase() == Well.WATER:
                  water_well_mask.setTaggedValue(I.name, 1)
         
         p_f_wells=interpolate(p_f, DiracDeltaFunctions(self.domain))
         S_fg_wells=interpolate(S_fg, DiracDeltaFunctions(self.domain))
         A_fw_wells= self.k_w(1-S_fg_wells)/self.mu_w(p_f_wells)*self.rho_w(p_f_wells)
         A_fg_wells= self.k_g(S_fg_wells)/self.mu_g(p_f_wells)*self.rho_g(p_f_wells)
         
         BHP=clip( p_f_wells - q/prod_index * ( m * self.rho_w.rho_surf/A_fw_wells + (1-m)* self.rho_g.rho_surf/A_fg_wells ), minval = BHP_limit)
         BHP=clip( p_f_wells - q/prod_index * self.rho_w.rho_surf/A_fw_wells, minval = BHP_limit)
         q_gas    = prod_index * A_fg_wells * (p_f_wells - BHP)/self.rho_g.rho_surf
         q_water = prod_index * A_fw_wells * (p_f_wells - BHP)/self.rho_w.rho_surf

         return (p_f,S_fg, c_mg, BHP, q_gas, q_water )
