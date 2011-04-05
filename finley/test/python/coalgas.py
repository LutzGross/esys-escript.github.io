######################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################
"""
Gas in Coal Seam (fully coupled version)
"""
__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from esys.escript.linearPDEs import LinearPDE
from esys.escript import unitsSI as U
from esys.escript.pdetools import Locator
from esys.escript import sqrt, log, whereNegative, sup, inf, whereNonPositive, integrate, Function, kronecker, grad, Lsup, clip
from esys.weipa import saveVTK
import math

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
       
          rho = rho_ref (1 + X + X*X/2) with X= C * ( p - p_ref )
          
    with rho_ref =rho_s/B_ref * gravity
    """
    def __init__(self, B_ref=1., p_ref = 1.*U.atm, gravity=1.,  C = 0./U.bar, rho_s= 998.2 * U.kg/U.m**3):
      self.rho_0 = rho_s * gravity
      self.rho_ref = self.rho_0/B_ref
      self.C=C
      self.p_ref=p_ref
    
    def getValue(self, p):
      """
      returns the density for given pressure p
      """
      X= self.C * ( p - self.p_ref )
      return self.rho_ref * (1+ X * (1+ X/2) )  
      
    def getValueDifferential(self,  p):
      """
      """
      X= self.C * ( p - self.p_ref )
      return self.rho_ref * self.C * (1+ X) 

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
      self.rho_ref =rho_air_s * gravity
      self.rho_0 =rho_air_s * gravity
      self.tab=InterpolationTable(x=p, y=B)
    
    def getValue(self, p):
      """
      returns the density for given pressure p
      """
      return self.rho_ref/self.tab.getValue(p)
      
    def getValueDifferential(self,  p):
      """
      """
      B    = self.tab.getValue(p)
      dBdp = self.tab.getValueDifferential(p)
      return  -self.rho_ref * dBdp/(B * B)

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
	 raise ValueError,"length of interpolation nodes and value lists need to be identical."
      if len(x) < 1 :
	 raise ValueError,"length of interpolation nodes a list needs to at least one."
      
      x_ref=x[0]
      for i in range(1,len(x)):
	 if x_ref >= x[i]:
	    raise ValueError,"interpolation nodes need to be given in increasing order."
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
	    raise ValueError,"interpolation argument out of range [%e, %e]"%(X[0],X[-1])
      out=self.__x[0]
      for i in range(1,len(X)):
	  z=(Y[i]-Y[i-1])/(X[i]-X[i-1]) * (x-X[i-1]) + Y[i-1]
	  out = out * m0 + z * (1-m0)
	  m0=whereNonPositive(x-X[i])
          
      if self.__obeyBounds:
	    if inf(m0) < 1 :
	       raise ValueError,"interpolation argument out of range [%e, %e]"%(X[0],X[-1])
      else:
	    out = out * m0 + V[-1] * (1-m0)
      return out

   def getValueDifferential(self, x):
      X=self.__x
      Y=self.__y

      x0=X[0]
      m0=whereNegative(x-x0)
      if self.__obeyBounds:
	 if sup(m0) > 0:
	    raise ValueError,"interpolation argument out of range [%e, %e]"%(X[0],X[-1])
      out=0.
      for i in range(1,len(X)):
	  z=(Y[i]-Y[i-1])/(X[i]-X[i-1])
	  out = out * m0 + z * (1-m0)
	  m0=whereNonPositive(x-X[i])
     
      if self.__obeyBounds:
	    if inf(m0) < 1:
	       raise ValueError,"interpolation argument out of range [%e, %e]"%(X[0],X[-1])
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
   def __init__(self, name, domain, Q=0., schedule=[0.], phase="Water", BHP_limit=1.*U.atm, *args, **kwargs):
       """
       set up well 
       """
       self.__schedule=schedule
       self.__phase=phase
       self.__Q = Q
       self.__BHP_limit=BHP_limit
       self.name=name
       self.domain=domain
       
   def getGeometry(self):
       """
       returns a function representing the geometry of the well.
       typically a Gaussian profile around the location of the well.
       
       :note: needs to be overwritten
       """
       raise NotImplementedError

   def getProductivityIndex(self):
       """
       returns the productivity index of the well.
       typically a Gaussian profile around the location of the well.
       
       :note: needs to be overwritten
       """
       raise NotImplementedError
   
   def getFlowRate(self):
      """
      returns the flow rate
      """
      return self.__Q
 	
   def getBHPLimit(self):
      """
      return bottom-hole pressure limit
      
      :note: needs to be overwritten
      """
      return self.__BHP_limit
      
   def getBHP(self):
      """
      return bottom-hole pressure
      
      :note: needs to be overwritten
      """
      return self.__BHP
   
   def setBHP(self, BHP):
      """
      sets bottom-hole pressure
      """
      self.__BHP= BHP
      
   def getPhase(self):
      """
      returns the pahse the well works on
      """
      return self.__phase
      
   def isOpen(self, t):
     """
     returns True is the well is open at time t
     """
     out = False
     t0=min(t, self.__schedule[0]) 
     i=0
     while t > t0:
       if out:
	 out=False
       else:
	 out=True  
       i+=1
       if i < len(self.__schedule):
	  t0=self.__schedule[i]
       else:
	  t0=t
     return out

class VerticalPeacemanWell(Well):
   """
   defines a well using Peaceman's formula
   """
   def __init__(self,name, domain, schedule = [0.], BHP_limit=1.*U.atm, Q=0, r=10*U.cm, X0=[0.,0.,0.], D=[1.*U.km,1.*U.km, 1.*U.m], 
		     perm=[1.*U.cPoise,1.*U.cPoise, 1.*U.cPoise], phase=Well.WATER, s=0):
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
       Well.__init__(self, name, domain, Q=Q, schedule=schedule, phase=phase,BHP_limit=BHP_limit)
       r_el=0.28 * sqrt( sqrt(perm[1]/perm[0]) * D[0]**2 +  sqrt(perm[0]/perm[1]) * D[1]**2 )\
                         / ( (perm[1]/perm[0])**0.25 + (perm[1]/perm[0])**0.25 )
       self.__PI = 2 * math.pi * D[2] * sqrt(perm[1]*perm[0]) / (log(r_el/r) + s)
       self.__r = r 

       self.__D=D
       self.r_el=r_el
       self.X0=X0[:self.domain.getDim()]
       
       x=Function(domain).getX()
       self.chi = 1.
       for i in range(domain.getDim()):
	    self.chi = self.chi * whereNegative(abs(x[i]-X0[i])-D[i]/2)

       self.chi*=1./(D[0]*D[1]*D[2])
       
       
       #self.chi=0.00000001
   def getLocation(self):
       return self.X0
       
   def getGeometry(self):
       """
       returns a function representing the geometry of the well.
       typically a Gaussian profile around the location of the well.
       
       :note: needs to be overwritten
       """
       return self.chi

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
                      wells=[]):
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
         self.u, self.u_old = self.u.copy(), self.u
         n=0
         rerr=1.
         while n < self.__iter_max and rerr > self.__rtol:
            u=self.solvePDE(dt)
	    norm_u=Lsup(u)
	    norm_e=Lsup(u-self.u)
	    
	    if norm_u > 0:
	        rerr=norm_e/norm_u
	    else:
	        rerr=norm_e
	    if self.verbose: print "iteration step %d: relative change = %e"%(n, rerr)
	    n+=1
	    self.u=u
         print "iteration completed."
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
			   wells=[]):
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
			      wells=wells)
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
	 return self.u[0], self.u[1],  self.u[2]

      def getOldState(self): 
	 return self.u_old[0], self.u_old[1],  self.u_old[2]

      def setInitialState(self, p=1.*U.atm, S_fg=0,  C_mg=None):
	    """
	    sets the initial state
	    
	    :param p: pressure
	    :param S_fg: gas saturation in fractured rock 
	    :param C_mg: gas concentration in coal matrix. if not given it is calculated 
			using the  gas adsorption curve.
	    """    
	    self.u=self.__pde.createSolution()
	    self.u[0]=p
	    self.u[1]=S_fg
	    if C_mg == None:
	      self.u[2]= self.L_g(p)
	    else:
	      self.u[2]=C_mg
	  
      def solvePDE(self, dt):
	 
	 C_couple=0
	 
	 p_f, S_fg, C_mg=self.getState() 
	 p_f_old, S_fg_old, C_mg_old=self.getOldState()

         S_fw=1-S_fg

	 if self.verbose: 
	      print "p_f range = ",inf(p_f),sup(p_f) 
	      print "S_fg range = ",inf(S_fg),sup(S_fg)
	      print "S_fw range = ",inf(S_fw),sup(S_fw)
	      print "C_mg range = ",inf(C_mg),sup(C_mg)

         k_fw=self.k_w(S_fw)
       	 if self.verbose: print "k_fw range = ",inf(k_fw),sup(k_fw) 

         k_fg=self.k_g(S_fg)
       	 if self.verbose: print "k_fg range = ",inf(k_fg),sup(k_fg) 

	 mu_fw=self.mu_w(p_f)
       	 if self.verbose: print "mu_fw range = ",inf(mu_fw),sup(mu_fw) 

	 mu_fg=self.mu_g(p_f)
       	 if self.verbose: print "mu_fg range = ",inf(mu_fg),sup(mu_fg) 
         

	 phi_f   =self.phi_f.getValue(p_f)
	 dphi_fdp=self.phi_f.getValueDifferential(p_f)
       	 if self.verbose: print "phi_f range = ",inf(phi_f),sup(phi_f)," (slope %e,%e)"%(inf(dphi_fdp), sup(dphi_fdp)) 
	 
	 rho_fw 	= self.rho_w.getValue(p_f)
	 drho_fwdp	= self.rho_w.getValueDifferential(p_f)
      	 if self.verbose: print "rho_fw range = ",inf(rho_fw),sup(rho_fw)," (slope %e,%e)"%(inf(drho_fwdp), sup(drho_fwdp)) 

         rho_fg = self.rho_g.getValue(p_f)
	 drho_fgdp	= self.rho_g.getValueDifferential(p_f)
      	 if self.verbose: print "rho_fg range = ",inf(rho_fg),sup(rho_fg)," (slope %e,%e)"%(inf(drho_fgdp), sup(drho_fgdp)) 
	 
	 L_g       = self.L_g.getValue(p_f)
	 dL_gdp =  self.rho_w.getValueDifferential(p_f)
      	 if self.verbose: print "L_g range = ",inf(L_g),sup(L_g)," (slope %e,%e)"%(inf(dL_gdp), sup(dL_gdp)) 
	 
	 
	 A_fw = rho_fw * k_fw/mu_fw 
	 A_fg = rho_fg * k_fg/mu_fg
	 
	 m = whereNegative(L_g - C_mg)  
	 B = self.sigma * self.A_mg * (m + (1-m) * self.f_rg * S_fg )
	 
	 
	 E_fpp= S_fw * (  rho_fw * dphi_fdp + phi_f  * drho_fwdp )
	 E_fps=  -  phi_f * rho_fw 
	 E_fsp= S_fg *( rho_fg * dphi_fdp - phi_f *  drho_fgdp )
	 E_fss= phi_f * rho_fg 
	 
	 F_fw=0.
	 F_fg=0.
	 D_fw=0.
	 D_fg=0.
	 
	 for I in self.wells:

	      chi_I= I.getGeometry()
	      loc=Locator(Function(self.domain),I.getLocation())
	      Pi_I = I.getProductivityIndex()
	      A_fw_I= loc(A_fw)
	      A_fg_I= loc(A_fg)
	      BHP_limit_I=I.getBHPLimit()
	      
	      if I.isOpen(self.t+dt):
		if self.verbose: print "well %s is open."%I.name
		if I.getPhase() == Well.WATER:
		    q_I = self.rho_w.rho_0 * I.getFlowRate()
		    p_I_ref=loc(p_f)-q_I/(A_fw_I * Pi_I) 
		else:
		    q_I = self.rho_g.rho_0 * I.getFlowRate()
		    p_I_ref=loc(p_f)-q_I/(A_fg_I * Pi_I) 
		    
		print "ZZZ =",loc(p_f), q_I, self.rho_w.rho_0, I.getFlowRate(), A_fw_I, Pi_I, q_I/(A_fw_I * Pi_I)
		1/0
		if BHP_limit_I > p_I_ref:
		  BHP_I=BHP_limit_I
		  D_fw = D_fw + Pi_I * A_fw_I *              chi_I
		  F_fw = F_fw - Pi_I * A_fw_I * BHP_limit_I *chi_I 
		  D_fg = D_fg + Pi_I * A_fg_I *              chi_I
		  F_fg = F_fg - Pi_I * A_fg_I * BHP_limit_I *chi_I 
		else:
		  if I.getPhase() == Well.WATER:
		      F_fw = F_fw +  q_I  * chi_I 
		      F_fg = F_fg +  A_fg_I/A_fw_I *  q_I *chi_I 
		  else:
		      F_fg = F_fg +  q_I  * chi_I 
		      F_fw = F_fw +  A_fw_I/A_fg_I *  q_I *chi_I 
		  BHP_I=p_I_ref
	      else:
		  if self.verbose: print "well %s is closed."%I.name
		  BHP_I=loc(p_f)
		  
	      if self.verbose: print "well %s:BHP = %e"%(I.name, BHP_I)
	      I.setBHP(BHP_I)

	 F_fp_Y = - F_fw
	 F_fs_Y = - F_fg
	 D_fpp =  D_fw
	 D_fsp =  D_fg

	 
	 if self.domain.getDim() == 3:
	    F_fp_X = ( A_fw * self.g * rho_fw * self.perm_f[2] ) * kronecker(self.domain)[2]
	    F_fs_X = ( A_fg * self.g * rho_fg * self.perm_f[2] ) * kronecker(self.domain)[2]
	 else:
	    F_fp_X = 0 * kronecker(self.domain)[1]
	    F_fs_X = 0 * kronecker(self.domain)[1]
	    
	 F_mg_Y = B * ( L_g - dL_gdp * p_f )



	 D=self.__pde.createCoefficient("D")
	 A=self.__pde.createCoefficient("A")
	 Y=self.__pde.createCoefficient("Y")
	 X=self.__pde.createCoefficient("X") 
	 
	 dtXI = dt*self.XI
	 dtcXI = dt*(1-self.XI)

	 D[0,0]=E_fpp + dtXI * D_fpp
	 D[0,1]=E_fps
	 D[1,0]=E_fsp + dtXI * D_fsp
	 D[1,1]=E_fss 
	 D[1,2]=rho_fg * C_couple
	 D[2,0]= - dtXI * B * dL_gdp
	 D[2,2]= 1 + dtXI * B
	 
	 
	 for i in range(self.domain.getDim()):
	    A[0,i,0,i] = dtXI * A_fw * self.perm_f[i]
	    A[1,i,1,i] = dtXI * A_fg * self.perm_f[i]

	 g= grad(p_f_old) * dtcXI * self.perm_f[0:self.domain.getDim()]
         X[0,:]=  - A_fw * g  + dt * F_fp_X
	 X[1,:]=  - A_fg * g  + dt * F_fs_X

	 Y[0] = E_fpp *  p_f_old + E_fps * S_fg_old +                                dt * F_fp_Y - dtcXI * D_fpp * p_f_old
	 Y[1] = E_fsp *  p_f_old + E_fss * S_fg_old + C_couple * rho_fg * C_mg_old + dt * F_fs_Y - dtcXI * D_fsp * p_f_old
	 Y[2] = C_mg_old                                                           + dt * F_mg_Y - dtcXI * B * ( C_mg_old - dL_gdp * p_f_old)
	 
 	 print "HHH Y[0] =", Y[0]
 	 print "HHH A_fw = ",A_fw
 	 print "HHH A_fg= ",A_fg
 	 print "HHH F_fg = ",F_fg
 	 print "HHH F_fw = ",F_fw
 	 print "HHH perm_f = ",self.perm_f
 	 print "HHH k = ",(A_fw*self.perm_f[0])/D[0,0]
 	 print "HHH k = ",(A_fw*self.perm_f[1])/D[0,0]
 	 print "HHH X[0,:]= ",X[0,:]
	 print "HHH D[0,0] = ",D[0,0]
	 print "HHH D[0,1] = ",D[0,1]
	 print "HHH D[0,2] = ",D[0,2]
	 

       
         
	 #print "HHH Y[1] = ",Y[1]
	 #print "HHH X[1,:] = ",X[1,:]
	 #print "HHH D[1,0] = ",D[1,0]
	 #print "HHH D[1,1] = ",D[1,1]
	 #print "HHH D[1,2] = ",D[1,2]
	 #print "HHH A_fg =",  A_fg
	 



	 #1/0
	 self.__pde.setValue(A=A, D=D, X=X, Y=Y)
	 
	 u2 = self.__pde.getSolution()
	 # this is deals with round-off errors to maintain physical meaningful values
	 # we should do this in a better way to detect values that are totally wrong
	 u2[1]=clip(u2[1],minval=0, maxval=1)
	 u2[2]=clip(u2[2],minval=0)
	 print "p _new =", u2[0]
	 #1/0
	 return u2
