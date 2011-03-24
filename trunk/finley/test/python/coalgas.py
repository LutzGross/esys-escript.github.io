#######################################################
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
from esys.escript import sqrt, log
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
      self.rho_ref = rho_s/B_ref * gravity
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
      x0=self.__x[0]
      m0=whereNegative(x-x0)
      if self.__obeyBounds:
	 if sup(m0) > 0:
	    raise ValueError,"interpolation argument out of range [%e, %e]"%(x[0],x[-1])
      out=self.__x[0]
      for i in range(1,len(self.__x)):
	  x1=self.__x[i]
	  z=(y[i]-y[i-1])/(x[i]-x[i-1]) * (x-x[i-1]) + y[i-1]
	  out = out * m0 + z * (1-m0)
	  m0=whereNegative(x-x[i])
     
      if self.__obeyBounds:
	    if inf(m0) < 1:
	       raise ValueError,"interpolation argument out of range [%e, %e]"%(x[0],x[-1])
      else:
	    out = out * m0 + y[-1] * (1-m0)
      return out

   def getValueDifferential(self, x):
      x0=self.__x[0]
      m0=whereNegative(x-x0)
      if self.__obeyBounds:
	 if sup(m0) > 0:
	    raise ValueError,"interpolation argument out of range [%e, %e]"%(x[0],x[-1])
      out=0.
      for i in range(1,len(self.__x)):
	  x1=self.__x[i]
	  z=(y[i]-y[i-1])/(x[i]-x[i-1])
	  out = out * m0 + z * (1-m0)
	  m0=whereNegative(x-x[i])
     
      if self.__obeyBounds:
	    if inf(m0) < 1:
	       raise ValueError,"interpolation argument out of range [%e, %e]"%(x[0],x[-1])
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
   def __init__(self,  Q=0., schedule=[0.], phase="Water", *args, **kwargs):
       """
       set up well 
       """
       self.__schedule=schedule
       self.__phase=phase
       self.__Q = Q
       
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
   
   def getFlowRate(self, t):
      """
      returns the flow rate
      """
      if self.isOpen(t):
	return self.__Q
      else:
	return 0.
	
   def getBHP(self):
      """
      return bottom-hole pressure
      
      :note: needs to be overwritten
      """
      raise NotImplementedError
      
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
     while t < t0:
       if out:
	 out=False
       else:
	 out=True  
       i+=1
       if i < len(self.__schedule[i]):
	  t0=self.__schedule[i]
       else:
	  t0=t
     return out

class VerticalPeacemanWell(Well):
   """
   defines a well using Peaceman's formula
   """
   def __init__(self,name, schedule = [0.], BHP=1.*U.atm, Q=0, r=10*U.cm, X0=[0.,0.,0.], D=[1.*U.km,1.*U.km, 1.*U.m], 
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
       Well.__init__(self, Q=Q, schedule=schedule, phase=phase,)
       r_el=0.28 * ( (perm[1]/perm[0])**0.5 * D[0]**2 +  (perm[0]/perm[1])**0.5 * D[1]**2 ) \
                         / ( (perm[1]/perm[0])**0.25 + (perm[1]/perm[0])**0.25 )
       self.__PI = 2 * math.pi * D[2] * sqrt(perm[1]*perm[0]) / (log(r_el/r) + s)
       self.__r = r 

       self.__BHP = BHP
       self.__D=D
       self.r_el=r_el
       self.X0=X0
       
   def getGeometry(self):
       """
       returns a function representing the geometry of the well.
       typically a Gaussian profile around the location of the well.
       
       :note: needs to be overwritten
       """
       raise NotImplementedError

   def getProductivityIndex(self,t):
       """
       returns the productivity index of the well.
       typically a Gaussian profile around the location of the well.
       """
       return self.__PI
      
   def getBHP(self):
      """
      return bottom-hole pressure
      """
      return self.__BHP

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
   




class PorosityOneHalfModel(DualPorosity):
      """
      Model for gas and water in double prosity model tracking water and gas 
      pressure in the fractured  rock and gas concentration in the coal matrix.
      This corresponds to the coal bed methan model in the ECLIPSE code.
      """ 

      def __init__(self, domain, phi_f=None, phi=None, L_g=None, 
			 perm_f_0=None, perm_f_1=None, perm_f_2=None,
			 k_w =None, k_g=None, mu_w =None, mu_g =None,
			 rho_w =None, rho_g=None, 
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
	 self.__pde=LinearPDE(self.domain, numEquations=3, numSolutions =3)
	 
    
      def getPDEOptions(self):
	 """
	 returns the `SolverOptions` of the underlying PDE
	 """
	 return self.__pde.getSolverOptions()
	 
      def setInitialState(p=1.*U.atm, S_fg=0,  C_mg=None):
	    """
	    sets the initial state
	    
	    :param p: pressure
	    :param S_fg: gas saturation in fractured rock 
	    :param C_mg: gas concentration in coal matrix. if not given it is calculated 
			using the  gas adsorption curve.
	    """
	    
	    self.u=self.__pde.getNewCoefficient("u")
	    u[0]=p
	    u[1]=S_fg
	    if C_mg == None:
	      u[2]= self.L_g(p)
	    else:
	      u[2]=C_mg
	  
      def solvePDE(self):
	 
	 p_f=self.u[0]
	 S_fg=self.u[1]
	 C_mg=self.u[3] 
      
      
	 p_f_old=self.u_old[0]
	 S_fg_old=self.u_old[1]
	 C_mg_old=self.u_old[3]


	 phi_f   =self.phi_f.getValue(S_fg)
	 dphi_fdp=self.phi_f.getValueDifferential(S_fg)
	 
	 S_fw=	1-S_fg
	 
	 rho_fw 	= self.rho_w.getValue(p_f)
	 drho_fwdp	= self.rho_w.getValueDifferential( p_f)
	 rho_fg = self.rho_g.getValue(p_f)
	 drho_fgdp	= self.rho_g.getValueDifferential(p_f)
	 
	 L_g       = self.getValue(p_f)
	 dL_gdp =  self.rho_w.getValueDifferential(p_f)
	 
	 A_fw = rho_fw * k_fw/mu_fw 
	 A_fg = rho_fg * k_fg/mu_fg
	 
	 m = whereNegative(L_g - C_mg)  
	 B = self.sigma * self.a_mg * (m + (1-m) * self.f_rg * S_fg )
	 
	 
	 E_fpp= S_fw * (  rho_fw * dphi_fdp + phi_f  * drho_fwdp )
	 E_fps=  -  phi_f * rho_fw 
	 E_fsp= S_fg *( rho_fg * dphi_fdp - phi_f *  drho_fgdp )
	 E_fss= phi_f * rho_fg 


	 
	 Q_w=0.
	 Q_g=0.
	 R_w=0.
	 R_g=0.
	 
	 for I in self.wells:
	    chi_I= I.getGeometry()
	    Pi_I = I.getProductivityIndex()
	    P_I = I.getBHP()
	    Q_I = I.getFlowRate()
	    
	    if self.getPhase == Well.WATER:
	       Q_w = Q_w + Q_I        * chi_I 
	       Q_g = Q_g + Pi_I * P_I * chi_I
	       R_g = R_g + Pi_I       * chi_I
	    else:
	       Q_g = Q_g + Q_I        * chi_I 
	       Q_w = Q_w + Pi_I * P_I * chi_I
	       R_w = R_w + Pi_I       * chi_I
	       
	 F_fp_Y = A_fw * Q_w 
	 F_fs_Y = A_fg * Q_g
	 D_fss = - A_fg * R_g 
	 D_fpp = - A_fw * R_w
	 
	 if self.domain.getDim() == 3:
	    F_fp_X = A_fw * self.g * rho_fw * self.perm_f[2] * kronecker(self.domain)[2]
	    F_fs_X = A_fg * self.g * rho_fg * self.perm_f[2] * kronecker(self.domain)[2]
	 else:
	    F_fp_X = 0 * kronecker(self.domain)[1]
	    F_fs_X = 0 * kronecker(self.domain)[1]
	    
	 F_mg = B * ( L_g - dL_gdp * S_fg )

	 
	 D=self.__pde.getNewCoefficient("D")
	 A=self.__pde.getNewCoefficient("A")
	 Y=self.__pde.getNewCoefficient("Y")
	 X=self.__pde.getNewCoefficient("X") 
	 
	 dtXI = dt*self.XI
	 dtcXI = dt*(1-self.XI)
	 
	 D[0,0]=E_fpp + dtXI * D_fpp
	 D[0,1]=E_fps
	 D[1,0]=E_fsp
	 D[1,1]=E_fss + dtXI * D_fss
	 D[1,2]=rho_fg
	 D[2,1]= dtXI * B * dL_gdp
	 D[2,2]= 1 + dtXI * B
	 
	 H_fw = dtcXI * A_fw * grad( p_f_old)
	 H_fg = dtcXI * A_fg * grad(S_fg_old)
	 
	 for i in range(self.domain.getDim()):
	    A[0,i,0,i] = dtXI * A_fw * self.perm_f[i]
	    A[1,i,1,i] = dtXI * A_fg * self.perm_f[i]
	    
	    X[0,i]= dt * F_fp_X - H_fw * self.perm_f[i]
	    X[1,i]= dt * F_fs_X - H_fg * self.perm_f[i]
	    
	 Y[0] = E_fpp *  p_f_old + E_fps * S_fg_old +                     dt * F_fp_Y - dtcXI * D_fpp *  p_f_old
	 Y[1] = E_fsp *  p_f_old + E_fss * S_fg_old + rho_fg * C_mg_old + dt * F_fs_Y - dtcXI * D_fss * S_fg_old
	 Y[2] = C_mg_old                                                + dt * F_mg_Y - dtcXI * ( dL_gdp * S_fg_old + C_mg_old)
	 
	 self.__pde.setValue(A=A, D=D, X=X, Y=Y)
	 
	 u2 = self.__pde.getSolution()
	 return u2
