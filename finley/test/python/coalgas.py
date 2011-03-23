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

from esys.finley import Rectangle, Brick
from esys.escript import unitsSI as U

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
      
class InterpolationTable(object):
   """
   a simple 1D interpolation table for escript Data object with non-equidistant nodes
   """
   def __init__(self,x=[], y=[], obeyBounds=True):
      """
      set up interpolation table. obeyBounds is set an exception is thrown if
      the interpolation argument is below min(x) or larger than max(x). Otherwise
      the value for x is set to y[0] for 
      """
      if len(x) != len(y):
	 raise ValueError,"length of interpolation nodes and value lists need to be identical."
      if len(x) > 0:
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

class Well(object):
   """
   generic well
   
   :var WATER: phase identifier for water well
   :var GAS: phase identifier for gas well
   """
   WATER="Water"
   GAS="Gas"
   def __init__(self,  *args, **kwargs):
       """
       set up well 
       """
       pass
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
      return self.WATER

class PeacemanWell(Well):
   """
   defines a well using Peaceman's formula
   """
   pass

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
      :param S_fg: gas saturation in the fractured rock as function of the capillary pressure p_fc=p_fg-p_fw
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

      def __init__(self, domain, phi_f=None, phi=None, L_g=None, S_fg=None, 
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
	 :param S_fg: gas saturation in the fractured rock as function of the capillary pressure p_fc=p_fg-p_fw
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
			      S_fg=S_fg, S_mg=None, 
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
	 
      def solvePDE(self):
	 
	 p_fw=self.u[0]
	 p_fg=self.u[1]
	 p_fc=p_fg-p_fw
	 C_mg=self.u[3] 
      
      
	 p_fw_old=self.u_old[0]
	 p_fg_old=self.u_old[1]
	 C_mg_old=self.u_old[3]


	 phi_f   =self.phi_f.getValue(p_fg)
	 dphi_fdp=self.phi_f.getValueDifferential(p_fg)
	 
	 S_fg=	self.S_fg.getValue(p_fc)
	 dS_fgdp=	self.S_fg.getValueDifferential(p_fc)
	 S_fw=	1-S_fg
	 
	 rho_fw 	= self.rho_w.getValue(p_fw)
	 drho_fwdp	= self.rho_w.getValueDifferential(p_fw)
	 rho_fg = self.rho_g.getValue(p_fg)
	 drho_fgdp	= self.rho_g.getValueDifferential(p_fg)
	 
	 L_g       = self.getValue(p_fg)
	 dL_gdp_fg =  self.rho_w.getValueDifferential(p_fg)
	 
	 A_fw = rho_fw * k_fw/mu_fw 
	 A_fg = rho_fg * k_fg/mu_fg
	 
	 m = whereNegative(L_g - C_mg)  
	 B = self.sigma * self.a_mg * (m + (1-m) * self.f_rg * S_fg )
	 
	 
	 E_fww= phi_f * ( rho_fw * dS_fgdp + S_fw * drho_fwdp )
	 E_fwg= rho_fw * ( S_fw * dphi_fdp - phi_f * dS_fgdp )
	 E_fgw= - phi_f * rho_fg * dS_fgdp
	 E_fgg= S_fg * rho_fg  * dphi_fdp + phi_f * rho_fg * dS_fgdp + phi_f  * S_fg * drho_fgdp


	 
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
	       
	 F_fw_Y = A_fw * Q_w 
	 F_fg_Y = A_fg * Q_g
	 D_fgg = - A_fg * R_g 
	 D_fww = - A_fw * R_w
	 
	 if self.domain.getDim() == 3:
	    F_fw_X = A_fw * self.g * rho_fw * self.perm_f[2] * kronecker(self.domain)[2]
	    F_fg_X = A_fg * self.g * rho_fg * self.perm_f[2] * kronecker(self.domain)[2]
	 else:
	    F_fw_X = 0 * kronecker(self.domain)[1]
	    F_fg_X = 0 * kronecker(self.domain)[1]
	    
	 F_mg = B * ( L_g - dL_gdp_fg * p_fg )

	 
	 D=self.__pde.getNewCoefficient("D")
	 A=self.__pde.getNewCoefficient("A")
	 Y=self.__pde.getNewCoefficient("Y")
	 X=self.__pde.getNewCoefficient("X") 
	 
	 dtXI = dt*self.XI
	 dtcXI = dt*(1-self.XI)
	 
	 D[0,0]=E_fww + dtXI * D_fww
	 D[0,1]=E_fwg
	 D[1,0]=E_fgw
	 D[1,1]=E_fgg + dtXI * D_fgg
	 D[1,2]=rho_fg
	 D[2,1]= dtXI * B * dL_gdp_fg
	 D[2,2]= 1 + dtXI * B
	 
	 H_fw = dtcXI * A_fw * grad(p_fw_old)
	 H_fg = dtcXI * A_fg * grad(p_fg_old)
	 
	 for i in range(self.domain.getDim()):
	    A[0,i,0,i] = dtXI * A_fw * self.perm_f[i]
	    A[1,i,1,i] = dtXI * A_fg * self.perm_f[i]
	    
	    X[0,i]= dt * F_fw_X - H_fw * self.perm_f[i]
	    X[1,i]= dt * F_fg_X - H_fg * self.perm_f[i]
	    
	 Y[0] = E_fww * p_fw_old + E_fwg * p_fg_old +                     dt * F_fw_Y - dtcXI * D_fww * p_fw_old
	 Y[1] = E_fgw * p_fw_old + E_fgg * p_fg_old + rho_fg * C_mg_old + dt * F_fg_Y - dtcXI * D_fgg * p_fg_old
	 Y[2] = C_mg_old                                                + dt * F_mg_Y - dtcXI * ( dL_gdp_fg * p_fg_old + C_mg_old)
	 
	 self.__pde.setValue(A=A, D=D, X=X, Y=Y)
	 
	 u2 = self.__pde.getSolution()
	 return u2
