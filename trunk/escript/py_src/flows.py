# -*- coding: utf-8 -*-
########################################################
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

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
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

import escript
import util
from linearPDEs import LinearPDE, LinearPDESystem, LinearSinglePDE, SolverOptions
from pdetools import HomogeneousSaddlePointProblem,Projector, ArithmeticTuple, PCG, NegativeNorm, GMRES

class DarcyFlow(object):
   """
   solves the problem
   
   *u_i+k_{ij}*p_{,j} = g_i*
   *u_{i,i} = f*
   
   where *p* represents the pressure and *u* the Darcy flux. *k* represents the permeability,
   
   :cvar SIMPLE: simple solver
   :cvar POST: solver using global postprocessing of flux
   :cvar STAB: solver uses (non-symmetric) stabilization
   :cvar SYMSTAB: solver uses symmetric stabilization
   """
   SIMPLE="SIMPLE"
   POST="POST"
   STAB="STAB"
   SYMSTAB="SYMSTAB"
   def __init__(self, domain, useReduced=False, solver="SYMSTAB", verbose=False, w=1.):
      """
      initializes the Darcy flux problem
      :param domain: domain of the problem
      :type domain: `Domain`
      :param useReduced: uses reduced oreder on flux and pressure
      :type useReduced: ``bool``
      :param solver: solver method 
      :type solver: in [`DarcyFlow.SIMPLE`, `DarcyFlow.POST', `DarcyFlow.STAB`, `DarcyFlow.SYMSTAB` ] 
      :param verbose: if ``True`` some information on the iteration progress are printed.
      :type verbose: ``bool``
      :param w: weighting factor for `DarcyFlow.POST` solver
      :type w: ``float``
      
      """
      self.domain=domain
      self.solver=solver
      self.useReduced=useReduced
      self.verbose=verbose
      self.scale=1.
      
      
      self.__pde_v=LinearPDESystem(domain)
      self.__pde_v.setSymmetryOn()
      if self.useReduced: self.__pde_v.setReducedOrderOn()
      self.__pde_p=LinearSinglePDE(domain)
      self.__pde_p.setSymmetryOn()
      if self.useReduced: self.__pde_p.setReducedOrderOn()
      
      if self.solver  == self.SIMPLE:
	 if self.verbose: print "DarcyFlow: simple solver is used."
         self.__pde_v.setValue(D=util.kronecker(self.domain.getDim()))
      elif self.solver  == self.POST:
	 self.w=w
	 if util.inf(w)<0.:
	    raise ValueError,"Weighting factor must be non-negative." 
	 if self.verbose: print "DarcyFlow: global postprocessing of flux is used."
      elif self.solver  == self.STAB:
	  if self.verbose: print "DarcyFlow: (non-symmetric) stabilization is used."
      elif  self.solver  == self.SYMSTAB:
	  if self.verbose: print "DarcyFlow: symmetric stabilization is used."
      else:
	raise ValueError,"unknown solver %s"%self.solver
      self.__f=escript.Scalar(0,self.__pde_p.getFunctionSpaceForCoefficient("X"))
      self.__g=escript.Vector(0,self.__pde_v.getFunctionSpaceForCoefficient("Y"))
      self.location_of_fixed_pressure = escript.Scalar(0, self.__pde_p.getFunctionSpaceForCoefficient("q"))
      self.location_of_fixed_flux = escript.Vector(0, self.__pde_v.getFunctionSpaceForCoefficient("q"))
      self.setTolerance()
     
        
   def setValue(self,f=None, g=None, location_of_fixed_pressure=None, location_of_fixed_flux=None, permeability=None):
      """
      assigns values to model parameters

      :param f: volumetic sources/sinks
      :type f: scalar value on the domain (e.g. `escript.Data`)
      :param g: flux sources/sinks
      :type g: vector values on the domain (e.g. `escript.Data`)
      :param location_of_fixed_pressure: mask for locations where pressure is fixed
      :type location_of_fixed_pressure: scalar value on the domain (e.g. `escript.Data`)
      :param location_of_fixed_flux:  mask for locations where flux is fixed.
      :type location_of_fixed_flux: vector values on the domain (e.g. `escript.Data`)
      :param permeability: permeability tensor. If scalar ``s`` is given the tensor with ``s`` on the main diagonal is used. 
      :type permeability: scalar or symmetric tensor values on the domain (e.g. `escript.Data`)

      :note: the values of parameters which are not set by calling ``setValue`` are not altered.
      :note: at any point on the boundary of the domain the pressure
             (``location_of_fixed_pressure`` >0) or the normal component of the
             flux (``location_of_fixed_flux[i]>0``) if direction of the normal
             is along the *x_i* axis.

      """
      if location_of_fixed_pressure!=None: 
           self.location_of_fixed_pressure=util.wherePositive(location_of_fixed_pressure)
           self.__pde_p.setValue(q=self.location_of_fixed_pressure)
      if location_of_fixed_flux!=None: 
          self.location_of_fixed_flux=util.wherePositive(location_of_fixed_flux)
          self.__pde_v.setValue(q=self.location_of_fixed_flux)
      
			
      # pressure  is rescaled by the factor 1/self.scale
      if permeability!=None:
	
	 perm=util.interpolate(permeability,self.__pde_v.getFunctionSpaceForCoefficient("A"))
         V=util.vol(self.domain)
         l=V**(1./self.domain.getDim())
         
	 if perm.getRank()==0:
	    perm_inv=(1./perm)
            self.scale=util.integrate(perm_inv)/V*l
	    perm_inv=perm_inv*((1./self.scale)*util.kronecker(self.domain.getDim()))
	    perm=perm*(self.scale*util.kronecker(self.domain.getDim()))
	    
	    
	 elif perm.getRank()==2:
	    perm_inv=util.inverse(perm)
            self.scale=util.sqrt(util.integrate(util.length(perm_inv)**2)/V)*l
	    perm_inv*=(1./self.scale)
	    perm=perm*self.scale
	 else:
	    raise ValueError,"illegal rank of permeability."
         
	 self.__permeability=perm
	 self.__permeability_inv=perm_inv
	 if self.verbose: print "DarcyFlow: scaling factor for pressure is %e."%self.scale
	 
	 if self.solver  == self.SIMPLE:
	    self.__pde_p.setValue(A=self.__permeability)
	 elif self.solver  == self.POST:
	    self.__pde_p.setValue(A=self.__permeability)
	    k=util.kronecker(self.domain.getDim())
	    self.lamb = self.w*util.length(perm_inv)*l 
	    self.__pde_v.setValue(D=self.__permeability_inv, A=self.lamb*self.domain.getSize()*util.outer(k,k))
	 elif self.solver  == self.STAB:
	    self.__pde_p.setValue(A=0.5*self.__permeability)
	    self.__pde_v.setValue(D=0.5*self.__permeability_inv)
	 elif  self.solver  == self.SYMSTAB:
	    self.__pde_p.setValue(A=0.5*self.__permeability)
	    self.__pde_v.setValue(D=0.5*self.__permeability_inv)

      if g != None:
	g=util.interpolate(g, self.__pde_v.getFunctionSpaceForCoefficient("Y"))
	if g.isEmpty():
	      g=Vector(0,self.__pde_v.getFunctionSpaceForCoefficient("Y"))
	else:
	    if not g.getShape()==(self.domain.getDim(),): raise ValueError,"illegal shape of g"
	self.__g=g
      if f !=None:
	 f=util.interpolate(f, self.__pde_p.getFunctionSpaceForCoefficient("Y"))
	 if f.isEmpty():	   
	      f=Scalar(0,self.__pde_p.getFunctionSpaceForCoefficient("Y"))
	 else:
	     if f.getRank()>0: raise ValueError,"illegal rank of f."
	 self.__f=f
   def getSolverOptionsFlux(self):
      """
      Returns the solver options used to solve the flux problems
      :return: `SolverOptions`
      """
      return self.__pde_v.getSolverOptions()
      
   def setSolverOptionsFlux(self, options=None):
      """
      Sets the solver options used to solve the flux problems
      If ``options`` is not present, the options are reset to default
      :param options: `SolverOptions`
      """
      return self.__pde_v.setSolverOptions(options)
	 
   def getSolverOptionsPressure(self):
      """
      Returns the solver options used to solve the pressure problems
      :return: `SolverOptions`
      """
      return self.__pde_p.getSolverOptions()
      
   def setSolverOptionsPressure(self, options=None):
      """
      Sets the solver options used to solve the pressure problems 
      If ``options`` is not present, the options are reset to default
      
      :param options: `SolverOptions`
      :note: if the adaption of subtolerance is choosen, the tolerance set by ``options`` will be overwritten before the solver is called.
      """
      return self.__pde_p.setSolverOptions(options)
      
   def setTolerance(self,rtol=1e-4):
      """
      sets the relative tolerance ``rtol`` for the pressure for the stabelized solvers.
      
      :param rtol: relative tolerance for the pressure
      :type rtol: non-negative ``float``
      """
      if rtol<0:
	 raise ValueError,"Relative tolerance needs to be non-negative."
      self.__rtol=rtol
      
   def getTolerance(self):
      """
      returns the relative tolerance
      :return: current relative tolerance
      :rtype: ``float``
      """
      return self.__rtol
      
   def solve(self,u0,p0, max_iter=100, iter_restart=20):
      """
      solves the problem.
      
      The iteration is terminated if the residual norm is less then self.getTolerance().

      :param u0: initial guess for the flux. At locations in the domain marked by ``location_of_fixed_flux`` the value of ``u0`` is kept unchanged.
      :type u0: vector value on the domain (e.g. `escript.Data`).
      :param p0: initial guess for the pressure. At locations in the domain marked by ``location_of_fixed_pressure`` the value of ``p0`` is kept unchanged.
      :type p0: scalar value on the domain (e.g. `escript.Data`).
      :param max_iter: maximum number of (outer) iteration steps for the stabilization solvers,
      :type max_iter: ``int``
      :param iter_restart: number of steps after which the iteration is restarted. The larger ``iter_restart`` the larger the required memory. 
                           A small value for ``iter_restart`` may require a large number of iteration steps or may even lead to a failure
                           of the iteration. ``iter_restart`` is relevant for the stabilization solvers only.
      :type iter_restart: ``int``
      :return: flux and pressure
      :rtype: ``tuple`` of `escript.Data`.

      """
      # rescale initial guess:
      p0=p0/self.scale
      if self.solver  == self.SIMPLE or self.solver  == self.POST :
	    self.__pde_p.setValue(X=self.__g , 
	                          Y=self.__f, 
	                          y= - util.inner(self.domain.getNormal(),u0 * self.location_of_fixed_flux), 
	                          r=p0)
	    p=self.__pde_p.getSolution()
	    u = self.getFlux(p, u0)
      elif  self.solver  == self.STAB:
	u,p = self.__solve_STAB(u0,p0, max_iter, iter_restart)
      elif  self.solver  == self.SYMSTAB:
	u,p = self.__solve_SYMSTAB(u0,p0, max_iter, iter_restart)
	
      if self.verbose:
	    KGp=util.tensor_mult(self.__permeability,util.grad(p))
	    def_p=self.__g-(u+KGp)
	    def_v=self.__f-util.div(u, self.__pde_v.getFunctionSpaceForCoefficient("X"))
	    print "DarcyFlux: |g-u-K*grad(p)|_2 = %e (|u|_2 = %e)."%(self.__L2(def_p),self.__L2(u))
	    print "DarcyFlux: |f-div(u)|_2 = %e (|grad(u)|_2 = %e)."%(self.__L2(def_v),self.__L2(util.grad(u)))
      #rescale result
      p=p*self.scale
      return u,p
      
   def getFlux(self,p, u0=None):
        """
        returns the flux for a given pressure ``p`` where the flux is equal to ``u0``
        on locations where ``location_of_fixed_flux`` is positive (see `setValue`).
        Notice that ``g`` and ``f`` are used, see `setValue`.

        :param p: pressure.
        :type p: scalar value on the domain (e.g. `escript.Data`).
        :param u0: flux on the locations of the domain marked be ``location_of_fixed_flux``.
        :type u0: vector values on the domain (e.g. `escript.Data`) or ``None``
        :return: flux
        :rtype: `escript.Data`
        """
        if self.solver  == self.SIMPLE or self.solver  == self.POST  :
            KGp=util.tensor_mult(self.__permeability,util.grad(p))
            self.__pde_v.setValue(Y=self.__g-KGp, X=escript.Data())
            if u0 == None:
	       self.__pde_v.setValue(r=escript.Data())
	    else:
	       self.__pde_v.setValue(r=u0)
            u= self.__pde_v.getSolution()
	elif self.solver  == self.POST:
            self.__pde_v.setValue(Y=util.tensor_mult(self.__permeability_inv,self.__g)-util.grad(p),
                                  X=self.lamb * self.__f * util.kronecker(self.domain.getDim()))
            if u0 == None:
	       self.__pde_v.setValue(r=escript.Data())
	    else:
	       self.__pde_v.setValue(r=u0)
            u= self.__pde_v.getSolution()
	elif self.solver  == self.STAB:
	     gp=util.grad(p)
	     self.__pde_v.setValue(Y=0.5*(util.tensor_mult(self.__permeability_inv,self.__g)+gp),
	                           X= p * util.kronecker(self.domain.getDim()),
	                           y= - p * self.domain.getNormal())                           
	     if u0 == None:
	       self.__pde_v.setValue(r=escript.Data())
	     else:
	       self.__pde_v.setValue(r=u0)
	     u= self.__pde_v.getSolution()
	elif  self.solver  == self.SYMSTAB:
	     gp=util.grad(p)
	     self.__pde_v.setValue(Y=0.5*(util.tensor_mult(self.__permeability_inv,self.__g)-gp),
	                           X= escript.Data() ,
	                           y= escript.Data() )                           
	     if u0 == None:
	       self.__pde_v.setValue(r=escript.Data())
	     else:
	       self.__pde_v.setValue(r=u0)
	     u= self.__pde_v.getSolution()
	return u
	  
     
   def __solve_STAB(self, u0, p0, max_iter, iter_restart):
          # p0 is used as an initial guess
	  u=self.getFlux(p0, u0)  
          self.__pde_p.setValue( Y=self.__f-util.div(u), 
                                 X=0.5*(self.__g - u - util.tensor_mult(self.__permeability,util.grad(p0)) ), 
                                 y= escript.Data(), 
                                 r=escript.Data())

	  dp=self.__pde_p.getSolution()
	  p=GMRES(dp, 
	          self.__STAB_Aprod, 
		  p0, 
		  self.__inner, 
		  atol=self.__norm(p0+dp)*self.getTolerance() ,
		  rtol=0.,
		  iter_max=max_iter,
		  iter_restart=iter_restart, 
		  verbose=self.verbose,P_R=None)
            
          u=self.getFlux(p, u0)
          return u,p

   def __solve_SYMSTAB(self, u0, p0, max_iter, iter_restart):
          # p0 is used as an initial guess
	  u=self.getFlux(p0, u0)  
          self.__pde_p.setValue( Y= self.__f, 
                                 X=  0.5*(self.__g + u - util.tensor_mult(self.__permeability,util.grad(p0)) ), 
                                 y=  -  util.inner(self.domain.getNormal(), u), 
                                 r=escript.Data())
	  dp=self.__pde_p.getSolution()
	  
	  print dp
          print p0+dp
          
	  p=GMRES(dp, 
	          self.__SYMSTAB_Aprod, 
		  p0, 
		  self.__inner, 
		  atol=self.__norm(p0+dp)*self.getTolerance() ,
		  rtol=0.,
		  iter_max=max_iter,
		  iter_restart=iter_restart, 
		  verbose=self.verbose,P_R=None)
            
          u=self.getFlux(p, u0)
          return u,p

   def __L2(self,v):
         return util.sqrt(util.integrate(util.length(util.interpolate(v,escript.Function(self.domain)))**2))      
   
   def __norm(self,r):
         return util.sqrt(self.__inner(r,r))
         
   def __inner(self,r,s):
         return util.integrate(util.inner(r,s), escript.Function(self.domain))
         
   def __STAB_Aprod(self,p):
      gp=util.grad(p)
      self.__pde_v.setValue(Y=-0.5*gp,
                            X=-p*util.kronecker(self.__pde_v.getDomain()), 
                            y= p * self.domain.getNormal(),  
                            r=escript.Data())
      u = -self.__pde_v.getSolution()
      self.__pde_p.setValue(Y=util.div(u), 
                            X=0.5*(u+util.tensor_mult(self.__permeability,gp)),
                            y=escript.Data(), 
                            r=escript.Data())
     
      return  self.__pde_p.getSolution()
   
   def __SYMSTAB_Aprod(self,p):
      gp=util.grad(p)
      self.__pde_v.setValue(Y=0.5*gp ,
                            X=escript.Data(), 
                            y=escript.Data(),  
                            r=escript.Data())
      u = -self.__pde_v.getSolution()
      self.__pde_p.setValue(Y=escript.Data(), 
                            X=0.5*(-u+util.tensor_mult(self.__permeability,gp)),
                            y=escript.Data(), 
                            r=escript.Data())
     
      return  self.__pde_p.getSolution()
      

class StokesProblemCartesian(HomogeneousSaddlePointProblem):
     """
     solves

          -(eta*(u_{i,j}+u_{j,i}))_j + p_i = f_i-stress_{ij,j}
                u_{i,i}=0

          u=0 where  fixed_u_mask>0
          eta*(u_{i,j}+u_{j,i})*n_j-p*n_i=surface_stress +stress_{ij}n_j

     if surface_stress is not given 0 is assumed.

     typical usage:

            sp=StokesProblemCartesian(domain)
            sp.setTolerance()
            sp.initialize(...)
            v,p=sp.solve(v0,p0)
     """
     def __init__(self,domain,**kwargs):
         """
         initialize the Stokes Problem

         The approximation spaces used for velocity (=Solution(domain)) and pressure (=ReducedSolution(domain)) must be
         LBB complient, for instance using quadratic and linear approximation on the same element or using linear approximation
         with macro elements for the pressure. 

         :param domain: domain of the problem.
         :type domain: `Domain`
         """
         HomogeneousSaddlePointProblem.__init__(self,**kwargs)
         self.domain=domain
         self.__pde_v=LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_v.setSymmetryOn()
	 
         self.__pde_prec=LinearPDE(domain)
         self.__pde_prec.setReducedOrderOn()
         self.__pde_prec.setSymmetryOn()

         self.__pde_proj=LinearPDE(domain)
         self.__pde_proj.setReducedOrderOn()
	 self.__pde_proj.setValue(D=1)
         self.__pde_proj.setSymmetryOn()

     def getSolverOptionsVelocity(self):
         """
	 returns the solver options used  solve the equation for velocity.
	 
	 :rtype: `SolverOptions`
	 """
	 return self.__pde_v.getSolverOptions()
     def setSolverOptionsVelocity(self, options=None):
         """
	 set the solver options for solving the equation for velocity.
	 
	 :param options: new solver  options
	 :type options: `SolverOptions`
	 """
         self.__pde_v.setSolverOptions(options)
     def getSolverOptionsPressure(self):
         """
	 returns the solver options used  solve the equation for pressure.
	 :rtype: `SolverOptions`
	 """
	 return self.__pde_prec.getSolverOptions()
     def setSolverOptionsPressure(self, options=None):
         """
	 set the solver options for solving the equation for pressure.
	 :param options: new solver  options
	 :type options: `SolverOptions`
	 """
	 self.__pde_prec.setSolverOptions(options)

     def setSolverOptionsDiv(self, options=None):
         """
	 set the solver options for solving the equation to project the divergence of
	 the velocity onto the function space of presure.
	 
	 :param options: new solver options
	 :type options: `SolverOptions`
	 """
	 self.__pde_proj.setSolverOptions(options)
     def getSolverOptionsDiv(self):
         """
	 returns the solver options for solving the equation to project the divergence of
	 the velocity onto the function space of presure.
	 
	 :rtype: `SolverOptions`
	 """
	 return self.__pde_proj.getSolverOptions()

     def updateStokesEquation(self, v, p):
         """
         updates the Stokes equation to consider dependencies from ``v`` and ``p``
         :note: This method can be overwritten by a subclass. Use `setStokesEquation` to set new values.
         """
         pass
     def setStokesEquation(self, f=None,fixed_u_mask=None,eta=None,surface_stress=None,stress=None, restoration_factor=None):
        """
        assigns new values to the model parameters. 

        :param f: external force
        :type f: `Vector` object in `FunctionSpace` `Function` or similar
        :param fixed_u_mask: mask of locations with fixed velocity.
        :type fixed_u_mask: `Vector` object on `FunctionSpace` `Solution` or similar
        :param eta: viscosity
        :type eta: `Scalar` object on `FunctionSpace` `Function` or similar
        :param surface_stress: normal surface stress
        :type surface_stress: `Vector` object on `FunctionSpace` `FunctionOnBoundary` or similar
        :param stress: initial stress
	:type stress: `Tensor` object on `FunctionSpace` `Function` or similar
        """
        if eta !=None:
            k=util.kronecker(self.domain.getDim())
            kk=util.outer(k,k)
            self.eta=util.interpolate(eta, escript.Function(self.domain))
	    self.__pde_prec.setValue(D=1/self.eta)
            self.__pde_v.setValue(A=self.eta*(util.swap_axes(kk,0,3)+util.swap_axes(kk,1,3)))
        if restoration_factor!=None:
            n=self.domain.getNormal()
            self.__pde_v.setValue(d=restoration_factor*util.outer(n,n))
        if fixed_u_mask!=None:
            self.__pde_v.setValue(q=fixed_u_mask)
        if f!=None: self.__f=f
        if surface_stress!=None: self.__surface_stress=surface_stress
        if stress!=None: self.__stress=stress

     def initialize(self,f=escript.Data(),fixed_u_mask=escript.Data(),eta=1,surface_stress=escript.Data(),stress=escript.Data(), restoration_factor=0):
        """
        assigns values to the model parameters

        :param f: external force
        :type f: `Vector` object in `FunctionSpace` `Function` or similar
        :param fixed_u_mask: mask of locations with fixed velocity.
        :type fixed_u_mask: `Vector` object on `FunctionSpace` `Solution` or similar
        :param eta: viscosity
        :type eta: `Scalar` object on `FunctionSpace` `Function` or similar
        :param surface_stress: normal surface stress
        :type surface_stress: `Vector` object on `FunctionSpace` `FunctionOnBoundary` or similar
        :param stress: initial stress
	:type stress: `Tensor` object on `FunctionSpace` `Function` or similar
        """
        self.setStokesEquation(f,fixed_u_mask, eta, surface_stress, stress, restoration_factor)

     def Bv(self,v,tol):
         """
         returns inner product of element p and div(v)

         :param v: a residual
         :return: inner product of element p and div(v)
         :rtype: ``float``
         """
         self.__pde_proj.setValue(Y=-util.div(v)) 
	 self.getSolverOptionsDiv().setTolerance(tol)
	 self.getSolverOptionsDiv().setAbsoluteTolerance(0.)
         out=self.__pde_proj.getSolution()
         return out

     def inner_pBv(self,p,Bv):
         """
         returns inner product of element p and Bv=-div(v)

         :param p: a pressure increment
         :param Bv: a residual
         :return: inner product of element p and Bv=-div(v)
         :rtype: ``float``
         """
         return util.integrate(util.interpolate(p,escript.Function(self.domain))*util.interpolate(Bv, escript.Function(self.domain)))

     def inner_p(self,p0,p1):
         """
         Returns inner product of p0 and p1

         :param p0: a pressure
         :param p1: a pressure
         :return: inner product of p0 and p1
         :rtype: ``float``
         """
         s0=util.interpolate(p0, escript.Function(self.domain))
         s1=util.interpolate(p1, escript.Function(self.domain))
         return util.integrate(s0*s1)

     def norm_v(self,v):
         """
         returns the norm of v

         :param v: a velovity
         :return: norm of v
         :rtype: non-negative ``float``
         """
         return util.sqrt(util.integrate(util.length(util.grad(v))**2))


     def getDV(self, p, v, tol):
         """
         return the value for v for a given p (overwrite)

         :param p: a pressure
         :param v: a initial guess for the value v to return.
         :return: dv given as *Adv=(f-Av-B^*p)*
         """
         self.updateStokesEquation(v,p)
         self.__pde_v.setValue(Y=self.__f, y=self.__surface_stress)
	 self.getSolverOptionsVelocity().setTolerance(tol)
	 self.getSolverOptionsVelocity().setAbsoluteTolerance(0.)
         if self.__stress.isEmpty():
            self.__pde_v.setValue(X=p*util.kronecker(self.domain)-2*self.eta*util.symmetric(util.grad(v)))
         else:
            self.__pde_v.setValue(X=self.__stress+p*util.kronecker(self.domain)-2*self.eta*util.symmetric(util.grad(v)))
         out=self.__pde_v.getSolution()
         return  out

     def norm_Bv(self,Bv):
        """
        Returns Bv (overwrite).

        :rtype: equal to the type of p
        :note: boundary conditions on p should be zero!
        """
        return util.sqrt(util.integrate(util.interpolate(Bv, escript.Function(self.domain))**2))

     def solve_AinvBt(self,p, tol):
         """
         Solves *Av=B^*p* with accuracy `tol`

         :param p: a pressure increment
         :return: the solution of *Av=B^*p*
         :note: boundary conditions on v should be zero!
         """
         self.__pde_v.setValue(Y=escript.Data(), y=escript.Data(), X=-p*util.kronecker(self.domain))
         out=self.__pde_v.getSolution()
         return  out

     def solve_prec(self,Bv, tol):
         """
         applies preconditioner for for *BA^{-1}B^** to *Bv*
         with accuracy `self.getSubProblemTolerance()` 

         :param Bv: velocity increment
         :return: *p=P(Bv)* where *P^{-1}* is an approximation of *BA^{-1}B^ * )*
         :note: boundary conditions on p are zero.
         """
         self.__pde_prec.setValue(Y=Bv)
	 self.getSolverOptionsPressure().setTolerance(tol)
	 self.getSolverOptionsPressure().setAbsoluteTolerance(0.)
         out=self.__pde_prec.getSolution()
         return out
