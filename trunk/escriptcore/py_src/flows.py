# -*- coding: utf-8 -*-
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
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
from . import linearPDEs as lpe
from . import pdetools as pdt

class DarcyFlow(object):
   """
   solves the problem
   
   *u_i+k_{ij}*p_{,j} = g_i*
   *u_{i,i} = f*
   
   where *p* represents the pressure and *u* the Darcy flux. *k* represents the permeability,
   
   :cvar EVAL: direct pressure gradient evaluation for flux 
   :cvar POST: global postprocessing of flux by solving the PDE *K_{ij} u_j + (w * K * l u_{k,k})_{,i}= - p_{,j} + K_{ij} g_j*
               where *l* is the length scale, *K* is the inverse of the permeability tensor, and *w* is a positive weighting factor.
   :cvar SMOOTH: global smoothing by solving the PDE *K_{ij} u_j= - p_{,j} + K_{ij} g_j*
   """
   EVAL="EVAL"
   SIMPLE="EVAL"
   POST="POST"
   SMOOTH="SMOOTH"
   def __init__(self, domain, useReduced=False, solver="POST", verbose=False, w=1.):
      """
      initializes the Darcy flux problem.

      :param domain: domain of the problem
      :type domain: `Domain`
      :param useReduced: uses reduced oreder on flux and pressure
      :type useReduced: ``bool``
      :param solver: solver method 
      :type solver: in [`DarcyFlow.EVAL`, `DarcyFlow.POST`, `DarcyFlow.SMOOTH` ]
      :param verbose: if ``True`` some information on the iteration progress are printed.
      :type verbose: ``bool``
      :param w: weighting factor for `DarcyFlow.POST` solver
      :type w: ``float``
      
      """
      if not solver in [DarcyFlow.EVAL, DarcyFlow.POST,  DarcyFlow.SMOOTH ] :
          raise ValueError("unknown solver %d."%solver)

      self.domain=domain
      self.solver=solver
      self.useReduced=useReduced
      self.verbose=verbose
      self.l=None
      self.w=None
     
      self.__pde_p=lpe.LinearSinglePDE(domain)
      self.__pde_p.setSymmetryOn()
      if self.useReduced: self.__pde_p.setReducedOrderOn()

      if self.solver  == self.EVAL:
         self.__pde_v=None
         if self.verbose: print("DarcyFlow: simple solver is used.")

      elif self.solver  == self.POST:
         if util.inf(w)<0.:
            raise ValueError("Weighting factor must be non-negative.") 
         if self.verbose: print("DarcyFlow: global postprocessing of flux is used.")
         self.__pde_v=lpe.LinearPDESystem(domain)
         self.__pde_v.setSymmetryOn()
         if self.useReduced: self.__pde_v.setReducedOrderOn()
         self.w=w
         x=self.domain.getX()
         self.l=min( [util.sup(x[i])-util.inf(x[i]) for i in range(self.domain.getDim()) ] )
         #self.l=util.vol(self.domain)**(1./self.domain.getDim()) # length scale

      elif self.solver  == self.SMOOTH:
         self.__pde_v=lpe.LinearPDESystem(domain)
         self.__pde_v.setSymmetryOn()
         if self.useReduced: self.__pde_v.setReducedOrderOn()
         if self.verbose: print("DarcyFlow: flux smoothing is used.")
         self.w=0

      self.__f=escore.Data(0,self.__pde_p.getFunctionSpaceForCoefficient("X"))
      self.__g=escore.Vector(0,self.__pde_p.getFunctionSpaceForCoefficient("Y"))
      self.__permeability_invXg=escore.Vector(0,self.__pde_p.getFunctionSpaceForCoefficient("Y"))
      self.__permeability_invXg_ref=util.numpy.zeros((self.domain.getDim()),util.numpy.float64)
      self.ref_point_id=None
      self.ref_point=util.numpy.zeros((self.domain.getDim()),util.numpy.float64)
      self.location_of_fixed_pressure = escore.Data(0, self.__pde_p.getFunctionSpaceForCoefficient("q"))
      self.location_of_fixed_flux = escore.Vector(0, self.__pde_p.getFunctionSpaceForCoefficient("q"))
      self.perm_scale=1.
     
        
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
      if location_of_fixed_pressure is not None: 
           self.location_of_fixed_pressure=util.wherePositive(util.interpolate(location_of_fixed_pressure, self.__pde_p.getFunctionSpaceForCoefficient("q")))
           self.ref_point_id=self.location_of_fixed_pressure.internal_maxGlobalDataPoint()
           if not self.location_of_fixed_pressure.getTupleForGlobalDataPoint(*self.ref_point_id)[0] > 0: raise ValueError("pressure needs to be fixed at least one point.")
           self.ref_point=self.__pde_p.getFunctionSpaceForCoefficient("q").getX().getTupleForGlobalDataPoint(*self.ref_point_id)
           if self.verbose: print(("DarcyFlow: reference point at %s."%(self.ref_point,)))
           self.__pde_p.setValue(q=self.location_of_fixed_pressure)
      if location_of_fixed_flux is not None: 
          self.location_of_fixed_flux=util.wherePositive(location_of_fixed_flux)
          if not self.__pde_v is None: 
              self.__pde_v.setValue(q=self.location_of_fixed_flux)
      
      if permeability is not None:
         perm=util.interpolate(permeability,self.__pde_p.getFunctionSpaceForCoefficient("A"))
         self.perm_scale=util.Lsup(util.length(perm))
         if self.verbose: print(("DarcyFlow: permeability scaling factor = %e."%self.perm_scale))
         perm=perm*(1./self.perm_scale)
         
         if perm.getRank()==0:

            perm_inv=(1./perm)
            perm_inv=perm_inv*util.kronecker(self.domain.getDim())
            perm=perm*util.kronecker(self.domain.getDim())
         elif perm.getRank()==2:
            perm_inv=util.inverse(perm)
         else:
            raise ValueError("illegal rank of permeability.")
         
         self.__permeability=perm
         self.__permeability_inv=perm_inv
 
         #====================
         self.__pde_p.setValue(A=self.__permeability)
         if self.solver  == self.EVAL:
              pass # no extra work required
         elif self.solver  == self.POST:
              k=util.kronecker(self.domain.getDim())
              self.omega = self.w*util.length(perm_inv)*self.l*self.domain.getSize()
              #self.__pde_v.setValue(D=self.__permeability_inv, A=self.omega*util.outer(k,k))
              self.__pde_v.setValue(D=self.__permeability_inv, A_reduced=self.omega*util.outer(k,k))
         elif self.solver  == self.SMOOTH:
            self.__pde_v.setValue(D=self.__permeability_inv)

      if g is not None:
        g=util.interpolate(g, self.__pde_p.getFunctionSpaceForCoefficient("Y"))
        if g.isEmpty():
             g=Vector(0,self.__pde_p.getFunctionSpaceForCoefficient("Y"))
        else:
             if not g.getShape()==(self.domain.getDim(),): raise ValueError("illegal shape of g")
        self.__g=g 
        self.__permeability_invXg=util.tensor_mult(self.__permeability_inv,self.__g * (1./self.perm_scale )) 
        self.__permeability_invXg_ref=util.integrate(self.__permeability_invXg)/util.vol(self.domain) 
      if f  is not None:
         f=util.interpolate(f, self.__pde_p.getFunctionSpaceForCoefficient("Y"))
         if f.isEmpty():   
             f=Scalar(0,self.__pde_p.getFunctionSpaceForCoefficient("Y"))
         else:
             if f.getRank()>0: raise ValueError("illegal rank of f.")
         self.__f=f

   def getSolverOptionsFlux(self):
      """
      Returns the solver options used to solve the flux problems
      :return: `SolverOptions`
      """
      if self.__pde_v is None:
          return None
      else:
          return self.__pde_v.getSolverOptions()
      
   def setSolverOptionsFlux(self, options=None):
      """
      Sets the solver options used to solve the flux problems
      If ``options`` is not present, the options are reset to default
      :param options: `SolverOptions`
      """
      if not self.__pde_v is None:
          self.__pde_v.setSolverOptions(options)
 
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
      
   def solve(self, u0, p0):
      """
      solves the problem.
      
      :param u0: initial guess for the flux. At locations in the domain marked by ``location_of_fixed_flux`` the value of ``u0`` is kept unchanged.
      :type u0: vector value on the domain (e.g. `escript.Data`).
      :param p0: initial guess for the pressure. At locations in the domain marked by ``location_of_fixed_pressure`` the value of ``p0`` is kept unchanged.
      :type p0: scalar value on the domain (e.g. `escript.Data`).
      :return: flux and pressure
      :rtype: ``tuple`` of `escript.Data`.

      """
      p0=util.interpolate(p0, self.__pde_p.getFunctionSpaceForCoefficient("q"))
      if self.ref_point_id is None:
          p_ref=0
      else:
          p_ref=p0.getTupleForGlobalDataPoint(*self.ref_point_id)[0]
      p0_hydrostatic=p_ref+util.inner(self.__permeability_invXg_ref, self.__pde_p.getFunctionSpaceForCoefficient("q").getX() - self.ref_point)
      g_2=self.__g - util.tensor_mult(self.__permeability, self.__permeability_invXg_ref * self.perm_scale)
      self.__pde_p.setValue(X=g_2 * 1./self.perm_scale, 
                            Y=self.__f * 1./self.perm_scale,
                            y= - util.inner(self.domain.getNormal(),u0 * self.location_of_fixed_flux * 1./self.perm_scale ), 
                            r=p0 - p0_hydrostatic)
      pp=self.__pde_p.getSolution()
      u = self._getFlux(pp, u0)
      return u,pp + p0_hydrostatic
      
   def getFlux(self,p, u0=None):
        """
        returns the flux for a given pressure ``p`` where the flux is equal to ``u0``
        on locations where ``location_of_fixed_flux`` is positive (see `setValue`).
        Notice that ``g`` is used, see `setValue`.

        :param p: pressure.
        :type p: scalar value on the domain (e.g. `escript.Data`).
        :param u0: flux on the locations of the domain marked be ``location_of_fixed_flux``.
        :type u0: vector values on the domain (e.g. `escript.Data`) or ``None``
        :return: flux
        :rtype: `escript.Data`
        """
        p=util.interpolate(p, self.__pde_p.getFunctionSpaceForCoefficient("q"))
        if self.ref_point_id is None:
            p_ref=0
        else:
            p_ref=p.getTupleForGlobalDataPoint(*self.ref_point_id)[0]
        p_hydrostatic=p_ref+util.inner(self.__permeability_invXg_ref, self.__pde_p.getFunctionSpaceForCoefficient("q").getX() - self.ref_point)
        return self._getFlux(p-p_hydrostatic, u0)

   def _getFlux(self, pp, u0=None):
        """
        returns the flux for a given pressure ``pp`` where the flux is equal to
        ``u0`` on locations where ``location_of_fixed_flux`` is positive (see
        `setValue`). Notice that ``g`` is used, see `setValue`.

        :param pp: pressure.
        :type pp: scalar value on the domain (i.e. `escript.Data`).
        :param u0: flux on the locations of the domain marked in ``location_of_fixed_flux``.
        :type u0: vector values on the domain (i.e. `escript.Data`) or ``None``
        :return: flux
        :rtype: `escript.Data`
        """
        if self.solver  == self.EVAL:
           u = self.__g - util.tensor_mult(self.__permeability, self.perm_scale * (util.grad(pp) + self.__permeability_invXg_ref))
        elif self.solver  == self.POST or self.solver  == self.SMOOTH:
            self.__pde_v.setValue(Y= self.__permeability_invXg - (util.grad(pp) + self.__permeability_invXg_ref))

            if u0 is None:
               self.__pde_v.setValue(r=escore.Data())
            else:
               if not isinstance(u0, escore.Data) : u0 = escore.Vector(u0, escore.Solution(self.domain))
               self.__pde_v.setValue(r=1./self.perm_scale * u0)
            u= self.__pde_v.getSolution() * self.perm_scale
        return u
  
class StokesProblemCartesian(pdt.HomogeneousSaddlePointProblem):
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
            sp.setStokesEquation(...) # new values for some parameters
            v1,p1=sp.solve(v,p)
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
         pdt.HomogeneousSaddlePointProblem.__init__(self,**kwargs)
         self.domain=domain
         self.__pde_v=lpe.LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_v.setSymmetryOn()
 
         self.__pde_prec=lpe.LinearPDE(domain)
         self.__pde_prec.setReducedOrderOn()
         self.__pde_prec.setSymmetryOn()

         self.__pde_proj=lpe.LinearPDE(domain)
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
         :note: This method can be overwritten by a subclass. Use `setStokesEquation` to set new values to model parameters.
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
        if eta  is not None:
            k=util.kronecker(self.domain.getDim())
            kk=util.outer(k,k)
            self.eta=util.interpolate(eta, escore.Function(self.domain))
            self.__pde_prec.setValue(D=1/self.eta)
            self.__pde_v.setValue(A=self.eta*(util.swap_axes(kk,0,3)+util.swap_axes(kk,1,3)))
        if restoration_factor is not None:
            n=self.domain.getNormal()
            self.__pde_v.setValue(d=restoration_factor*util.outer(n,n))
        if fixed_u_mask is not None:
            self.__pde_v.setValue(q=fixed_u_mask)
        if f is not None: self.__f=f
        if surface_stress is not None: self.__surface_stress=surface_stress
        if stress is not None: self.__stress=stress

     def initialize(self,f=escore.Data(),fixed_u_mask=escore.Data(),eta=1,surface_stress=escore.Data(),stress=escore.Data(), restoration_factor=0):
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
         return util.integrate(util.interpolate(p,escore.Function(self.domain))*util.interpolate(Bv, escore.Function(self.domain)))

     def inner_p(self,p0,p1):
         """
         Returns inner product of p0 and p1

         :param p0: a pressure
         :param p1: a pressure
         :return: inner product of p0 and p1
         :rtype: ``float``
         """
         s0=util.interpolate(p0, escore.Function(self.domain))
         s1=util.interpolate(p1, escore.Function(self.domain))
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
         return the value for v for a given p 

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
        return util.sqrt(util.integrate(util.interpolate(Bv, escore.Function(self.domain))**2))

     def solve_AinvBt(self,p, tol):
         """
         Solves *Av=B^*p* with accuracy `tol`

         :param p: a pressure increment
         :return: the solution of *Av=B^*p*
         :note: boundary conditions on v should be zero!
         """
         self.__pde_v.setValue(Y=escore.Data(), y=escore.Data(), X=-p*util.kronecker(self.domain))
         out=self.__pde_v.getSolution()
         return  out

     def solve_prec(self,Bv, tol):
         """
         applies preconditioner for for *BA^{-1}B^** to *Bv*
         with accuracy ``self.getSubProblemTolerance()``

         :param Bv: velocity increment
         :return: *p=P(Bv)* where *P^{-1}* is an approximation of *BA^{-1}B^ * )*
         :note: boundary conditions on p are zero.
         """
         self.__pde_prec.setValue(Y=Bv)
         self.getSolverOptionsPressure().setTolerance(tol)
         self.getSolverOptionsPressure().setAbsoluteTolerance(0.)
         out=self.__pde_prec.getSolution()
         return out
