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

from escript import *
import util
from linearPDEs import LinearPDE, LinearPDESystem, LinearSinglePDE, SolverOptions
from pdetools import HomogeneousSaddlePointProblem,Projector, ArithmeticTuple, PCG, NegativeNorm, GMRES

class DarcyFlow(object):
    """
    solves the problem

    *u_i+k_{ij}*p_{,j} = g_i*
    *u_{i,i} = f*

    where *p* represents the pressure and *u* the Darcy flux. *k* represents the permeability,

    :note: The problem is solved in a least squares formulation.
    """

    def __init__(self, domain, weight=None, useReduced=False, adaptSubTolerance=True):
        """
        initializes the Darcy flux problem
        :param domain: domain of the problem
        :type domain: `Domain`
	:param useReduced: uses reduced oreder on flux and pressure
	:type useReduced: ``bool``
	:param adaptSubTolerance: switches on automatic subtolerance selection
	:type adaptSubTolerance: ``bool``	
        """
        self.domain=domain
        if weight == None:
           s=self.domain.getSize()
           self.__l=(3.*util.longestEdge(self.domain)*s/util.sup(s))**2
           # self.__l=(3.*util.longestEdge(self.domain))**2
           #self.__l=(0.1*util.longestEdge(self.domain)*s/util.sup(s))**2
        else:
           self.__l=weight
        self.__pde_v=LinearPDESystem(domain)
        if useReduced: self.__pde_v.setReducedOrderOn()
        self.__pde_v.setSymmetryOn()
        self.__pde_v.setValue(D=util.kronecker(domain), A=self.__l*util.outer(util.kronecker(domain),util.kronecker(domain)))
        self.__pde_p=LinearSinglePDE(domain)
        self.__pde_p.setSymmetryOn()
        if useReduced: self.__pde_p.setReducedOrderOn()
        self.__f=Scalar(0,self.__pde_v.getFunctionSpaceForCoefficient("X"))
        self.__g=Vector(0,self.__pde_v.getFunctionSpaceForCoefficient("Y"))
        self.setTolerance()
        self.setAbsoluteTolerance()
	self.__adaptSubTolerance=adaptSubTolerance
	self.verbose=False
    def getSolverOptionsFlux(self):
	"""
	Returns the solver options used to solve the flux problems
	
	*(I+D^*D)u=F*
	
	:return: `SolverOptions`
	"""
	return self.__pde_v.getSolverOptions()
    def setSolverOptionsFlux(self, options=None):
	"""
	Sets the solver options used to solve the flux problems
	
	*(I+D^*D)u=F*
	
	If ``options`` is not present, the options are reset to default
	:param options: `SolverOptions`
	:note: if the adaption of subtolerance is choosen, the tolerance set by ``options`` will be overwritten before the solver is called.
	"""
	return self.__pde_v.setSolverOptions(options)
    def getSolverOptionsPressure(self):
	"""
	Returns the solver options used to solve the pressure problems
	
	*(Q^*Q)p=Q^*G*
	
	:return: `SolverOptions`
	"""
	return self.__pde_p.getSolverOptions()
    def setSolverOptionsPressure(self, options=None):
	"""
	Sets the solver options used to solve the pressure problems
	
	*(Q^*Q)p=Q^*G*
	
	If ``options`` is not present, the options are reset to default
	:param options: `SolverOptions`
	:note: if the adaption of subtolerance is choosen, the tolerance set by ``options`` will be overwritten before the solver is called.
	"""
	return self.__pde_p.setSolverOptions(options)

    def setValue(self,f=None, g=None, location_of_fixed_pressure=None, location_of_fixed_flux=None, permeability=None):
        """
        assigns values to model parameters

        :param f: volumetic sources/sinks
        :type f: scalar value on the domain (e.g. `Data`)
        :param g: flux sources/sinks
        :type g: vector values on the domain (e.g. `Data`)
        :param location_of_fixed_pressure: mask for locations where pressure is fixed
        :type location_of_fixed_pressure: scalar value on the domain (e.g. `Data`)
        :param location_of_fixed_flux:  mask for locations where flux is fixed.
        :type location_of_fixed_flux: vector values on the domain (e.g. `Data`)
        :param permeability: permeability tensor. If scalar ``s`` is given the tensor with
                             ``s`` on the main diagonal is used. If vector ``v`` is given the tensor with
                             ``v`` on the main diagonal is used.
        :type permeability: scalar, vector or tensor values on the domain (e.g. `Data`)

        :note: the values of parameters which are not set by calling ``setValue`` are not altered.
        :note: at any point on the boundary of the domain the pressure (``location_of_fixed_pressure`` >0)
               or the normal component of the flux (``location_of_fixed_flux[i]>0`` if direction of the normal
               is along the *x_i* axis.
        """
        if f !=None:
           f=util.interpolate(f, self.__pde_v.getFunctionSpaceForCoefficient("X"))
           if f.isEmpty():
               f=Scalar(0,self.__pde_v.getFunctionSpaceForCoefficient("X"))
           else:
               if f.getRank()>0: raise ValueError,"illegal rank of f."
           self.__f=f
        if g !=None:
           g=util.interpolate(g, self.__pde_p.getFunctionSpaceForCoefficient("Y"))
           if g.isEmpty():
             g=Vector(0,self.__pde_v.getFunctionSpaceForCoefficient("Y"))
           else:
             if not g.getShape()==(self.domain.getDim(),):
               raise ValueError,"illegal shape of g"
           self.__g=g

        if location_of_fixed_pressure!=None: self.__pde_p.setValue(q=location_of_fixed_pressure)
        if location_of_fixed_flux!=None: self.__pde_v.setValue(q=location_of_fixed_flux)

        if permeability!=None:
           perm=util.interpolate(permeability,self.__pde_p.getFunctionSpaceForCoefficient("A"))
           if perm.getRank()==0:
               perm=perm*util.kronecker(self.domain.getDim())
           elif perm.getRank()==1:
               perm, perm2=Tensor(0.,self.__pde_p.getFunctionSpaceForCoefficient("A")), perm
               for i in range(self.domain.getDim()): perm[i,i]=perm2[i]
           elif perm.getRank()==2:
              pass
           else:
              raise ValueError,"illegal rank of permeability."
           self.__permeability=perm
           self.__pde_p.setValue(A=util.transposed_tensor_mult(self.__permeability,self.__permeability))

    def setTolerance(self,rtol=1e-4):
        """
        sets the relative tolerance ``rtol`` used to terminate the solution process. The iteration is terminated if

        *|g-v-Qp| <= atol + rtol * min( max( |g-v|, |Qp| ), max( |v|, |g-Qp| ) )*

        where ``atol`` is an absolut tolerance (see `setAbsoluteTolerance`), *|f|^2 = integrate(length(f)^2)* and *(Qp)_i=k_{ij}p_{,j}* for the permeability *k_{ij}*.

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

    def setAbsoluteTolerance(self,atol=0.):
        """
        sets the absolute tolerance ``atol`` used to terminate the solution process. The iteration is terminated if

        *|g-v-Qp| <= atol + rtol * min( max( |g-v|, |Qp| ), max( |v|, |g-Qp| ) )*

        where ``rtol`` is an absolut tolerance (see `setTolerance`), *|f|^2 = integrate(length(f)^2)* and *(Qp)_i=k_{ij}p_{,j}* for the permeability *k_{ij}*.

        :param atol: absolute tolerance for the pressure
        :type atol: non-negative ``float``
        """
        if atol<0:
            raise ValueError,"Absolute tolerance needs to be non-negative."
        self.__atol=atol
    def getAbsoluteTolerance(self):
       """
       returns the absolute tolerance 
       
       :return: current absolute tolerance
       :rtype: ``float``
       """
       return self.__atol
    def getSubProblemTolerance(self):
	"""
	Returns a suitable subtolerance
	@type: ``float``
	"""
	return max(util.EPSILON**(0.75),self.getTolerance()**2)
    def setSubProblemTolerance(self):
         """
         Sets the relative tolerance to solve the subproblem(s) if subtolerance adaption is selected.
         """
	 if self.__adaptSubTolerance:
		 sub_tol=self.getSubProblemTolerance()
	         self.getSolverOptionsFlux().setTolerance(sub_tol)
		 self.getSolverOptionsFlux().setAbsoluteTolerance(0.)
		 self.getSolverOptionsPressure().setTolerance(sub_tol)
		 self.getSolverOptionsPressure().setAbsoluteTolerance(0.)
		 if self.verbose: print "DarcyFlux: relative subtolerance is set to %e."%sub_tol

    def solve(self,u0,p0, max_iter=100, verbose=False, max_num_corrections=10):
         """
         solves the problem.

         The iteration is terminated if the residual norm is less then self.getTolerance().

         :param u0: initial guess for the flux. At locations in the domain marked by ``location_of_fixed_flux`` the value of ``u0`` is kept unchanged.
         :type u0: vector value on the domain (e.g. `Data`).
         :param p0: initial guess for the pressure. At locations in the domain marked by ``location_of_fixed_pressure`` the value of ``p0`` is kept unchanged.
         :type p0: scalar value on the domain (e.g. `Data`).
         :param verbose: if set some information on iteration progress are printed
         :type verbose: ``bool``
         :return: flux and pressure
         :rtype: ``tuple`` of `Data`.

         :note: The problem is solved as a least squares form

         *(I+D^*D)u+Qp=D^*f+g*
         *Q^*u+Q^*Qp=Q^*g*

         where *D* is the *div* operator and *(Qp)_i=k_{ij}p_{,j}* for the permeability *k_{ij}*.
         We eliminate the flux form the problem by setting

         *u=(I+D^*D)^{-1}(D^*f-g-Qp)* with u=u0 on location_of_fixed_flux

         form the first equation. Inserted into the second equation we get

         *Q^*(I-(I+D^*D)^{-1})Qp= Q^*(g-(I+D^*D)^{-1}(D^*f+g))* with p=p0  on location_of_fixed_pressure

         which is solved using the PCG method (precondition is *Q^*Q*). In each iteration step
         PDEs with operator *I+D^*D* and with *Q^*Q* needs to be solved using a sub iteration scheme.
         """
         self.verbose=verbose
         rtol=self.getTolerance()
         atol=self.getAbsoluteTolerance()
	 self.setSubProblemTolerance()
         num_corrections=0
         converged=False
         p=p0
         norm_r=None
         while not converged:
               v=self.getFlux(p, fixed_flux=u0)
               Qp=self.__Q(p)
               norm_v=self.__L2(v)
               norm_Qp=self.__L2(Qp)
               if norm_v == 0.: 
                  if norm_Qp == 0.: 
                     return v,p
                  else:
                    fac=norm_Qp
               else:
                  if norm_Qp == 0.: 
                    fac=norm_v
                  else:
                    fac=2./(1./norm_v+1./norm_Qp)
               ATOL=(atol+rtol*fac)
               if self.verbose: 
                    print "DarcyFlux: L2 norm of v = %e."%norm_v
                    print "DarcyFlux: L2 norm of k.grad(p) = %e."%norm_Qp
                    print "DarcyFlux: L2 defect u = %e."%(util.integrate(util.length(self.__g-util.interpolate(v,Function(self.domain))-Qp)**2)**(0.5),)
                    print "DarcyFlux: L2 defect div(v) = %e."%(util.integrate((self.__f-util.div(v))**2)**(0.5),)
                    print "DarcyFlux: absolute tolerance ATOL = %e."%ATOL
               if norm_r == None or norm_r>ATOL:
                   if num_corrections>max_num_corrections:
                         raise ValueError,"maximum number of correction steps reached."
                   p,r, norm_r=PCG(self.__g-util.interpolate(v,Function(self.domain))-Qp,self.__Aprod,p,self.__Msolve_PCG,self.__inner_PCG,atol=0.5*ATOL, rtol=0.,iter_max=max_iter, verbose=self.verbose)
                   num_corrections+=1
               else: 
                   converged=True
         return v,p 
    def __L2(self,v):
         return util.sqrt(util.integrate(util.length(util.interpolate(v,Function(self.domain)))**2))

    def __Q(self,p):
          return util.tensor_mult(self.__permeability,util.grad(p))

    def __Aprod(self,dp):
          if self.getSolverOptionsFlux().isVerbose(): print "DarcyFlux: Applying operator"
          Qdp=self.__Q(dp)
          self.__pde_v.setValue(Y=-Qdp,X=Data(), r=Data())
          du=self.__pde_v.getSolution()
          # self.__pde_v.getOperator().saveMM("proj.mm")
          return Qdp+du
    def __inner_GMRES(self,r,s):
         return util.integrate(util.inner(r,s))
 
    def __inner_PCG(self,p,r):
         return util.integrate(util.inner(self.__Q(p), r))

    def __Msolve_PCG(self,r):
	  if self.getSolverOptionsPressure().isVerbose(): print "DarcyFlux: Applying preconditioner"
          self.__pde_p.setValue(X=util.transposed_tensor_mult(self.__permeability,r), Y=Data(), r=Data())
          # self.__pde_p.getOperator().saveMM("prec.mm")
          return self.__pde_p.getSolution()

    def getFlux(self,p=None, fixed_flux=Data()):
        """
        returns the flux for a given pressure ``p`` where the flux is equal to ``fixed_flux``
        on locations where ``location_of_fixed_flux`` is positive (see `setValue`).
        Note that ``g`` and ``f`` are used, see `setValue`.

        :param p: pressure.
        :type p: scalar value on the domain (e.g. `Data`).
        :param fixed_flux: flux on the locations of the domain marked be ``location_of_fixed_flux``.
        :type fixed_flux: vector values on the domain (e.g. `Data`).
        :return: flux
        :rtype: `Data`
        :note: the method uses the least squares solution *u=(I+D^*D)^{-1}(D^*f-g-Qp)* where *D* is the *div* operator and *(Qp)_i=k_{ij}p_{,j}*
               for the permeability *k_{ij}*
        """
	self.setSubProblemTolerance()
        g=self.__g
        f=self.__f
        self.__pde_v.setValue(X=self.__l*f*util.kronecker(self.domain), r=fixed_flux)
        if p == None:
           self.__pde_v.setValue(Y=g)
        else:
           self.__pde_v.setValue(Y=g-self.__Q(p))
        return self.__pde_v.getSolution()

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
         self.__pde_u=LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_u.setSymmetryOn()
	 
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
	 return self.__pde_u.getSolverOptions()
     def setSolverOptionsVelocity(self, options=None):
         """
	 set the solver options for solving the equation for velocity.
	 
	 :param options: new solver  options
	 :type options: `SolverOptions`
	 """
         self.__pde_u.setSolverOptions(options)
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
            self.eta=util.interpolate(eta, Function(self.domain))
	    self.__pde_prec.setValue(D=1/self.eta)
            self.__pde_u.setValue(A=self.eta*(util.swap_axes(kk,0,3)+util.swap_axes(kk,1,3)))
        if restoration_factor!=None:
            n=self.domain.getNormal()
            self.__pde_u.setValue(d=restoration_factor*util.outer(n,n))
        if fixed_u_mask!=None:
            self.__pde_u.setValue(q=fixed_u_mask)
        if f!=None: self.__f=f
        if surface_stress!=None: self.__surface_stress=surface_stress
        if stress!=None: self.__stress=stress

     def initialize(self,f=Data(),fixed_u_mask=Data(),eta=1,surface_stress=Data(),stress=Data(), restoration_factor=0):
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
         return util.integrate(util.interpolate(p,Function(self.domain))*util.interpolate(Bv,Function(self.domain)))

     def inner_p(self,p0,p1):
         """
         Returns inner product of p0 and p1

         :param p0: a pressure
         :param p1: a pressure
         :return: inner product of p0 and p1
         :rtype: ``float``
         """
         s0=util.interpolate(p0,Function(self.domain))
         s1=util.interpolate(p1,Function(self.domain))
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
         self.__pde_u.setValue(Y=self.__f, y=self.__surface_stress)
	 self.getSolverOptionsVelocity().setTolerance(tol)
	 self.getSolverOptionsVelocity().setAbsoluteTolerance(0.)
         if self.__stress.isEmpty():
            self.__pde_u.setValue(X=p*util.kronecker(self.domain)-2*self.eta*util.symmetric(util.grad(v)))
         else:
            self.__pde_u.setValue(X=self.__stress+p*util.kronecker(self.domain)-2*self.eta*util.symmetric(util.grad(v)))
         out=self.__pde_u.getSolution()
         return  out

     def norm_Bv(self,Bv):
        """
        Returns Bv (overwrite).

        :rtype: equal to the type of p
        :note: boundary conditions on p should be zero!
        """
        return util.sqrt(util.integrate(util.interpolate(Bv,Function(self.domain))**2))

     def solve_AinvBt(self,p, tol):
         """
         Solves *Av=B^*p* with accuracy `tol`

         :param p: a pressure increment
         :return: the solution of *Av=B^*p*
         :note: boundary conditions on v should be zero!
         """
         self.__pde_u.setValue(Y=Data(), y=Data(), X=-p*util.kronecker(self.domain))
         out=self.__pde_u.getSolution()
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
