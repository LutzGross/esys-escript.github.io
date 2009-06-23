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
from linearPDEs import LinearPDE, LinearPDESystem, LinearSinglePDE, SolverOptions
from pdetools import HomogeneousSaddlePointProblem,Projector, ArithmeticTuple, PCG, NegativeNorm, GMRES

class DarcyFlow(object):
    """
    solves the problem

    M{u_i+k_{ij}*p_{,j} = g_i}
    M{u_{i,i} = f}

    where M{p} represents the pressure and M{u} the Darcy flux. M{k} represents the permeability,

    @note: The problem is solved in a least squares formulation.
    """

    def __init__(self, domain, weight=None, useReduced=False, adaptSubTolerance=True):
        """
        initializes the Darcy flux problem
        @param domain: domain of the problem
        @type domain: L{Domain}
	@param useReduced: uses reduced oreder on flux and pressure
	@type useReduced: C{bool}
	@param adaptSubTolerance: switches on automatic subtolerance selection
	@type adaptSubTolerance: C{bool}	
        """
        self.domain=domain
        if weight == None:
           s=self.domain.getSize()
           self.__l=(3.*util.longestEdge(self.domain)*s/util.sup(s))**2
           # self.__l=(3.*util.longestEdge(self.domain))**2
           # self.__l=(0.1*util.longestEdge(self.domain)*s/util.sup(s))**2
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
	
	M{(I+D^*D)u=F}
	
	@return: L{SolverOptions}
	"""
	return self.__pde_v.getSolverOptions()
    def setSolverOptionsFlux(self, options=None):
	"""
	Sets the solver options used to solve the flux problems
	
	M{(I+D^*D)u=F}
	
	If C{options} is not present, the options are reset to default
	@param options: L{SolverOptions}
	@note: if the adaption of subtolerance is choosen, the tolerance set by C{options} will be overwritten before the solver is called.
	"""
	return self.__pde_v.setSolverOptions(options)
    def getSolverOptionsPressure(self):
	"""
	Returns the solver options used to solve the pressure problems
	
	M{(Q^*Q)p=Q^*G}
	
	@return: L{SolverOptions}
	"""
	return self.__pde_p.getSolverOptions()
    def setSolverOptionsPressure(self, options=None):
	"""
	Sets the solver options used to solve the pressure problems
	
	M{(Q^*Q)p=Q^*G}
	
	If C{options} is not present, the options are reset to default
	@param options: L{SolverOptions}
	@note: if the adaption of subtolerance is choosen, the tolerance set by C{options} will be overwritten before the solver is called.
	"""
	return self.__pde_p.setSolverOptions(options)

    def setValue(self,f=None, g=None, location_of_fixed_pressure=None, location_of_fixed_flux=None, permeability=None):
        """
        assigns values to model parameters

        @param f: volumetic sources/sinks
        @type f: scalar value on the domain (e.g. L{Data})
        @param g: flux sources/sinks
        @type g: vector values on the domain (e.g. L{Data})
        @param location_of_fixed_pressure: mask for locations where pressure is fixed
        @type location_of_fixed_pressure: scalar value on the domain (e.g. L{Data})
        @param location_of_fixed_flux:  mask for locations where flux is fixed.
        @type location_of_fixed_flux: vector values on the domain (e.g. L{Data})
        @param permeability: permeability tensor. If scalar C{s} is given the tensor with
                             C{s} on the main diagonal is used. If vector C{v} is given the tensor with
                             C{v} on the main diagonal is used.
        @type permeability: scalar, vector or tensor values on the domain (e.g. L{Data})

        @note: the values of parameters which are not set by calling C{setValue} are not altered.
        @note: at any point on the boundary of the domain the pressure (C{location_of_fixed_pressure} >0)
               or the normal component of the flux (C{location_of_fixed_flux[i]>0} if direction of the normal
               is along the M{x_i} axis.
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
        sets the relative tolerance C{rtol} used to terminate the solution process. The iteration is terminated if

        M{|g-v-Qp| <= atol + rtol * min( max( |g-v|, |Qp| ), max( |v|, |g-Qp| ) ) }

        where C{atol} is an absolut tolerance (see L{setAbsoluteTolerance}), M{|f|^2 = integrate(length(f)^2)} and M{(Qp)_i=k_{ij}p_{,j}} for the permeability M{k_{ij}}.

        @param rtol: relative tolerance for the pressure
        @type rtol: non-negative C{float}
        """
        if rtol<0:
            raise ValueError,"Relative tolerance needs to be non-negative."
        self.__rtol=rtol
    def getTolerance(self):
        """
        returns the relative tolerance

        @return: current relative tolerance
        @rtype: C{float}
        """
        return self.__rtol

    def setAbsoluteTolerance(self,atol=0.):
        """
        sets the absolute tolerance C{atol} used to terminate the solution process. The iteration is terminated if

        M{|g-v-Qp| <= atol + rtol * min( max( |g-v|, |Qp| ), max( |v|, |g-Qp| ) ) }

        where C{rtol} is an absolut tolerance (see L{setTolerance}), M{|f|^2 = integrate(length(f)^2)} and M{(Qp)_i=k_{ij}p_{,j}} for the permeability M{k_{ij}}.

        @param atol: absolute tolerance for the pressure
        @type atol: non-negative C{float}
        """
        if atol<0:
            raise ValueError,"Absolute tolerance needs to be non-negative."
        self.__atol=atol
    def getAbsoluteTolerance(self):
       """
       returns the absolute tolerance 
       
       @return: current absolute tolerance
       @rtype: C{float}
       """
       return self.__atol
    def getSubProblemTolerance(self):
	"""
	Returns a suitable subtolerance
	@type: C{float}
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

         @param u0: initial guess for the flux. At locations in the domain marked by C{location_of_fixed_flux} the value of C{u0} is kept unchanged.
         @type u0: vector value on the domain (e.g. L{Data}).
         @param p0: initial guess for the pressure. At locations in the domain marked by C{location_of_fixed_pressure} the value of C{p0} is kept unchanged.
         @type p0: scalar value on the domain (e.g. L{Data}).
         @param verbose: if set some information on iteration progress are printed
         @type verbose: C{bool}
         @return: flux and pressure
         @rtype: C{tuple} of L{Data}.

         @note: The problem is solved as a least squares form

         M{(I+D^*D)u+Qp=D^*f+g}
         M{Q^*u+Q^*Qp=Q^*g}

         where M{D} is the M{div} operator and M{(Qp)_i=k_{ij}p_{,j}} for the permeability M{k_{ij}}.
         We eliminate the flux form the problem by setting

         M{u=(I+D^*D)^{-1}(D^*f-g-Qp)} with u=u0 on location_of_fixed_flux

         form the first equation. Inserted into the second equation we get

         M{Q^*(I-(I+D^*D)^{-1})Qp= Q^*(g-(I+D^*D)^{-1}(D^*f+g))} with p=p0  on location_of_fixed_pressure

         which is solved using the PCG method (precondition is M{Q^*Q}). In each iteration step
         PDEs with operator M{I+D^*D} and with M{Q^*Q} needs to be solved using a sub iteration scheme.
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
        returns the flux for a given pressure C{p} where the flux is equal to C{fixed_flux}
        on locations where C{location_of_fixed_flux} is positive (see L{setValue}).
        Note that C{g} and C{f} are used, see L{setValue}.

        @param p: pressure.
        @type p: scalar value on the domain (e.g. L{Data}).
        @param fixed_flux: flux on the locations of the domain marked be C{location_of_fixed_flux}.
        @type fixed_flux: vector values on the domain (e.g. L{Data}).
        @param tol: relative tolerance to be used.
        @type tol: positive C{float}.
        @return: flux
        @rtype: L{Data}
        @note: the method uses the least squares solution M{u=(I+D^*D)^{-1}(D^*f-g-Qp)} where M{D} is the M{div} operator and M{(Qp)_i=k_{ij}p_{,j}}
               for the permeability M{k_{ij}}
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
     def __init__(self,domain,adaptSubTolerance=True, **kwargs):
         """
         initialize the Stokes Problem

         @param domain: domain of the problem. The approximation order needs to be two.
         @type domain: L{Domain}
	 @param adaptSubTolerance: If True the tolerance for subproblem is set automatically.
	 @type adaptSubTolerance: C{bool}
         @warning: The apprximation order needs to be two otherwise you may see oscilations in the pressure.
         """
         HomogeneousSaddlePointProblem.__init__(self,adaptSubTolerance=adaptSubTolerance,**kwargs)
         self.domain=domain
         self.vol=util.integrate(1.,Function(self.domain))
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
	 
	 @rtype: L{SolverOptions}
	 """
	 return self.__pde_u.getSolverOptions()
     def setSolverOptionsVelocity(self, options=None):
         """
	 set the solver options for solving the equation for velocity.
	 
	 @param options: new solver  options
	 @type options: L{SolverOptions}
	 """
         self.__pde_u.setSolverOptions(options)
     def getSolverOptionsPressure(self):
         """
	 returns the solver options used  solve the equation for pressure.
	 @rtype: L{SolverOptions}
	 """
	 return self.__pde_prec.getSolverOptions()
     def setSolverOptionsPressure(self, options=None):
         """
	 set the solver options for solving the equation for pressure.
	 @param options: new solver  options
	 @type options: L{SolverOptions}
	 """
	 self.__pde_prec.setSolverOptions(options)

     def setSolverOptionsDiv(self, options=None):
         """
	 set the solver options for solving the equation to project the divergence of
	 the velocity onto the function space of presure.
	 
	 @param options: new solver options
	 @type options: L{SolverOptions}
	 """
	 self.__pde_prec.setSolverOptions(options)
     def getSolverOptionsDiv(self):
         """
	 returns the solver options for solving the equation to project the divergence of
	 the velocity onto the function space of presure.
	 
	 @rtype: L{SolverOptions}
	 """
	 return self.__pde_prec.getSolverOptions()
     def setSubProblemTolerance(self):
         """
	 Updates the tolerance for subproblems
         """
	 if self.adaptSubTolerance():
             sub_tol=self.getSubProblemTolerance()
	     self.getSolverOptionsDiv().setTolerance(sub_tol)
	     self.getSolverOptionsDiv().setAbsoluteTolerance(0.)
	     self.getSolverOptionsPressure().setTolerance(sub_tol)
	     self.getSolverOptionsPressure().setAbsoluteTolerance(0.)
	     self.getSolverOptionsVelocity().setTolerance(sub_tol)
	     self.getSolverOptionsVelocity().setAbsoluteTolerance(0.)
	     

     def initialize(self,f=Data(),fixed_u_mask=Data(),eta=1,surface_stress=Data(),stress=Data()):
        """
        assigns values to the model parameters

        @param f: external force
        @type f: L{Vector} object in L{FunctionSpace} L{Function} or similar
        @param fixed_u_mask: mask of locations with fixed velocity.
        @type fixed_u_mask: L{Vector} object on L{FunctionSpace} L{Solution} or similar
        @param eta: viscosity
        @type eta: L{Scalar} object on L{FunctionSpace} L{Function} or similar
        @param surface_stress: normal surface stress
        @type eta: L{Vector} object on L{FunctionSpace} L{FunctionOnBoundary} or similar
        @param stress: initial stress
	@type stress: L{Tensor} object on L{FunctionSpace} L{Function} or similar
        @note: All values needs to be set.
        """
        self.eta=eta
        A =self.__pde_u.createCoefficient("A")
	self.__pde_u.setValue(A=Data())
        for i in range(self.domain.getDim()):
		for j in range(self.domain.getDim()):
			A[i,j,j,i] += 1.
			A[i,j,i,j] += 1.
	self.__pde_prec.setValue(D=1/self.eta)
        self.__pde_u.setValue(A=A*self.eta,q=fixed_u_mask)
        self.__f=f
        self.__surface_stress=surface_stress
        self.__stress=stress

     def Bv(self,v):
         """
         returns inner product of element p and div(v)

         @param p: a pressure increment
         @param v: a residual
         @return: inner product of element p and div(v)
         @rtype: C{float}
         """
         self.__pde_proj.setValue(Y=-util.div(v))
         return self.__pde_proj.getSolution()

     def inner_pBv(self,p,Bv):
         """
         returns inner product of element p and Bv=-div(v)

         @param p: a pressure increment
         @param v: a residual
         @return: inner product of element p and Bv=-div(v)
         @rtype: C{float}
         """
         return util.integrate(util.interpolate(p,Function(self.domain))*util.interpolate(Bv,Function(self.domain)))

     def inner_p(self,p0,p1):
         """
         Returns inner product of p0 and p1

         @param p0: a pressure
         @param p1: a pressure
         @return: inner product of p0 and p1
         @rtype: C{float}
         """
         s0=util.interpolate(p0/self.eta,Function(self.domain))
         s1=util.interpolate(p1/self.eta,Function(self.domain))
         return util.integrate(s0*s1)

     def norm_v(self,v):
         """
         returns the norm of v

         @param v: a velovity
         @return: norm of v
         @rtype: non-negative C{float}
         """
         return util.sqrt(util.integrate(util.length(util.grad(v))))

     def getV(self, p, v0):
         """
         return the value for v for a given p (overwrite)

         @param p: a pressure
         @param v0: a initial guess for the value v to return.
         @return: v given as M{v= A^{-1} (f-B^*p)}
         """
         self.__pde_u.setValue(Y=self.__f, y=self.__surface_stress, r=v0)
         if self.__stress.isEmpty():
            self.__pde_u.setValue(X=p*util.kronecker(self.domain))
         else:
            self.__pde_u.setValue(X=self.__stress+p*util.kronecker(self.domain))
         out=self.__pde_u.getSolution()
         return  out

     def norm_Bv(self,Bv):
        """
        Returns Bv (overwrite).

        @rtype: equal to the type of p
        @note: boundary conditions on p should be zero!
        """
        return util.sqrt(util.integrate(util.interpolate(Bv,Function(self.domain))**2))

     def solve_AinvBt(self,p):
         """
         Solves M{Av=B^*p} with accuracy L{self.getSubProblemTolerance()}

         @param p: a pressure increment
         @return: the solution of M{Av=B^*p}
         @note: boundary conditions on v should be zero!
         """
         self.__pde_u.setValue(Y=Data(), y=Data(), r=Data(),X=-p*util.kronecker(self.domain))
         out=self.__pde_u.getSolution()
         return  out

     def solve_prec(self,Bv):
         """
         applies preconditioner for for M{BA^{-1}B^*} to M{Bv}
         with accuracy L{self.getSubProblemTolerance()} 

         @param v: velocity increment
         @return: M{p=P(Bv)} where M{P^{-1}} is an approximation of M{BA^{-1}B^*}
         @note: boundary conditions on p are zero.
         """
         self.__pde_prec.setValue(Y=Bv)
         return self.__pde_prec.getSolution()
