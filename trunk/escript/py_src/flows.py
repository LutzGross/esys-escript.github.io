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
__url__="http://www.uq.edu.au/esscc/escript-finley"

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
from linearPDEs import LinearPDE, LinearPDESystem, LinearSinglePDE
from pdetools import HomogeneousSaddlePointProblem,Projector, ArithmeticTuple, PCG, NegativeNorm

class DarcyFlow(object):
    """
    solves the problem 

    M{u_i+k_{ij}*p_{,j} = g_i}
    M{u_{i,i} = f}

    where M{p} represents the pressure and M{u} the Darcy flux. M{k} represents the permeability, 

    @note: The problem is solved in a least squares formulation.
    """

    def __init__(self, domain,useReduced=False):
        """
        initializes the Darcy flux problem
        @param domain: domain of the problem
        @type domain: L{Domain}
        """
        self.domain=domain
        self.__pde_v=LinearPDESystem(domain)
        if useReduced: self.__pde_v.setReducedOrderOn()
        self.__pde_v.setSymmetryOn()
        self.__pde_v.setValue(D=util.kronecker(domain), A=util.outer(util.kronecker(domain),util.kronecker(domain)))
        self.__pde_p=LinearSinglePDE(domain)
        self.__pde_p.setSymmetryOn()
        if useReduced: self.__pde_p.setReducedOrderOn()
        self.__f=Scalar(0,self.__pde_v.getFunctionSpaceForCoefficient("X"))
        self.__g=Vector(0,self.__pde_v.getFunctionSpaceForCoefficient("Y"))
        self.__ATOL= None

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
           self.f=f
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


    def getFlux(self,p=None, fixed_flux=Data(),tol=1.e-8, show_details=False):
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
        self.__pde_v.setTolerance(tol)
        g=self.__g
        f=self.__f
        self.__pde_v.setValue(X=f*util.kronecker(self.domain), r=fixed_flux)
        if p == None:
           self.__pde_v.setValue(Y=g)
        else:
           self.__pde_v.setValue(Y=g-util.tensor_mult(self.__permeability,util.grad(p)))
        return self.__pde_v.getSolution(verbose=show_details)

    def getPressure(self,v=None, fixed_pressure=Data(),tol=1.e-8, show_details=False):
        """
        returns the pressure for a given flux C{v} where the pressure is equal to C{fixed_pressure}
        on locations where C{location_of_fixed_pressure} is positive (see L{setValue}). 
        Note that C{g} is used, see L{setValue}.
        
        @param v: flux.
        @type v: vector-valued on the domain (e.g. L{Data}).
        @param fixed_pressure: pressure on the locations of the domain marked be C{location_of_fixed_pressure}.
        @type fixed_pressure: vector values on the domain (e.g. L{Data}).
        @param tol: relative tolerance to be used.
        @type tol: positive C{float}.
        @return: pressure 
        @rtype: L{Data}
        @note: the method uses the least squares solution M{p=(Q^*Q)^{-1}Q^*(g-u)} where and M{(Qp)_i=k_{ij}p_{,j}}
               for the permeability M{k_{ij}}
        """
        self.__pde_v.setTolerance(tol)
        g=self.__g
        self.__pde_p.setValue(r=fixed_pressure)
        if v == None:
           self.__pde_p.setValue(X=util.transposed_tensor_mult(self.__permeability,g-v))
        else:
           self.__pde_p.setValue(X=util.transposed_tensor_mult(self.__permeability,g))
        return self.__pde_p.getSolution(verbose=show_details)

    def setTolerance(self,atol=0,rtol=1e-8,p_ref=None,v_ref=None):
        """
        set the tolerance C{ATOL} used to terminate the solution process. It is used

        M{ATOL = atol + rtol * max( |g-v_ref|, |Qp_ref| )}

        where M{|f|^2 = integrate(length(f)^2)} and M{(Qp)_i=k_{ij}p_{,j}} for the permeability M{k_{ij}}. If C{v_ref} or C{p_ref} is not present zero is assumed.

        The iteration is terminated if for the current approximation C{p}, flux C{v=(I+D^*D)^{-1}(D^*f-g-Qp)} and their residual

        M{r=Q^*(g-Qp-v)}

        the condition 

        M{<(Q^*Q)^{-1} r,r> <= ATOL} 

        holds. M{D} is the M{div} operator and M{(Qp)_i=k_{ij}p_{,j}} for the permeability M{k_{ij}}

        @param atol: absolute tolerance for the pressure
        @type atol: non-negative C{float}
        @param rtol: relative tolerance for the pressure
        @type rtol: non-negative C{float}
        @param p_ref: reference pressure. If not present zero is used. You may use physical arguments to set a resonable value for C{p_ref}, use the 
        L{getPressure} method or use  the value from a previous time step.
        @type p_ref: scalar value on the domain (e.g. L{Data}).
        @param v_ref: reference velocity.  If not present zero is used. You may use physical arguments to set a resonable value for C{v_ref}, use the 
        L{getFlux} method or use  the value from a previous time step.
        @type v_ref: vector-valued on the domain (e.g. L{Data}).
        @return: used absolute tolerance.
        @rtype: positive C{float}
        """
        g=self.__g
        if not v_ref == None:
           f1=util.integrate(util.length(util.interpolate(g-v_ref,Function(self.domain)))**2)
        else:
           f1=util.integrate(util.length(util.interpolate(g))**2)
        if not p_ref == None:
           f2=util.integrate(util.length(util.tensor_mult(self.__permeability,util.grad(p_ref)))**2)
        else:
           f2=0
        self.__ATOL= atol + rtol * util.sqrt(max(f1,f2))
        if self.__ATOL<=0:
           raise ValueError,"Positive tolerance (=%e) is expected."%self.__ATOL
        return self.__ATOL
        
    def getTolerance(self):
        """
        returns the current tolerance.
   
        @return: used absolute tolerance.
        @rtype: positive C{float}
        """
        if self.__ATOL==None:
           raise ValueError,"no tolerance is defined."
        return self.__ATOL

    def solve(self,u0,p0, max_iter=100, verbose=False, show_details=False, sub_rtol=1.e-8):
         """ 
         solves the problem.
 
         The iteration is terminated if the residual norm is less then self.getTolerance().

         @param u0: initial guess for the flux. At locations in the domain marked by C{location_of_fixed_flux} the value of C{u0} is kept unchanged.
         @type u0: vector value on the domain (e.g. L{Data}).
         @param p0: initial guess for the pressure. At locations in the domain marked by C{location_of_fixed_pressure} the value of C{p0} is kept unchanged.
         @type p0: scalar value on the domain (e.g. L{Data}).
         @param sub_rtol: tolerance to be used in the sub iteration. It is recommended that M{sub_rtol<rtol*5.e-3}
         @type sub_rtol: positive-negative C{float}
         @param verbose: if set some information on iteration progress are printed
         @type verbose: C{bool}
         @param show_details:  if set information on the subiteration process are printed.
         @type show_details: C{bool}
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
         self.show_details= show_details and self.verbose
         self.__pde_v.setTolerance(sub_rtol) 
         self.__pde_p.setTolerance(sub_rtol)
         ATOL=self.getTolerance()
         if self.verbose: print "DarcyFlux: absolute tolerance = %e"%ATOL
         #########################################################################################################################
         # 
         #   we solve:
         #  
         #      Q^*(I-(I+D^*D)^{-1})Q dp =  Q^* (g-u0-Qp0 - (I+D^*D)^{-1} ( D^*(f-Du0)+g-u0-Qp0) ) 
         #
         #   residual is 
         #
         #    r=  Q^* (g-u0-Qp0 - (I+D^*D)^{-1} ( D^*(f-Du0)+g-u0-Qp0) - Q dp +(I+D^*D)^{-1})Q dp ) = Q^* (g - Qp - v)
         #
         #        with v = (I+D^*D)^{-1} (D^*f+g-Qp) including BC 
         #
         #    we use (g - Qp, v) to represent the residual. not that 
         #
         #    dr(dp)=( -Q(dp), dv) with dv = - (I+D^*D)^{-1} Q(dp)
         #
         #   while the initial residual is 
         #
         #      r0=( g - Qp0, v00) with v00=(I+D^*D)^{-1} (D^*f+g-Qp0) including BC
         #   
         d0=self.__g-util.tensor_mult(self.__permeability,util.grad(p0))
         self.__pde_v.setValue(Y=d0, X=self.__f*util.kronecker(self.domain), r=u0)
         v00=self.__pde_v.getSolution(verbose=show_details)
         if self.verbose: print "DarcyFlux: range of initial flux = ",util.inf(v00), util.sup(v00)
         self.__pde_v.setValue(r=Data())
         # start CG
         r=ArithmeticTuple(d0, v00)
         p,r=PCG(r,self.__Aprod_PCG,p0,self.__Msolve_PCG,self.__inner_PCG,atol=ATOL, rtol=0.,iter_max=max_iter, verbose=self.verbose)
         return r[1],p

    def __Aprod_PCG(self,dp):
          if self.show_details: print "DarcyFlux: Applying operator"
          #  -dr(dp) = (Qdp,du) where du = (I+D^*D)^{-1} (Qdp)
          mQdp=util.tensor_mult(self.__permeability,util.grad(dp))
          self.__pde_v.setValue(Y=mQdp,X=Data(), r=Data())
          du=self.__pde_v.getSolution(verbose=self.show_details)
          return ArithmeticTuple(mQdp,du)

    def __inner_PCG(self,p,r):
         a=util.tensor_mult(self.__permeability,util.grad(p))
         f0=util.integrate(util.inner(a,r[0]))
         f1=util.integrate(util.inner(a,r[1]))
         # print "__inner_PCG:",f0,f1,"->",f0-f1
         return f0-f1

    def __Msolve_PCG(self,r):
          if self.show_details: print "DarcyFlux: Applying preconditioner"
          self.__pde_p.setValue(X=util.transposed_tensor_mult(self.__permeability,r[0]-r[1]), r=Data())
          return self.__pde_p.getSolution(verbose=self.show_details)

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

         @param domain: domain of the problem. The approximation order needs to be two.
         @type domain: L{Domain}
         @warning: The apprximation order needs to be two otherwise you may see oscilations in the pressure.
         """
         HomogeneousSaddlePointProblem.__init__(self,**kwargs)
         self.domain=domain
         self.vol=util.integrate(1.,Function(self.domain))
         self.__pde_u=LinearPDE(domain,numEquations=self.domain.getDim(),numSolutions=self.domain.getDim())
         self.__pde_u.setSymmetryOn()
         # self.__pde_u.setSolverMethod(self.__pde_u.DIRECT)
         # self.__pde_u.setSolverMethod(preconditioner=LinearPDE.RILU)
            
         self.__pde_prec=LinearPDE(domain)
         self.__pde_prec.setReducedOrderOn()
         # self.__pde_prec.setSolverMethod(self.__pde_prec.LUMPING)
         self.__pde_prec.setSymmetryOn()

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

#===============================================================================================================
     def inner_pBv(self,p,v):
         """
         returns inner product of element p and div(v)

         @param p: a pressure increment
         @param v: a residual
         @return: inner product of element p and div(v)
         @rtype: C{float}
         """
         return util.integrate(-p*util.div(v))

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
         self.__pde_u.setTolerance(self.getSubProblemTolerance())
         self.__pde_u.setValue(Y=self.__f, y=self.__surface_stress, r=v0)
         if self.__stress.isEmpty():
            self.__pde_u.setValue(X=p*util.kronecker(self.domain))
         else:
            self.__pde_u.setValue(X=self.__stress+p*util.kronecker(self.domain))
         out=self.__pde_u.getSolution(verbose=self.show_details)
         return  out


         raise NotImplementedError,"no v calculation implemented."

 
     def norm_Bv(self,v):
        """
        Returns Bv (overwrite).

        @rtype: equal to the type of p
        @note: boundary conditions on p should be zero!
        """
        return util.sqrt(util.integrate(util.div(v)**2))

     def solve_AinvBt(self,p):
         """
         Solves M{Av=B^*p} with accuracy L{self.getSubProblemTolerance()}

         @param p: a pressure increment
         @return: the solution of M{Av=B^*p} 
         @note: boundary conditions on v should be zero!
         """
         self.__pde_u.setTolerance(self.getSubProblemTolerance())
         self.__pde_u.setValue(Y=Data(), y=Data(), r=Data(),X=-p*util.kronecker(self.domain))
         out=self.__pde_u.getSolution(verbose=self.show_details)
         return  out

     def solve_precB(self,v):
         """
         applies preconditioner for for M{BA^{-1}B^*} to M{Bv}
         with accuracy L{self.getSubProblemTolerance()} (overwrite).

         @param v: velocity increment
         @return: M{p=P(Bv)} where M{P^{-1}} is an approximation of M{BA^{-1}B^*}
         @note: boundary conditions on p are zero.
         """
         self.__pde_prec.setValue(Y=-util.div(v))
         self.__pde_prec.setTolerance(self.getSubProblemTolerance())
         return self.__pde_prec.getSolution(verbose=self.show_details)
