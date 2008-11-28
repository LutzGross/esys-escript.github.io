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
from pdetools import HomogeneousSaddlePointProblem,Projector, ArithmeticTuple, PCG

class DarcyFlow(object):
    """
    solves the problem 

    M{u_i+k_{ij}*p_{,j} = g_i}
    M{u_{i,i} = f}

    where M{p} represents the pressure and M{u} the Darcy flux. M{k} represents the permeability, 

    @note: The problem is solved in a least squares formulation.
    """

    def __init__(self, domain):
        """
        initializes the Darcy flux problem
        @param domain: domain of the problem
        @type domain: L{Domain}
        """
        self.domain=domain
        self.__pde_v=LinearPDESystem(domain)
        self.__pde_v.setValue(D=util.kronecker(domain), A=util.outer(util.kronecker(domain),util.kronecker(domain)))
        self.__pde_v.setSymmetryOn()
        self.__pde_p=LinearSinglePDE(domain)
        self.__pde_p.setSymmetryOn()
        self.__f=Scalar(0,self.__pde_v.getFunctionSpaceForCoefficient("X"))
        self.__g=Vector(0,self.__pde_v.getFunctionSpaceForCoefficient("Y"))

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


    def getFlux(self,p, fixed_flux=Data(),tol=1.e-8):
        """
        returns the flux for a given pressure C{p} where the flux is equal to C{fixed_flux}
        on locations where C{location_of_fixed_flux} is positive (see L{setValue}). 
        Note that C{g} and C{f} are used, L{setValue}.
        
        @param p: pressure.
        @type p: scalar value on the domain (e.g. L{Data}).
        @param fixed_flux: flux on the locations of the domain marked be C{location_of_fixed_flux}.
        @type fixed_flux: vector values on the domain (e.g. L{Data}).
        @param tol: relative tolerance to be used.
        @type tol: positive float.
        @return: flux 
        @rtype: L{Data}
        @note: the method uses the least squares solution M{u=(I+D^*D)^{-1}(D^*f-g-Qp)} where M{D} is the M{div} operator and M{(Qp)_i=k_{ij}p_{,j}}
               for the permeability M{k_{ij}}
        """
        self.__pde_v.setTolerance(tol)
        self.__pde_v.setValue(Y=self.__g, X=self.__f*util.kronecker(self.domain), r=boundary_flux)
        return self.__pde_v.getSolution()

    def solve(self,u0,p0,atol=0,rtol=1e-8, max_iter=100, verbose=False, show_details=False, sub_rtol=1.e-8):
         """ 
         solves the problem.
 
         The iteration is terminated if the error in the pressure is less then C{rtol * |q| + atol} where 
         C{|q|} denotes the norm of the right hand side (see escript user's guide for details).

         @param u0: initial guess for the flux. At locations in the domain marked by C{location_of_fixed_flux} the value of C{u0} is kept unchanged.
         @type u0: vector value on the domain (e.g. L{Data}).
         @param p0: initial guess for the pressure. At locations in the domain marked by C{location_of_fixed_pressure} the value of C{p0} is kept unchanged.
         @type p0: scalar value on the domain (e.g. L{Data}).
         @param atol: absolute tolerance for the pressure
         @type atol: non-negative C{float}
         @param rtol: relative tolerance for the pressure
         @type rtol: non-negative C{float}
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
         p2=p0*self.__pde_p.getCoefficient("q")
         u2=u0*self.__pde_v.getCoefficient("q")
         g=self.__g-u2-util.tensor_mult(self.__permeability,util.grad(p2))
         f=self.__f-util.div(u2)
         self.__pde_v.setValue(Y=g, X=f*util.kronecker(self.domain), r=Data())
         dv=self.__pde_v.getSolution(verbose=show_details)
         self.__pde_p.setValue(X=util.transposed_tensor_mult(self.__permeability,g-dv))
         self.__pde_p.setValue(r=Data())
         dp=self.__pde_p.getSolution(verbose=self.show_details)
         norm_rhs=self.__inner_PCG(dp,ArithmeticTuple(g,dv))
         if norm_rhs<0:
             raise NegativeNorm,"negative norm. Maybe the sub-tolerance is too large."
         ATOL=util.sqrt(norm_rhs)*rtol +atol
         if not ATOL>0:
             raise ValueError,"Negative absolute tolerance (rtol = %e, norm right hand side =%, atol =%e)."%(rtol, util.sqrt(norm_rhs), atol)
         rhs=ArithmeticTuple(g,dv)
         dp,r=PCG(rhs,self.__Aprod_PCG,self.__Msolve_PCG,self.__inner_PCG,atol=ATOL, rtol=0.,iter_max=max_iter, x=p0-p2, verbose=self.verbose, initial_guess=True)
         return u2+r[1],p2+dp
        
    def __Aprod_PCG(self,p):
          if self.show_details: print "DarcyFlux: Applying operator"
          Qp=util.tensor_mult(self.__permeability,util.grad(p))
          self.__pde_v.setValue(Y=Qp,X=Data())
          w=self.__pde_v.getSolution(verbose=self.show_details)
          return ArithmeticTuple(Qp,w)

    def __inner_PCG(self,p,r):
         a=util.tensor_mult(self.__permeability,util.grad(p))
         return util.integrate(util.inner(a,r[0]-r[1]))

    def __Msolve_PCG(self,r):
          if self.show_details: print "DarcyFlux: Applying preconditioner"
          self.__pde_p.setValue(X=util.transposed_tensor_mult(self.__permeability,r[0]-r[1]))
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
         self.__pde_prec.setSolverMethod(self.__pde_prec.LUMPING)
         self.__pde_prec.setSymmetryOn()

         self.__pde_proj=LinearPDE(domain)
         self.__pde_proj.setReducedOrderOn()
         self.__pde_proj.setSymmetryOn()
         self.__pde_proj.setValue(D=1.)

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
        self.__pde_u.setValue(A=A*self.eta,q=fixed_u_mask,Y=f,y=surface_stress)
        self.__stress=stress

      def B(self,v):
        """
        returns div(v)
        @rtype: equal to the type of p

        @note: boundary conditions on p should be zero!
        """
        if self.show_details: print "apply divergence:"
        self.__pde_proj.setValue(Y=-util.div(v))
        self.__pde_proj.setTolerance(self.getSubProblemTolerance())
        return self.__pde_proj.getSolution(verbose=self.show_details)

      def inner_pBv(self,p,Bv):
         """
         returns inner product of element p and Bv  (overwrite)
         
         @type p: equal to the type of p
         @type Bv: equal to the type of result of operator B
         @rtype: C{float}

         @rtype: equal to the type of p
         """
         s0=util.interpolate(p,Function(self.domain))
         s1=util.interpolate(Bv,Function(self.domain))
         return util.integrate(s0*s1)

      def inner_p(self,p0,p1):
         """
         returns inner product of element p0 and p1  (overwrite)
         
         @type p0: equal to the type of p
         @type p1: equal to the type of p
         @rtype: C{float}

         @rtype: equal to the type of p
         """
         s0=util.interpolate(p0/self.eta,Function(self.domain))
         s1=util.interpolate(p1/self.eta,Function(self.domain))
         return util.integrate(s0*s1)

      def inner_v(self,v0,v1):
         """
         returns inner product of two element v0 and v1  (overwrite)
         
         @type v0: equal to the type of v
         @type v1: equal to the type of v
         @rtype: C{float}

         @rtype: equal to the type of v
         """
	 gv0=util.grad(v0)
	 gv1=util.grad(v1)
         return util.integrate(util.inner(gv0,gv1))

      def solve_A(self,u,p):
         """
         solves Av=f-Au-B^*p (v=0 on fixed_u_mask)
         """
         if self.show_details: print "solve for velocity:"
         self.__pde_u.setTolerance(self.getSubProblemTolerance())
         if self.__stress.isEmpty():
            self.__pde_u.setValue(X=-2*self.eta*util.symmetric(util.grad(u))+p*util.kronecker(self.domain))
         else:
            self.__pde_u.setValue(X=self.__stress-2*self.eta*util.symmetric(util.grad(u))+p*util.kronecker(self.domain))
         out=self.__pde_u.getSolution(verbose=self.show_details)
         return  out

      def solve_prec(self,p):
         if self.show_details: print "apply preconditioner:"
         self.__pde_prec.setTolerance(self.getSubProblemTolerance())
         self.__pde_prec.setValue(Y=p)
         q=self.__pde_prec.getSolution(verbose=self.show_details)
         return q
