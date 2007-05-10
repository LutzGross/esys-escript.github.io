# $Id$
"""
The module provides an interface to define and solve linear partial
differential equations (PDEs) within L{escript}. L{linearPDEs} does not provide any
solver capabilities in itself but hands the PDE over to
the PDE solver library defined through the L{Domain<escript.Domain>} of the PDE.
The general interface is provided through the L{LinearPDE} class. The
L{AdvectivePDE} which is derived from the L{LinearPDE} class
provides an interface to PDE dominated by its advective terms. The L{Poisson},
L{Helmholtz}, L{LameEquation}, L{AdvectivePDE}
classs which are also derived form the L{LinearPDE} class should be used
to define of solve these sepecial PDEs.

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

import escript
import util
import numarray

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


class IllegalCoefficient(ValueError):
   """
   raised if an illegal coefficient of the general ar particular PDE is requested.
   """
   pass

class IllegalCoefficientValue(ValueError):
   """
   raised if an incorrect value for a coefficient is used.
   """
   pass

class IllegalCoefficientFunctionSpace(ValueError):
   """
   raised if an incorrect function space for a coefficient is used.
   """

class UndefinedPDEError(ValueError):
   """
   raised if a PDE is not fully defined yet.
   """
   pass

class PDECoefficient(object):
    """
    A class for describing a PDE coefficient

    @cvar INTERIOR: indicator that coefficient is defined on the interior of the PDE domain
    @cvar BOUNDARY: indicator that coefficient is defined on the boundary of the PDE domain
    @cvar CONTACT: indicator that coefficient is defined on the contact region within the PDE domain
    @cvar INTERIOR_REDUCED: indicator that coefficient is defined on the interior of the PDE domain using a reduced integration order
    @cvar BOUNDARY_REDUCED: indicator that coefficient is defined on the boundary of the PDE domain using a reduced integration order
    @cvar CONTACT_REDUCED: indicator that coefficient is defined on the contact region within the PDE domain using a reduced integration order
    @cvar SOLUTION: indicator that coefficient is defined trough a solution of the PDE
    @cvar REDUCED: indicator that coefficient is defined trough a reduced solution of the PDE
    @cvar BY_EQUATION: indicator that the dimension of the coefficient shape is defined by the number PDE equations
    @cvar BY_SOLUTION: indicator that the dimension of the coefficient shape is defined by the number PDE solutions
    @cvar BY_DIM: indicator that the dimension of the coefficient shape is defined by the spatial dimension
    @cvar OPERATOR: indicator that the the coefficient alters the operator of the PDE
    @cvar RIGHTHANDSIDE: indicator that the the coefficient alters the right hand side of the PDE
    @cvar BOTH: indicator that the the coefficient alters the operator as well as the right hand side of the PDE

    """
    INTERIOR=0
    BOUNDARY=1
    CONTACT=2
    SOLUTION=3
    REDUCED=4
    BY_EQUATION=5
    BY_SOLUTION=6
    BY_DIM=7
    OPERATOR=10
    RIGHTHANDSIDE=11
    BOTH=12
    INTERIOR_REDUCED=13
    BOUNDARY_REDUCED=14
    CONTACT_REDUCED=15

    def __init__(self, where, pattern, altering):
       """
       Initialise a PDE Coefficient type

       @param where: describes where the coefficient lives
       @type where: one of L{INTERIOR}, L{BOUNDARY}, L{CONTACT}, L{SOLUTION}, L{REDUCED}, 
                           L{INTERIOR_REDUCED}, L{BOUNDARY_REDUCED}, L{CONTACT_REDUCED}.
       @param pattern: describes the shape of the coefficient and how the shape is build for a given
              spatial dimension and numbers of equation and solution in then PDE. For instance,
              (L{BY_EQUATION},L{BY_SOLUTION},L{BY_DIM}) descrbes a rank 3 coefficient which
              is instanciated as shape (3,2,2) in case of a three equations and two solution components
              on a 2-dimensional domain. In the case of single equation and a single solution component
              the shape compoments marked by L{BY_EQUATION} or L{BY_SOLUTION} are dropped. In this case
              the example would be read as (2,).
       @type pattern: C{tuple} of L{BY_EQUATION}, L{BY_SOLUTION}, L{BY_DIM}
       @param altering: indicates what part of the PDE is altered if the coefficiennt is altered
       @type altering: one of L{OPERATOR}, L{RIGHTHANDSIDE}, L{BOTH}
       @param reduced: indicates if reduced 
       @type reduced: C{bool}
       """
       
       super(PDECoefficient, self).__init__()
       self.what=where
       self.pattern=pattern
       self.altering=altering
       self.resetValue()

    def resetValue(self):
       """
       resets coefficient value to default
       """
       self.value=escript.Data()

    def getFunctionSpace(self,domain,reducedEquationOrder=False,reducedSolutionOrder=False):
       """
       defines the L{FunctionSpace<escript.FunctionSpace>} of the coefficient

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param reducedEquationOrder: True to indicate that reduced order is used to represent the equation
       @type reducedEquationOrder: C{bool}
       @param reducedSolutionOrder: True to indicate that reduced order is used to represent the solution
       @type reducedSolutionOrder: C{bool}
       @return:  L{FunctionSpace<escript.FunctionSpace>} of the coefficient
       @rtype:  L{FunctionSpace<escript.FunctionSpace>}
       """
       if self.what==self.INTERIOR:
            return escript.Function(domain)
       elif self.what==self.INTERIOR_REDUCED:
            return escript.ReducedFunction(domain)
       elif self.what==self.BOUNDARY:
            return escript.FunctionOnBoundary(domain)
       elif self.what==self.BOUNDARY_REDUCED:
            return escript.ReducedFunctionOnBoundary(domain)
       elif self.what==self.CONTACT:
            return escript.FunctionOnContactZero(domain)
       elif self.what==self.CONTACT_REDUCED:
            return escript.ReducedFunctionOnContactZero(domain)
       elif self.what==self.SOLUTION:
            if reducedEquationOrder and reducedSolutionOrder:
                return escript.ReducedSolution(domain)
            else:
                return escript.Solution(domain)
       elif self.what==self.REDUCED:
            return escript.ReducedSolution(domain)

    def getValue(self):
       """
       returns the value of the coefficient

       @return:  value of the coefficient
       @rtype:  L{Data<escript.Data>}
       """
       return self.value

    def setValue(self,domain,numEquations=1,numSolutions=1,reducedEquationOrder=False,reducedSolutionOrder=False,newValue=None):
       """
       set the value of the coefficient to a new value

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param numEquations: number of equations of the PDE
       @type numEquations: C{int}
       @param numSolutions: number of components of the PDE solution
       @type numSolutions: C{int}
       @param reducedEquationOrder: True to indicate that reduced order is used to represent the equation
       @type reducedEquationOrder: C{bool}
       @param reducedSolutionOrder: True to indicate that reduced order is used to represent the solution
       @type reducedSolutionOrder: C{bool}
       @param newValue: number of components of the PDE solution
       @type newValue: any object that can be converted into a L{Data<escript.Data>} object with the appropriate shape and L{FunctionSpace<escript.FunctionSpace>}
       @raise IllegalCoefficientValue: if the shape of the assigned value does not match the shape of the coefficient
       @raise IllegalCoefficientFunctionSpace: if unable to interploate value to appropriate function space
       """
       if newValue==None:
           newValue=escript.Data()
       elif isinstance(newValue,escript.Data):
           if not newValue.isEmpty():
              if not newValue.getFunctionSpace() == self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder):
                try:
                  newValue=escript.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
                except:
                  raise IllegalCoefficientFunctionSpace,"Unable to interpolate coefficient to function space %s"%self.getFunctionSpace(domain)
       else:
           newValue=escript.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
       if not newValue.isEmpty():
           if not self.getShape(domain,numEquations,numSolutions)==newValue.getShape():
               raise IllegalCoefficientValue,"Expected shape of coefficient is %s but actual shape is %s."%(self.getShape(domain,numEquations,numSolutions),newValue.getShape())
       self.value=newValue

    def isAlteringOperator(self):
        """
        checks if the coefficient alters the operator of the PDE

        @return:  True if the operator of the PDE is changed when the coefficient is changed
        @rtype:  C{bool}
	"""
        if self.altering==self.OPERATOR or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isAlteringRightHandSide(self):
        """
        checks if the coefficeint alters the right hand side of the PDE

	@rtype:  C{bool}
        @return:  True if the right hand side of the PDE is changed when the coefficient is changed
	"""
        if self.altering==self.RIGHTHANDSIDE or self.altering==self.BOTH:
            return not None
        else:
            return None

    def estimateNumEquationsAndNumSolutions(self,domain,shape=()):
       """
       tries to estimate the number of equations and number of solutions if the coefficient has the given shape

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param shape: suggested shape of the coefficient
       @type shape: C{tuple} of C{int} values
       @return: the number of equations and number of solutions of the PDE is the coefficient has shape s.
                 If no appropriate numbers could be identified, C{None} is returned
       @rtype: C{tuple} of two C{int} values or C{None}
       """
       dim=domain.getDim()
       if len(shape)>0:
           num=max(shape)+1
       else:
           num=1
       search=[]
       if self.definesNumEquation() and self.definesNumSolutions():
          for u in range(num):
             for e in range(num):
                search.append((e,u))
          search.sort(self.__CompTuple2)
          for item in search:
             s=self.getShape(domain,item[0],item[1])
             if len(s)==0 and len(shape)==0:
                 return (1,1)
             else:
                 if s==shape: return item
       elif self.definesNumEquation():
          for e in range(num,0,-1):
             s=self.getShape(domain,e,0)
             if len(s)==0 and len(shape)==0:
                 return (1,None)
             else:
                 if s==shape: return (e,None)

       elif self.definesNumSolutions():
          for u in range(num,0,-1):
             s=self.getShape(domain,0,u)
             if len(s)==0 and len(shape)==0:
                 return (None,1)
             else:
                 if s==shape: return (None,u)
       return None
    def definesNumSolutions(self):
       """
       checks if the coefficient allows to estimate the number of solution components

       @return: True if the coefficient allows an estimate of the number of solution components
       @rtype: C{bool}
       """
       for i in self.pattern:
             if i==self.BY_SOLUTION: return True
       return False

    def definesNumEquation(self):
       """
       checks if the coefficient allows to estimate the number of equations

       @return: True if the coefficient allows an estimate of the number of equations
       @rtype: C{bool}
       """
       for i in self.pattern:
             if i==self.BY_EQUATION: return True
       return False

    def __CompTuple2(self,t1,t2):
      """
      Compare two tuples of possible number of equations and number of solutions

      @param t1: The first tuple
      @param t2: The second tuple

      """

      dif=t1[0]+t1[1]-(t2[0]+t2[1])
      if dif<0: return 1
      elif dif>0: return -1
      else: return 0

    def getShape(self,domain,numEquations=1,numSolutions=1):
       """
       builds the required shape of the coefficient

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param numEquations: number of equations of the PDE
       @type numEquations: C{int}
       @param numSolutions: number of components of the PDE solution
       @type numSolutions: C{int}
       @return: shape of the coefficient
       @rtype: C{tuple} of C{int} values
       """
       dim=domain.getDim()
       s=()
       for i in self.pattern:
             if i==self.BY_EQUATION:
                if numEquations>1: s=s+(numEquations,)
             elif i==self.BY_SOLUTION:
                if numSolutions>1: s=s+(numSolutions,)
             else:
                s=s+(dim,)
       return s

class LinearPDE(object):
   """
   This class is used to define a general linear, steady, second order PDE
   for an unknown function M{u} on a given domain defined through a L{Domain<escript.Domain>} object.

   For a single PDE with a solution with a single component the linear PDE is defined in the following form:

   M{-(grad(A[j,l]+A_reduced[j,l])*grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)=-grad(X+X_reduced)[j,j]+(Y+Y_reduced)}


   where M{grad(F)} denotes the spatial derivative of M{F}. Einstein's summation convention,
   ie. summation over indexes appearing twice in a term of a sum is performed, is used.
   The coefficients M{A}, M{B}, M{C}, M{D}, M{X} and M{Y} have to be specified through L{Data<escript.Data>} objects in the
   L{Function<escript.Function>} and the coefficients M{A_reduced}, M{B_reduced}, M{C_reduced}, M{D_reduced}, M{X_reduced} and M{Y_reduced} have to be specified through L{Data<escript.Data>} objects in the
   L{ReducedFunction<escript.ReducedFunction>}. It is also allowd to use objects that can be converted into 
   such L{Data<escript.Data>} objects. M{A} and M{A_reduced} are rank two, M{B_reduced}, M{C_reduced}, M{X_reduced} 
   M{B_reduced}, M{C_reduced} and M{X_reduced} are rank one and M{D}, M{D_reduced} and M{Y_reduced} are scalar.

   The following natural boundary conditions are considered:

   M{n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u)+(d+d_reduced)*u=n[j]*(X[j]+X_reduced[j])+y}

   where M{n} is the outer normal field. Notice that the coefficients M{A}, M{A_reduced}, M{B}, M{B_reduced}, M{X} and M{X_reduced} are defined in the PDE. The coefficients M{d} and M{y} and are each a scalar in the L{FunctionOnBoundary<escript.FunctionOnBoundary>} and the coefficients M{d_reduced} and M{y_reduced} and are each a scalar in the L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}.


   Constraints for the solution prescribing the value of the solution at certain locations in the domain. They have the form

   M{u=r}  where M{q>0}

   M{r} and M{q} are each scalar where M{q} is the characteristic function defining where the constraint is applied.
   The constraints override any other condition set by the PDE or the boundary condition.

   The PDE is symmetrical if

   M{A[i,j]=A[j,i]}  and M{B[j]=C[j]} and M{A_reduced[i,j]=A_reduced[j,i]}  and M{B_reduced[j]=C_reduced[j] 

   For a system of PDEs and a solution with several components the PDE has the form

   M{-grad((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])[j]+(C[i,k,l]+C_reduced[i,k,l])*grad(u[k])[l]+(D[i,k]+D_reduced[i,k]*u[k] =-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i] }

   M{A} and M{A_reduced} are of rank four, M{B}, M{B_reduced}, M{C} and M{C_reduced} are each of rank three, M{D}, M{D_reduced}, M{X_reduced} and M{X} are each a rank two and M{Y} and M{Y_reduced} are of rank one.
   The natural boundary conditions take the form:

   M{n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])+(d[i,k]+d_reduced[i,k])*u[k]=n[j]*(X[i,j]+X_reduced[i,j])+y[i]+y_reduced[i]}


   The coefficient M{d} is a rank two and M{y} is a rank one both in the L{FunctionOnBoundary<escript.FunctionOnBoundary>}. Constraints take the form and the coefficients M{d_reduced} is a rank two and M{y_reduced} is a rank one both in the L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}. 

   Constraints take the form 

   M{u[i]=r[i]}  where  M{q[i]>0}

   M{r} and M{q} are each rank one. Notice that at some locations not necessarily all components must have a constraint.

   The system of PDEs is symmetrical if

        - M{A[i,j,k,l]=A[k,l,i,j]}
        - M{A_reduced[i,j,k,l]=A_reduced[k,l,i,j]}
        - M{B[i,j,k]=C[k,i,j]}
        - M{B_reduced[i,j,k]=C_reduced[k,i,j]}
        - M{D[i,k]=D[i,k]}
        - M{D_reduced[i,k]=D_reduced[i,k]}
        - M{d[i,k]=d[k,i]}
        - M{d_reduced[i,k]=d_reduced[k,i]}

   L{LinearPDE} also supports solution discontinuities over a contact region in the domain. To specify the conditions across the
   discontinuity we are using the generalised flux M{J} which is in the case of a systems of PDEs and several components of the solution
   defined as

   M{J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]-X[i,j]-X_reduced[i,j]}

   For the case of single solution component and single PDE M{J} is defined

   M{J_{j}=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u-X[i]-X_reduced[i]}

   In the context of discontinuities M{n} denotes the normal on the discontinuity pointing from side 0 towards side 1
   calculated from L{getNormal<escript.FunctionSpace.getNormal>} of L{FunctionOnContactZero<escript.FunctionOnContactZero>}. For a system of PDEs
   the contact condition takes the form

   M{n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]}

   where M{J0} and M{J1} are the fluxes on side 0 and side 1 of the discontinuity, respectively. M{jump(u)}, which is the difference
   of the solution at side 1 and at side 0, denotes the jump of M{u} across discontinuity along the normal calcualted by
   L{jump<util.jump>}.
   The coefficient M{d_contact} is a rank two and M{y_contact} is a rank one both in the L{FunctionOnContactZero<escript.FunctionOnContactZero>} or L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
    The coefficient M{d_contact_reduced} is a rank two and M{y_contact_reduced} is a rank one both in the L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>} or L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}.
   In case of a single PDE and a single component solution the contact condition takes the form

   M{n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)}

   In this case the coefficient M{d_contact} and M{y_contact} are each scalar both in the L{FunctionOnContactZero<escript.FunctionOnContactZero>} or L{FunctionOnContactOne<escript.FunctionOnContactOne>} and the coefficient M{d_contact_reduced} and M{y_contact_reduced} are each scalar both in the L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>} or L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}

   @cvar DEFAULT: The default method used to solve the system of linear equations
   @cvar DIRECT: The direct solver based on LDU factorization
   @cvar CHOLEVSKY: The direct solver based on LDLt factorization (can only be applied for symmetric PDEs)
   @cvar PCG: The preconditioned conjugate gradient method (can only be applied for symmetric PDEs)
   @cvar CR: The conjugate residual method
   @cvar CGS: The conjugate gardient square method
   @cvar BICGSTAB: The stabilized BiConjugate Gradient method.
   @cvar SSOR: The symmetric overrealaxtion method
   @cvar ILU0: The incomplete LU factorization preconditioner  with no fill in
   @cvar ILUT: The incomplete LU factorization preconditioner with will in
   @cvar JACOBI: The Jacobi preconditioner
   @cvar GMRES: The Gram-Schmidt minimum residual method
   @cvar PRES20: Special GMRES with restart after 20 steps and truncation after 5 residuals
   @cvar LUMPING: Matrix lumping.
   @cvar NO_REORDERING: No matrix reordering allowed
   @cvar MINIMUM_FILL_IN: Reorder matrix to reduce fill-in during factorization
   @cvar NESTED_DISSECTION: Reorder matrix to improve load balancing during factorization
   @cvar PASO: PASO solver package
   @cvar SCSL: SGI SCSL solver library
   @cvar MKL: Intel's MKL solver library
   @cvar UMFPACK: the UMFPACK library
   @cvar ITERATIVE: The default iterative solver
   @cvar AMG: algebraic multi grid
   @cvar RILU: recursive ILU

   """
   DEFAULT= 0
   DIRECT= 1
   CHOLEVSKY= 2
   PCG= 3
   CR= 4
   CGS= 5
   BICGSTAB= 6
   SSOR= 7
   ILU0= 8
   ILUT= 9
   JACOBI= 10
   GMRES= 11
   PRES20= 12
   LUMPING= 13
   NO_REORDERING= 17
   MINIMUM_FILL_IN= 18
   NESTED_DISSECTION= 19
   SCSL= 14
   MKL= 15
   UMFPACK= 16
   ITERATIVE= 20
   PASO= 21
   AMG= 22
   RILU = 23

   SMALL_TOLERANCE=1.e-13
   __PACKAGE_KEY="package"
   __METHOD_KEY="method"
   __SYMMETRY_KEY="symmetric"
   __TOLERANCE_KEY="tolerance"
   __PRECONDITIONER_KEY="preconditioner"


   def __init__(self,domain,numEquations=None,numSolutions=None,debug=False):
     """
     initializes a new linear PDE

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param numEquations: number of equations. If numEquations==None the number of equations
                          is exracted from the PDE coefficients.
     @param numSolutions: number of solution components. If  numSolutions==None the number of solution components
                          is exracted from the PDE coefficients.
     @param debug: if True debug informations are printed.

     """
     super(LinearPDE, self).__init__()
     #
     #   the coefficients of the general PDE:
     #
     self.__COEFFICIENTS_OF_GENEARL_PDE={
       "A"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM,PDECoefficient.BY_SOLUTION,PDECoefficient.BY_DIM),PDECoefficient.OPERATOR),
       "B"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "C"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION,PDECoefficient.BY_DIM),PDECoefficient.OPERATOR),
       "D"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "X"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM),PDECoefficient.RIGHTHANDSIDE),
       "Y"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "d"         : PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "y"         : PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "d_contact" : PDECoefficient(PDECoefficient.CONTACT,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "y_contact" : PDECoefficient(PDECoefficient.CONTACT,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "A_reduced"         : PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM,PDECoefficient.BY_SOLUTION,PDECoefficient.BY_DIM),PDECoefficient.OPERATOR),
       "B_reduced"         : PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "C_reduced"         : PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION,PDECoefficient.BY_DIM),PDECoefficient.OPERATOR),
       "D_reduced"         : PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "X_reduced"         : PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM),PDECoefficient.RIGHTHANDSIDE),
       "Y_reduced"         : PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "d_reduced"         : PDECoefficient(PDECoefficient.BOUNDARY_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "y_reduced"         : PDECoefficient(PDECoefficient.BOUNDARY_REDUCED,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "d_contact_reduced" : PDECoefficient(PDECoefficient.CONTACT_REDUCED,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_SOLUTION),PDECoefficient.OPERATOR),
       "y_contact_reduced" : PDECoefficient(PDECoefficient.CONTACT_REDUCED,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "r"         : PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_SOLUTION,),PDECoefficient.RIGHTHANDSIDE),
       "q"         : PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_SOLUTION,),PDECoefficient.BOTH)}

     # COEFFICIENTS can be overwritten by subclasses:
     self.COEFFICIENTS=self.__COEFFICIENTS_OF_GENEARL_PDE
     self.__altered_coefficients=False
     # initialize attributes
     self.__debug=debug
     self.__domain=domain
     self.__numEquations=numEquations
     self.__numSolutions=numSolutions
     self.__resetSystem()

     # set some default values:
     self.__reduce_equation_order=False
     self.__reduce_solution_order=False
     self.__tolerance=1.e-8
     self.__solver_method=self.DEFAULT
     self.__solver_package=self.DEFAULT
     self.__preconditioner=self.DEFAULT
     self.__matrix_type=self.__domain.getSystemMatrixTypeId(self.DEFAULT,self.DEFAULT,False)
     self.__sym=False

     self.resetCoefficients()
     self.trace("PDE Coeffients are %s"%str(self.COEFFICIENTS.keys()))
   # =============================================================================
   #    general stuff:
   # =============================================================================
   def __str__(self):
     """
     returns string representation of the PDE

     @return: a simple representation of the PDE
     @rtype: C{str}
     """
     return "<LinearPDE %d>"%id(self)
   # =============================================================================
   #    debug :
   # =============================================================================
   def setDebugOn(self):
     """
     switches on debugging
     """
     self.__debug=not None

   def setDebugOff(self):
     """
     switches off debugging
     """
     self.__debug=None

   def trace(self,text):
     """
     print the text message if debugging is swiched on.
     @param text: message
     @type text: C{string}
     """
     if self.__debug: print "%s: %s"%(str(self),text)

   # =============================================================================
   # some service functions:
   # =============================================================================
   def getDomain(self):
     """
     returns the domain of the PDE

     @return: the domain of the PDE
     @rtype: L{Domain<escript.Domain>}
     """
     return self.__domain

   def getDim(self):
     """
     returns the spatial dimension of the PDE

     @return: the spatial dimension of the PDE domain
     @rtype: C{int}
     """
     return self.getDomain().getDim()

   def getNumEquations(self):
     """
     returns the number of equations

     @return: the number of equations
     @rtype: C{int}
     @raise UndefinedPDEError: if the number of equations is not be specified yet.
     """
     if self.__numEquations==None:
         raise UndefinedPDEError,"Number of equations is undefined. Please specify argument numEquations."
     else:
         return self.__numEquations

   def getNumSolutions(self):
     """
     returns the number of unknowns

     @return: the number of unknowns
     @rtype: C{int}
     @raise UndefinedPDEError: if the number of unknowns is not be specified yet.
     """
     if self.__numSolutions==None:
        raise UndefinedPDEError,"Number of solution is undefined. Please specify argument numSolutions."
     else:
        return self.__numSolutions

   def reduceEquationOrder(self):
     """
     return status for order reduction for equation

     @return: return True is reduced interpolation order is used for the represenation of the equation
     @rtype: L{bool}
     """
     return self.__reduce_equation_order

   def reduceSolutionOrder(self):
     """
     return status for order reduction for the solution

     @return: return True is reduced interpolation order is used for the represenation of the solution
     @rtype: L{bool}
     """
     return self.__reduce_solution_order

   def getFunctionSpaceForEquation(self):
     """
     returns the L{FunctionSpace<escript.FunctionSpace>} used to discretize the equation

     @return: representation space of equation
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     """
     if self.reduceEquationOrder():
         return escript.ReducedSolution(self.getDomain())
     else:
         return escript.Solution(self.getDomain())

   def getFunctionSpaceForSolution(self):
     """
     returns the L{FunctionSpace<escript.FunctionSpace>} used to represent the solution

     @return: representation space of solution
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     """
     if self.reduceSolutionOrder():
         return escript.ReducedSolution(self.getDomain())
     else:
         return escript.Solution(self.getDomain())


   def getOperator(self):
     """
     provides access to the operator of the PDE

     @return: the operator of the PDE
     @rtype: L{Operator<escript.Operator>}
     """
     m=self.getSystem()[0]
     if self.isUsingLumping():
         return self.copyConstraint(1./m)
     else:
         return m

   def getRightHandSide(self):
     """
     provides access to the right hand side of the PDE
     @return: the right hand side of the PDE
     @rtype: L{Data<escript.Data>}
     """
     r=self.getSystem()[1]
     if self.isUsingLumping():
         return self.copyConstraint(r)
     else:
         return r

   def applyOperator(self,u=None):
     """
     applies the operator of the PDE to a given u or the solution of PDE if u is not present.

     @param u: argument of the operator. It must be representable in C{elf.getFunctionSpaceForSolution()}. If u is not present or equals L{None}
               the current solution is used.
     @type u: L{Data<escript.Data>} or None
     @return: image of u
     @rtype: L{Data<escript.Data>}
     """
     if u==None:
        return self.getOperator()*self.getSolution()
     else:
        return self.getOperator()*escript.Data(u,self.getFunctionSpaceForSolution())

   def getResidual(self,u=None):
     """
     return the residual of u or the current solution if u is not present.

     @param u: argument in the residual calculation. It must be representable in C{elf.getFunctionSpaceForSolution()}. If u is not present or equals L{None}
               the current solution is used.
     @type u: L{Data<escript.Data>} or None
     @return: residual of u
     @rtype: L{Data<escript.Data>}
     """
     return self.applyOperator(u)-self.getRightHandSide()

   def checkSymmetry(self,verbose=True):
      """
      test the PDE for symmetry.

      @param verbose: if equal to True or not present a report on coefficients which are breaking the symmetry is printed.
      @type verbose: C{bool}
      @return:  True if the PDE is symmetric.
      @rtype: L{Data<escript.Data>}
      @note: This is a very expensive operation. It should be used for degugging only! The symmetry flag is not altered.
      """
      verbose=verbose or self.__debug
      out=True
      if self.getNumSolutions()!=self.getNumEquations():
         if verbose: print "non-symmetric PDE because of different number of equations and solutions"
         out=False
      else:
         A=self.getCoefficientOfGeneralPDE("A")
         if not A.isEmpty():
            tol=util.Lsup(A)*self.SMALL_TOLERANCE
            if self.getNumSolutions()>1:
               for i in range(self.getNumEquations()):
                  for j in range(self.getDim()):
                     for k in range(self.getNumSolutions()):
                        for l in range(self.getDim()):
                            if util.Lsup(A[i,j,k,l]-A[k,l,i,j])>tol:
                               if verbose: print "non-symmetric PDE because A[%d,%d,%d,%d]!=A[%d,%d,%d,%d]"%(i,j,k,l,k,l,i,j)
                               out=False
            else:
               for j in range(self.getDim()):
                  for l in range(self.getDim()):
                     if util.Lsup(A[j,l]-A[l,j])>tol:
                        if verbose: print "non-symmetric PDE because A[%d,%d]!=A[%d,%d]"%(j,l,l,j)
                        out=False
         B=self.getCoefficientOfGeneralPDE("B")
         C=self.getCoefficientOfGeneralPDE("C")
         if B.isEmpty() and not C.isEmpty():
            if verbose: print "non-symmetric PDE because B is not present but C is"
            out=False
         elif not B.isEmpty() and C.isEmpty():
            if verbose: print "non-symmetric PDE because C is not present but B is"
            out=False
         elif not B.isEmpty() and not C.isEmpty():
            tol=(util.Lsup(B)+util.Lsup(C))*self.SMALL_TOLERANCE/2.
            if self.getNumSolutions()>1:
               for i in range(self.getNumEquations()):
                   for j in range(self.getDim()):
                      for k in range(self.getNumSolutions()):
                         if util.Lsup(B[i,j,k]-C[k,i,j])>tol:
                              if verbose: print "non-symmetric PDE because B[%d,%d,%d]!=C[%d,%d,%d]"%(i,j,k,k,i,j)
                              out=False
            else:
               for j in range(self.getDim()):
                  if util.Lsup(B[j]-C[j])>tol:
                     if verbose: print "non-symmetric PDE because B[%d]!=C[%d]"%(j,j)
                     out=False
         if self.getNumSolutions()>1:
           D=self.getCoefficientOfGeneralPDE("D")
           if not D.isEmpty():
             tol=util.Lsup(D)*self.SMALL_TOLERANCE
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(D[i,k]-D[k,i])>tol:
                      if verbose: print "non-symmetric PDE because D[%d,%d]!=D[%d,%d]"%(i,k,k,i)
                      out=False
           d=self.getCoefficientOfGeneralPDE("d")
           if not d.isEmpty():
             tol=util.Lsup(d)*self.SMALL_TOLERANCE
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(d[i,k]-d[k,i])>tol:
                      if verbose: print "non-symmetric PDE because d[%d,%d]!=d[%d,%d]"%(i,k,k,i)
                      out=False
           d_contact=self.getCoefficientOfGeneralPDE("d_contact")
           if not d_contact.isEmpty():
             tol=util.Lsup(d_contact)*self.SMALL_TOLERANCE
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(d_contact[i,k]-d_contact[k,i])>tol:
                      if verbose: print "non-symmetric PDE because d_contact[%d,%d]!=d_contact[%d,%d]"%(i,k,k,i)
                      out=False
         # and now the reduced coefficients
         A_reduced=self.getCoefficientOfGeneralPDE("A_reduced")
         if not A_reduced.isEmpty():
            tol=util.Lsup(A_reduced)*self.SMALL_TOLERANCE
            if self.getNumSolutions()>1:
               for i in range(self.getNumEquations()):
                  for j in range(self.getDim()):
                     for k in range(self.getNumSolutions()):
                        for l in range(self.getDim()):
                            if util.Lsup(A_reduced[i,j,k,l]-A_reduced[k,l,i,j])>tol:
                               if verbose: print "non-symmetric PDE because A_reduced[%d,%d,%d,%d]!=A_reduced[%d,%d,%d,%d]"%(i,j,k,l,k,l,i,j)
                               out=False
            else:
               for j in range(self.getDim()):
                  for l in range(self.getDim()):
                     if util.Lsup(A_reduced[j,l]-A_reduced[l,j])>tol:
                        if verbose: print "non-symmetric PDE because A_reduced[%d,%d]!=A_reduced[%d,%d]"%(j,l,l,j)
                        out=False
         B_reduced=self.getCoefficientOfGeneralPDE("B_reduced")
         C_reduced=self.getCoefficientOfGeneralPDE("C_reduced")
         if B_reduced.isEmpty() and not C_reduced.isEmpty():
            if verbose: print "non-symmetric PDE because B_reduced is not present but C_reduced is"
            out=False
         elif not B_reduced.isEmpty() and C_reduced.isEmpty():
            if verbose: print "non-symmetric PDE because C_reduced is not present but B_reduced is"
            out=False
         elif not B_reduced.isEmpty() and not C_reduced.isEmpty():
            tol=(util.Lsup(B_reduced)+util.Lsup(C_reduced))*self.SMALL_TOLERANCE/2.
            if self.getNumSolutions()>1:
               for i in range(self.getNumEquations()):
                   for j in range(self.getDim()):
                      for k in range(self.getNumSolutions()):
                         if util.Lsup(B_reduced[i,j,k]-C_reduced[k,i,j])>tol:
                              if verbose: print "non-symmetric PDE because B_reduced[%d,%d,%d]!=C_reduced[%d,%d,%d]"%(i,j,k,k,i,j)
                              out=False
            else:
               for j in range(self.getDim()):
                  if util.Lsup(B_reduced[j]-C_reduced[j])>tol:
                     if verbose: print "non-symmetric PDE because B_reduced[%d]!=C_reduced[%d]"%(j,j)
                     out=False
         if self.getNumSolutions()>1:
           D_reduced=self.getCoefficientOfGeneralPDE("D_reduced")
           if not D_reduced.isEmpty():
             tol=util.Lsup(D_reduced)*self.SMALL_TOLERANCE
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(D_reduced[i,k]-D_reduced[k,i])>tol:
                      if verbose: print "non-symmetric PDE because D_reduced[%d,%d]!=D_reduced[%d,%d]"%(i,k,k,i)
                      out=False
           d_reduced=self.getCoefficientOfGeneralPDE("d_reduced")
           if not d_reduced.isEmpty():
             tol=util.Lsup(d_reduced)*self.SMALL_TOLERANCE
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(d_reduced[i,k]-d_reduced[k,i])>tol:
                      if verbose: print "non-symmetric PDE because d_reduced[%d,%d]!=d_reduced[%d,%d]"%(i,k,k,i)
                      out=False
           d_contact_reduced=self.getCoefficientOfGeneralPDE("d_contact_reduced")
           if not d_contact_reduced.isEmpty():
             tol=util.Lsup(d_contact_reduced)*self.SMALL_TOLERANCE
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(d_contact_reduced[i,k]-d_contact_reduced[k,i])>tol:
                      if verbose: print "non-symmetric PDE because d_contact_reduced[%d,%d]!=d_contact_reduced[%d,%d]"%(i,k,k,i)
                      out=False
      return out

   def getSolution(self,**options):
       """
       returns the solution of the PDE. If the solution is not valid the PDE is solved.

       @return: the solution
       @rtype: L{Data<escript.Data>}
       @param options: solver options
       @keyword verbose: True to get some information during PDE solution
       @type verbose: C{bool}
       @keyword reordering: reordering scheme to be used during elimination. Allowed values are
                            L{NO_REORDERING}, L{MINIMUM_FILL_IN}, L{NESTED_DISSECTION}
       @keyword iter_max: maximum number of iteration steps allowed.
       @keyword drop_tolerance: threshold for drupping in L{ILUT}
       @keyword drop_storage: maximum of allowed memory in L{ILUT}
       @keyword truncation: maximum number of residuals in L{GMRES}
       @keyword restart: restart cycle length in L{GMRES}
       """
       if not self.__solution_isValid:
          mat,f=self.getSystem()
          if self.isUsingLumping():
             self.__solution=self.copyConstraint(f*mat)
          else:
             options[self.__TOLERANCE_KEY]=self.getTolerance()
             options[self.__METHOD_KEY]=self.getSolverMethod()[0]
             options[self.__PRECONDITIONER_KEY]=self.getSolverMethod()[1]
             options[self.__PACKAGE_KEY]=self.getSolverPackage()
             options[self.__SYMMETRY_KEY]=self.isSymmetric()
             self.trace("PDE is resolved.")
             self.trace("solver options: %s"%str(options))
             self.__solution=mat.solve(f,options)
          self.__solution_isValid=True
       return self.__solution

   def getFlux(self,u=None):
     """
     returns the flux M{J} for a given M{u}

     M{J[i,j]=(A[i,j,k,l]+A_reduced[A[i,j,k,l]]*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])u[k]-X[i,j]-X_reduced[i,j]}

     or

     M{J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[l]+(B[j]+B_reduced[j])u-X[j]-X_reduced[j]}

     @param u: argument in the flux. If u is not present or equals L{None} the current solution is used.
     @type u: L{Data<escript.Data>} or None
     @return: flux
     @rtype: L{Data<escript.Data>}
     """
     if u==None: u=self.getSolution()
     return util.tensormult(self.getCoefficientOfGeneralPDE("A"),util.grad(u,Funtion(self.getDomain))) \
           +util.matrixmult(self.getCoefficientOfGeneralPDE("B"),u) \
           -util.self.getCoefficientOfGeneralPDE("X") \
           +util.tensormult(self.getCoefficientOfGeneralPDE("A_reduced"),util.grad(u,ReducedFuntion(self.getDomain))) \
           +util.matrixmult(self.getCoefficientOfGeneralPDE("B_reduced"),u) \
           -util.self.getCoefficientOfGeneralPDE("X_reduced")
   # =============================================================================
   #   solver settings:
   # =============================================================================
   def setSolverMethod(self,solver=None,preconditioner=None):
       """
       sets a new solver

       @param solver: sets a new solver method.
       @type solver: one of L{DEFAULT}, L{ITERATIVE} L{DIRECT}, L{CHOLEVSKY}, L{PCG}, L{CR}, L{CGS}, L{BICGSTAB}, L{SSOR}, L{GMRES}, L{PRES20}, L{LUMPING}, L{AMG}
       @param preconditioner: sets a new solver method.
       @type preconditioner: one of L{DEFAULT}, L{JACOBI} L{ILU0}, L{ILUT},L{SSOR}, L{RILU}
       """
       if solver==None: solve=self.DEFAULT
       if preconditioner==None: preconditioner=self.DEFAULT
       if not (solver,preconditioner)==self.getSolverMethod():
           self.__solver_method=solver
           self.__preconditioner=preconditioner
           self.__checkMatrixType()
           self.trace("New solver is %s"%self.getSolverMethodName())

   def getSolverMethodName(self):
       """
       returns the name of the solver currently used

       @return: the name of the solver currently used.
       @rtype: C{string}
       """

       m=self.getSolverMethod()
       p=self.getSolverPackage()
       method=""
       if m[0]==self.DEFAULT: method="DEFAULT"
       elif m[0]==self.DIRECT: method= "DIRECT"
       elif m[0]==self.ITERATIVE: method= "ITERATIVE"
       elif m[0]==self.CHOLEVSKY: method= "CHOLEVSKY"
       elif m[0]==self.PCG: method= "PCG"
       elif m[0]==self.CR: method= "CR"
       elif m[0]==self.CGS: method= "CGS"
       elif m[0]==self.BICGSTAB: method= "BICGSTAB"
       elif m[0]==self.SSOR: method= "SSOR"
       elif m[0]==self.GMRES: method= "GMRES"
       elif m[0]==self.PRES20: method= "PRES20"
       elif m[0]==self.LUMPING: method= "LUMPING"
       elif m[0]==self.AMG: method= "AMG"
       if m[1]==self.DEFAULT: method+="+DEFAULT"
       elif m[1]==self.JACOBI: method+= "+JACOBI"
       elif m[1]==self.ILU0: method+= "+ILU0"
       elif m[1]==self.ILUT: method+= "+ILUT"
       elif m[1]==self.SSOR: method+= "+SSOR"
       elif m[1]==self.AMG: method+= "+AMG"
       elif m[1]==self.RILU: method+= "+RILU"
       if p==self.DEFAULT: package="DEFAULT"
       elif p==self.PASO: package= "PASO"
       elif p==self.MKL: package= "MKL"
       elif p==self.SCSL: package= "SCSL"
       elif p==self.UMFPACK: package= "UMFPACK"
       else : method="unknown"
       return "%s solver of %s package"%(method,package)


   def getSolverMethod(self):
       """
       returns the solver method

       @return: the solver method currently be used.
       @rtype: C{int}
       """
       return self.__solver_method,self.__preconditioner

   def setSolverPackage(self,package=None):
       """
       sets a new solver package

       @param package: sets a new solver method.
       @type package: one of L{DEFAULT}, L{PASO} L{SCSL}, L{MKL}, L{UMFPACK}
       """
       if package==None: package=self.DEFAULT
       if not package==self.getSolverPackage():
           self.__solver_package=package
           self.__checkMatrixType()
           self.trace("New solver is %s"%self.getSolverMethodName())

   def getSolverPackage(self):
       """
       returns the package of the solver

       @return: the solver package currently being used.
       @rtype: C{int}
       """
       return self.__solver_package

   def isUsingLumping(self):
      """
      checks if matrix lumping is used a solver method

      @return: True is lumping is currently used a solver method.
      @rtype: C{bool}
      """
      return self.getSolverMethod()[0]==self.LUMPING

   def setTolerance(self,tol=1.e-8):
       """
       resets the tolerance for the solver method to tol where for an appropriate norm M{|.|}

       M{|L{getResidual}()|<tol*|L{getRightHandSide}()|}

       defines the stopping criterion.

       @param tol: new tolerance for the solver. If the tol is lower then the current tolerence
                   the system will be resolved.
       @type tol: positive C{float}
       @raise ValueError: if tolerance is not positive.
       """
       if not tol>0:
           raise ValueError,"Tolerance as to be positive"
       if tol<self.getTolerance(): self.__invalidateSolution()
       self.trace("New tolerance %e"%tol)
       self.__tolerance=tol
       return

   def getTolerance(self):
       """
       returns the tolerance set for the solution

       @return: tolerance currently used.
       @rtype: C{float}
       """
       return self.__tolerance

   # =============================================================================
   #    symmetry  flag:
   # =============================================================================
   def isSymmetric(self):
      """
      checks if symmetry is indicated.

      @return: True is a symmetric PDE is indicated, otherwise False is returned
      @rtype: C{bool}
      """
      return self.__sym

   def setSymmetryOn(self):
      """
      sets the symmetry flag.
      """
      if not self.isSymmetric():
         self.trace("PDE is set to be symmetric")
         self.__sym=True
         self.__checkMatrixType()

   def setSymmetryOff(self):
      """
      removes the symmetry flag.
      """
      if self.isSymmetric():
         self.trace("PDE is set to be unsymmetric")
         self.__sym=False
         self.__checkMatrixType()

   def setSymmetryTo(self,flag=False):
      """
      sets the symmetry flag to flag

      @param flag: If flag, the symmetry flag is set otherwise the symmetry flag is released.
      @type flag: C{bool}
      """
      if flag:
         self.setSymmetryOn()
      else:
         self.setSymmetryOff()

   # =============================================================================
   # function space handling for the equation as well as the solution
   # =============================================================================
   def setReducedOrderOn(self):
     """
     switches on reduced order for solution and equation representation

     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     switches off reduced order for solution and equation representation

     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     sets order reduction for both solution and equation representation according to flag.
     @param flag: if flag is True, the order reduction is switched on for both  solution and equation representation, otherwise or
                  if flag is not present order reduction is switched off
     @type flag: C{bool}
     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)


   def setReducedOrderForSolutionOn(self):
     """
     switches on reduced order for solution representation

     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     if not self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Reduced order is used to solution representation.")
         self.__reduce_solution_order=True
         self.__resetSystem()

   def setReducedOrderForSolutionOff(self):
     """
     switches off reduced order for solution representation

     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     if self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Full order is used to interpolate solution.")
         self.__reduce_solution_order=False
         self.__resetSystem()

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     sets order for test functions according to flag

     @param flag: if flag is True, the order reduction is switched on for solution representation, otherwise or
                  if flag is not present order reduction is switched off
     @type flag: C{bool}
     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()

   def setReducedOrderForEquationOn(self):
     """
     switches on reduced order for equation representation

     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     if not self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Reduced order is used for test functions.")
         self.__reduce_equation_order=True
         self.__resetSystem()

   def setReducedOrderForEquationOff(self):
     """
     switches off reduced order for equation representation

     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     if self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Full order is used for test functions.")
         self.__reduce_equation_order=False
         self.__resetSystem()

   def setReducedOrderForEquationTo(self,flag=False):
     """
     sets order for test functions according to flag

     @param flag: if flag is True, the order reduction is switched on for equation representation, otherwise or
                  if flag is not present order reduction is switched off
     @type flag: C{bool}
     @raise RuntimeError: if order reduction is altered after a coefficient has been set.
     """
     if flag:
        self.setReducedOrderForEquationOn()
     else:
        self.setReducedOrderForEquationOff()

   # =============================================================================
   # private method:
   # =============================================================================
   def __checkMatrixType(self):
     """
     reassess the matrix type and, if a new matrix is needed, resets the system.
     """
     new_matrix_type=self.getDomain().getSystemMatrixTypeId(self.getSolverMethod()[0],self.getSolverPackage(),self.isSymmetric())
     if not new_matrix_type==self.__matrix_type:
         self.trace("Matrix type is now %d."%new_matrix_type)
         self.__matrix_type=new_matrix_type
         self.__resetSystem()
   #
   #   rebuild switches :
   #
   def __invalidateSolution(self):
       """
       indicates the PDE has to be resolved if the solution is requested
       """
       if self.__solution_isValid: self.trace("PDE has to be resolved.")
       self.__solution_isValid=False

   def __invalidateOperator(self):
       """
       indicates the operator has to be rebuilt next time it is used
       """
       if self.__operator_is_Valid: self.trace("Operator has to be rebuilt.")
       self.__invalidateSolution()
       self.__operator_is_Valid=False

   def __invalidateRightHandSide(self):
       """
       indicates the right hand side has to be rebuild next time it is used
       """
       if self.__righthandside_isValid: self.trace("Right hand side has to be rebuilt.")
       self.__invalidateSolution()
       self.__righthandside_isValid=False

   def __invalidateSystem(self):
       """
       annonced that everthing has to be rebuild:
       """
       if self.__righthandside_isValid: self.trace("System has to be rebuilt.")
       self.__invalidateSolution()
       self.__invalidateOperator()
       self.__invalidateRightHandSide()

   def __resetSystem(self):
       """
       annonced that everthing has to be rebuild:
       """
       self.trace("New System is built from scratch.")
       self.__operator=escript.Operator()
       self.__operator_is_Valid=False
       self.__righthandside=escript.Data()
       self.__righthandside_isValid=False
       self.__solution=escript.Data()
       self.__solution_isValid=False
   #
   #    system initialization:
   #
   def __getNewOperator(self):
       """
       returns an instance of a new operator
       """
       self.trace("New operator is allocated.")
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           self.__matrix_type)

   def __getNewRightHandSide(self):
       """
       returns an instance of a new right hand side
       """
       self.trace("New right hand side is allocated.")
       if self.getNumEquations()>1:
           return escript.Data(0.,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),True)
       else:
           return escript.Data(0.,(),self.getFunctionSpaceForEquation(),True)

   def __getNewSolution(self):
       """
       returns an instance of a new solution
       """
       self.trace("New solution is allocated.")
       if self.getNumSolutions()>1:
           return escript.Data(0.,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)
       else:
           return escript.Data(0.,(),self.getFunctionSpaceForSolution(),True)

   def __makeFreshSolution(self):
       """
       makes sure that the solution is instantiated and returns it initialized by zeros
       """
       if self.__solution.isEmpty():
           self.__solution=self.__getNewSolution()
       else:
           self.__solution*=0
           self.trace("Solution is reset to zero.")
       return self.__solution

   def __makeFreshRightHandSide(self):
       """
       makes sure that the right hand side is instantiated and returns it initialized by zeros
       """
       if self.__righthandside.isEmpty():
           self.__righthandside=self.__getNewRightHandSide()
       else:
           self.__righthandside.setToZero()
           self.trace("Right hand side is reset to zero.")
       return self.__righthandside

   def __makeFreshOperator(self):
       """
       makes sure that the operator is instantiated and returns it initialized by zeros
       """
       if self.__operator.isEmpty():
           self.__operator=self.__getNewOperator()
       else:
           self.__operator.resetValues()
           self.trace("Operator reset to zero")
       return self.__operator

   def __applyConstraint(self):
       """
       applies the constraints defined by q and r to the system
       """
       if not self.isUsingLumping():
          q=self.getCoefficientOfGeneralPDE("q")
          r=self.getCoefficientOfGeneralPDE("r")
          if not q.isEmpty() and not self.__operator.isEmpty():
             # q is the row and column mask to indicate where constraints are set:
             row_q=escript.Data(q,self.getFunctionSpaceForEquation())
             col_q=escript.Data(q,self.getFunctionSpaceForSolution())
             u=self.__getNewSolution()
             if r.isEmpty():
                r_s=self.__getNewSolution()
             else:
                r_s=escript.Data(r,self.getFunctionSpaceForSolution())
             u.copyWithMask(r_s,col_q)
             if not self.__righthandside.isEmpty():
                self.__righthandside-=self.__operator*u
                self.__righthandside=self.copyConstraint(self.__righthandside)
             self.__operator.nullifyRowsAndCols(row_q,col_q,1.)
   # =============================================================================
   # function giving access to coefficients of the general PDE:
   # =============================================================================
   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE.

     @note: This method is called by the assembling routine it can be overwritten
           to map coefficients of a particular PDE to the general PDE.
     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: the value of the coefficient  name
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, 
                  M{A_reduced}, M{B_reduced}, M{C_reduced}, M{D_reduced}, M{X_reduced}, M{Y_reduced}, M{d_reduced}, M{y_reduced}, M{d_contact_reduced}, M{y_contact_reduced}, M{r} or M{q}.
     """
     if self.hasCoefficientOfGeneralPDE(name):
        return self.getCoefficient(name)
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def hasCoefficientOfGeneralPDE(self,name):
     """
     checks if name is a the name of a coefficient of the general PDE.

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return: True if name is the name of a coefficient of the general PDE. Otherwise False.
     @rtype: C{bool}

     """
     return self.__COEFFICIENTS_OF_GENEARL_PDE.has_key(name)

   def createCoefficientOfGeneralPDE(self,name):
     """
     returns a new instance of a coefficient for coefficient name of the general PDE

     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: a coefficient name initialized to 0.
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, 
                  M{A_reduced}, M{B_reduced}, M{C_reduced}, M{D_reduced}, M{X_reduced}, M{Y_reduced}, M{d_reduced}, M{y_reduced}, M{d_contact_reduced}, M{y_contact_reduced}, M{r} or M{q}.
     """
     if self.hasCoefficientOfGeneralPDE(name):
        return escript.Data(0,self.getShapeOfCoefficientOfGeneralPDE(name),self.getFunctionSpaceForCoefficientOfGeneralPDE(name))
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def getFunctionSpaceForCoefficientOfGeneralPDE(self,name):
     """
     return the L{FunctionSpace<escript.FunctionSpace>} to be used for coefficient name of the general PDE

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return: the function space to be used for coefficient name
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, 
                  M{A_reduced}, M{B_reduced}, M{C_reduced}, M{D_reduced}, M{X_reduced}, M{Y_reduced}, M{d_reduced}, M{y_reduced}, M{d_contact_reduced}, M{y_contact_reduced}, M{r} or M{q}.
     """
     if self.hasCoefficientOfGeneralPDE(name):
        return self.__COEFFICIENTS_OF_GENEARL_PDE[name].getFunctionSpace(self.getDomain())
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def getShapeOfCoefficientOfGeneralPDE(self,name):
     """
     return the shape of the coefficient name of the general PDE

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return: the shape of the coefficient name
     @rtype: C{tuple} of C{int}
     @raise IllegalCoefficient: if name is not one of coefficients
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, 
                  M{A_reduced}, M{B_reduced}, M{C_reduced}, M{D_reduced}, M{X_reduced}, M{Y_reduced}, M{d_reduced}, M{y_reduced}, M{d_contact_reduced}, M{y_contact_reduced}, M{r} or M{q}.
     """
     if self.hasCoefficientOfGeneralPDE(name):
        return self.__COEFFICIENTS_OF_GENEARL_PDE[name].getShape(self.getDomain(),self.getNumEquations(),self.getNumSolutions())
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   # =============================================================================
   # functions giving access to coefficients of a particular PDE implementation:
   # =============================================================================
   def getCoefficient(self,name):
     """
     returns the value of the coefficient name

     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: the value of the coefficient name
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not a coefficient of the PDE.
     """
     if self.hasCoefficient(name):
         return self.COEFFICIENTS[name].getValue()
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def hasCoefficient(self,name):
     """
     return True if name is the name of a coefficient

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return: True if name is the name of a coefficient of the general PDE. Otherwise False.
     @rtype: C{bool}
     """
     return self.COEFFICIENTS.has_key(name)

   def createCoefficient(self, name):
     """
     create a L{Data<escript.Data>} object corresponding to coefficient name

     @return: a coefficient name initialized to 0.
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not a coefficient of the PDE.
     """
     if self.hasCoefficient(name):
        return escript.Data(0.,self.getShapeOfCoefficient(name),self.getFunctionSpaceForCoefficient(name))
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def getFunctionSpaceForCoefficient(self,name):
     """
     return the L{FunctionSpace<escript.FunctionSpace>} to be used for coefficient name

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return: the function space to be used for coefficient name
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     @raise IllegalCoefficient: if name is not a coefficient of the PDE.
     """
     if self.hasCoefficient(name):
        return self.COEFFICIENTS[name].getFunctionSpace(self.getDomain())
     else:
        raise ValueError,"unknown coefficient %s requested"%name
   def getShapeOfCoefficient(self,name):
     """
     return the shape of the coefficient name

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return: the shape of the coefficient name
     @rtype: C{tuple} of C{int}
     @raise IllegalCoefficient: if name is not a coefficient of the PDE.
     """
     if self.hasCoefficient(name):
        return self.COEFFICIENTS[name].getShape(self.getDomain(),self.getNumEquations(),self.getNumSolutions())
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def resetCoefficients(self):
     """
     resets all coefficients to there default values.
     """
     for i in self.COEFFICIENTS.iterkeys():
         self.COEFFICIENTS[i].resetValue()

   def alteredCoefficient(self,name):
     """
     announce that coefficient name has been changed

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @raise IllegalCoefficient: if name is not a coefficient of the PDE.
     @note: if name is q or r, the method will not trigger a rebuilt of the system as constraints are applied to the solved system.
     """
     if self.hasCoefficient(name):
        self.trace("Coefficient %s has been altered."%name)
        if not ((name=="q" or name=="r") and self.isUsingLumping()):
           if self.COEFFICIENTS[name].isAlteringOperator(): self.__invalidateOperator()
           if self.COEFFICIENTS[name].isAlteringRightHandSide(): self.__invalidateRightHandSide()
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def copyConstraint(self,u):
      """
      copies the constraint into u and returns u.

      @param u: a function of rank 0 is a single PDE is solved and of shape (numSolution,) for a system of PDEs
      @type u: L{Data<escript.Data>}
      @return: the input u modified by the constraints.
      @rtype: L{Data<escript.Data>}
      @warning: u is altered if it has the appropriate L{FunctionSpace<escript.FunctionSpace>}
      """
      q=self.getCoefficientOfGeneralPDE("q")
      r=self.getCoefficientOfGeneralPDE("r")
      if not q.isEmpty():
         if u.isEmpty(): u=escript.Data(0.,q.getShape(),q.getFunctionSpace())
         if r.isEmpty():
             r=escript.Data(0,u.getShape(),u.getFunctionSpace())
         else:
             r=escript.Data(r,u.getFunctionSpace())
         u.copyWithMask(r,escript.Data(q,u.getFunctionSpace()))
      return u

   def setValue(self,**coefficients):
      """
      sets new values to coefficients

      @param coefficients: new values assigned to coefficients
      @keyword A: value for coefficient A.
      @type A: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword A_reduced: value for coefficient A_reduced.
      @type A_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunction<escript.ReducedFunction>}.
      @keyword B: value for coefficient B
      @type B: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword B_reduced: value for coefficient B_reduced
      @type B_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunction<escript.ReducedFunction>}.
      @keyword C: value for coefficient C
      @type C: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword C_reduced: value for coefficient C_reduced
      @type C_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunction<escript.ReducedFunction>}.
      @keyword D: value for coefficient D
      @type D: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword D_reduced: value for coefficient D_reduced
      @type D_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunction<escript.ReducedFunction>}.
      @keyword X: value for coefficient X
      @type X: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword X_reduced: value for coefficient X_reduced
      @type X_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunction<escript.ReducedFunction>}.
      @keyword Y: value for coefficient Y
      @type Y: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword Y_reduced: value for coefficient Y_reduced
      @type Y_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunction<escript.Function>}.
      @keyword d: value for coefficient d
      @type d: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
      @keyword d_reduced: value for coefficient d_reduced
      @type d_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}.
      @keyword y: value for coefficient y
      @type y: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
      @keyword d_contact: value for coefficient d_contact
      @type d_contact: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnContactOne<escript.FunctionOnContactOne>} or  L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
      @keyword d_contact_reduced: value for coefficient d_contact_reduced
      @type d_contact_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>} or  L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>}.
      @keyword y_contact: value for coefficient y_contact
      @type y_contact: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnContactOne<escript.FunctionOnContactOne>} or  L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
      @keyword y_contact_reduced: value for coefficient y_contact_reduced
      @type y_contact_reduced: any type that can be casted to L{Data<escript.Data>} object on L{ReducedFunctionOnContactOne<escript.FunctionOnContactOne>} or L{ReducedFunctionOnContactZero<escript.FunctionOnContactZero>}.
      @keyword r: values prescribed to the solution at the locations of constraints
      @type r: any type that can be casted to L{Data<escript.Data>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the solution.
      @keyword q: mask for location of constraints
      @type q: any type that can be casted to L{Data<escript.Data>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
      @raise IllegalCoefficient: if an unknown coefficient keyword is used.
      """
      # check if the coefficients are  legal:
      for i in coefficients.iterkeys():
         if not self.hasCoefficient(i):
            raise IllegalCoefficient,"Attempt to set unknown coefficient %s"%i
      # if the number of unknowns or equations is still unknown we try to estimate them:
      if self.__numEquations==None or self.__numSolutions==None:
         for i,d in coefficients.iteritems():
            if hasattr(d,"shape"):
                s=d.shape
            elif hasattr(d,"getShape"):
                s=d.getShape()
            else:
                s=numarray.array(d).shape
            if s!=None:
                # get number of equations and number of unknowns:
                res=self.COEFFICIENTS[i].estimateNumEquationsAndNumSolutions(self.getDomain(),s)
                if res==None:
                    raise IllegalCoefficientValue,"Illegal shape %s of coefficient %s"%(s,i)
                else:
                    if self.__numEquations==None: self.__numEquations=res[0]
                    if self.__numSolutions==None: self.__numSolutions=res[1]
      if self.__numEquations==None: raise UndefinedPDEError,"unidententified number of equations"
      if self.__numSolutions==None: raise UndefinedPDEError,"unidententified number of solutions"
      # now we check the shape of the coefficient if numEquations and numSolutions are set:
      for i,d in coefficients.iteritems():
        try:
           self.COEFFICIENTS[i].setValue(self.getDomain(),
                                         self.getNumEquations(),self.getNumSolutions(),
                                         self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
           self.alteredCoefficient(i)
        except IllegalCoefficientFunctionSpace,m:
            # if the function space is wrong then we try the reduced version:
            i_red=i+"_reduced"
            if (not i_red in coefficients.keys()) and i_red in self.COEFFICIENTS.keys():
                try:
                    self.COEFFICIENTS[i_red].setValue(self.getDomain(),
                                                      self.getNumEquations(),self.getNumSolutions(),
                                                      self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
                    self.alteredCoefficient(i_red)
                except IllegalCoefficientValue,m:
                    raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
                except IllegalCoefficientFunctionSpace,m:
                    raise IllegalCoefficientFunctionSpace("Coefficient %s:%s"%(i,m))
            else:
                raise IllegalCoefficientFunctionSpace("Coefficient %s:%s"%(i,m))
        except IllegalCoefficientValue,m:
           raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
      self.__altered_coefficients=True
      # check if the systrem is inhomogeneous:
      if len(coefficients)>0 and not self.isUsingLumping():
         q=self.getCoefficientOfGeneralPDE("q")
         r=self.getCoefficientOfGeneralPDE("r")
         homogeneous_constraint=True
         if not q.isEmpty() and not r.isEmpty():
             if util.Lsup(q*r)>0.:
               self.trace("Inhomogeneous constraint detected.")
               self.__invalidateSystem()

   def getSystem(self):
       """
       return the operator and right hand side of the PDE

       @return: the discrete version of the PDE
       @rtype: C{tuple} of L{Operator,<escript.Operator>} and L{Data<escript.Data>}.
       """
       if not self.__operator_is_Valid or not self.__righthandside_isValid:
          if self.isUsingLumping():
              if not self.__operator_is_Valid:
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution(): 
                      raise TypeError,"Lumped matrix requires same order for equations and unknowns"
                 if not self.getCoefficientOfGeneralPDE("A").isEmpty(): 
                      raise ValueError,"coefficient A in lumped matrix may not be present."
                 if not self.getCoefficientOfGeneralPDE("B").isEmpty():
                      raise ValueError,"coefficient B in lumped matrix may not be present."
                 if not self.getCoefficientOfGeneralPDE("C").isEmpty():
                      raise ValueError,"coefficient C in lumped matrix may not be present."
                 if not self.getCoefficientOfGeneralPDE("A_reduced").isEmpty(): 
                      raise ValueError,"coefficient A_reduced in lumped matrix may not be present."
                 if not self.getCoefficientOfGeneralPDE("B_reduced").isEmpty():
                      raise ValueError,"coefficient B_reduced in lumped matrix may not be present."
                 if not self.getCoefficientOfGeneralPDE("C_reduced").isEmpty():
                      raise ValueError,"coefficient C_reduced in lumped matrix may not be present."
                 D=self.getCoefficientOfGeneralPDE("D")
                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numarray.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escript.Data()
                 d=self.getCoefficientOfGeneralPDE("d")
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numarray.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escript.Data()
                 d_contact=self.getCoefficientOfGeneralPDE("d_contact")
                 if not d_contact.isEmpty():
                     if self.getNumSolutions()>1:
                        d_contact_times_e=util.matrixmult(d_contact,numarray.ones((self.getNumSolutions(),)))
                     else:
                        d_contact_times_e=d_contact
                 else:
                    d_contact_times_e=escript.Data()
    
                 self.__operator=self.__getNewRightHandSide()
                 self.getDomain().addPDEToRHS(self.__operator, \
                                              escript.Data(), \
                                              D_times_e, \
                                              d_times_e,\
                                              d_contact_times_e)
                 D_reduced=self.getCoefficientOfGeneralPDE("D_reduced")
                 if not D_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        D_reduced_times_e=util.matrix_mult(D_reduced,numarray.ones((self.getNumSolutions(),)))
                     else:
                        D_reduced_times_e=D_reduced
                 else:
                    D_reduced_times_e=escript.Data()
                 d_reduced=self.getCoefficientOfGeneralPDE("d_reduced")
                 if not d_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        d_reduced_times_e=util.matrix_mult(d_reduced,numarray.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_times_e=d_reduced
                 else:
                    d_reduced_times_e=escript.Data()
                 d_contact_reduced=self.getCoefficientOfGeneralPDE("d_contact_reduced")
                 if not d_contact_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        d_contact_reduced_times_e=util.matrixmult(d_contact_reduced,numarray.ones((self.getNumSolutions(),)))
                     else:
                        d_contact_reduced_times_e=d_contact_reduced
                 else:
                    d_contact_reduced_times_e=escript.Data()
    
                 self.__operator=self.__getNewRightHandSide()
                 self.getDomain().addPDEToRHS(self.__operator, \
                                              escript.Data(), \
                                              D_times_e, \
                                              d_times_e,\
                                              d_contact_times_e)
                 self.getDomain().addPDEToRHS(self.__operator, \
                                              escript.Data(), \
                                              D_reduced_times_e, \
                                              d_reduced_times_e,\
                                              d_contact_reduced_times_e)
                 self.__operator=1./self.__operator
                 self.trace("New lumped operator has been built.")
                 self.__operator_is_Valid=True
              if not self.__righthandside_isValid:
                 self.getDomain().addPDEToRHS(self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfGeneralPDE("X"), \
                               self.getCoefficientOfGeneralPDE("Y"),\
                               self.getCoefficientOfGeneralPDE("y"),\
                               self.getCoefficientOfGeneralPDE("y_contact"))
                 self.getDomain().addPDEToRHS(self.__righthandside, \
                               self.getCoefficientOfGeneralPDE("X_reduced"), \
                               self.getCoefficientOfGeneralPDE("Y_reduced"),\
                               self.getCoefficientOfGeneralPDE("y_reduced"),\
                               self.getCoefficientOfGeneralPDE("y_contact_reduced"))
                 self.trace("New right hand side as been built.")
                 self.__righthandside_isValid=True
          else:
             if not self.__operator_is_Valid and not self.__righthandside_isValid:
                 self.getDomain().addPDEToSystem(self.__makeFreshOperator(),self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfGeneralPDE("A"), \
                               self.getCoefficientOfGeneralPDE("B"), \
                               self.getCoefficientOfGeneralPDE("C"), \
                               self.getCoefficientOfGeneralPDE("D"), \
                               self.getCoefficientOfGeneralPDE("X"), \
                               self.getCoefficientOfGeneralPDE("Y"), \
                               self.getCoefficientOfGeneralPDE("d"), \
                               self.getCoefficientOfGeneralPDE("y"), \
                               self.getCoefficientOfGeneralPDE("d_contact"), \
                               self.getCoefficientOfGeneralPDE("y_contact"))
                 self.getDomain().addPDEToSystem(self.__operator,self.__righthandside, \
                               self.getCoefficientOfGeneralPDE("A_reduced"), \
                               self.getCoefficientOfGeneralPDE("B_reduced"), \
                               self.getCoefficientOfGeneralPDE("C_reduced"), \
                               self.getCoefficientOfGeneralPDE("D_reduced"), \
                               self.getCoefficientOfGeneralPDE("X_reduced"), \
                               self.getCoefficientOfGeneralPDE("Y_reduced"), \
                               self.getCoefficientOfGeneralPDE("d_reduced"), \
                               self.getCoefficientOfGeneralPDE("y_reduced"), \
                               self.getCoefficientOfGeneralPDE("d_contact_reduced"), \
                               self.getCoefficientOfGeneralPDE("y_contact_reduced"))
                 self.__applyConstraint()
                 self.__righthandside=self.copyConstraint(self.__righthandside)
                 self.trace("New system has been built.")
                 self.__operator_is_Valid=True
                 self.__righthandside_isValid=True
             elif not self.__righthandside_isValid:
                 self.getDomain().addPDEToRHS(self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfGeneralPDE("X"), \
                               self.getCoefficientOfGeneralPDE("Y"),\
                               self.getCoefficientOfGeneralPDE("y"),\
                               self.getCoefficientOfGeneralPDE("y_contact"))
                 self.getDomain().addPDEToRHS(self.__righthandside, \
                               self.getCoefficientOfGeneralPDE("X_reduced"), \
                               self.getCoefficientOfGeneralPDE("Y_reduced"),\
                               self.getCoefficientOfGeneralPDE("y_reduced"),\
                               self.getCoefficientOfGeneralPDE("y_contact_reduced"))
                 self.__righthandside=self.copyConstraint(self.__righthandside)
                 self.trace("New right hand side has been built.")
                 self.__righthandside_isValid=True
             elif not self.__operator_is_Valid:
                 self.getDomain().addPDEToSystem(self.__makeFreshOperator(),escript.Data(), \
                            self.getCoefficientOfGeneralPDE("A"), \
                            self.getCoefficientOfGeneralPDE("B"), \
                            self.getCoefficientOfGeneralPDE("C"), \
                            self.getCoefficientOfGeneralPDE("D"), \
                            escript.Data(), \
                            escript.Data(), \
                            self.getCoefficientOfGeneralPDE("d"), \
                            escript.Data(),\
                            self.getCoefficientOfGeneralPDE("d_contact"), \
                            escript.Data())
                 self.getDomain().addPDEToSystem(self.__operator,escript.Data(), \
                            self.getCoefficientOfGeneralPDE("A_reduced"), \
                            self.getCoefficientOfGeneralPDE("B_reduced"), \
                            self.getCoefficientOfGeneralPDE("C_reduced"), \
                            self.getCoefficientOfGeneralPDE("D_reduced"), \
                            escript.Data(), \
                            escript.Data(), \
                            self.getCoefficientOfGeneralPDE("d_reduced"), \
                            escript.Data(),\
                            self.getCoefficientOfGeneralPDE("d_contact_reduced"), \
                            escript.Data())
                 self.__applyConstraint()
                 self.trace("New operator has been built.")
                 self.__operator_is_Valid=True
       return (self.__operator,self.__righthandside)


class Poisson(LinearPDE):
   """
   Class to define a Poisson equation problem, which is genear L{LinearPDE} of the form

   M{-grad(grad(u)[j])[j] = f}

   with natural boundary conditons

   M{n[j]*grad(u)[j] = 0 }

   and constraints:

   M{u=0} where M{q>0}

   """

   def __init__(self,domain,debug=False):
     """
     initializes a new Poisson equation

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param debug: if True debug informations are printed.

     """
     super(Poisson, self).__init__(domain,1,1,debug)
     self.COEFFICIENTS={"f": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "f_reduced": PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "q": PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH)}
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     sets new values to coefficients

     @param coefficients: new values assigned to coefficients
     @keyword f: value for right hand side M{f}
     @type f: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword q: mask for location of constraints
     @type q: any type that can be casted to rank zeo L{Data<escript.Data>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
     @raise IllegalCoefficient: if an unknown coefficient keyword is used.
     """
     super(Poisson, self).setValue(**coefficients)

   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE
     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: the value of the coefficient  name
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
     @note: This method is called by the assembling routine to map the Possion equation onto the general PDE.
     """
     if name == "A" :
         return escript.Data(util.kronecker(self.getDim()),escript.Function(self.getDomain()))
     elif name == "B" :
         return escript.Data()
     elif name == "C" :
         return escript.Data()
     elif name == "D" :
         return escript.Data()
     elif name == "X" :
         return escript.Data()
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "d" :
         return escript.Data()
     elif name == "y" :
         return escript.Data()
     elif name == "d_contact" :
         return escript.Data()
     elif name == "y_contact" :
         return escript.Data()
     elif name == "A_reduced" :
         return escript.Data()
     elif name == "B_reduced" :
         return escript.Data()
     elif name == "C_reduced" :
         return escript.Data()
     elif name == "D_reduced" :
         return escript.Data()
     elif name == "X_reduced" :
         return escript.Data()
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     elif name == "d_reduced" :
         return escript.Data()
     elif name == "y_reduced" :
         return escript.Data()
     elif name == "d_contact_reduced" :
         return escript.Data()
     elif name == "y_contact_reduced" :
         return escript.Data()
     elif name == "r" :
         return escript.Data()
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

class Helmholtz(LinearPDE):
   """
   Class to define a Helmhotz equation problem, which is genear L{LinearPDE} of the form

   M{S{omega}*u - grad(k*grad(u)[j])[j] = f}

   with natural boundary conditons

   M{k*n[j]*grad(u)[j] = g- S{alpha}u }

   and constraints:

   M{u=r} where M{q>0}

   """

   def __init__(self,domain,debug=False):
     """
     initializes a new Poisson equation

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param debug: if True debug informations are printed.

     """
     super(Helmholtz, self).__init__(domain,1,1,debug)
     self.COEFFICIENTS={"omega": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.OPERATOR),
                        "k": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.OPERATOR),
                        "f": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "f_reduced": PDECoefficient(PDECoefficient.INTERIOR_REDUCED,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "alpha": PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.OPERATOR),
                        "g": PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "g_reduced": PDECoefficient(PDECoefficient.BOUNDARY_REDUCED,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "r": PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH),
                        "q": PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH)}
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     sets new values to coefficients

     @param coefficients: new values assigned to coefficients
     @keyword omega: value for coefficient M{S{omega}}
     @type omega: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword k: value for coefficeint M{k}
     @type k: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword f: value for right hand side M{f}
     @type f: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword alpha: value for right hand side M{S{alpha}}
     @type alpha: any type that can be casted to L{Scalar<escript.Scalar>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
     @keyword g: value for right hand side M{g}
     @type g: any type that can be casted to L{Scalar<escript.Scalar>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
     @keyword r: prescribed values M{r} for the solution in constraints.
     @type r: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
     @keyword q: mask for location of constraints
     @type q: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
     @raise IllegalCoefficient: if an unknown coefficient keyword is used.
     """
     super(Helmholtz, self).setValue(**coefficients)

   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: the value of the coefficient  name
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  "A", M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
     @note: This method is called by the assembling routine to map the Possion equation onto the general PDE.
     """
     if name == "A" :
         return escript.Data(numarray.identity(self.getDim()),escript.Function(self.getDomain()))*self.getCoefficient("k")
     elif name == "B" :
         return escript.Data()
     elif name == "C" :
         return escript.Data()
     elif name == "D" :
         return self.getCoefficient("omega")
     elif name == "X" :
         return escript.Data()
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "d" :
         return self.getCoefficient("alpha")
     elif name == "y" :
         return self.getCoefficient("g")
     elif name == "d_contact" :
         return escript.Data()
     elif name == "y_contact" :
         return escript.Data()
     elif name == "A_reduced" :
         return escript.Data()
     elif name == "B_reduced" :
         return escript.Data()
     elif name == "C_reduced" :
         return escript.Data()
     elif name == "D_reduced" :
         return escript.Data()
     elif name == "X_reduced" :
         return escript.Data()
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     elif name == "d_reduced" :
         return escript.Data()
     elif name == "y_reduced" :
        return self.getCoefficient("g_reduced")
     elif name == "d_contact_reduced" :
         return escript.Data()
     elif name == "y_contact_reduced" :
         return escript.Data()
     elif name == "r" :
         return self.getCoefficient("r")
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

class LameEquation(LinearPDE):
   """
   Class to define a Lame equation problem:

   M{-grad(S{mu}*(grad(u[i])[j]+grad(u[j])[i]))[j] - grad(S{lambda}*grad(u[k])[k])[j] = F_i -grad(S{sigma}[ij])[j] }

   with natural boundary conditons:

   M{n[j]*(S{mu}*(grad(u[i])[j]+grad(u[j])[i]) + S{lambda}*grad(u[k])[k]) = f_i +n[j]*S{sigma}[ij] }

   and constraints:

   M{u[i]=r[i]} where M{q[i]>0}

   """

   def __init__(self,domain,debug=False):
      super(LameEquation, self).__init__(domain,\
                                         domain.getDim(),domain.getDim(),debug)
      self.COEFFICIENTS={ "lame_lambda"  : PDECoefficient(PDECoefficient.INTERIOR,(),PDECoefficient.OPERATOR),
                          "lame_mu"      : PDECoefficient(PDECoefficient.INTERIOR,(),PDECoefficient.OPERATOR),
                          "F"            : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                          "sigma"        : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,PDECoefficient.BY_DIM),PDECoefficient.RIGHTHANDSIDE),
                          "f"            : PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                          "r"            : PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH),
                          "q"            : PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH)}
      self.setSymmetryOn()

   def setValues(self,**coefficients):
     """
     sets new values to coefficients

     @param coefficients: new values assigned to coefficients
     @keyword lame_mu: value for coefficient M{S{mu}}
     @type lame_mu: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword lame_lambda: value for coefficient M{S{lambda}}
     @type lame_lambda: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword F: value for internal force M{F}
     @type F: any type that can be casted to L{Vector<escript.Vector>} object on L{Function<escript.Function>}
     @keyword sigma: value for initial stress M{S{sigma}}
     @type sigma: any type that can be casted to L{Tensor<escript.Tensor>} object on L{Function<escript.Function>}
     @keyword f: value for extrenal force M{f}
     @type f: any type that can be casted to L{Vector<escript.Vector>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}
     @keyword r: prescribed values M{r} for the solution in constraints.
     @type r: any type that can be casted to L{Vector<escript.Vector>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
     @keyword q: mask for location of constraints
     @type q: any type that can be casted to L{Vector<escript.Vector>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
     @raise IllegalCoefficient: if an unknown coefficient keyword is used.
     """
     super(LameEquation, self).setValues(**coefficients)

   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: the value of the coefficient  name
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  "A", M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
     @note: This method is called by the assembling routine to map the Possion equation onto the general PDE.
     """
     if name == "A" :
         out =self.createCoefficientOfGeneralPDE("A")
         for i in range(self.getDim()):
           for j in range(self.getDim()):
             out[i,i,j,j] += self.getCoefficient("lame_lambda")
             out[i,j,j,i] += self.getCoefficient("lame_mu")
             out[i,j,i,j] += self.getCoefficient("lame_mu")
         return out
     elif name == "B" :
         return escript.Data()
     elif name == "C" :
         return escript.Data()
     elif name == "D" :
         return escript.Data()
     elif name == "X" :
         return self.getCoefficient("sigma")
     elif name == "Y" :
         return self.getCoefficient("F")
     elif name == "d" :
         return escript.Data()
     elif name == "y" :
         return self.getCoefficient("f")
     elif name == "d_contact" :
         return escript.Data()
     elif name == "y_contact" :
         return escript.Data()
     elif name == "A_reduced" :
         return escript.Data()
     elif name == "B_reduced" :
         return escript.Data()
     elif name == "C_reduced" :
         return escript.Data()
     elif name == "D_reduced" :
         return escript.Data()
     elif name == "X_reduced" :
         return escript.Data()
     elif name == "Y_reduced" :
         return escript.Data()
     elif name == "d_reduced" :
         return escript.Data()
     elif name == "y_reduced" :
         return escript.Data()
     elif name == "d_contact_reduced" :
         return escript.Data()
     elif name == "y_contact_reduced" :
         return escript.Data()
     elif name == "r" :
         return self.getCoefficient("r")
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

