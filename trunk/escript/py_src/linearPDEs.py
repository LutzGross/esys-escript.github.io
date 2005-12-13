# $Id$

#
#      COPYRIGHT ACcESS 2004 -  All Rights Reserved
#
#   This software is the property of ACcESS.  No part of this code
#   may be copied in any form or by any means without the expressed written
#   consent of ACcESS.  Copying, use or modification of this software
#   by any unauthorised person is illegal unless that
#   person has a software license agreement with ACcESS.
#
"""
The module provides an interface to define and solve linear partial
differential equations (PDEs) within L{escript}. L{linearPDEs} does not provide any
solver capabilities in itself but hands the PDE over to
the PDE solver library defined through the L{Domain<escript.Domain>} of the PDE.
The general interface is provided through the L{LinearPDE} class. The
L{AdvectivePDE} which is derived from the L{LinearPDE} class
provides an interface to PDE dominated by its advective terms. The L{Poisson},
L{Helmholtz}, L{LameEquation}, L{AdvectionDiffusion}
classs which are also derived form the L{LinearPDE} class should be used
to define of solve these sepecial PDEs.

@var __author__: name of author
@var __licence__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

import escript
import util
import numarray

__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"


class IllegalCoefficient(ValueError):
   """
   raised if an illegal coefficient of the general ar particular PDE is requested.
   """

class IllegalCoefficientValue(ValueError):
   """
   raised if an incorrect value for a coefficient is used.
   """

class UndefinedPDEError(ValueError):
   """
   raised if a PDE is not fully defined yet.
   """

class PDECoefficient(object):
    """
    A class for describing a PDE coefficient

    @cvar INTERIOR: indicator that coefficient is defined on the interior of the PDE domain
    @cvar BOUNDARY: indicator that coefficient is defined on the boundary of the PDE domain
    @cvar CONTACT: indicator that coefficient is defined on the contact region within the PDE domain
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

    def __init__(self,where,pattern,altering):
       """
       Initialise a PDE Coefficient type

       @param where: describes where the coefficient lives
       @type where: one of L{INTERIOR}, L{BOUNDARY}, L{CONTACT}, L{SOLUTION}, L{REDUCED}
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
       @type domain: C{bool}
       @param reducedSolutionOrder: True to indicate that reduced order is used to represent the solution
       @type domain: C{bool}
       @return:  L{FunctionSpace<escript.FunctionSpace>} of the coefficient
       @rtype:  L{FunctionSpace<escript.FunctionSpace>}
       """
       if self.what==self.INTERIOR:
            return escript.Function(domain)
       elif self.what==self.BOUNDARY:
            return escript.FunctionOnBoundary(domain)
       elif self.what==self.CONTACT:
            return escript.FunctionOnContactZero(domain)
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
       @type domain: C{bool}
       @param reducedSolutionOrder: True to indicate that reduced order is used to represent the solution
       @type domain: C{bool}
       @param newValue: number of components of the PDE solution
       @type newValue: any object that can be converted into a L{Data<escript.Data>} object with the appropriate shape and L{FunctionSpace<escript.FunctionSpace>}
       @raise IllegalCoefficientValue: if the shape of the assigned value does not match the shape of the coefficient
       """
       if newValue==None:
           newValue=escript.Data()
       elif isinstance(newValue,escript.Data):
           if not newValue.isEmpty():
              try:
                 newValue=escript.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
              except:
                 raise IllegalCoefficientValue,"Unable to interpolate coefficient to function space %s"%self.getFunctionSpace(domain)
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

   M{-grad(A[j,l]*grad(u)[l]+B[j]u)[j]+C[l]*grad(u)[l]+D*u =-grad(X)[j,j]+Y}

   where M{grad(F)} denotes the spatial derivative of M{F}. Einstein's summation convention,
   ie. summation over indexes appearing twice in a term of a sum is performed, is used.
   The coefficients M{A}, M{B}, M{C}, M{D}, M{X} and M{Y} have to be specified through L{Data<escript.Data>} objects in the
   L{Function<escript.Function>} on the PDE or objects that can be converted into such L{Data<escript.Data>} objects.
   M{A} is a rank two, M{B}, M{C} and M{X} are rank one and M{D} and M{Y} are scalar.

   The following natural boundary conditions are considered:

   M{n[j]*(A[i,j]*grad(u)[l]+B[j]*u)+d*u=n[j]*X[j]+y}

   where M{n} is the outer normal field calculated by L{getNormal<escript.FunctionSpace.getNormal>} of L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
   Notice that the coefficients M{A}, M{B} and M{X} are defined in the PDE. The coefficients M{d} and M{y} are
   each a scalar in the L{FunctionOnBoundary<escript.FunctionOnBoundary>}.


   Constraints for the solution prescribing the value of the solution at certain locations in the domain. They have the form

   M{u=r}  where M{q>0}

   M{r} and M{q} are each scalar where M{q} is the characteristic function defining where the constraint is applied.
   The constraints override any other condition set by the PDE or the boundary condition.

   The PDE is symmetrical if

   M{A[i,j]=A[j,i]}  and M{B[j]=C[j]}

   For a system of PDEs and a solution with several components the PDE has the form

   M{-grad(A[i,j,k,l]*grad(u[k])[l]+B[i,j,k]*u[k])[j]+C[i,k,l]*grad(u[k])[l]+D[i,k]*u[k] =-grad(X[i,j])[j]+Y[i] }

   M{A} is a ramk four, M{B} and M{C} are each a rank three, M{D} and M{X} are each a rank two and M{Y} is a rank one.
   The natural boundary conditions take the form:

   M{n[j]*(A[i,j,k,l]*grad(u[k])[l]+B[i,j,k]*u[k])+d[i,k]*u[k]=n[j]*X[i,j]+y[i]}


   The coefficient M{d} is a rank two and M{y} is a  rank one both in the L{FunctionOnBoundary<escript.FunctionOnBoundary>}. Constraints take the form


   M{u[i]=r[i]}  where  M{q[i]>0}

   M{r} and M{q} are each rank one. Notice that at some locations not necessarily all components must have a constraint.

   The system of PDEs is symmetrical if

        - M{A[i,j,k,l]=A[k,l,i,j]}
        - M{B[i,j,k]=C[k,i,j]}
        - M{D[i,k]=D[i,k]}
        - M{d[i,k]=d[k,i]}

   L{LinearPDE} also supports solution discontinuities over a contact region in the domain. To specify the conditions across the
   discontinuity we are using the generalised flux M{J} which is in the case of a systems of PDEs and several components of the solution
   defined as

   M{J[i,j]=A[i,j,k,l]*grad(u[k])[l]+B[i,j,k]*u[k]-X[i,j]}

   For the case of single solution component and single PDE M{J} is defined

   M{J_{j}=A[i,j]*grad(u)[j]+B[i]*u-X[i]}

   In the context of discontinuities M{n} denotes the normal on the discontinuity pointing from side 0 towards side 1
   calculated from L{getNormal<escript.FunctionSpace.getNormal>} of L{FunctionOnContactZero<escript.FunctionOnContactZero>}. For a system of PDEs
   the contact condition takes the form

   M{n[j]*J0[i,j]=n[j]*J1[i,j]=y_contact[i]- d_contact[i,k]*jump(u)[k]}

   where M{J0} and M{J1} are the fluxes on side 0 and side 1 of the discontinuity, respectively. M{jump(u)}, which is the difference
   of the solution at side 1 and at side 0, denotes the jump of M{u} across discontinuity along the normal calcualted by
   L{jump<util.jump>}.
   The coefficient M{d_contact} is a rank two and M{y_contact} is a rank one both in the L{FunctionOnContactZero<escript.FunctionOnContactZero>} or L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
   In case of a single PDE and a single component solution the contact condition takes the form

   M{n[j]*J0_{j}=n[j]*J1_{j}=y_contact-d_contact*jump(u)}

   In this case the the coefficient M{d_contact} and M{y_contact} are eaach scalar
   both in the L{FunctionOnContactZero<escript.FunctionOnContactZero>} or L{FunctionOnContactOne<escript.FunctionOnContactOne>}.

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

   __TOL=1.e-13
   __PACKAGE_KEY="package"
   __METHOD_KEY="method"
   __SYMMETRY_KEY="symmetric"
   __TOLERANCE_KEY="tolerance"


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
        self.getOperator()*escript.Data(u,self.getFunctionSpaceForSolution())

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
            tol=util.Lsup(A)*self.__TOL
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
            tol=(util.Lsup(B)+util.Lsup(C))*self.__TOL/2.
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
             tol=util.Lsup(D)*self.__TOL
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(D[i,k]-D[k,i])>tol:
                      if verbose: print "non-symmetric PDE because D[%d,%d]!=D[%d,%d]"%(i,k,k,i)
                      out=False
           d=self.getCoefficientOfGeneralPDE("d")
           if not d.isEmpty():
             tol=util.Lsup(d)*self.__TOL
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(d[i,k]-d[k,i])>tol:
                      if verbose: print "non-symmetric PDE because d[%d,%d]!=d[%d,%d]"%(i,k,k,i)
                      out=False
           d_contact=self.getCoefficientOfGeneralPDE("d_contact")
           if not d_contact.isEmpty():
             tol=util.Lsup(d_contact)*self.__TOL
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(d_contact[i,k]-d_contact[k,i])>tol:
                      if verbose: print "non-symmetric PDE because d_contact[%d,%d]!=d_contact[%d,%d]"%(i,k,k,i)
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
       @keyword preconditioner: preconditioner method to be used. Allowed values are
                                L{SSOR}, L{ILU0}, L{ILUT}, L{JACOBI}
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
             options[self.__METHOD_KEY]=self.getSolverMethod()
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

     M{J[i,j]=A[i,j,k,l]*grad(u[k])[l]+B[i,j,k]u[k]-X[i,j]}

     or

     M{J[j]=A[i,j]*grad(u)[l]+B[j]u-X[j]}

     @param u: argument in the flux. If u is not present or equals L{None} the current solution is used.
     @type u: L{Data<escript.Data>} or None
     @return: flux
     @rtype: L{Data<escript.Data>}
     """
     if u==None: u=self.getSolution()
     return util.tensormult(self.getCoefficientOfGeneralPDE("A"),util.grad(u))+util.matrixmult(self.getCoefficientOfGeneralPDE("B"),u)-util.self.getCoefficientOfGeneralPDE("X")
   # =============================================================================
   #   solver settings:
   # =============================================================================
   def setSolverMethod(self,solver=None):
       """
       sets a new solver

       @param solver: sets a new solver method.
       @type solver: one of L{DEFAULT}, L{ITERATIVE} L{DIRECT}, L{CHOLEVSKY}, L{PCG}, L{CR}, L{CGS}, L{BICGSTAB}, L{SSOR}, L{GMRES}, L{PRES20}, L{LUMPING}.
       """
       if solver==None: solve=self.DEFAULT
       if not solver==self.getSolverMethod():
           self.__solver_method=solver
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
       if m==self.DEFAULT: method="DEFAULT"
       elif m==self.DIRECT: method= "DIRECT"
       elif m==self.ITERATIVE: method= "ITERATIVE"
       elif m==self.CHOLEVSKY: method= "CHOLEVSKY"
       elif m==self.PCG: method= "PCG"
       elif m==self.CR: method= "CR"
       elif m==self.CGS: method= "CGS"
       elif m==self.BICGSTAB: method= "BICGSTAB"
       elif m==self.SSOR: method= "SSOR"
       elif m==self.GMRES: method= "GMRES"
       elif m==self.PRES20: method= "PRES20"
       elif m==self.LUMPING: method= "LUMPING"
       else : method="unknown"
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
       return self.__solver_method

   def setSolverPackage(self,package=None):
       """
       sets a new solver package

       @param solver: sets a new solver method.
       @type solver: one of L{DEFAULT}, L{PASO} L{SCSL}, L{MKL}, L{UMLPACK}
       """
       if package==None: package=self.DEFAULT
       if not package==self.getSolverPackage():
           self.__solver_method=solver
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
      return self.getSolverMethod()==self.LUMPING

   def setTolerance(self,tol=1.e-8):
       """
       resets the tolerance for the solver method to tol where for an appropriate norm M{|.|}

       M{|L{getResidual}()|<tol*|L{getRightHandSide}()|}

       defines the stopping criterion.

       @param tol: new tolerance for the solver. If the tol is lower then the current tolerence
                   the system will be resolved.
       @type tol: positive C{float}
       @raise ValueException: if tolerance is not positive.
       """
       if not tol>0:
           raise ValueException,"Tolerance as to be positive"
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
     new_matrix_type=self.getDomain().getSystemMatrixTypeId(self.getSolverMethod(),self.getSolverPackage(),self.isSymmetric())
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
           self.__righthandside*=0
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
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
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
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
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
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
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
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
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
      @keyword B: value for coefficient B
      @type B: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword C: value for coefficient C
      @type C: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword D: value for coefficient D
      @type D: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword X: value for coefficient X
      @type X: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword Y: value for coefficient Y
      @type Y: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword d: value for coefficient d
      @type d: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
      @keyword y: value for coefficient y
      @type y: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
      @keyword d_contact: value for coefficient d_contact
      @type d_contact: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
                       or  L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
      @keyword y_contact: value for coefficient y_contact
      @type y_contact: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
                       or  L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
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
           self.COEFFICIENTS[i].setValue(self.getDomain(),self.getNumEquations(),self.getNumSolutions(),self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
        except IllegalCoefficientValue,m:
           raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
        self.alteredCoefficient(i)

      self.__altered_coefficients=True
      # check if the systrem is inhomogeneous:
      if len(coefficients)>0 and not self.isUsingLumping():
         q=self.getCoefficientOfGeneralPDE("q")
         r=self.getCoefficientOfGeneralPDE("r")
         homogeneous_constraint=True
         if not q.isEmpty() and not r.isEmpty():
             if util.Lsup(q*r)>=1.e-13*util.Lsup(r):
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
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution(): raise TypeError,"Lumped matrix requires same order for equations and unknowns"
                 if not self.getCoefficientOfGeneralPDE("A").isEmpty(): raise Warning,"Using coefficient A in lumped matrix can produce wrong results"
                 if not self.getCoefficientOfGeneralPDE("B").isEmpty(): raise Warning,"Using coefficient B in lumped matrix can produce wrong results"
                 if not self.getCoefficientOfGeneralPDE("C").isEmpty(): raise Warning,"Using coefficient C in lumped matrix can produce wrong results"
                 mat=self.__getNewOperator()
                 self.getDomain().addPDEToSystem(mat,escript.Data(), \
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
                 self.__operator=1./(mat*escript.Data(1,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True))
                 del mat
                 self.trace("New lumped operator has been built.")
                 self.__operator_is_Valid=True
              if not self.__righthandside_isValid:
                 self.getDomain().addPDEToRHS(self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfGeneralPDE("X"), \
                               self.getCoefficientOfGeneralPDE("Y"),\
                               self.getCoefficientOfGeneralPDE("y"),\
                               self.getCoefficientOfGeneralPDE("y_contact"))
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
                        "alpha": PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.OPERATOR),
                        "g": PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
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
     elif name == "r" :
         return self.getCoefficient("r")
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

class LameEquation(LinearPDE):
   """
   Class to define a Lame equation problem:

   M{-grad(S{mu}*(grad(u[i])[j]+grad(u[j])[i]))[j] - grad(S{lambda}*grad(u[j])[i])[j] = F_i -grad(S{sigma}[i,j])[j] }

   with natural boundary conditons:

   M{n[j]*(S{mu}*(grad(u[i])[j]+grad(u[j])[i]) - S{lambda}*grad(u[j])[i]) = f_i -n[j]*S{sigma}[i,j] }

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

   def setValue(self,**coefficients):
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
     super(LameEquation, self).setValue(**coefficients)

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
     elif name == "r" :
         return self.getCoefficient("r")
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

class AdvectivePDE(LinearPDE):
   """
   In cases of PDEs dominated by the advection terms M{B} and M{C} against the adevctive terms M{A}
   up-winding has been used.  The L{AdvectivePDE} class applies SUPG upwinding to the advective terms.

   In the following we set

   M{Z[j]=C[j]-B[j]}

   or

   M{Z[i,k,l]=C[i,k,l]-B[i,l,k]}

   To measure the dominance of the advective terms over the diffusive term M{A} the
   X{Pelclet number} M{P} is used. It is defined as

   M{P=h|Z|/(2|A|)}

   where M{|.|} denotes the L{length<util.length>} of the arument and M{h} is the local cell size
   from L{getSize<escript.Domain.getSize>}. Where M{|A|==0} M{P} is M{S{infinity}}.

   From the X{Pelclet number} the stabilization parameters M{S{Xi}} and M{S{Xi}} are calculated:

   M{S{Xi}=S{xi}(P) h/|Z|}

   where M{S{xi}} is a suitable function of the Peclet number.

   In the case of a single PDE the coefficient are up-dated in the following way:
         - M{A[i,j] S{<-} A[i,j] + S{Xi} * Z[j] * Z[l]}
         - M{B[j] S{<-} B[j] + S{Xi} * C[j] * D}
         - M{C[j] S{<-} C[j] + S{Xi} * B[j] * D}
         - M{X[j] S{<-} X[j] + S{Xi} * Z[j] * Y}

   Similar for the case of a systems of PDEs:
         - M{A[i,j,k,l] S{<-} A[i,j,k,l]+ S{delta}[p,m] * S{Xi} * Z[p,i,j] * Z[m,k,l]}
         - M{B[i,j,k] S{<-} B[i,j,k] +  S{delta}[p,m] * S{Xi} * D[p,k] * C[m,i,j]}
         - M{C[i,k,l] S{<-} C[i,k,l] +  S{delta}[p,m] * S{Xi} * D[p,k] * B[m,l,i]}
         - M{X[i,j] S{<-} X[i,j] + S{delta}[p,m] * S{Xi}  * Y[p] * Z[m,i,j]}

   where M{S{delta}} is L{kronecker}.
   Using upwinding in this form, introduces an additonal error which is proprtional to the cell size M{h}
   but with the intension to stabilize the solution.

   """
   def __init__(self,domain,numEquations=None,numSolutions=None,xi=None,debug=False):
      """
      creates a linear, steady, second order PDE on a L{Domain<escript.Domain>}

      @param domain: domain of the PDE
      @type domain: L{Domain<escript.Domain>}
      @param numEquations: number of equations. If numEquations==None the number of equations
                           is exracted from the PDE coefficients.
      @param numSolutions: number of solution components. If  numSolutions==None the number of solution components
                           is exracted from the PDE coefficients.
      @param xi: defines a function which returns for any given Preclet number as L{Scalar<escript.Scalar>} object the
                 M{S{xi}}-value used to define the stabilization parameters. If equal to None, L{ELMAN_RAMAGE} is used.
      @type xi: callable object which returns a L{Scalar<escript.Scalar>} object.
      @param debug: if True debug informations are printed.
      """
      super(AdvectivePDE, self).__init__(domain,\
                                         numEquations,numSolutions,debug)
      if xi==None:
         self.__xi=AdvectivePDE.ELMAN_RAMAGE
      else:
         self.__xi=xi
      self.__Xi=escript.Data()

   def setValue(**coefficients):
      """
      sets new values to coefficients

      @param coefficients: new values assigned to coefficients
      @keyword A: value for coefficient A.
      @type A: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword B: value for coefficient B
      @type B: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword C: value for coefficient C
      @type C: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword D: value for coefficient D
      @type D: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword X: value for coefficient X
      @type X: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword Y: value for coefficient Y
      @type Y: any type that can be casted to L{Data<escript.Data>} object on L{Function<escript.Function>}.
      @keyword d: value for coefficient d
      @type d: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
      @keyword y: value for coefficient y
      @type y: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}.
      @keyword d_contact: value for coefficient d_contact
      @type d_contact: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
                       or  L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
      @keyword y_contact: value for coefficient y_contact
      @type y_contact: any type that can be casted to L{Data<escript.Data>} object on L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
                       or  L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
      @keyword r: values prescribed to the solution at the locations of constraints
      @type r: any type that can be casted to L{Data<escript.Data>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the solution.
      @keyword q: mask for location of constraints
      @type q: any type that can be casted to L{Data<escript.Data>} object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending of reduced order is used for the representation of the equation.
      @raise IllegalCoefficient: if an unknown coefficient keyword is used.

      """
      if "A" in coefficients.keys()   or "B" in coefficients.keys() or "C" in coefficients.keys(): self.__Xi=escript.Data()
      super(AdvectivePDE, self).setValue(**coefficients)

   def ELMAN_RAMAGE(self,P):
     """
     Predefined function to set a values for M{S{xi}} from a Preclet number M{P}.
     This function uses the method suggested by H.C. Elman and A. Ramage, I{SIAM J. Numer. Anal.}, B{40} (2002)
          - M{S{xi}(P)=0} for M{P<1}
          - M{S{xi}(P)=(1-1/P)/2} otherwise

     @param P: Preclet number
     @type P: L{Scalar<escript.Scalar>}
     @return: up-wind weightimg factor
     @rtype: L{Scalar<escript.Scalar>}
     """
     return util.wherePositive(P-1.)*0.5*(1.-1./(P+1.e-15))

   def SIMPLIFIED_BROOK_HUGHES(self,P):
     """
     Predefined function to set a values for M{S{xi}} from a Preclet number M{P}.
     The original methods is

     M{S{xi}(P)=coth(P)-1/P}

     As the evaluation of M{coth} is expensive we are using the approximation:

         - M{S{xi}(P)=P/3} where M{P<3}
         - M{S{xi}(P)=1/2} otherwise

     @param P: Preclet number
     @type P: L{Scalar<escript.Scalar>}
     @return: up-wind weightimg factor
     @rtype: L{Scalar<escript.Scalar>}
     """
     c=util.whereNegative(P-3.)
     return P/6.*c+1./2.*(1.-c)

   def HALF(self,P):
     """
     Predefined function to set value M{1/2} for M{S{xi}}

     @param P: Preclet number
     @type P: L{Scalar<escript.Scalar>}
     @return: up-wind weightimg factor
     @rtype: L{Scalar<escript.Scalar>}
     """
     return escript.Scalar(0.5,P.getFunctionSpace())

   def __calculateXi(self,peclet_factor,flux,h):
       flux=util.Lsup(flux)
       if flux_max>0.:
          return h*self.__xi(flux*peclet_factor)/(flux+flux_max*self.__TOL)
       else:
          return 0.

   def __getXi(self):
      if self.__Xi.isEmpty():
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         A=self.getCoefficient("A")
         h=self.getDomain().getSize()
         self.__Xi=escript.Scalar(0.,self.getFunctionSpaceForCoefficient("A"))
         if not C.isEmpty() or not B.isEmpty():
            if not C.isEmpty() and not B.isEmpty():
                flux2=escript.Scalar(0,self.getFunctionSpaceForCoefficient("A"))
                if self.getNumEquations()>1:
                   if self.getNumSolutions()>1:
                      for i in range(self.getNumEquations()):
                         for k in range(self.getNumSolutions()):
                            for l in range(self.getDim()): flux2+=(C[i,k,l]-B[i,l,k])**2
                      # flux=C-util.reorderComponents(B,[0,2,1])
                   else:
                      for i in range(self.getNumEquations()):
                         for l in range(self.getDim()): flux2+=(C[i,l]-B[i,l])**2
                      # flux=C-B
                else:
                   if self.getNumSolutions()>1:
                      for k in range(self.getNumSolutions()):
                         for l in range(self.getDim()): flux2+=(C[k,l]-B[l,k])**2
                      # flux=C-util.reorderComponents(B,[1,0])
                   else:
                      for l in range(self.getDim()): flux2+=(C[l]-B[l])**2
                      #flux=C-B
                length_of_flux=util.sqrt(flux2)
            elif C.isEmpty():
              length_of_flux=util.length(B)
              #flux=B
            else:
              length_of_flux=util.length(C)
              #flux=C

            #length_of_flux=util.length(flux)
            flux_max=util.Lsup(length_of_flux)
            if flux_max>0.:
               # length_of_A=util.inner(flux,util.tensormutiply(A,flux))
               length_of_A=util.length(A)
               A_max=util.Lsup(length_of_A)
               if A_max>0:
                    inv_A=1./(length_of_A+A_max*self.__TOL)
               else:
                    inv_A=1./self.__TOL
               peclet_number=length_of_flux*h/2*inv_A
               xi=self.__xi(peclet_number)
               self.__Xi=h*xi/(length_of_flux+flux_max*self.__TOL)
               self.trace("preclet number = %e"%util.Lsup(peclet_number))
      return self.__Xi


   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name: name of the coefficient requested.
     @type name: C{string}
     @return: the value of the coefficient name
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if name is not one of coefficients
                  M{A}, M{B}, M{C}, M{D}, M{X}, M{Y}, M{d}, M{y}, M{d_contact}, M{y_contact}, M{r} or M{q}.
     @note: This method is called by the assembling routine to map the Possion equation onto the general PDE.
     """
     if not self.getNumEquations() == self.getNumSolutions():
          raise ValueError,"AdvectivePDE expects the number of solution componets and the number of equations to be equal."

     if name == "A" :
         A=self.getCoefficient("A")
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         if B.isEmpty() and C.isEmpty():
            Aout=A
         else:
            if A.isEmpty():
               Aout=self.createNewCoefficient("A")
            else:
               Aout=A[:]
            Xi=self.__getXi()
            if self.getNumEquations()>1:
                for i in range(self.getNumEquations()):
                   for j in range(self.getDim()):
                      for k in range(self.getNumSolutions()):
                         for l in range(self.getDim()):
                            if not C.isEmpty() and not B.isEmpty():
                               # tmp=C-util.reorderComponents(B,[0,2,1])
                               # Aout=Aout+Xi*util.generalTensorProduct(util.reorder(tmp,[1,2,0]),tmp,offset=1)
                               for p in range(self.getNumEquations()): Aout[i,j,k,l]+=Xi*(C[p,i,j]-B[p,j,i])*(C[p,k,l]-B[p,l,k])
                            elif C.isEmpty():
                               for p in range(self.getNumEquations()): Aout[i,j,k,l]+=Xi*B[p,j,i]*B[p,l,k]
                               # Aout=Aout+Xi*util.generalTensorProduct(util.reorder(B,[2,1,0]),util.reorder(B,[0,2,1]),offset=1)
                            else:
                               for p in range(self.getNumEquations()): Aout[i,j,k,l]+=Xi*C[p,i,j]*C[p,k,l]
                               # Aout=Aout+Xi*util.generalTensorProduct(util.reorder(C,[1,2,0]),C,offset=1)
            else:
                for j in range(self.getDim()):
                   for l in range(self.getDim()):
                      if not C.isEmpty() and not B.isEmpty():
                          Aout[j,l]+=Xi*(C[j]-B[j])*(C[l]-B[l])
                      elif C.isEmpty():
                          Aout[j,l]+=Xi*B[j]*B[l]
                      else:
                          Aout[j,l]+=Xi*C[j]*C[l]
                 # if not C.isEmpty() and not B.isEmpty():
                 #    tmp=C-B
                 #    Aout=Aout+Xi*util.outer(tmp,tmp)
                 # elif C.isEmpty():
                 #    Aout=Aout+Xi*util.outer(B,B)
                 # else:
                 # Aout=Aout+Xi*util.outer(C,C)
         return Aout
     elif name == "B" :
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         D=self.getCoefficient("D")
         if C.isEmpty() or D.isEmpty():
            Bout=B
         else:
            Xi=self.__getXi()
            if B.isEmpty():
                Bout=self.createNewCoefficient("B")
            else:
                Bout=B[:]
            if self.getNumEquations()>1:
               for k in range(self.getNumSolutions()):
                  for p in range(self.getNumEquations()):
                     tmp=Xi*D[p,k]
                     for i in range(self.getNumEquations()):
                        for j in range(self.getDim()):
                           Bout[i,j,k]+=tmp*C[p,i,j]
                           # Bout=Bout+Xi*util.generalTensorProduct(util.reorder(C,[1,2,0]),D,offset=1)
            else:
               tmp=Xi*D
               for j in range(self.getDim()): Bout[j]+=tmp*C[j]
               # Bout=Bout+Xi*D*C
         return Bout
     elif name == "C" :
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         D=self.getCoefficient("D")
         if B.isEmpty() or D.isEmpty():
            Cout=C
         else:
            Xi=self.__getXi()
            if C.isEmpty():
                Cout=self.createNewCoefficient("C")
            else:
                Cout=C[:]
            if self.getNumEquations()>1:
               for k in range(self.getNumSolutions()):
                   for p in range(self.getNumEquations()):
                      tmp=Xi*D[p,k]
                      for i in range(self.getNumEquations()):
                        for l in range(self.getDim()):
                                 Cout[i,k,l]+=tmp*B[p,l,i]
                                 # Cout=Cout+Xi*B[p,l,i]*D[p,k]
            else:
               tmp=Xi*D
               for j in range(self.getDim()): Cout[j]+=tmp*B[j]
               # Cout=Cout+tmp*D*B
         return Cout
     elif name == "D" :
         return self.getCoefficient("D")
     elif name == "X" :
         X=self.getCoefficient("X")
         Y=self.getCoefficient("Y")
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         if Y.isEmpty() or (B.isEmpty() and C.isEmpty()):
            Xout=X
         else:
            if X.isEmpty():
                Xout=self.createNewCoefficient("X")
            else:
                Xout=X[:]
            Xi=self.__getXi()
            if self.getNumEquations()>1:
                 for p in range(self.getNumEquations()):
                    tmp=Xi*Y[p]
                    for i in range(self.getNumEquations()):
                       for j in range(self.getDim()):
                          if not C.isEmpty() and not B.isEmpty():
                             Xout[i,j]+=tmp*(C[p,i,j]-B[p,j,i])
                             # Xout=X_out+Xi*util.inner(Y,C-util.reorderComponents(B,[0,2,1]),offset=1)
                          elif C.isEmpty():
                             Xout[i,j]-=tmp*B[p,j,i]
                             # Xout=X_out-Xi*util.inner(Y,util.reorderComponents(B,[0,2,1]),offset=1)
                          else:
                             Xout[i,j]+=tmp*C[p,i,j]
                             # Xout=X_out+Xi*util.inner(Y,C,offset=1)
            else:
                 tmp=Xi*Y
                 for j in range(self.getDim()):
                    if not C.isEmpty() and not B.isEmpty():
                       Xout[j]+=tmp*(C[j]-B[j])
                       # Xout=Xout+Xi*Y*(C-B)
                    elif C.isEmpty():
                       Xout[j]-=tmp*B[j]
                       # Xout=Xout-Xi*Y*B
                    else:
                       Xout[j]+=tmp*C[j]
                       # Xout=Xout+Xi*Y*C
         return Xout
     elif name == "Y" :
         return self.getCoefficient("Y")
     elif name == "d" :
         return self.getCoefficient("d")
     elif name == "y" :
         return self.getCoefficient("y")
     elif name == "d_contact" :
         return self.getCoefficient("d_contact")
     elif name == "y_contact" :
         return self.getCoefficient("y_contact")
     elif name == "r" :
         return self.getCoefficient("r")
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

class AdvectionDiffusion(LinearPDE):
   """
   Class to define PDE equation of the unisotropic advection-diffusion problem, which is genear L{LinearPDE} of the form

   M{S{omega}*u + inner(v,grad(u))- grad(matrixmult(k,grad(u))[j])[j] = f}

   with natural boundary conditons

   M{n[j]*matrixmult(k,grad(u))[j] = g- S{alpha}u }

   and constraints:

   M{u=r} where M{q>0}

   @remark: no upwinding is applied yet.

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
                        "k": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_DIM,PDECoefficient.BY_DIM),PDECoefficient.OPERATOR),
                        "f": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "v": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.BY_DIM,),PDECoefficient.OPERATOR),
                        "alpha": PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.OPERATOR),
                        "g": PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.BY_EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                        "r": PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH),
                        "q": PDECoefficient(PDECoefficient.SOLUTION,(PDECoefficient.BY_EQUATION,),PDECoefficient.BOTH)}

   def setValue(self,**coefficients):
     """
     sets new values to coefficients

     @param coefficients: new values assigned to coefficients
     @keyword omega: value for coefficient M{S{omega}}
     @type omega: any type that can be casted to L{Scalar<escript.Scalar>} object on L{Function<escript.Function>}.
     @keyword k: value for coefficient M{k}
     @type k: any type that can be casted to L{Tensor<escript.Tensor>} object on L{Function<escript.Function>}.
     @keyword v: value for coefficient M{v}
     @type v: any type that can be casted to L{Vector<escript.Vector>} object on L{Function<escript.Function>}.
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
         return escript.getCoefficient("v")
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
     elif name == "r" :
         return self.getCoefficient("r")
     elif name == "q" :
         return self.getCoefficient("q")
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name


# $Log$
# Revision 1.14  2005/09/22 01:54:57  jgs
# Merge of development branch dev-02 back to main trunk on 2005-09-22
#
# Revision 1.13  2005/09/15 03:44:19  jgs
# Merge of development branch dev-02 back to main trunk on 2005-09-15
#
# Revision 1.12  2005/09/01 03:31:28  jgs
# Merge of development branch dev-02 back to main trunk on 2005-09-01
#
# Revision 1.11  2005/08/23 01:24:28  jgs
# Merge of development branch dev-02 back to main trunk on 2005-08-23
#
# Revision 1.10  2005/08/12 01:45:36  jgs
# erge of development branch dev-02 back to main trunk on 2005-08-12
#
# Revision 1.9.2.17  2005/09/21 07:03:33  matt
# PDECoefficient and LinearPDE are now new style classes (introduced in Python
# 2.2). Classes Poisson, Helmholtz, LameEquation and AdvectivePDE have been
# modified to instead use portable/cooperative "super" calls to extend base
# class methods.
#
# Revision 1.9.2.16  2005/09/16 01:54:37  matt
# Removed redundant if-loop.
#
# Revision 1.9.2.15  2005/09/14 08:09:18  matt
# Added "REDUCED" solution PDECoefficient descriptors for LinearPDEs.
#
# Revision 1.9.2.14  2005/09/07 06:26:16  gross
# the solver from finley are put into the standalone package paso now
#
# Revision 1.9.2.13  2005/08/31 08:45:03  gross
# in the case of lumping no new system is allocated if the constraint is changed.
#
# Revision 1.9.2.12  2005/08/31 07:10:23  gross
# test for Lumping added
#
# Revision 1.9.2.11  2005/08/30 01:53:45  gross
# bug in format fixed.
#
# Revision 1.9.2.10  2005/08/26 07:14:17  gross
# a few more bugs in linearPDE fixed. remaining problem are finley problems
#
# Revision 1.9.2.9  2005/08/26 06:30:45  gross
# fix for reported bug  0000004. test_linearPDE passes a few more tests
#
# Revision 1.9.2.8  2005/08/26 04:30:13  gross
# gneric unit testing for linearPDE
#
# Revision 1.9.2.7  2005/08/25 07:06:50  gross
# linearPDE documentation is parsed now by epydoc. there is still a problem with links into escriptcpp.so
#
# Revision 1.9.2.6  2005/08/24 05:01:24  gross
# problem with resetting the matrix in case of resetting its values to 0 fixed.
#
# Revision 1.9.2.5  2005/08/24 02:03:28  gross
# epydoc mark up partially fixed
#
# Revision 1.9.2.4  2005/08/22 07:11:09  gross
# some problems with LinearPDEs fixed.
#
# Revision 1.9.2.3  2005/08/18 04:48:48  gross
# the methods SetLumping*() are removed. Lumping is set trough setSolverMethod(LinearPDE.LUMPING)
#
# Revision 1.9.2.2  2005/08/18 04:39:32  gross
# the constants have been removed from util.py as they not needed anymore. PDE related constants are accessed through LinearPDE attributes now
#
# Revision 1.9.2.1  2005/07/29 07:10:27  gross
# new functions in util and a new pde type in linearPDEs
#
# Revision 1.1.2.25  2005/07/28 04:21:09  gross
# Lame equation: (linear elastic, isotropic) added
#
# Revision 1.1.2.24  2005/07/22 06:37:11  gross
# some extensions to modellib and linearPDEs
#
# Revision 1.1.2.23  2005/05/13 00:55:20  cochrane
# Fixed up some docstrings.  Moved module-level functions to top of file so
# that epydoc and doxygen can pick them up properly.
#
# Revision 1.1.2.22  2005/05/12 11:41:30  gross
# some basic Models have been added
#
# Revision 1.1.2.21  2005/05/12 07:16:12  cochrane
# Moved ELMAN_RAMAGE, SIMPLIFIED_BROOK_HUGHES, and HALF functions to bottom of
# file so that the AdvectivePDE class is picked up by doxygen.  Some
# reformatting of docstrings.  Addition of code to make equations come out
# as proper LaTeX.
#
# Revision 1.1.2.20  2005/04/15 07:09:08  gross
# some problems with functionspace and linearPDEs fixed.
#
# Revision 1.1.2.19  2005/03/04 05:27:07  gross
# bug in SystemPattern fixed.
#
# Revision 1.1.2.18  2005/02/08 06:16:45  gross
# Bugs in AdvectivePDE fixed, AdvectiveTest is stable but more testing is needed
#
# Revision 1.1.2.17  2005/02/08 05:56:19  gross
# Reference Number handling added
#
# Revision 1.1.2.16  2005/02/07 04:41:28  gross
# some function exposed to python to make mesh merging running
#
# Revision 1.1.2.15  2005/02/03 00:14:44  gross
# timeseries add and ESySParameter.py renames esysXML.py for consistence
#
# Revision 1.1.2.14  2005/02/01 06:44:10  gross
# new implementation of AdvectivePDE which now also updates right hand side. systems of PDEs are still not working
#
# Revision 1.1.2.13  2005/01/25 00:47:07  gross
# updates in the documentation
#
# Revision 1.1.2.12  2005/01/12 01:28:04  matt
# Added createCoefficient method for linearPDEs.
#
# Revision 1.1.2.11  2005/01/11 01:55:34  gross
# a problem in linearPDE class fixed
#
# Revision 1.1.2.10  2005/01/07 01:13:29  gross
# some bugs in linearPDE fixed
#
# Revision 1.1.2.9  2005/01/06 06:24:58  gross
# some bugs in slicing fixed
#
# Revision 1.1.2.8  2005/01/05 04:21:40  gross
# FunctionSpace checking/matchig in slicing added
#
# Revision 1.1.2.7  2004/12/29 10:03:41  gross
# bug in setValue fixed
#
# Revision 1.1.2.6  2004/12/29 05:29:59  gross
# AdvectivePDE successfully tested for Peclet number 1000000. there is still a problem with setValue and Data()
#
# Revision 1.1.2.5  2004/12/29 00:18:41  gross
# AdvectivePDE added
#
# Revision 1.1.2.4  2004/12/24 06:05:41  gross
# some changes in linearPDEs to add AdevectivePDE
#
# Revision 1.1.2.3  2004/12/16 00:12:34  gross
# __init__ of LinearPDE does not accept any coefficient anymore
#
# Revision 1.1.2.2  2004/12/14 03:55:01  jgs
# *** empty log message ***
#
# Revision 1.1.2.1  2004/12/12 22:53:47  gross
# linearPDE has been renamed LinearPDE
#
# Revision 1.1.1.1.2.7  2004/12/07 10:13:08  gross
# GMRES added
#
# Revision 1.1.1.1.2.6  2004/12/07 03:19:50  gross
# options for GMRES and PRES20 added
#
# Revision 1.1.1.1.2.5  2004/12/01 06:25:15  gross
# some small changes
#
# Revision 1.1.1.1.2.4  2004/11/24 01:50:21  gross
# Finley solves 4M unknowns now
#
# Revision 1.1.1.1.2.3  2004/11/15 06:05:26  gross
# poisson solver added
#
# Revision 1.1.1.1.2.2  2004/11/12 06:58:15  gross
# a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
#
# Revision 1.1.1.1.2.1  2004/10/28 22:59:22  gross
# finley's RecTest.py is running now: problem in SystemMatrixAdapater fixed
#
# Revision 1.1.1.1  2004/10/26 06:53:56  jgs
# initial import of project esys2
#
# Revision 1.3.2.3  2004/10/26 06:43:48  jgs
# committing Lutz's and Paul's changes to brach jgs
#
# Revision 1.3.4.1  2004/10/20 05:32:51  cochrane
# Added incomplete Doxygen comments to files, or merely put the docstrings that already exist into Doxygen form.
#
# Revision 1.3  2004/09/23 00:53:23  jgs
# minor fixes
#
# Revision 1.1  2004/08/28 12:58:06  gross
# SimpleSolve is not running yet: problem with == of functionsspace
#
