
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
The module provides an interface to define and solve linear partial
differential equations (PDEs) and Transport problems within L{escript}.
L{linearPDEs} does not provide any solver capabilities in itself but hands the
PDE over to the PDE solver library defined through the L{Domain<escript.Domain>}
of the PDE. The general interface is provided through the L{LinearPDE} class.
L{TransportProblem} provides an interface to initial value problems dominated
by its advective terms.

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

import math
import escript
import util
import numpy

__author__="Lutz Gross, l.gross@uq.edu.au"


class IllegalCoefficient(ValueError):
   """
   Exception that is raised if an illegal coefficient of the general or
   particular PDE is requested.
   """
   pass

class IllegalCoefficientValue(ValueError):
   """
   Exception that is raised if an incorrect value for a coefficient is used.
   """
   pass

class IllegalCoefficientFunctionSpace(ValueError):
   """
   Exception that is raised if an incorrect function space for a coefficient
   is used.
   """

class UndefinedPDEError(ValueError):
   """
   Exception that is raised if a PDE is not fully defined yet.
   """
   pass

class PDECoef(object):
    """
    A class for describing a PDE coefficient.

    @cvar INTERIOR: indicator that coefficient is defined on the interior of
                    the PDE domain
    @cvar BOUNDARY: indicator that coefficient is defined on the boundary of
                    the PDE domain
    @cvar CONTACT: indicator that coefficient is defined on the contact region
                   within the PDE domain
    @cvar INTERIOR_REDUCED: indicator that coefficient is defined on the
                            interior of the PDE domain using a reduced
                            integration order
    @cvar BOUNDARY_REDUCED: indicator that coefficient is defined on the
                            boundary of the PDE domain using a reduced
                            integration order
    @cvar CONTACT_REDUCED: indicator that coefficient is defined on the contact
                           region within the PDE domain using a reduced
                           integration order
    @cvar SOLUTION: indicator that coefficient is defined through a solution of
                    the PDE
    @cvar REDUCED: indicator that coefficient is defined through a reduced
                   solution of the PDE
    @cvar BY_EQUATION: indicator that the dimension of the coefficient shape is
                       defined by the number of PDE equations
    @cvar BY_SOLUTION: indicator that the dimension of the coefficient shape is
                       defined by the number of PDE solutions
    @cvar BY_DIM: indicator that the dimension of the coefficient shape is
                  defined by the spatial dimension
    @cvar OPERATOR: indicator that the the coefficient alters the operator of
                    the PDE
    @cvar RIGHTHANDSIDE: indicator that the the coefficient alters the right
                         hand side of the PDE
    @cvar BOTH: indicator that the the coefficient alters the operator as well
                as the right hand side of the PDE

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
       Initialises a PDE coefficient type.

       @param where: describes where the coefficient lives
       @type where: one of L{INTERIOR}, L{BOUNDARY}, L{CONTACT}, L{SOLUTION},
                    L{REDUCED}, L{INTERIOR_REDUCED}, L{BOUNDARY_REDUCED},
                    L{CONTACT_REDUCED}
       @param pattern: describes the shape of the coefficient and how the shape
                       is built for a given spatial dimension and numbers of
                       equations and solutions in then PDE. For instance,
                       (L{BY_EQUATION},L{BY_SOLUTION},L{BY_DIM}) describes a
                       rank 3 coefficient which is instantiated as shape (3,2,2)
                       in case of three equations and two solution components
                       on a 2-dimensional domain. In the case of single equation
                       and a single solution component the shape components
                       marked by L{BY_EQUATION} or L{BY_SOLUTION} are dropped.
                       In this case the example would be read as (2,).
       @type pattern: C{tuple} of L{BY_EQUATION}, L{BY_SOLUTION}, L{BY_DIM}
       @param altering: indicates what part of the PDE is altered if the
                        coefficient is altered
       @type altering: one of L{OPERATOR}, L{RIGHTHANDSIDE}, L{BOTH}
       """
       super(PDECoef, self).__init__()
       self.what=where
       self.pattern=pattern
       self.altering=altering
       self.resetValue()

    def resetValue(self):
       """
       Resets the coefficient value to the default.
       """
       self.value=escript.Data()

    def getFunctionSpace(self,domain,reducedEquationOrder=False,reducedSolutionOrder=False):
       """
       Returns the L{FunctionSpace<escript.FunctionSpace>} of the coefficient.

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param reducedEquationOrder: True to indicate that reduced order is used
                                    to represent the equation
       @type reducedEquationOrder: C{bool}
       @param reducedSolutionOrder: True to indicate that reduced order is used
                                    to represent the solution
       @type reducedSolutionOrder: C{bool}
       @return: L{FunctionSpace<escript.FunctionSpace>} of the coefficient
       @rtype: L{FunctionSpace<escript.FunctionSpace>}
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
       Returns the value of the coefficient.

       @return: value of the coefficient
       @rtype: L{Data<escript.Data>}
       """
       return self.value

    def setValue(self,domain,numEquations=1,numSolutions=1,reducedEquationOrder=False,reducedSolutionOrder=False,newValue=None):
       """
       Sets the value of the coefficient to a new value.

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param numEquations: number of equations of the PDE
       @type numEquations: C{int}
       @param numSolutions: number of components of the PDE solution
       @type numSolutions: C{int}
       @param reducedEquationOrder: True to indicate that reduced order is used
                                    to represent the equation
       @type reducedEquationOrder: C{bool}
       @param reducedSolutionOrder: True to indicate that reduced order is used
                                    to represent the solution
       @type reducedSolutionOrder: C{bool}
       @param newValue: number of components of the PDE solution
       @type newValue: any object that can be converted into a
                       L{Data<escript.Data>} object with the appropriate shape
                       and L{FunctionSpace<escript.FunctionSpace>}
       @raise IllegalCoefficientValue: if the shape of the assigned value does
                                       not match the shape of the coefficient
       @raise IllegalCoefficientFunctionSpace: if unable to interpolate value
                                               to appropriate function space
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
        Checks if the coefficient alters the operator of the PDE.

        @return: True if the operator of the PDE is changed when the
                 coefficient is changed
        @rtype: C{bool}
        """
        if self.altering==self.OPERATOR or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isAlteringRightHandSide(self):
        """
        Checks if the coefficient alters the right hand side of the PDE.

        @rtype: C{bool}
        @return: True if the right hand side of the PDE is changed when the
                 coefficient is changed, C{None} otherwise.
        """
        if self.altering==self.RIGHTHANDSIDE or self.altering==self.BOTH:
            return not None
        else:
            return None

    def estimateNumEquationsAndNumSolutions(self,domain,shape=()):
       """
       Tries to estimate the number of equations and number of solutions if
       the coefficient has the given shape.

       @param domain: domain on which the PDE uses the coefficient
       @type domain: L{Domain<escript.Domain>}
       @param shape: suggested shape of the coefficient
       @type shape: C{tuple} of C{int} values
       @return: the number of equations and number of solutions of the PDE if
                the coefficient has given shape. If no appropriate numbers
                could be identified, C{None} is returned
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
       Checks if the coefficient allows to estimate the number of solution
       components.

       @return: True if the coefficient allows an estimate of the number of
                solution components, False otherwise
       @rtype: C{bool}
       """
       for i in self.pattern:
             if i==self.BY_SOLUTION: return True
       return False

    def definesNumEquation(self):
       """
       Checks if the coefficient allows to estimate the number of equations.

       @return: True if the coefficient allows an estimate of the number of
                equations, False otherwise
       @rtype: C{bool}
       """
       for i in self.pattern:
             if i==self.BY_EQUATION: return True
       return False

    def __CompTuple2(self,t1,t2):
      """
      Compares two tuples of possible number of equations and number of
      solutions.

      @param t1: the first tuple
      @param t2: the second tuple
      @return: 0, 1, or -1
      """

      dif=t1[0]+t1[1]-(t2[0]+t2[1])
      if dif<0: return 1
      elif dif>0: return -1
      else: return 0

    def getShape(self,domain,numEquations=1,numSolutions=1):
       """
       Builds the required shape of the coefficient.

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

class LinearProblem(object):
   """
   This is the base class to define a general linear PDE-type problem for
   for an unknown function M{u} on a given domain defined through a
   L{Domain<escript.Domain>} object. The problem can be given as a single
   equation or as a system of equations.

   The class assumes that some sort of assembling process is required to form
   a problem of the form

   M{L u=f}

   where M{L} is an operator and M{f} is the right hand side. This operator
   problem will be solved to get the unknown M{u}.

   @cvar DEFAULT: The default method used to solve the system of linear
                  equations
   @cvar DIRECT: The direct solver based on LDU factorization
   @cvar CHOLEVSKY: The direct solver based on LDLt factorization (can only be
                    applied for symmetric PDEs)
   @cvar PCG: The preconditioned conjugate gradient method (can only be applied
              for symmetric PDEs)
   @cvar CR: The conjugate residual method
   @cvar CGS: The conjugate gradient square method
   @cvar BICGSTAB: The stabilized BiConjugate Gradient method
   @cvar TFQMR: Transport Free Quasi Minimal Residual method
   @cvar MINRES: Minimum residual method
   @cvar SSOR: The symmetric overrelaxation method
   @cvar ILU0: The incomplete LU factorization preconditioner with no fill-in
   @cvar ILUT: The incomplete LU factorization preconditioner with fill-in
   @cvar JACOBI: The Jacobi preconditioner
   @cvar GMRES: The Gram-Schmidt minimum residual method
   @cvar PRES20: Special GMRES with restart after 20 steps and truncation after
                 5 residuals
   @cvar LUMPING: Matrix lumping
   @cvar NO_REORDERING: No matrix reordering allowed
   @cvar MINIMUM_FILL_IN: Reorder matrix to reduce fill-in during factorization
   @cvar NESTED_DISSECTION: Reorder matrix to improve load balancing during
                            factorization
   @cvar PASO: PASO solver package
   @cvar SCSL: SGI SCSL solver library
   @cvar MKL: Intel's MKL solver library
   @cvar UMFPACK: The UMFPACK library
   @cvar TRILINOS: The TRILINOS parallel solver class library from Sandia Natl
                   Labs
   @cvar ITERATIVE: The default iterative solver
   @cvar AMG: Algebraic Multi Grid
   @cvar RILU: Recursive ILU
   @cvar GS: Gauss-Seidel solver

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
   TRILINOS = 24
   NONLINEAR_GMRES = 25
   TFQMR = 26
   MINRES = 27
   GS=28

   SMALL_TOLERANCE=1.e-13
   PACKAGE_KEY="package"
   METHOD_KEY="method"
   SYMMETRY_KEY="symmetric"
   TOLERANCE_KEY="tolerance"
   PRECONDITIONER_KEY="preconditioner"


   def __init__(self,domain,numEquations=None,numSolutions=None,debug=False):
     """
     Initializes a linear problem.

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param numEquations: number of equations. If C{None} the number of
                          equations is extracted from the coefficients.
     @param numSolutions: number of solution components. If C{None} the number
                          of solution components is extracted from the
                          coefficients.
     @param debug: if True debug information is printed

     """
     super(LinearProblem, self).__init__()

     self.__debug=debug
     self.__domain=domain
     self.__numEquations=numEquations
     self.__numSolutions=numSolutions
     self.__altered_coefficients=False
     self.__reduce_equation_order=False
     self.__reduce_solution_order=False
     self.__tolerance=1.e-8
     self.__solver_method=self.DEFAULT
     self.__solver_package=self.DEFAULT
     self.__preconditioner=self.DEFAULT
     self.__is_RHS_valid=False
     self.__is_operator_valid=False
     self.__sym=False
     self.__COEFFICIENTS={}
     # initialize things:
     self.resetAllCoefficients()
     self.__system_type=None
     self.updateSystemType()
   # ==========================================================================
   #    general stuff:
   # ==========================================================================
   def __str__(self):
     """
     Returns a string representation of the PDE.

     @return: a simple representation of the PDE
     @rtype: C{str}
     """
     return "<LinearProblem %d>"%id(self)
   # ==========================================================================
   #    debug :
   # ==========================================================================
   def setDebugOn(self):
     """
     Switches debug output on.
     """
     self.__debug=not None

   def setDebugOff(self):
     """
     Switches debug output off.
     """
     self.__debug=None

   def setDebug(self, flag):
     """
     Switches debug output on if C{flag} is True otherwise it is switched off.

     @param flag: desired debug status
     @type flag: C{bool}
     """
     if flag:
         self.setDebugOn()
     else:
         self.setDebugOff()

   def trace(self,text):
     """
     Prints the text message if debug mode is switched on.

     @param text: message to be printed
     @type text: C{string}
     """
     if self.__debug: print "%s: %s"%(str(self),text)

   # ==========================================================================
   # some service functions:
   # ==========================================================================
   def introduceCoefficients(self,**coeff):
       """
       Introduces new coefficients into the problem.

       Use::

       p.introduceCoefficients(A=PDECoef(...), B=PDECoef(...))

       to introduce the coefficients M{A} ans M{B}.
       """
       for name, type in coeff.items():
           if not isinstance(type,PDECoef):
              raise ValueError,"coefficient %s has no type."%name
           self.__COEFFICIENTS[name]=type
           self.__COEFFICIENTS[name].resetValue()
           self.trace("coefficient %s has been introduced."%name)

   def getDomain(self):
     """
     Returns the domain of the PDE.

     @return: the domain of the PDE
     @rtype: L{Domain<escript.Domain>}
     """
     return self.__domain

   def getDim(self):
     """
     Returns the spatial dimension of the PDE.

     @return: the spatial dimension of the PDE domain
     @rtype: C{int}
     """
     return self.getDomain().getDim()

   def getNumEquations(self):
     """
     Returns the number of equations.

     @return: the number of equations
     @rtype: C{int}
     @raise UndefinedPDEError: if the number of equations is not specified yet
     """
     if self.__numEquations==None:
         if self.__numSolutions==None:
            raise UndefinedPDEError,"Number of equations is undefined. Please specify argument numEquations."
         else:
            self.__numEquations=self.__numSolutions
     return self.__numEquations

   def getNumSolutions(self):
     """
     Returns the number of unknowns.

     @return: the number of unknowns
     @rtype: C{int}
     @raise UndefinedPDEError: if the number of unknowns is not specified yet
     """
     if self.__numSolutions==None:
        if self.__numEquations==None:
            raise UndefinedPDEError,"Number of solution is undefined. Please specify argument numSolutions."
        else:
            self.__numSolutions=self.__numEquations
     return self.__numSolutions

   def reduceEquationOrder(self):
     """
     Returns the status of order reduction for the equation.

     @return: True if reduced interpolation order is used for the
              representation of the equation, False otherwise
     @rtype: L{bool}
     """
     return self.__reduce_equation_order

   def reduceSolutionOrder(self):
     """
     Returns the status of order reduction for the solution.

     @return: True if reduced interpolation order is used for the
              representation of the solution, False otherwise
     @rtype: L{bool}
     """
     return self.__reduce_solution_order

   def getFunctionSpaceForEquation(self):
     """
     Returns the L{FunctionSpace<escript.FunctionSpace>} used to discretize
     the equation.

     @return: representation space of equation
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     """
     if self.reduceEquationOrder():
         return escript.ReducedSolution(self.getDomain())
     else:
         return escript.Solution(self.getDomain())

   def getFunctionSpaceForSolution(self):
     """
     Returns the L{FunctionSpace<escript.FunctionSpace>} used to represent
     the solution.

     @return: representation space of solution
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     """
     if self.reduceSolutionOrder():
         return escript.ReducedSolution(self.getDomain())
     else:
         return escript.Solution(self.getDomain())

   # ==========================================================================
   #   solver settings:
   # ==========================================================================
   def setSolverMethod(self,solver=None,preconditioner=None):
       """
       Sets a new solver method and/or preconditioner.

       @param solver: new solver method to use
       @type solver: one of L{DEFAULT}, L{ITERATIVE} L{DIRECT}, L{CHOLEVSKY},
                     L{PCG}, L{CR}, L{CGS}, L{BICGSTAB}, L{SSOR}, L{GMRES},
                     L{TFQMR}, L{MINRES}, L{PRES20}, L{LUMPING}, L{AMG}
       @param preconditioner: new preconditioner to use
       @type preconditioner: one of L{DEFAULT}, L{JACOBI} L{ILU0}, L{ILUT},
                             L{SSOR}, L{RILU}, L{GS}

       @note: the solver method chosen may not be available by the selected
              package.
       """
       if solver==None: solver=self.__solver_method
       if preconditioner==None: preconditioner=self.__preconditioner
       if solver==None: solver=self.DEFAULT
       if preconditioner==None: preconditioner=self.DEFAULT
       if not (solver,preconditioner)==self.getSolverMethod():
           self.__solver_method=solver
           self.__preconditioner=preconditioner
           self.updateSystemType()
           self.trace("New solver is %s"%self.getSolverMethodName())

   def getSolverMethodName(self):
       """
       Returns the name of the solver currently used.

       @return: the name of the solver currently used
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
       elif m[0]==self.TFQMR: method= "TFQMR"
       elif m[0]==self.MINRES: method= "MINRES"
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
       elif m[1]==self.GS: method+= "+GS"
       if p==self.DEFAULT: package="DEFAULT"
       elif p==self.PASO: package= "PASO"
       elif p==self.MKL: package= "MKL"
       elif p==self.SCSL: package= "SCSL"
       elif p==self.UMFPACK: package= "UMFPACK"
       elif p==self.TRILINOS: package= "TRILINOS"
       else : method="unknown"
       return "%s solver of %s package"%(method,package)

   def getSolverMethod(self):
       """
       Returns the solver method and preconditioner used.

       @return: the solver method currently used.
       @rtype: C{tuple} of C{int}
       """
       return self.__solver_method,self.__preconditioner

   def setSolverPackage(self,package=None):
       """
       Sets a new solver package.

       @param package: the new solver package
       @type package: one of L{DEFAULT}, L{PASO} L{SCSL}, L{MKL}, L{UMFPACK},
                      L{TRILINOS}
       """
       if package==None: package=self.DEFAULT
       if not package==self.getSolverPackage():
           self.__solver_package=package
           self.updateSystemType()
           self.trace("New solver is %s"%self.getSolverMethodName())

   def getSolverPackage(self):
       """
       Returns the package of the solver.

       @return: the solver package currently being used
       @rtype: C{int}
       """
       return self.__solver_package

   def isUsingLumping(self):
      """
      Checks if matrix lumping is the current solver method.

      @return: True if the current solver method is lumping
      @rtype: C{bool}
      """
      return self.getSolverMethod()[0]==self.LUMPING

   def setTolerance(self,tol=1.e-8):
       """
       Resets the tolerance for the solver method to C{tol}.

       @param tol: new tolerance for the solver. If C{tol} is lower than the
                   current tolerence the system will be resolved.
       @type tol: positive C{float}
       @raise ValueError: if tolerance is not positive
       """
       if not tol>0:
           raise ValueError,"Tolerance has to be positive"
       if tol<self.getTolerance(): self.invalidateSolution()
       self.trace("New tolerance %e"%tol)
       self.__tolerance=tol
       return

   def getTolerance(self):
       """
       Returns the tolerance currently set for the solution.

       @return: tolerance currently used
       @rtype: C{float}
       """
       return self.__tolerance

   # ==========================================================================
   #    symmetry  flag:
   # ==========================================================================
   def isSymmetric(self):
      """
      Checks if symmetry is indicated.

      @return: True if a symmetric PDE is indicated, False otherwise
      @rtype: C{bool}
      """
      return self.__sym

   def setSymmetryOn(self):
      """
      Sets the symmetry flag.
      """
      if not self.isSymmetric():
         self.trace("PDE is set to be symmetric")
         self.__sym=True
         self.updateSystemType()

   def setSymmetryOff(self):
      """
      Clears the symmetry flag.
      """
      if self.isSymmetric():
         self.trace("PDE is set to be nonsymmetric")
         self.__sym=False
         self.updateSystemType()

   def setSymmetryTo(self,flag=False):
      """
      Sets the symmetry flag to C{flag}.

      @param flag: If True, the symmetry flag is set otherwise reset.
      @type flag: C{bool}
      """
      if flag:
         self.setSymmetryOn()
      else:
         self.setSymmetryOff()

   # ==========================================================================
   # function space handling for the equation as well as the solution
   # ==========================================================================
   def setReducedOrderOn(self):
     """
     Switches reduced order on for solution and equation representation.

     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     Switches reduced order off for solution and equation representation

     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     Sets order reduction state for both solution and equation representation
     according to flag.

     @param flag: if True, the order reduction is switched on for both solution
                  and equation representation, otherwise or if flag is not
                  present order reduction is switched off
     @type flag: C{bool}
     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)


   def setReducedOrderForSolutionOn(self):
     """
     Switches reduced order on for solution representation.

     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if not self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Reduced order is used for solution representation.")
         self.__reduce_solution_order=True
         self.initializeSystem()

   def setReducedOrderForSolutionOff(self):
     """
     Switches reduced order off for solution representation

     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set.
     """
     if self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Full order is used to interpolate solution.")
         self.__reduce_solution_order=False
         self.initializeSystem()

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     Sets order reduction state for solution representation according to flag.

     @param flag: if flag is True, the order reduction is switched on for
                  solution representation, otherwise or if flag is not present
                  order reduction is switched off
     @type flag: C{bool}
     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()

   def setReducedOrderForEquationOn(self):
     """
     Switches reduced order on for equation representation.

     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if not self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Reduced order is used for test functions.")
         self.__reduce_equation_order=True
         self.initializeSystem()

   def setReducedOrderForEquationOff(self):
     """
     Switches reduced order off for equation representation.

     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError,"order cannot be altered after coefficients have been defined."
         self.trace("Full order is used for test functions.")
         self.__reduce_equation_order=False
         self.initializeSystem()

   def setReducedOrderForEquationTo(self,flag=False):
     """
     Sets order reduction state for equation representation according to flag.

     @param flag: if flag is True, the order reduction is switched on for
                  equation representation, otherwise or if flag is not present
                  order reduction is switched off
     @type flag: C{bool}
     @raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if flag:
        self.setReducedOrderForEquationOn()
     else:
        self.setReducedOrderForEquationOff()

   def updateSystemType(self):
     """
     Reassesses the system type and, if a new matrix is needed, resets the
     system.
     """
     new_system_type=self.getRequiredSystemType()
     if not new_system_type==self.__system_type:
         self.trace("Matrix type is now %d."%new_system_type)
         self.__system_type=new_system_type
         self.initializeSystem()

   def getSystemType(self):
      """
      Returns the current system type.
      """
      return self.__system_type

   def checkSymmetricTensor(self,name,verbose=True):
      """
      Tests a coefficient for symmetry.

      @param name: name of the coefficient
      @type name: C{str}
      @param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      @type verbose: C{bool}
      @return: True if coefficient C{name} is symmetric
      @rtype: C{bool}
      """
      A=self.getCoefficient(name)
      verbose=verbose or self.__debug
      out=True
      if not A.isEmpty():
         tol=util.Lsup(A)*self.SMALL_TOLERANCE
         s=A.getShape()
         if A.getRank() == 4:
            if s[0]==s[2] and s[1] == s[3]:
               for i in range(s[0]):
                  for j in range(s[1]):
                     for k in range(s[2]):
                        for l in range(s[3]):
                            if util.Lsup(A[i,j,k,l]-A[k,l,i,j])>tol:
                               if verbose: print "non-symmetric problem as %s[%d,%d,%d,%d]!=%s[%d,%d,%d,%d]"%(name,i,j,k,l,name,k,l,i,j)
                               out=False
            else:
                 if verbose: print "non-symmetric problem because of inappropriate shape %s of coefficient %s."%(s,name)
                 out=False
         elif A.getRank() == 2:
            if s[0]==s[1]:
               for j in range(s[0]):
                  for l in range(s[1]):
                     if util.Lsup(A[j,l]-A[l,j])>tol:
                        if verbose: print "non-symmetric problem because %s[%d,%d]!=%s[%d,%d]"%(name,j,l,name,l,j)
                        out=False
            else:
                 if verbose: print "non-symmetric problem because of inappropriate shape %s of coefficient %s."%(s,name)
                 out=False
         elif A.getRank() == 0:
            pass
         else:
             raise ValueError,"Cannot check rank %s of %s."%(A.getRank(),name)
      return out

   def checkReciprocalSymmetry(self,name0,name1,verbose=True):
      """
      Tests two coefficients for reciprocal symmetry.

      @param name0: name of the first coefficient
      @type name0: C{str}
      @param name1: name of the second coefficient
      @type name1: C{str}
      @param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed
      @type verbose: C{bool}
      @return: True if coefficients C{name0} and C{name1} are reciprocally
               symmetric.
      @rtype: C{bool}
      """
      B=self.getCoefficient(name0)
      C=self.getCoefficient(name1)
      verbose=verbose or self.__debug
      out=True
      if B.isEmpty() and not C.isEmpty():
         if verbose: print "non-symmetric problem because %s is not present but %s is"%(name0,name1)
         out=False
      elif not B.isEmpty() and C.isEmpty():
         if verbose: print "non-symmetric problem because %s is not present but %s is"%(name0,name1)
         out=False
      elif not B.isEmpty() and not C.isEmpty():
         sB=B.getShape()
         sC=C.getShape()
         tol=(util.Lsup(B)+util.Lsup(C))*self.SMALL_TOLERANCE/2.
         if len(sB) != len(sC):
             if verbose: print "non-symmetric problem because ranks of %s (=%s) and %s (=%s) are different."%(name0,len(sB),name1,len(sC))
             out=False
         else:
             if len(sB)==0:
               if util.Lsup(B-C)>tol:
                  if verbose: print "non-symmetric problem because %s!=%s"%(name0,name1)
                  out=False
             elif len(sB)==1:
               if sB[0]==sC[0]:
                  for j in range(sB[0]):
                     if util.Lsup(B[j]-C[j])>tol:
                        if verbose: print "non-symmetric PDE because %s[%d]!=%s[%d]"%(name0,j,name1,j)
                        out=False
               else:
                 if verbose: print "non-symmetric problem because of inappropriate shapes %s and %s of coefficients %s and %s, respectively."%(sB,sC,name0,name1)
             elif len(sB)==3:
               if sB[0]==sC[1] and sB[1]==sC[2] and sB[2]==sC[0]:
                   for i in range(sB[0]):
                      for j in range(sB[1]):
                         for k in range(sB[2]):
                            if util.Lsup(B[i,j,k]-C[k,i,j])>tol:
                                 if verbose: print "non-symmetric problem because %s[%d,%d,%d]!=%s[%d,%d,%d]"%(name0,i,j,k,name1,k,i,j)
                                 out=False
               else:
                 if verbose: print "non-symmetric problem because of inappropriate shapes %s and %s of coefficients %s and %s, respectively."%(sB,sC,name0,name1)
             else:
                 raise ValueError,"Cannot check rank %s of %s and %s."%(len(sB),name0,name1)
      return out

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient C{name}.

     @param name: name of the coefficient requested
     @type name: C{string}
     @return: the value of the coefficient
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if C{name} is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
         return self.__COEFFICIENTS[name].getValue()
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def hasCoefficient(self,name):
     """
     Returns True if C{name} is the name of a coefficient.

     @param name: name of the coefficient enquired
     @type name: C{string}
     @return: True if C{name} is the name of a coefficient of the general PDE,
              False otherwise
     @rtype: C{bool}
     """
     return self.__COEFFICIENTS.has_key(name)

   def createCoefficient(self, name):
     """
     Creates a L{Data<escript.Data>} object corresponding to coefficient
     C{name}.

     @return: the coefficient C{name} initialized to 0
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: if C{name} is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return escript.Data(0.,self.getShapeOfCoefficient(name),self.getFunctionSpaceForCoefficient(name))
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def getFunctionSpaceForCoefficient(self,name):
     """
     Returns the L{FunctionSpace<escript.FunctionSpace>} to be used for
     coefficient C{name}.

     @param name: name of the coefficient enquired
     @type name: C{string}
     @return: the function space to be used for coefficient C{name}
     @rtype: L{FunctionSpace<escript.FunctionSpace>}
     @raise IllegalCoefficient: if C{name} is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return self.__COEFFICIENTS[name].getFunctionSpace(self.getDomain())
     else:
        raise ValueError,"unknown coefficient %s requested"%name

   def getShapeOfCoefficient(self,name):
     """
     Returns the shape of the coefficient C{name}.

     @param name: name of the coefficient enquired
     @type name: C{string}
     @return: the shape of the coefficient C{name}
     @rtype: C{tuple} of C{int}
     @raise IllegalCoefficient: if C{name} is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return self.__COEFFICIENTS[name].getShape(self.getDomain(),self.getNumEquations(),self.getNumSolutions())
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def resetAllCoefficients(self):
     """
     Resets all coefficients to their default values.
     """
     for i in self.__COEFFICIENTS.iterkeys():
         self.__COEFFICIENTS[i].resetValue()

   def alteredCoefficient(self,name):
     """
     Announces that coefficient C{name} has been changed.

     @param name: name of the coefficient affected
     @type name: C{string}
     @raise IllegalCoefficient: if C{name} is not a coefficient of the PDE
     @note: if C{name} is q or r, the method will not trigger a rebuild of the
            system as constraints are applied to the solved system.
     """
     if self.hasCoefficient(name):
        self.trace("Coefficient %s has been altered."%name)
        if not ((name=="q" or name=="r") and self.isUsingLumping()):
           if self.__COEFFICIENTS[name].isAlteringOperator(): self.invalidateOperator()
           if self.__COEFFICIENTS[name].isAlteringRightHandSide(): self.invalidateRightHandSide()
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def validSolution(self):
       """
       Marks the solution as valid.
       """
       self.__is_solution_valid=True

   def invalidateSolution(self):
       """
       Indicates the PDE has to be resolved if the solution is requested.
       """
       self.trace("System will be resolved.")
       self.__is_solution_valid=False

   def isSolutionValid(self):
       """
       Returns True if the solution is still valid.
       """
       return self.__is_solution_valid

   def validOperator(self):
       """
       Marks the operator as valid.
       """
       self.__is_operator_valid=True

   def invalidateOperator(self):
       """
       Indicates the operator has to be rebuilt next time it is used.
       """
       self.trace("Operator will be rebuilt.")
       self.invalidateSolution()
       self.__is_operator_valid=False

   def isOperatorValid(self):
       """
       Returns True if the operator is still valid.
       """
       return self.__is_operator_valid

   def validRightHandSide(self):
       """
       Marks the right hand side as valid.
       """
       self.__is_RHS_valid=True

   def invalidateRightHandSide(self):
       """
       Indicates the right hand side has to be rebuilt next time it is used.
       """
       if self.isRightHandSideValid(): self.trace("Right hand side has to be rebuilt.")
       self.invalidateSolution()
       self.__is_RHS_valid=False

   def isRightHandSideValid(self):
       """
       Returns True if the operator is still valid.
       """
       return self.__is_RHS_valid

   def invalidateSystem(self):
       """
       Announces that everything has to be rebuilt.
       """
       if self.isRightHandSideValid(): self.trace("System has to be rebuilt.")
       self.invalidateSolution()
       self.invalidateOperator()
       self.invalidateRightHandSide()

   def isSystemValid(self):
       """
       Returns True if the system (including solution) is still vaild.
       """
       return self.isSolutionValid() and self.isOperatorValid() and self.isRightHandSideValid()

   def initializeSystem(self):
       """
       Resets the system clearing the operator, right hand side and solution.
       """
       self.trace("New System has been created.")
       self.__operator=escript.Operator()
       self.__righthandside=escript.Data()
       self.__solution=escript.Data()
       self.invalidateSystem()

   def getOperator(self):
     """
     Returns the operator of the linear problem.

     @return: the operator of the problem
     """
     return self.getSystem()[0]

   def getRightHandSide(self):
     """
     Returns the right hand side of the linear problem.

     @return: the right hand side of the problem
     @rtype: L{Data<escript.Data>}
     """
     return self.getSystem()[1]

   def createRightHandSide(self):
       """
       Returns an instance of a new right hand side.
       """
       self.trace("New right hand side is allocated.")
       if self.getNumEquations()>1:
           return escript.Data(0.,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),True)
       else:
           return escript.Data(0.,(),self.getFunctionSpaceForEquation(),True)

   def createSolution(self):
       """
       Returns an instance of a new solution.
       """
       self.trace("New solution is allocated.")
       if self.getNumSolutions()>1:
           return escript.Data(0.,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)
       else:
           return escript.Data(0.,(),self.getFunctionSpaceForSolution(),True)

   def resetSolution(self):
       """
       Sets the solution to zero.
       """
       if self.__solution.isEmpty():
           self.__solution=self.createSolution()
       else:
           self.__solution.setToZero()
           self.trace("Solution is reset to zero.")

   def setSolution(self,u):
       """
       Sets the solution assuming that makes the system valid.
       """
       self.__solution=u
       self.validSolution()

   def getCurrentSolution(self):
       """
       Returns the solution in its current state.
       """
       if self.__solution.isEmpty(): self.__solution=self.createSolution()
       return self.__solution

   def resetRightHandSide(self):
       """
       Sets the right hand side to zero.
       """
       if self.__righthandside.isEmpty():
           self.__righthandside=self.createRightHandSide()
       else:
           self.__righthandside.setToZero()
           self.trace("Right hand side is reset to zero.")

   def getCurrentRightHandSide(self):
       """
       Returns the right hand side in its current state.
       """
       if self.__righthandside.isEmpty(): self.__righthandside=self.createRightHandSide()
       return self.__righthandside

   def resetOperator(self):
       """
       Makes sure that the operator is instantiated and returns it initialized
       with zeros.
       """
       if self.__operator.isEmpty():
           if self.isUsingLumping():
               self.__operator=self.createSolution()
           else:
               self.__operator=self.createOperator()
       else:
           if self.isUsingLumping():
               self.__operator.setToZero()
           else:
               self.__operator.resetValues()
           self.trace("Operator reset to zero")

   def getCurrentOperator(self):
       """
       Returns the operator in its current state.
       """
       if self.__operator.isEmpty(): self.resetOperator()
       return self.__operator

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      @raise IllegalCoefficient: if an unknown coefficient keyword is used
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
                s=numpy.array(d).shape
            if s!=None:
                # get number of equations and number of unknowns:
                res=self.__COEFFICIENTS[i].estimateNumEquationsAndNumSolutions(self.getDomain(),s)
                if res==None:
                    raise IllegalCoefficientValue,"Illegal shape %s of coefficient %s"%(s,i)
                else:
                    if self.__numEquations==None: self.__numEquations=res[0]
                    if self.__numSolutions==None: self.__numSolutions=res[1]
      if self.__numEquations==None: raise UndefinedPDEError,"unidentified number of equations"
      if self.__numSolutions==None: raise UndefinedPDEError,"unidentified number of solutions"
      # now we check the shape of the coefficient if numEquations and numSolutions are set:
      for i,d in coefficients.iteritems():
        try:
           self.__COEFFICIENTS[i].setValue(self.getDomain(),
                     self.getNumEquations(),self.getNumSolutions(),
                     self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
           self.alteredCoefficient(i)
        except IllegalCoefficientFunctionSpace,m:
            # if the function space is wrong then we try the reduced version:
            i_red=i+"_reduced"
            if (not i_red in coefficients.keys()) and i_red in self.__COEFFICIENTS.keys():
                try:
                    self.__COEFFICIENTS[i_red].setValue(self.getDomain(),
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

   # ==========================================================================
   # methods that are typically overwritten when implementing a particular
   # linear problem
   # ==========================================================================
   def getRequiredSystemType(self):
      """
      Returns the system type which needs to be used by the current set up.

      @note: Typically this method is overwritten when implementing a
             particular linear problem.
      """
      return 0

   def createOperator(self):
       """
       Returns an instance of a new operator.

       @note: This method is overwritten when implementing a particular
              linear problem.
       """
       return escript.Operator()

   def checkSymmetry(self,verbose=True):
      """
      Tests the PDE for symmetry.

      @param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed
      @type verbose: C{bool}
      @return: True if the problem is symmetric
      @rtype: C{bool}
      @note: Typically this method is overwritten when implementing a
             particular linear problem.
      """
      out=True
      return out

   def getSolution(self,**options):
       """
       Returns the solution of the problem.

       @return: the solution
       @rtype: L{Data<escript.Data>}

       @note: This method is overwritten when implementing a particular
              linear problem.
       """
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       @return: the discrete version of the PDE
       @rtype: C{tuple} of L{Operator,<escript.Operator>} and L{Data<escript.Data>}.

       @note: This method is overwritten when implementing a particular
              linear problem.
       """
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())
#=====================

class LinearPDE(LinearProblem):
   """
   This class is used to define a general linear, steady, second order PDE
   for an unknown function M{u} on a given domain defined through a
   L{Domain<escript.Domain>} object.

   For a single PDE having a solution with a single component the linear PDE
   is defined in the following form:

   M{-(grad(A[j,l]+A_reduced[j,l])*grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)=-grad(X+X_reduced)[j,j]+(Y+Y_reduced)}

   where M{grad(F)} denotes the spatial derivative of M{F}. Einstein's
   summation convention, ie. summation over indexes appearing twice in a term
   of a sum performed, is used.
   The coefficients M{A}, M{B}, M{C}, M{D}, M{X} and M{Y} have to be specified
   through L{Data<escript.Data>} objects in L{Function<escript.Function>} and
   the coefficients M{A_reduced}, M{B_reduced}, M{C_reduced}, M{D_reduced},
   M{X_reduced} and M{Y_reduced} have to be specified through
   L{Data<escript.Data>} objects in L{ReducedFunction<escript.ReducedFunction>}.
   It is also allowed to use objects that can be converted into such
   L{Data<escript.Data>} objects. M{A} and M{A_reduced} are rank two, M{B},
   M{C}, M{X}, M{B_reduced}, M{C_reduced} and M{X_reduced} are rank one and
   M{D}, M{D_reduced}, M{Y} and M{Y_reduced} are scalar.

   The following natural boundary conditions are considered:

   M{n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u)+(d+d_reduced)*u=n[j]*(X[j]+X_reduced[j])+y}

   where M{n} is the outer normal field. Notice that the coefficients M{A},
   M{A_reduced}, M{B}, M{B_reduced}, M{X} and M{X_reduced} are defined in the
   PDE. The coefficients M{d} and M{y} are each a scalar in
   L{FunctionOnBoundary<escript.FunctionOnBoundary>} and the coefficients
   M{d_reduced} and M{y_reduced} are each a scalar in
   L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}.

   Constraints for the solution prescribe the value of the solution at certain
   locations in the domain. They have the form

   M{u=r} where M{q>0}

   M{r} and M{q} are each scalar where M{q} is the characteristic function
   defining where the constraint is applied. The constraints override any
   other condition set by the PDE or the boundary condition.

   The PDE is symmetrical if

   M{A[i,j]=A[j,i]}  and M{B[j]=C[j]} and M{A_reduced[i,j]=A_reduced[j,i]}
   and M{B_reduced[j]=C_reduced[j]}

   For a system of PDEs and a solution with several components the PDE has the
   form

   M{-grad((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])[j]+(C[i,k,l]+C_reduced[i,k,l])*grad(u[k])[l]+(D[i,k]+D_reduced[i,k]*u[k] =-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i] }

   M{A} and M{A_reduced} are of rank four, M{B}, M{B_reduced}, M{C} and
   M{C_reduced} are each of rank three, M{D}, M{D_reduced}, M{X_reduced} and
   M{X} are each of rank two and M{Y} and M{Y_reduced} are of rank one.
   The natural boundary conditions take the form:

   M{n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])+(d[i,k]+d_reduced[i,k])*u[k]=n[j]*(X[i,j]+X_reduced[i,j])+y[i]+y_reduced[i]}

   The coefficient M{d} is of rank two and M{y} is of rank one both in
   L{FunctionOnBoundary<escript.FunctionOnBoundary>}. The coefficients
   M{d_reduced} is of rank two and M{y_reduced} is of rank one both in
   L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}.

   Constraints take the form

   M{u[i]=r[i]}  where  M{q[i]>0}

   M{r} and M{q} are each rank one. Notice that at some locations not
   necessarily all components must have a constraint.

   The system of PDEs is symmetrical if

      - M{A[i,j,k,l]=A[k,l,i,j]}
      - M{A_reduced[i,j,k,l]=A_reduced[k,l,i,j]}
      - M{B[i,j,k]=C[k,i,j]}
      - M{B_reduced[i,j,k]=C_reduced[k,i,j]}
      - M{D[i,k]=D[i,k]}
      - M{D_reduced[i,k]=D_reduced[i,k]}
      - M{d[i,k]=d[k,i]}
      - M{d_reduced[i,k]=d_reduced[k,i]}

   L{LinearPDE} also supports solution discontinuities over a contact region
   in the domain. To specify the conditions across the discontinuity we are
   using the generalised flux M{J} which, in the case of a system of PDEs
   and several components of the solution, is defined as

   M{J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]-X[i,j]-X_reduced[i,j]}

   For the case of single solution component and single PDE M{J} is defined as

   M{J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u-X[i]-X_reduced[i]}

   In the context of discontinuities M{n} denotes the normal on the
   discontinuity pointing from side 0 towards side 1 calculated from
   L{getNormal<escript.FunctionSpace.getNormal>} of L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
   For a system of PDEs the contact condition takes the form

   M{n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]}

   where M{J0} and M{J1} are the fluxes on side 0 and side 1 of the
   discontinuity, respectively. M{jump(u)}, which is the difference of the
   solution at side 1 and at side 0, denotes the jump of M{u} across
   discontinuity along the normal calculated by L{jump<util.jump>}.
   The coefficient M{d_contact} is of rank two and M{y_contact} is of rank one
   both in L{FunctionOnContactZero<escript.FunctionOnContactZero>} or
   L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
   The coefficient M{d_contact_reduced} is of rank two and M{y_contact_reduced}
   is of rank one both in L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>}
   or L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}.
   In case of a single PDE and a single component solution the contact
   condition takes the form

   M{n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)}

   In this case the coefficient M{d_contact} and M{y_contact} are each scalar
   both in L{FunctionOnContactZero<escript.FunctionOnContactZero>} or
   L{FunctionOnContactOne<escript.FunctionOnContactOne>} and the coefficient
   M{d_contact_reduced} and M{y_contact_reduced} are each scalar both in
   L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>} or
   L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}.

   Typical usage::

       p = LinearPDE(dom)
       p.setValue(A=kronecker(dom), D=1, Y=0.5)
       u = p.getSolution()

   """

   def __init__(self,domain,numEquations=None,numSolutions=None,debug=False):
     """
     Initializes a new linear PDE.

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param numEquations: number of equations. If C{None} the number of
                          equations is extracted from the PDE coefficients.
     @param numSolutions: number of solution components. If C{None} the number
                          of solution components is extracted from the PDE
                          coefficients.
     @param debug: if True debug information is printed

     """
     super(LinearPDE, self).__init__(domain,numEquations,numSolutions,debug)
     #
     #   the coefficients of the PDE:
     #
     self.introduceCoefficients(
       A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       A_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),
       q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH) )

   def __str__(self):
     """
     Returns the string representation of the PDE.

     @return: a simple representation of the PDE
     @rtype: C{str}
     """
     return "<LinearPDE %d>"%id(self)

   def getRequiredSystemType(self):
      """
      Returns the system type which needs to be used by the current set up.
      """
      return self.getDomain().getSystemMatrixTypeId(self.getSolverMethod()[0],self.getSolverMethod()[1],self.getSolverPackage(),self.isSymmetric())

   def checkSymmetry(self,verbose=True):
      """
      Tests the PDE for symmetry.

      @param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      @type verbose: C{bool}
      @return: True if the PDE is symmetric
      @rtype: L{bool}
      @note: This is a very expensive operation. It should be used for
             degugging only! The symmetry flag is not altered.
      """
      out=True
      out=out and self.checkSymmetricTensor("A", verbose)
      out=out and self.checkSymmetricTensor("A_reduced", verbose)
      out=out and self.checkReciprocalSymmetry("B","C", verbose)
      out=out and self.checkReciprocalSymmetry("B_reduced","C_reduced", verbose)
      out=out and self.checkSymmetricTensor("D", verbose)
      out=out and self.checkSymmetricTensor("D_reduced", verbose)
      out=out and self.checkSymmetricTensor("d", verbose)
      out=out and self.checkSymmetricTensor("d_reduced", verbose)
      out=out and self.checkSymmetricTensor("d_contact", verbose)
      out=out and self.checkSymmetricTensor("d_contact_reduced", verbose)
      return out

   def createOperator(self):
       """
       Returns an instance of a new operator.
       """
       self.trace("New operator is allocated.")
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           self.getSystemType())

   def getSolution(self,**options):
       """
       Returns the solution of the PDE.

       @return: the solution
       @rtype: L{Data<escript.Data>}
       @param options: solver options
       @keyword verbose: True to get some information during PDE solution
       @type verbose: C{bool}
       @keyword reordering: reordering scheme to be used during elimination.
                            Allowed values are L{NO_REORDERING},
                            L{MINIMUM_FILL_IN} and L{NESTED_DISSECTION}
       @keyword iter_max: maximum number of iteration steps allowed
       @keyword drop_tolerance: threshold for dropping in L{ILUT}
       @keyword drop_storage: maximum of allowed memory in L{ILUT}
       @keyword truncation: maximum number of residuals in L{GMRES}
       @keyword restart: restart cycle length in L{GMRES}
       """
       if not self.isSolutionValid():
          mat,f=self.getSystem()
          if self.isUsingLumping():
             self.setSolution(f*1/mat)
          else:
             options[self.TOLERANCE_KEY]=self.getTolerance()
             options[self.METHOD_KEY]=self.getSolverMethod()[0]
             options[self.PRECONDITIONER_KEY]=self.getSolverMethod()[1]
             options[self.PACKAGE_KEY]=self.getSolverPackage()
             options[self.SYMMETRY_KEY]=self.isSymmetric()
             self.trace("PDE is resolved.")
             self.trace("solver options: %s"%str(options))
             self.setSolution(mat.solve(f,options))
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       @return: the discrete version of the PDE
       @rtype: C{tuple} of L{Operator,<escript.Operator>} and
               L{Data<escript.Data>}
       """
       if not self.isOperatorValid() or not self.isRightHandSideValid():
          if self.isUsingLumping():
              if not self.isOperatorValid():
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                      raise TypeError,"Lumped matrix requires same order for equations and unknowns"
                 if not self.getCoefficient("A").isEmpty():
                      raise ValueError,"coefficient A in lumped matrix may not be present."
                 if not self.getCoefficient("B").isEmpty():
                      raise ValueError,"coefficient B in lumped matrix may not be present."
                 if not self.getCoefficient("C").isEmpty():
                      raise ValueError,"coefficient C in lumped matrix may not be present."
                 if not self.getCoefficient("d_contact").isEmpty():
                      raise ValueError,"coefficient d_contact in lumped matrix may not be present."
                 if not self.getCoefficient("A_reduced").isEmpty():
                      raise ValueError,"coefficient A_reduced in lumped matrix may not be present."
                 if not self.getCoefficient("B_reduced").isEmpty():
                      raise ValueError,"coefficient B_reduced in lumped matrix may not be present."
                 if not self.getCoefficient("C_reduced").isEmpty():
                      raise ValueError,"coefficient C_reduced in lumped matrix may not be present."
                 if not self.getCoefficient("d_contact_reduced").isEmpty():
                      raise ValueError,"coefficient d_contact_reduced in lumped matrix may not be present."
                 D=self.getCoefficient("D")
                 d=self.getCoefficient("d")
                 D_reduced=self.getCoefficient("D_reduced")
                 d_reduced=self.getCoefficient("d_reduced")
                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escript.Data()
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escript.Data()

                 if not D_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        D_reduced_times_e=util.matrix_mult(D_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_reduced_times_e=D_reduced
                 else:
                    D_reduced_times_e=escript.Data()
                 if not d_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        d_reduced_times_e=util.matrix_mult(d_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_times_e=d_reduced
                 else:
                    d_reduced_times_e=escript.Data()

                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 if False and hasattr(self.getDomain(), "addPDEToLumpedSystem") :
                    self.getDomain().addPDEToLumpedSystem(operator, D_times_e, d_times_e)
                    self.getDomain().addPDEToLumpedSystem(operator, D_reduced_times_e, d_reduced_times_e)
                 else:
                    self.getDomain().addPDEToRHS(operator, \
                                                 escript.Data(), \
                                                 D_times_e, \
                                                 d_times_e,\
                                                 escript.Data())
                    self.getDomain().addPDEToRHS(operator, \
                                                 escript.Data(), \
                                                 D_reduced_times_e, \
                                                 d_reduced_times_e,\
                                                 escript.Data())
                 self.trace("New lumped operator has been built.")
              if not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.getDomain().addPDEToRHS(righthandside, \
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"),\
                               self.getCoefficient("y"),\
                               self.getCoefficient("y_contact"))
                 self.getDomain().addPDEToRHS(righthandside, \
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"),\
                               self.getCoefficient("y_reduced"),\
                               self.getCoefficient("y_contact_reduced"))
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
              self.insertConstraint(rhs_only=False)
              self.validOperator()
          else:
             if not self.isOperatorValid() and not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 self.getDomain().addPDEToSystem(operator,righthandside, \
                               self.getCoefficient("A"), \
                               self.getCoefficient("B"), \
                               self.getCoefficient("C"), \
                               self.getCoefficient("D"), \
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"), \
                               self.getCoefficient("d"), \
                               self.getCoefficient("y"), \
                               self.getCoefficient("d_contact"), \
                               self.getCoefficient("y_contact"))
                 self.getDomain().addPDEToSystem(operator,righthandside, \
                               self.getCoefficient("A_reduced"), \
                               self.getCoefficient("B_reduced"), \
                               self.getCoefficient("C_reduced"), \
                               self.getCoefficient("D_reduced"), \
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"), \
                               self.getCoefficient("d_reduced"), \
                               self.getCoefficient("y_reduced"), \
                               self.getCoefficient("d_contact_reduced"), \
                               self.getCoefficient("y_contact_reduced"))
                 self.insertConstraint(rhs_only=False)
                 self.trace("New system has been built.")
                 self.validOperator()
                 self.validRightHandSide()
             elif not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.getDomain().addPDEToRHS(righthandside,
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"),\
                               self.getCoefficient("y"),\
                               self.getCoefficient("y_contact"))
                 self.getDomain().addPDEToRHS(righthandside,
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"),\
                               self.getCoefficient("y_reduced"),\
                               self.getCoefficient("y_contact_reduced"))
                 self.insertConstraint(rhs_only=True)
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
             elif not self.isOperatorValid():
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 self.getDomain().addPDEToSystem(operator,escript.Data(), \
                            self.getCoefficient("A"), \
                            self.getCoefficient("B"), \
                            self.getCoefficient("C"), \
                            self.getCoefficient("D"), \
                            escript.Data(), \
                            escript.Data(), \
                            self.getCoefficient("d"), \
                            escript.Data(),\
                            self.getCoefficient("d_contact"), \
                            escript.Data())
                 self.getDomain().addPDEToSystem(operator,escript.Data(), \
                            self.getCoefficient("A_reduced"), \
                            self.getCoefficient("B_reduced"), \
                            self.getCoefficient("C_reduced"), \
                            self.getCoefficient("D_reduced"), \
                            escript.Data(), \
                            escript.Data(), \
                            self.getCoefficient("d_reduced"), \
                            escript.Data(),\
                            self.getCoefficient("d_contact_reduced"), \
                            escript.Data())
                 self.insertConstraint(rhs_only=False)
                 self.trace("New operator has been built.")
                 self.validOperator()
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

   def insertConstraint(self, rhs_only=False):
      """
      Applies the constraints defined by q and r to the PDE.

      @param rhs_only: if True only the right hand side is altered by the
                       constraint
      @type rhs_only: C{bool}
      """
      q=self.getCoefficient("q")
      r=self.getCoefficient("r")
      righthandside=self.getCurrentRightHandSide()
      operator=self.getCurrentOperator()

      if not q.isEmpty():
         if r.isEmpty():
            r_s=self.createSolution()
         else:
            r_s=r
         if not rhs_only and not operator.isEmpty():
             if self.isUsingLumping():
                 operator.copyWithMask(escript.Data(1.,q.getShape(),q.getFunctionSpace()),q)
             else:
                 row_q=escript.Data(q,self.getFunctionSpaceForEquation())
                 col_q=escript.Data(q,self.getFunctionSpaceForSolution())
                 u=self.createSolution()
                 u.copyWithMask(r_s,col_q)
                 righthandside-=operator*u
                 operator.nullifyRowsAndCols(row_q,col_q,1.)
         righthandside.copyWithMask(r_s,q)

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      @param coefficients: new values assigned to coefficients
      @keyword A: value for coefficient C{A}
      @type A: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword A_reduced: value for coefficient C{A_reduced}
      @type A_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword B: value for coefficient C{B}
      @type B: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword B_reduced: value for coefficient C{B_reduced}
      @type B_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword C: value for coefficient C{C}
      @type C: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword C_reduced: value for coefficient C{C_reduced}
      @type C_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword D: value for coefficient C{D}
      @type D: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword D_reduced: value for coefficient C{D_reduced}
      @type D_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword X: value for coefficient C{X}
      @type X: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword X_reduced: value for coefficient C{X_reduced}
      @type X_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword Y: value for coefficient C{Y}
      @type Y: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword Y_reduced: value for coefficient C{Y_reduced}
      @type Y_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.Function>}
      @keyword d: value for coefficient C{d}
      @type d: any type that can be cast to a L{Data<escript.Data>} object on
               L{FunctionOnBoundary<escript.FunctionOnBoundary>}
      @keyword d_reduced: value for coefficient C{d_reduced}
      @type d_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}
      @keyword y: value for coefficient C{y}
      @type y: any type that can be cast to a L{Data<escript.Data>} object on
               L{FunctionOnBoundary<escript.FunctionOnBoundary>}
      @keyword d_contact: value for coefficient C{d_contact}
      @type d_contact: any type that can be cast to a L{Data<escript.Data>}
                       object on L{FunctionOnContactOne<escript.FunctionOnContactOne>}
                       or L{FunctionOnContactZero<escript.FunctionOnContactZero>}
      @keyword d_contact_reduced: value for coefficient C{d_contact_reduced}
      @type d_contact_reduced: any type that can be cast to a L{Data<escript.Data>}
                               object on L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}
                               or L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>}
      @keyword y_contact: value for coefficient C{y_contact}
      @type y_contact: any type that can be cast to a L{Data<escript.Data>}
                       object on L{FunctionOnContactOne<escript.FunctionOnContactOne>}
                       or L{FunctionOnContactZero<escript.FunctionOnContactZero>}
      @keyword y_contact_reduced: value for coefficient C{y_contact_reduced}
      @type y_contact_reduced: any type that can be cast to a L{Data<escript.Data>}
                               object on L{ReducedFunctionOnContactOne<escript.FunctionOnContactOne>}
                               or L{ReducedFunctionOnContactZero<escript.FunctionOnContactZero>}
      @keyword r: values prescribed to the solution at the locations of
                  constraints
      @type r: any type that can be cast to a L{Data<escript.Data>} object on
               L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending on whether reduced order is used for the solution
      @keyword q: mask for location of constraints
      @type q: any type that can be cast to a L{Data<escript.Data>} object on
               L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending on whether reduced order is used for the
               representation of the equation
      @raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      super(LinearPDE,self).setValue(**coefficients)
      # check if the systrem is inhomogeneous:
      if len(coefficients)>0 and not self.isUsingLumping():
         q=self.getCoefficient("q")
         r=self.getCoefficient("r")
         if not q.isEmpty() and not r.isEmpty():
             if util.Lsup(q*r)>0.:
               self.trace("Inhomogeneous constraint detected.")
               self.invalidateSystem()


   def getResidual(self,u=None):
     """
     Returns the residual of u or the current solution if u is not present.

     @param u: argument in the residual calculation. It must be representable
               in L{self.getFunctionSpaceForSolution()}. If u is not present
               or equals C{None} the current solution is used.
     @type u: L{Data<escript.Data>} or None
     @return: residual of u
     @rtype: L{Data<escript.Data>}
     """
     if u==None:
        return self.getOperator()*self.getSolution()-self.getRightHandSide()
     else:
        return self.getOperator()*escript.Data(u,self.getFunctionSpaceForSolution())-self.getRightHandSide()

   def getFlux(self,u=None):
     """
     Returns the flux M{J} for a given M{u}.

     M{J[i,j]=(A[i,j,k,l]+A_reduced[A[i,j,k,l]]*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])u[k]-X[i,j]-X_reduced[i,j]}

     or

     M{J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[l]+(B[j]+B_reduced[j])u-X[j]-X_reduced[j]}

     @param u: argument in the flux. If u is not present or equals L{None} the
               current solution is used.
     @type u: L{Data<escript.Data>} or None
     @return: flux
     @rtype: L{Data<escript.Data>}
     """
     if u==None: u=self.getSolution()
     return util.tensormult(self.getCoefficient("A"),util.grad(u,Funtion(self.getDomain))) \
           +util.matrixmult(self.getCoefficient("B"),u) \
           -util.self.getCoefficient("X") \
           +util.tensormult(self.getCoefficient("A_reduced"),util.grad(u,ReducedFuntion(self.getDomain))) \
           +util.matrixmult(self.getCoefficient("B_reduced"),u) \
           -util.self.getCoefficient("X_reduced")


class Poisson(LinearPDE):
   """
   Class to define a Poisson equation problem. This is generally a
   L{LinearPDE} of the form

   M{-grad(grad(u)[j])[j] = f}

   with natural boundary conditions

   M{n[j]*grad(u)[j] = 0 }

   and constraints:

   M{u=0} where M{q>0}

   """

   def __init__(self,domain,debug=False):
     """
     Initializes a new Poisson equation.

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param debug: if True debug information is printed

     """
     super(Poisson, self).__init__(domain,1,1,debug)
     self.introduceCoefficients(
                        f=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        f_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     Sets new values to coefficients.

     @param coefficients: new values assigned to coefficients
     @keyword f: value for right hand side M{f}
     @type f: any type that can be cast to a L{Scalar<escript.Scalar>} object
              on L{Function<escript.Function>}
     @keyword q: mask for location of constraints
     @type q: any type that can be cast to a rank zero L{Data<escript.Data>}
              object on L{Solution<escript.Solution>} or
              L{ReducedSolution<escript.ReducedSolution>} depending on whether
              reduced order is used for the representation of the equation
     @raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(Poisson, self).setValue(**coefficients)


   def getCoefficient(self,name):
     """
     Returns the value of the coefficient C{name} of the general PDE.

     @param name: name of the coefficient requested
     @type name: C{string}
     @return: the value of the coefficient C{name}
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: invalid coefficient name
     @note: This method is called by the assembling routine to map the Poisson
            equation onto the general PDE.
     """
     if name == "A" :
         return escript.Data(util.kronecker(self.getDim()),escript.Function(self.getDomain()))
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     else:
         return super(Poisson, self).getCoefficient(name)

class Helmholtz(LinearPDE):
   """
   Class to define a Helmholtz equation problem. This is generally a
   L{LinearPDE} of the form

   M{S{omega}*u - grad(k*grad(u)[j])[j] = f}

   with natural boundary conditions

   M{k*n[j]*grad(u)[j] = g- S{alpha}u }

   and constraints:

   M{u=r} where M{q>0}

   """

   def __init__(self,domain,debug=False):
     """
     Initializes a new Helmholtz equation.

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param debug: if True debug information is printed

     """
     super(Helmholtz, self).__init__(domain,1,1,debug)
     self.introduceCoefficients(
                        omega=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.OPERATOR),
                        k=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.OPERATOR),
                        f=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        f_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        alpha=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.OPERATOR),
                        g=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        g_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     Sets new values to coefficients.

     @param coefficients: new values assigned to coefficients
     @keyword omega: value for coefficient M{S{omega}}
     @type omega: any type that can be cast to a L{Scalar<escript.Scalar>}
                  object on L{Function<escript.Function>}
     @keyword k: value for coefficient M{k}
     @type k: any type that can be cast to a L{Scalar<escript.Scalar>} object
              on L{Function<escript.Function>}
     @keyword f: value for right hand side M{f}
     @type f: any type that can be cast to a L{Scalar<escript.Scalar>} object
              on L{Function<escript.Function>}
     @keyword alpha: value for right hand side M{S{alpha}}
     @type alpha: any type that can be cast to a L{Scalar<escript.Scalar>}
                  object on L{FunctionOnBoundary<escript.FunctionOnBoundary>}
     @keyword g: value for right hand side M{g}
     @type g: any type that can be cast to a L{Scalar<escript.Scalar>} object
              on L{FunctionOnBoundary<escript.FunctionOnBoundary>}
     @keyword r: prescribed values M{r} for the solution in constraints
     @type r: any type that can be cast to a L{Scalar<escript.Scalar>} object
              on L{Solution<escript.Solution>} or
              L{ReducedSolution<escript.ReducedSolution>} depending on whether
              reduced order is used for the representation of the equation
     @keyword q: mask for the location of constraints
     @type q: any type that can be cast to a L{Scalar<escript.Scalar>} object
              on L{Solution<escript.Solution>} or
              L{ReducedSolution<escript.ReducedSolution>} depending on whether
              reduced order is used for the representation of the equation
     @raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(Helmholtz, self).setValue(**coefficients)

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient C{name} of the general PDE.

     @param name: name of the coefficient requested
     @type name: C{string}
     @return: the value of the coefficient C{name}
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: invalid name
     """
     if name == "A" :
         if self.getCoefficient("k").isEmpty():
              return escript.Data(numpy.identity(self.getDim()),escript.Function(self.getDomain()))
         else:
              return escript.Data(numpy.identity(self.getDim()),escript.Function(self.getDomain()))*self.getCoefficient("k")
     elif name == "D" :
         return self.getCoefficient("omega")
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "d" :
         return self.getCoefficient("alpha")
     elif name == "y" :
         return self.getCoefficient("g")
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     elif name == "y_reduced" :
        return self.getCoefficient("g_reduced")
     else:
        return super(Helmholtz, self).getCoefficient(name)

class LameEquation(LinearPDE):
   """
   Class to define a Lame equation problem. This problem is defined as:

   M{-grad(S{mu}*(grad(u[i])[j]+grad(u[j])[i]))[j] - grad(S{lambda}*grad(u[k])[k])[j] = F_i -grad(S{sigma}[ij])[j] }

   with natural boundary conditions:

   M{n[j]*(S{mu}*(grad(u[i])[j]+grad(u[j])[i]) + S{lambda}*grad(u[k])[k]) = f_i +n[j]*S{sigma}[ij] }

   and constraints:

   M{u[i]=r[i]} where M{q[i]>0}

   """

   def __init__(self,domain,debug=False):
      """
      Initializes a new Lame equation.

      @param domain: domain of the PDE
      @type domain: L{Domain<escript.Domain>}
      @param debug: if True debug information is printed

      """
      super(LameEquation, self).__init__(domain,\
                                         domain.getDim(),domain.getDim(),debug)
      self.introduceCoefficients(lame_lambda=PDECoef(PDECoef.INTERIOR,(),PDECoef.OPERATOR),
                                 lame_mu=PDECoef(PDECoef.INTERIOR,(),PDECoef.OPERATOR),
                                 F=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                                 sigma=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
                                 f=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
      self.setSymmetryOn()

   def setValues(self,**coefficients):
     """
     Sets new values to coefficients.

     @param coefficients: new values assigned to coefficients
     @keyword lame_mu: value for coefficient M{S{mu}}
     @type lame_mu: any type that can be cast to a L{Scalar<escript.Scalar>}
                    object on L{Function<escript.Function>}
     @keyword lame_lambda: value for coefficient M{S{lambda}}
     @type lame_lambda: any type that can be cast to a L{Scalar<escript.Scalar>}
                        object on L{Function<escript.Function>}
     @keyword F: value for internal force M{F}
     @type F: any type that can be cast to a L{Vector<escript.Vector>} object
              on L{Function<escript.Function>}
     @keyword sigma: value for initial stress M{S{sigma}}
     @type sigma: any type that can be cast to a L{Tensor<escript.Tensor>}
                  object on L{Function<escript.Function>}
     @keyword f: value for external force M{f}
     @type f: any type that can be cast to a L{Vector<escript.Vector>} object
              on L{FunctionOnBoundary<escript.FunctionOnBoundary>}
     @keyword r: prescribed values M{r} for the solution in constraints
     @type r: any type that can be cast to a L{Vector<escript.Vector>} object
              on L{Solution<escript.Solution>} or
              L{ReducedSolution<escript.ReducedSolution>} depending on whether
              reduced order is used for the representation of the equation
     @keyword q: mask for the location of constraints
     @type q: any type that can be cast to a L{Vector<escript.Vector>} object
              on L{Solution<escript.Solution>} or
              L{ReducedSolution<escript.ReducedSolution>} depending on whether
              reduced order is used for the representation of the equation
     @raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(LameEquation, self).setValues(**coefficients)

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient C{name} of the general PDE.

     @param name: name of the coefficient requested
     @type name: C{string}
     @return: the value of the coefficient C{name}
     @rtype: L{Data<escript.Data>}
     @raise IllegalCoefficient: invalid coefficient name
     """
     out =self.createCoefficient("A")
     if name == "A" :
         if self.getCoefficient("lame_lambda").isEmpty(): 
            if self.getCoefficient("lame_mu").isEmpty():
                pass
            else:
                for i in range(self.getDim()):
                  for j in range(self.getDim()):
                    out[i,j,j,i] += self.getCoefficient("lame_mu")
                    out[i,j,i,j] += self.getCoefficient("lame_mu")
         else: 
            if self.getCoefficient("lame_mu").isEmpty():
                for i in range(self.getDim()):
                  for j in range(self.getDim()):
                    out[i,i,j,j] += self.getCoefficient("lame_lambda")
            else:
                for i in range(self.getDim()):
                  for j in range(self.getDim()):
                    out[i,i,j,j] += self.getCoefficient("lame_lambda")
                    out[i,j,j,i] += self.getCoefficient("lame_mu")
                    out[i,j,i,j] += self.getCoefficient("lame_mu")
         return out
     elif name == "X" :
         return self.getCoefficient("sigma")
     elif name == "Y" :
         return self.getCoefficient("F")
     elif name == "y" :
         return self.getCoefficient("f")
     else:
        return super(LameEquation, self).getCoefficient(name)

def LinearSinglePDE(domain,debug=False):
   """
   Defines a single linear PDE.

   @param domain: domain of the PDE
   @type domain: L{Domain<escript.Domain>}
   @param debug: if True debug information is printed
   @rtype: L{LinearPDE}
   """
   return LinearPDE(domain,numEquations=1,numSolutions=1,debug=debug)

def LinearPDESystem(domain,debug=False):
   """
   Defines a system of linear PDEs.

   @param domain: domain of the PDEs
   @type domain: L{Domain<escript.Domain>}
   @param debug: if True debug information is printed
   @rtype: L{LinearPDE}
   """
   return LinearPDE(domain,numEquations=domain.getDim(),numSolutions=domain.getDim(),debug=debug)


class TransportPDE(LinearProblem):
   """
   This class is used to define a transport problem given by a general linear,
   time dependent, second order PDE for an unknown, non-negative,
   time-dependent function M{u} on a given domain defined through a
   L{Domain<escript.Domain>} object.

   For a single equation with a solution with a single component the transport
   problem is defined in the following form:

   M{(M+M_reduced)*u_t=-(grad(A[j,l]+A_reduced[j,l])*grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)-grad(X+X_reduced)[j,j]+(Y+Y_reduced)}

   where M{u_t} denotes the time derivative of M{u} and M{grad(F)} denotes the
   spatial derivative of M{F}. Einstein's summation convention,  ie. summation
   over indexes appearing twice in a term of a sum performed, is used.
   The coefficients M{M}, M{A}, M{B}, M{C}, M{D}, M{X} and M{Y} have to be
   specified through L{Data<escript.Data>} objects in L{Function<escript.Function>}
   and the coefficients M{M_reduced}, M{A_reduced}, M{B_reduced}, M{C_reduced},
   M{D_reduced}, M{X_reduced} and M{Y_reduced} have to be specified through
   L{Data<escript.Data>} objects in L{ReducedFunction<escript.ReducedFunction>}.
   It is also allowed to use objects that can be converted into such
   L{Data<escript.Data>} objects. M{M} and M{M_reduced} are scalar, M{A} and
   M{A_reduced} are rank two, M{B}, M{C}, M{X}, M{B_reduced}, M{C_reduced} and
   M{X_reduced} are rank one and M{D}, M{D_reduced}, M{Y} and M{Y_reduced} are
   scalar.

   The following natural boundary conditions are considered:

   M{n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u+X[j]+X_reduced[j])+(d+d_reduced)*u+y+y_reduced=(m+m_reduced)*u_t}

   where M{n} is the outer normal field. Notice that the coefficients M{A},
   M{A_reduced}, M{B}, M{B_reduced}, M{X} and M{X_reduced} are defined in the
   transport problem. The coefficients M{m}, M{d} and M{y} are each a scalar in
   L{FunctionOnBoundary<escript.FunctionOnBoundary>} and the coefficients
   M{m_reduced}, M{d_reduced} and M{y_reduced} are each a scalar in
   L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}.

   Constraints for the solution prescribing the value of the solution at
   certain locations in the domain have the form

   M{u_t=r} where M{q>0}

   M{r} and M{q} are each scalar where M{q} is the characteristic function
   defining where the constraint is applied. The constraints override any other
   condition set by the transport problem or the boundary condition.

   The transport problem is symmetrical if

   M{A[i,j]=A[j,i]} and M{B[j]=C[j]} and M{A_reduced[i,j]=A_reduced[j,i]} and
   M{B_reduced[j]=C_reduced[j]}

   For a system and a solution with several components the transport problem
   has the form

   M{(M[i,k]+M_reduced[i,k])*u[k]_t=-grad((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])[j]+(C[i,k,l]+C_reduced[i,k,l])*grad(u[k])[l]+(D[i,k]+D_reduced[i,k]*u[k]-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i] }

   M{A} and M{A_reduced} are of rank four, M{B}, M{B_reduced}, M{C} and
   M{C_reduced} are each of rank three, M{M}, M{M_reduced}, M{D}, M{D_reduced},
   M{X_reduced} and M{X} are each of rank two and M{Y} and M{Y_reduced} are of
   rank one. The natural boundary conditions take the form:

   M{n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]+X[i,j]+X_reduced[i,j])+(d[i,k]+d_reduced[i,k])*u[k]+y[i]+y_reduced[i]= (m[i,k]+m_reduced[i,k])*u[k]_t}

   The coefficient M{d} and M{m} are of rank two and M{y} is of rank one with
   all in L{FunctionOnBoundary<escript.FunctionOnBoundary>}. The coefficients
   M{d_reduced} and M{m_reduced} are of rank two and M{y_reduced} is of rank
   one all in L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}.

   Constraints take the form

   M{u[i]_t=r[i]} where M{q[i]>0}

   M{r} and M{q} are each rank one. Notice that at some locations not
   necessarily all components must have a constraint.

   The transport problem is symmetrical if

      - M{M[i,k]=M[i,k]}
      - M{M_reduced[i,k]=M_reduced[i,k]}
      - M{A[i,j,k,l]=A[k,l,i,j]}
      - M{A_reduced[i,j,k,l]=A_reduced[k,l,i,j]}
      - M{B[i,j,k]=C[k,i,j]}
      - M{B_reduced[i,j,k]=C_reduced[k,i,j]}
      - M{D[i,k]=D[i,k]}
      - M{D_reduced[i,k]=D_reduced[i,k]}
      - M{m[i,k]=m[k,i]}
      - M{m_reduced[i,k]=m_reduced[k,i]}
      - M{d[i,k]=d[k,i]}
      - M{d_reduced[i,k]=d_reduced[k,i]}

   L{TransportPDE} also supports solution discontinuities over a contact region
   in the domain. To specify the conditions across the discontinuity we are
   using the generalised flux M{J} which, in the case of a system of PDEs and
   several components of the solution, is defined as

   M{J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]+X[i,j]+X_reduced[i,j]}

   For the case of single solution component and single PDE M{J} is defined as

   M{J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u+X[i]+X_reduced[i]}

   In the context of discontinuities M{n} denotes the normal on the
   discontinuity pointing from side 0 towards side 1 calculated from
   L{getNormal<escript.FunctionSpace.getNormal>} of L{FunctionOnContactZero<escript.FunctionOnContactZero>}.
   For a system of transport problems the contact condition takes the form

   M{n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]}

   where M{J0} and M{J1} are the fluxes on side 0 and side 1 of the
   discontinuity, respectively. M{jump(u)}, which is the difference of the
   solution at side 1 and at side 0, denotes the jump of M{u} across
   discontinuity along the normal calculated by L{jump<util.jump>}.
   The coefficient M{d_contact} is of rank two and M{y_contact} is of rank one
   both in L{FunctionOnContactZero<escript.FunctionOnContactZero>} or L{FunctionOnContactOne<escript.FunctionOnContactOne>}.
   The coefficient M{d_contact_reduced} is of rank two and M{y_contact_reduced}
   is of rank one both in L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>} or L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}.
   In case of a single PDE and a single component solution the contact
   condition takes the form

   M{n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)}

   In this case the coefficient M{d_contact} and M{y_contact} are each scalar
   both in L{FunctionOnContactZero<escript.FunctionOnContactZero>} or
   L{FunctionOnContactOne<escript.FunctionOnContactOne>} and the coefficient
   M{d_contact_reduced} and M{y_contact_reduced} are each scalar both in
   L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>} or
   L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>}.

   Typical usage::

       p = TransportPDE(dom)
       p.setValue(M=1., C=[-1.,0.])
       p.setInitialSolution(u=exp(-length(dom.getX()-[0.1,0.1])**2)
       t = 0
       dt = 0.1
       while (t < 1.):
           u = p.solve(dt)

   """
   def __init__(self,domain,numEquations=None,numSolutions=None, useBackwardEuler=False, debug=False):
     """
     Initializes a transport problem.

     @param domain: domain of the PDE
     @type domain: L{Domain<escript.Domain>}
     @param numEquations: number of equations. If C{None} the number of
                          equations is extracted from the coefficients.
     @param numSolutions: number of solution components. If C{None} the number
                          of solution components is extracted from the
                          coefficients.
     @param debug: if True debug information is printed
     @param useBackwardEuler: if set the backward Euler scheme is used. Otherwise the Crank-Nicholson scheme is applied. Note that backward Euler scheme will return a safe time step size which is practically infinity as the scheme is unconditional unstable. The Crank-Nicholson scheme provides a higher accuracy but requires to limit the time step size to be stable.
     @type useBackwardEuler: C{bool}
     """
     if useBackwardEuler:
         self.__useBackwardEuler=True
     else:
         self.__useBackwardEuler=False
     super(TransportPDE, self).__init__(domain,numEquations,numSolutions,debug)

     self.setConstraintWeightingFactor()
     #
     #   the coefficients of the transport problem
     #
     self.introduceCoefficients(
       M=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       m=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       M_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       A_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       B_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       C_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
       D_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       X_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
       Y_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       m_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       d_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       d_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),
       q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH) )

   def __str__(self):
     """
     Returns the string representation of the transport problem.

     @return: a simple representation of the transport problem
     @rtype: C{str}
     """
     return "<TransportPDE %d>"%id(self)

   def useBackwardEuler(self):
      """
      Returns true if backward Euler is used. Otherwise false is returned.
      @rtype: bool
      """
      return self.__useBackwardEuler


   def checkSymmetry(self,verbose=True):
      """
      Tests the transport problem for symmetry.

      @param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      @type verbose: C{bool}
      @return:  True if the PDE is symmetric
      @rtype: C{bool}
      @note: This is a very expensive operation. It should be used for
             degugging only! The symmetry flag is not altered.
      """
      out=True
      out=out and self.checkSymmetricTensor("M", verbose)
      out=out and self.checkSymmetricTensor("M_reduced", verbose)
      out=out and self.checkSymmetricTensor("A", verbose)
      out=out and self.checkSymmetricTensor("A_reduced", verbose)
      out=out and self.checkReciprocalSymmetry("B","C", verbose)
      out=out and self.checkReciprocalSymmetry("B_reduced","C_reduced", verbose)
      out=out and self.checkSymmetricTensor("D", verbose)
      out=out and self.checkSymmetricTensor("D_reduced", verbose)
      out=out and self.checkSymmetricTensor("m", verbose)
      out=out and self.checkSymmetricTensor("m_reduced", verbose)
      out=out and self.checkSymmetricTensor("d", verbose)
      out=out and self.checkSymmetricTensor("d_reduced", verbose)
      out=out and self.checkSymmetricTensor("d_contact", verbose)
      out=out and self.checkSymmetricTensor("d_contact_reduced", verbose)
      return out

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      @param coefficients: new values assigned to coefficients
      @keyword M: value for coefficient C{M}
      @type M: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword M_reduced: value for coefficient C{M_reduced}
      @type M_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{Function<escript.ReducedFunction>}
      @keyword A: value for coefficient C{A}
      @type A: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword A_reduced: value for coefficient C{A_reduced}
      @type A_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword B: value for coefficient C{B}
      @type B: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword B_reduced: value for coefficient C{B_reduced}
      @type B_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword C: value for coefficient C{C}
      @type C: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword C_reduced: value for coefficient C{C_reduced}
      @type C_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword D: value for coefficient C{D}
      @type D: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword D_reduced: value for coefficient C{D_reduced}
      @type D_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword X: value for coefficient C{X}
      @type X: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword X_reduced: value for coefficient C{X_reduced}
      @type X_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.ReducedFunction>}
      @keyword Y: value for coefficient C{Y}
      @type Y: any type that can be cast to a L{Data<escript.Data>} object on
               L{Function<escript.Function>}
      @keyword Y_reduced: value for coefficient C{Y_reduced}
      @type Y_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunction<escript.Function>}
      @keyword m: value for coefficient C{m}
      @type m: any type that can be cast to a L{Data<escript.Data>} object on
               L{FunctionOnBoundary<escript.FunctionOnBoundary>}
      @keyword m_reduced: value for coefficient C{m_reduced}
      @type m_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{FunctionOnBoundary<escript.ReducedFunctionOnBoundary>}
      @keyword d: value for coefficient C{d}
      @type d: any type that can be cast to a L{Data<escript.Data>} object on
               L{FunctionOnBoundary<escript.FunctionOnBoundary>}
      @keyword d_reduced: value for coefficient C{d_reduced}
      @type d_reduced: any type that can be cast to a L{Data<escript.Data>}
                       object on L{ReducedFunctionOnBoundary<escript.ReducedFunctionOnBoundary>}
      @keyword y: value for coefficient C{y}
      @type y: any type that can be cast to a L{Data<escript.Data>} object on
               L{FunctionOnBoundary<escript.FunctionOnBoundary>}
      @keyword d_contact: value for coefficient C{d_contact}
      @type d_contact: any type that can be cast to a L{Data<escript.Data>}
                       object on L{FunctionOnContactOne<escript.FunctionOnContactOne>} or L{FunctionOnContactZero<escript.FunctionOnContactZero>}
      @keyword d_contact_reduced: value for coefficient C{d_contact_reduced}
      @type d_contact_reduced: any type that can be cast to a L{Data<escript.Data>} object on L{ReducedFunctionOnContactOne<escript.ReducedFunctionOnContactOne>} or L{ReducedFunctionOnContactZero<escript.ReducedFunctionOnContactZero>}
      @keyword y_contact: value for coefficient C{y_contact}
      @type y_contact: any type that can be cast to a L{Data<escript.Data>}
                       object on L{FunctionOnContactOne<escript.FunctionOnContactOne>} or L{FunctionOnContactZero<escript.FunctionOnContactZero>}
      @keyword y_contact_reduced: value for coefficient C{y_contact_reduced}
      @type y_contact_reduced: any type that can be cast to a L{Data<escript.Data>} object on L{ReducedFunctionOnContactOne<escript.FunctionOnContactOne>} or L{ReducedFunctionOnContactZero<escript.FunctionOnContactZero>}
      @keyword r: values prescribed to the solution at the locations of constraints
      @type r: any type that can be cast to a L{Data<escript.Data>} object on
               L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
               depending on whether reduced order is used for the solution
      @keyword q: mask for the location of constraints
      @type q: any type that can be cast to a L{Data<escript.Data>} object on
               L{Solution<escript.Solution>} or
               L{ReducedSolution<escript.ReducedSolution>} depending on whether
               reduced order is used for the representation of the equation
      @raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      super(TransportPDE,self).setValue(**coefficients)

   def createOperator(self):
       """
       Returns an instance of a new transport operator.
       """
       if self.useBackwardEuler():
         theta=1.
       else:
         theta=0.5

       self.trace("New Transport problem is allocated.")
       return self.getDomain().newTransportProblem( \
                               theta,
                               self.getNumEquations(), \
                               self.getFunctionSpaceForSolution(), \
                               self.getSystemType())

   def setInitialSolution(self,u):
       """
       Sets the initial solution.

       @param u: new initial solution
       @type u: any object that can be interpolated to a L{Data<escript.Data>}
                object on L{Solution<escript.Solution>} or L{ReducedSolution<escript.ReducedSolution>}
       @note: C{u} must be non-negative
       """
       u2=util.interpolate(u,self.getFunctionSpaceForSolution())
       if self.getNumSolutions() == 1:
          if u2.getShape()!=():
              raise ValueError,"Illegal shape %s of initial solution."%(u2.getShape(),)
       else:
          if u2.getShape()!=(self.getNumSolutions(),):
              raise ValueError,"Illegal shape %s of initial solution."%(u2.getShape(),)
       self.getOperator().setInitialValue(u2)

   def getRequiredSystemType(self):
      """
      Returns the system type which needs to be used by the current set up.

      @return: a code to indicate the type of transport problem scheme used
      @rtype: C{float}
      """
      return self.getDomain().getTransportTypeId(self.getSolverMethod()[0],self.getSolverMethod()[1],self.getSolverPackage(),self.isSymmetric())

   def getUnlimitedTimeStepSize(self):
      """
      Returns the value returned by the C{getSafeTimeStepSize} method to
      indicate no limit on the safe time step size.

       @return: the value used to indicate that no limit is set to the time
                step size
       @rtype: C{float}
       @note: Typically the biggest positive float is returned
      """
      return self.getOperator().getUnlimitedTimeStepSize()

   def getSafeTimeStepSize(self):
       """
       Returns a safe time step size to do the next time step.

       @return: safe time step size
       @rtype: C{float}
       @note: If not C{getSafeTimeStepSize()} < C{getUnlimitedTimeStepSize()}
              any time step size can be used.
       """
       return self.getOperator().getSafeTimeStepSize()

   def setConstraintWeightingFactor(self,value=1./util.sqrt(util.EPSILON)):
       """
       Sets the weighting factor used to insert the constraints into the problem

       @param value: value for the weighting factor
       @type value: large positive C{float}
       """
       if not value>0:
         raise ValueError,"weighting factor needs to be positive."
       self.__constraint_factor=value
       self.trace("Weighting factor for constraints is set to %e."%value)

   def getConstraintWeightingFactor(self):
       """
       returns the weighting factor used to insert the constraints into the problem
       @return: value for the weighting factor
       @rtype: C{float}
       """
       return self.__constraint_factor
   #====================================================================
   def getSolution(self,dt,**options):
       """
       Returns the solution of the problem.

       @return: the solution
       @rtype: L{Data<escript.Data>}

       """
       if dt<=0:
           raise ValueError,"step size needs to be positive."
       options[self.TOLERANCE_KEY]=self.getTolerance()
       self.setSolution(self.getOperator().solve(self.getRightHandSide(),dt,options))
       self.validSolution()
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       @return: the discrete version of the PDE
       @rtype: C{tuple} of L{Operator,<escript.Operator>} and
               L{Data<escript.Data>}

       """
       if not self.isOperatorValid() or not self.isRightHandSideValid():
          self.resetRightHandSide()
          righthandside=self.getCurrentRightHandSide()
          self.resetOperator()
          operator=self.getCurrentOperator()
          self.getDomain().addPDEToTransportProblem(
                            operator,
                            righthandside,
                            self.getCoefficient("M"),
                            self.getCoefficient("A"),
                            self.getCoefficient("B"),
                            self.getCoefficient("C"),
                            self.getCoefficient("D"),
                            self.getCoefficient("X"),
                            self.getCoefficient("Y"),
                            self.getCoefficient("d"),
                            self.getCoefficient("y"),
                            self.getCoefficient("d_contact"),
                            self.getCoefficient("y_contact"))
          self.getDomain().addPDEToTransportProblem(
                            operator,
                            righthandside,
                            self.getCoefficient("M_reduced"),
                            self.getCoefficient("A_reduced"),
                            self.getCoefficient("B_reduced"),
                            self.getCoefficient("C_reduced"),
                            self.getCoefficient("D_reduced"),
                            self.getCoefficient("X_reduced"),
                            self.getCoefficient("Y_reduced"),
                            self.getCoefficient("d_reduced"),
                            self.getCoefficient("y_reduced"),
                            self.getCoefficient("d_contact_reduced"),
                            self.getCoefficient("y_contact_reduced"))
          operator.insertConstraint(righthandside,self.getCoefficient("q"),self.getCoefficient("r"),self.getConstraintWeightingFactor())
          self.trace("New system has been built.")
          self.validOperator()
          self.validRightHandSide()
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

   def setDebug(self, flag):
     """
     Switches debug output on if C{flag} is True,
     otherwise it is switched off.

     @param flag: desired debug status
     @type flag: C{bool}
     """
     if flag:
         self.setDebugOn()
     else:
         self.setDebugOff()

   def setDebugOn(self):
     """
     Switches debug output on.
     """
     super(TransportPDE,self).setDebugOn()
     
   def setDebugOff(self):
     """
     Switches debug output off.
     """
     super(TransportPDE,self).setDebugOff()
     