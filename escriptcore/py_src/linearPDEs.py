# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################

"""
The module provides an interface to define and solve linear partial
differential equations (PDEs) and Transport problems within `escript`.
`linearPDEs` does not provide any solver capabilities in itself but hands the
PDE over to the PDE solver library defined through the `Domain`
of the PDE. The general interface is provided through the `LinearPDE` class.
`TransportProblem` provides an interface to initial value problems dominated
by their advective terms.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import *
from . import escriptcpp as escore
from . import util
import math
import numpy

__author__="Lutz Gross, l.gross@uq.edu.au"

SolverOptions = escore.SolverOptions
SolverBuddy = escore.SolverBuddy

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

    :cvar INTERIOR: indicator that coefficient is defined on the interior of
                    the PDE domain
    :cvar BOUNDARY: indicator that coefficient is defined on the boundary of
                    the PDE domain
    :cvar CONTACT: indicator that coefficient is defined on the contact region
                   within the PDE domain
    :cvar INTERIOR_REDUCED: indicator that coefficient is defined on the
                            interior of the PDE domain using a reduced
                            integration order
    :cvar BOUNDARY_REDUCED: indicator that coefficient is defined on the
                            boundary of the PDE domain using a reduced
                            integration order
    :cvar CONTACT_REDUCED: indicator that coefficient is defined on the contact
                           region within the PDE domain using a reduced
                           integration order
    :cvar SOLUTION: indicator that coefficient is defined through a solution of
                    the PDE
    :cvar REDUCED: indicator that coefficient is defined through a reduced
                   solution of the PDE
    :cvar DIRACDELTA: indicator that coefficient is defined as Dirac delta functions
    :cvar BY_EQUATION: indicator that the dimension of the coefficient shape is
                       defined by the number of PDE equations
    :cvar BY_SOLUTION: indicator that the dimension of the coefficient shape is
                       defined by the number of PDE solutions
    :cvar BY_DIM: indicator that the dimension of the coefficient shape is
                  defined by the spatial dimension
    :cvar OPERATOR: indicator that the coefficient alters the operator of
                    the PDE
    :cvar RIGHTHANDSIDE: indicator that the coefficient alters the right
                         hand side of the PDE
    :cvar BOTH: indicator that the coefficient alters the operator as well
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
    DIRACDELTA=16

    def __init__(self, where, pattern, altering, isComplex=False):
       """
       Initialises a PDE coefficient type.

       :param where: describes where the coefficient lives
       :type where: one of `INTERIOR`, `BOUNDARY`, `CONTACT`, `SOLUTION`,
                    `REDUCED`, `INTERIOR_REDUCED`, `BOUNDARY_REDUCED`,
                    `CONTACT_REDUCED`, 'DIRACDELTA'
       :param pattern: describes the shape of the coefficient and how the shape
                       is built for a given spatial dimension and numbers of
                       equations and solutions in then PDE. For instance,
                       (`BY_EQUATION`,`BY_SOLUTION`,`BY_DIM`) describes a
                       rank 3 coefficient which is instantiated as shape (3,2,2)
                       in case of three equations and two solution components
                       on a 2-dimensional domain. In the case of single equation
                       and a single solution component the shape components
                       marked by `BY_EQUATION` or `BY_SOLUTION` are dropped.
                       In this case the example would be read as (2,).
       :type pattern: ``tuple`` of `BY_EQUATION`, `BY_SOLUTION`, `BY_DIM`
       :param altering: indicates what part of the PDE is altered if the
                        coefficient is altered
       :type altering: one of `OPERATOR`, `RIGHTHANDSIDE`, `BOTH`
       :param isComplex: if true, this coefficient is part of a complex-valued
                       PDE and values will be converted to complex.
       :type isComplex: ``boolean``
       """
       super(PDECoef, self).__init__()
       self.what = where
       self.pattern = pattern
       self.altering = altering
       self.__complex = isComplex
       self.resetValue()

    def isComplex(self):
        """
        Returns true if the coefficient is complex.

        :return: ``True`` if the coefficient is complex-valued
        :rtype: ``bool``
        """
        return self.__complex

    def resetValue(self):
       """
       Resets the coefficient value to the default.
       """
       self.value=escore.Data()

    def getFunctionSpace(self,domain,reducedEquationOrder=False,reducedSolutionOrder=False):
       """
       Returns the `FunctionSpace` of the coefficient.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param reducedEquationOrder: True to indicate that reduced order is used
                                    to represent the equation
       :type reducedEquationOrder: ``bool``
       :param reducedSolutionOrder: True to indicate that reduced order is used
                                    to represent the solution
       :type reducedSolutionOrder: ``bool``
       :return: `FunctionSpace` of the coefficient
       :rtype: `FunctionSpace`
       """
       if self.what==self.INTERIOR:
            return escore.Function(domain)
       elif self.what==self.INTERIOR_REDUCED:
            return escore.ReducedFunction(domain)
       elif self.what==self.BOUNDARY:
            return escore.FunctionOnBoundary(domain)
       elif self.what==self.BOUNDARY_REDUCED:
            return escore.ReducedFunctionOnBoundary(domain)
       elif self.what==self.CONTACT:
            return escore.FunctionOnContactZero(domain)
       elif self.what==self.CONTACT_REDUCED:
            return escore.ReducedFunctionOnContactZero(domain)
       elif self.what==self.DIRACDELTA:
            return escore.DiracDeltaFunctions(domain)
       elif self.what==self.SOLUTION:
            if reducedEquationOrder and reducedSolutionOrder:
                return escore.ReducedSolution(domain)
            else:
                return escore.Solution(domain)
       elif self.what==self.REDUCED:
            return escore.ReducedSolution(domain)

    def getValue(self):
       """
       Returns the value of the coefficient.

       :return: value of the coefficient
       :rtype: `Data`
       """
       return self.value

    def setValue(self,domain,numEquations=1,numSolutions=1,reducedEquationOrder=False,reducedSolutionOrder=False,newValue=None):
       """
       Sets the value of the coefficient to a new value.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param numEquations: number of equations of the PDE
       :type numEquations: ``int``
       :param numSolutions: number of components of the PDE solution
       :type numSolutions: ``int``
       :param reducedEquationOrder: True to indicate that reduced order is used
                                    to represent the equation
       :type reducedEquationOrder: ``bool``
       :param reducedSolutionOrder: True to indicate that reduced order is used
                                    to represent the solution
       :type reducedSolutionOrder: ``bool``
       :param newValue: new value of coefficient
       :type newValue: any object that can be converted into a
                       `Data` object with the appropriate shape
                       and `FunctionSpace`
       :raise IllegalCoefficientValue: if the shape of the assigned value does
                                       not match the shape of the coefficient
       :raise IllegalCoefficientFunctionSpace: if unable to interpolate value
                                               to appropriate function space
       """
       if newValue is None:
           newValue=escore.Data()
       elif isinstance(newValue, escore.Data):
           if not newValue.isEmpty():
              if not newValue.getFunctionSpace() == self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder):
                try:
                    newValue=escore.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
                except RuntimeError as er:
                 msg="Attempting to interpolate coefficient to function space %s encountered the following error: %s"%(self.getFunctionSpace(domain),str(er))
                 raise IllegalCoefficientFunctionSpace(msg)
                except:
                  raise IllegalCoefficientFunctionSpace("Unable to interpolate coefficient to function space %s"%self.getFunctionSpace(domain))
       else:
           try:
                newValue=escore.Data(newValue,self.getFunctionSpace(domain,reducedEquationOrder,reducedSolutionOrder))
           except:
                if not (isinstance(newValue, Data) or \
                isinstance(newValue, float) or \
                isinstance(newValue, complex) or \
                isinstance(newValue, numpy.ndarray)):
                    raise IllegalCoefficient("Illegal coefficient type: Type should be castable to Data.")
                else:
                    raise RuntimeError("Check input data")
       if not newValue.isEmpty():
           if not self.getShape(domain,numEquations,numSolutions)==newValue.getShape():
               raise IllegalCoefficientValue("Expected shape of coefficient is %s but actual shape is %s."%(self.getShape(domain,numEquations,numSolutions),newValue.getShape()))
           if newValue.isComplex() and not self.isComplex():
               raise IllegalCoefficientValue("Cannot assign a complex value to a real-valued coefficient!")
           elif not newValue.isComplex() and self.isComplex():
               newValue.promote()
       self.value = newValue

    def isAlteringOperator(self):
        """
        Checks if the coefficient alters the operator of the PDE.

        :return: True if the operator of the PDE is changed when the
                 coefficient is changed
        :rtype: ``bool``
        """
        if self.altering==self.OPERATOR or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isAlteringRightHandSide(self):
        """
        Checks if the coefficient alters the right hand side of the PDE.

        :rtype: ``bool``
        :return: True if the right hand side of the PDE is changed when the
                 coefficient is changed, ``None`` otherwise.
        """
        if self.altering==self.RIGHTHANDSIDE or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isComplex(self):
        """
        Checks if the coefficient is complex-valued.

        :rtype: ``bool``
        :return: True if the coefficient is complex-valued, False otherwise.
        """
        return self.__complex

    def estimateNumEquationsAndNumSolutions(self,domain,shape=()):
       """
       Tries to estimate the number of equations and number of solutions if
       the coefficient has the given shape.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param shape: suggested shape of the coefficient
       :type shape: ``tuple`` of ``int`` values
       :return: the number of equations and number of solutions of the PDE if
                the coefficient has given shape. If no appropriate numbers
                could be identified, ``None`` is returned
       :rtype: ``tuple`` of two ``int`` values or ``None``
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
          search.sort(key=lambda x: -(x[0]+x[1]))
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

       :return: True if the coefficient allows an estimate of the number of
                solution components, False otherwise
       :rtype: ``bool``
       """
       for i in self.pattern:
             if i==self.BY_SOLUTION: return True
       return False

    def definesNumEquation(self):
       """
       Checks if the coefficient allows to estimate the number of equations.

       :return: True if the coefficient allows an estimate of the number of
                equations, False otherwise
       :rtype: ``bool``
       """
       for i in self.pattern:
             if i==self.BY_EQUATION: return True
       return False

    def __CompTuple2(self,t1,t2):
      """
      Compares two tuples of possible number of equations and number of
      solutions.

      :param t1: the first tuple
      :param t2: the second tuple
      :return: 0, 1, or -1
      """

      dif=t1[0]+t1[1]-(t2[0]+t2[1])
      if dif<0: return 1
      elif dif>0: return -1
      else: return 0

    def getShape(self,domain,numEquations=1,numSolutions=1):
       """
       Builds the required shape of the coefficient.

       :param domain: domain on which the PDE uses the coefficient
       :type domain: `Domain`
       :param numEquations: number of equations of the PDE
       :type numEquations: ``int``
       :param numSolutions: number of components of the PDE solution
       :type numSolutions: ``int``
       :return: shape of the coefficient
       :rtype: ``tuple`` of ``int`` values
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

#====================================================================================================================

class LinearProblem(object):
   """
   This is the base class to define a general linear PDE-type problem for
   for an unknown function *u* on a given domain defined through a
   `Domain` object. The problem can be given as a single
   equation or as a system of equations.

   The class assumes that some sort of assembling process is required to form
   a problem of the form

   *L u=f*

   where *L* is an operator and *f* is the right hand side. This operator
   problem will be solved to get the unknown *u*.

   """
   def __init__(self,domain,numEquations=None,numSolutions=None,isComplex=False,debug=False):
     """
     Initializes a linear problem.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the coefficients.
     :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the
                          coefficients.
     :param isComplex: if True this problem will have complex coefficients and
                     a complex-valued result.
     :param debug: if True debug information is printed

     """
     super(LinearProblem, self).__init__()

     self.__complex=isComplex
     self.__debug=debug
     self.__domain=domain
     self.domainSupportsAssemblers = hasattr(domain, "createAssembler")
     self.assembler = None
     if self.domainSupportsAssemblers:
        options=[]
        if isComplex:
            options=[('dummy', escore.Data(0.j))]
        self.assembler = domain.createAssembler("DefaultAssembler", options)
     self.__numEquations=numEquations
     self.__numSolutions=numSolutions
     self.__preservePreconditioner=False
     self.__altered_coefficients=False
     self.__reduce_equation_order=False
     self.__reduce_solution_order=False
     self.__sym=False
     self.__herm=False
     self.__is_RHS_valid=False
     self.__is_operator_valid=False
     self.__COEFFICIENTS={}
     self.__solution_rtol=1.e99
     self.__solution_atol=1.e99
     # Record if we are using oxley
     self.__have_oxley=False
     if domain.getDescription() == 'oxley::rectangle' or domain.getDescription() == 'oxley::brick':
        self.__have_oxley=True
     self.setSolverOptions()
     self.setSymmetryOff()
     # Set on lumping if we are using Speckley
     try:
          if domain.getDescription() == 'speckley::Rectangle' or domain.getDescription() == 'speckley::Brick':
               self.getSolverOptions().setSolverMethod(SolverOptions.LUMPING)
     except:
          pass
     # set number of equations in trilinos
     self.getSolverOptions().setTrilinosParameter("number of equations", numEquations)
     # initialize things:
     self.resetAllCoefficients()
     self.initializeSystem()

   # ==========================================================================
   #    general stuff:
   # ==========================================================================
   def __str__(self):
     """
     Returns a string representation of the PDE.

     :return: a simple representation of the PDE
     :rtype: ``str``
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
     Switches debug output on if ``flag`` is True otherwise it is switched off.

     :param flag: desired debug status
     :type flag: ``bool``
     """
     if flag:
         self.setDebugOn()
     else:
         self.setDebugOff()

   def trace(self,text):
     """
     Prints the text message if debug mode is switched on.

     :param text: message to be printed
     :type text: ``string``
     """
     if self.__debug: print(("%s: %s"%(str(self),text)))

   # ==========================================================================
   # some service functions:
   # ==========================================================================
   def introduceCoefficients(self,**coeff):
       """
       Introduces new coefficients into the problem.

       Use:

       p.introduceCoefficients(A=PDECoef(...), B=PDECoef(...))

       to introduce the coefficients *A* and *B*.
       """
       for name, type in sorted(coeff.items(), key=lambda x: x[0]):
           if not isinstance(type,PDECoef):
              raise ValueError("coefficient %s has no type."%name)
           self.__COEFFICIENTS[name]=type
           self.__COEFFICIENTS[name].resetValue()
           self.trace("coefficient %s has been introduced."%name)

   def resetRightHandSideCoefficients(self):
       """
       Resets all coefficients defining the right hand side
       """
       for name in self.__COEFFICIENTS:
         if self.__COEFFICIENTS[name].altering == PDECoef.RIGHTHANDSIDE  :
              self.__COEFFICIENTS[name].resetValue()
              self.trace("coefficient %s has been reset."%name)

   def getDomain(self):
     """
     Returns the domain of the PDE.

     :return: the domain of the PDE
     :rtype: `Domain`
     """
     return self.__domain

   def getDomainStatus(self):
     """
     Return the status indicator of the domain
     """
     return self.getDomain().getStatus()

   def getSystemStatus(self):
     """
     Return the domain status used to build the current system
     """
     return self.__system_status

   def setSystemStatus(self,status=None):
     """
     Sets the system status to ``status`` if ``status`` is not present the
     current status of the domain is used.
     """
     if status is None:
         self.__system_status=self.getDomainStatus()
     else:
         self.__system_status=status

   def getDim(self):
     """
     Returns the spatial dimension of the PDE.

     :return: the spatial dimension of the PDE domain
     :rtype: ``int``
     """
     return self.getDomain().getDim()

   def getNumEquations(self):
     """
     Returns the number of equations.

     :return: the number of equations
     :rtype: ``int``
     :raise UndefinedPDEError: if the number of equations is not specified yet
     """
     if self.__numEquations is None:
         if self.__numSolutions is None:
            raise UndefinedPDEError("Number of equations is undefined. Please specify argument numEquations.")
         else:
            self.__numEquations=self.__numSolutions
            self.getSolverOptions().setTrilinosParameter("number of equations", self.__numEquations)
     return self.__numEquations

   def getNumSolutions(self):
     """
     Returns the number of unknowns.

     :return: the number of unknowns
     :rtype: ``int``
     :raise UndefinedPDEError: if the number of unknowns is not specified yet
     """
     if self.__numSolutions is None:
        if self.__numEquations is None:
            raise UndefinedPDEError("Number of solution is undefined. Please specify argument numSolutions.")
        else:
            self.__numSolutions=self.__numEquations
     return self.__numSolutions

   def reduceEquationOrder(self):
     """
     Returns the status of order reduction for the equation.

     :return: True if reduced interpolation order is used for the
              representation of the equation, False otherwise
     :rtype: `bool`
     """
     return self.__reduce_equation_order

   def reduceSolutionOrder(self):
     """
     Returns the status of order reduction for the solution.

     :return: True if reduced interpolation order is used for the
              representation of the solution, False otherwise
     :rtype: `bool`
     """
     return self.__reduce_solution_order

   def getFunctionSpaceForEquation(self):
     """
     Returns the `FunctionSpace` used to discretize
     the equation.

     :return: representation space of equation
     :rtype: `FunctionSpace`
     """
     if self.reduceEquationOrder():
         return escore.ReducedSolution(self.getDomain())
     else:
         return escore.Solution(self.getDomain())

   def getFunctionSpaceForSolution(self):
     """
     Returns the `FunctionSpace` used to represent
     the solution.

     :return: representation space of solution
     :rtype: `FunctionSpace`
     """
     if self.reduceSolutionOrder():
         return escore.ReducedSolution(self.getDomain())
     else:
         return escore.Solution(self.getDomain())

   # ==========================================================================
   #   solver settings:
   # ==========================================================================
   def setSolverOptions(self,options=None):
       """
       Sets the solver options.

       :param options: the new solver options. If equal ``None``, the solver options are set to the default.
       :type options: `SolverOptions` or ``None``
       :note: The symmetry flag of options is overwritten by the symmetry flag of the `LinearProblem`.
       """
       if options is None:
          self.__solver_options=SolverBuddy()
       elif isinstance(options, SolverBuddy):
          self.__solver_options=options
       else:
          raise ValueError("options must be a SolverOptions object.")
       self.__solver_options.setComplex(self.isComplex())
       self.__solver_options.setSymmetry(self.__sym)
       self.__solver_options.setHermitian(self.__herm)
       self.__solver_options.setDim(self.getDim())
       self.__solver_options.setOxleyDomain(self.hasOxley())

   def getSolverOptions(self):
       """
       Returns the solver options

       :rtype: `SolverOptions`
       """
       self.__solver_options.setSymmetry(self.__sym)
       return self.__solver_options

   def isUsingLumping(self):
      """
      Checks if matrix lumping is the current solver method.

      :return: True if the current solver method is lumping
      :rtype: ``bool``
      """
      return self.getSolverOptions().getSolverMethod() in [ SolverOptions.ROWSUM_LUMPING, SolverOptions.HRZ_LUMPING ]

   def isComplex(self):
       """
       Returns true if this is a complex-valued LinearProblem, false if
       real-valued.

       :rtype: ``bool``
       """
       return self.__complex

   def shouldPreservePreconditioner(self):
       """
       Returns true if the preconditioner / factorisation should be kept even
       when resetting the operator.

       :rtype: ``bool``
       """
       return self.__preservePreconditioner

   def preservePreconditioner(self, preserve = True):
        """
        Notifies the PDE that the preconditioner should not be reset when
        making changes to the operator.

        Building the preconditioner data can be quite expensive (e.g. for
        multigrid methods) so if it is known that changes to the operator are
        going to be minor calling this method can speed up successive PDE
        solves.

        :note: Not all operator types support this.
        :param preserve: if True, preconditioner will be preserved, otherwise
                            it will be reset when making changes to the operator,
                            which is the default behaviour.
        :type preserve: ``bool``
        """
        self.__preservePreconditioner =  preserve

   def hasOxley(self):
        """
        return True if an ``esys.escript.oxley`` domain is used
        """
        return self.__have_oxley
   # ==========================================================================
   #    symmetry  flag:
   # ==========================================================================
   def isSymmetric(self):
      """
      Checks if symmetry is indicated.

      :return: True if a symmetric PDE is indicated, False otherwise
      :rtype: ``bool``
      :note: the method is equivalent to use getSolverOptions().isSymmetric()
      """
      self.getSolverOptions().isSymmetric()

   def setSymmetryOn(self):
      """
      Sets the symmetry flag.
      :note: The method overwrites the symmetry flag set by the solver options
      """
      self.__sym=True
      self.getSolverOptions().setSymmetryOn()

   def setSymmetryOff(self):
      """
      Clears the symmetry flag.
      :note: The method overwrites the symmetry flag set by the solver options
      """
      self.__sym=False
      self.getSolverOptions().setSymmetryOff()

   def setSymmetry(self,flag=False):
      """
      Sets the symmetry flag to ``flag``.

      :param flag: If True, the symmetry flag is set otherwise reset.
      :type flag: ``bool``
      :note: The method overwrites the symmetry flag set by the solver options
      """
      self.getSolverOptions().setSymmetry(flag)

   # ==========================================================================
   #    hermitian  flag:
   # ==========================================================================
   def isHermitian(self):
      """
      Checks if the pde is indicated to be Hermitian.

      :return: True if a Hermitian PDE is indicated, False otherwise
      :rtype: ``bool``
      :note: the method is equivalent to use getSolverOptions().isHermitian()
      """
      self.getSolverOptions().isHermitian()

   def setHermitianOn(self):
      """
      Sets the Hermitian flag.
      :note: The method overwrites the Hermitian flag set by the solver options
      """
      self.__herm=True
      self.getSolverOptions().setHermitianOn()

   def setHermitianOff(self):
      """
      Clears the Hermitian flag.
      :note: The method overwrites the Hermitian flag set by the solver options
      """
      self.__herm=False
      self.getSolverOptions().setHermitianOff()

   def setHermitian(self,flag=False):
      """
      Sets the Hermitian flag to ``flag``.

      :param flag: If True, the Hermitian flag is set otherwise reset.
      :type flag: ``bool``
      :note: The method overwrites the Hermitian flag set by the solver options
      """
      self.getSolverOptions().setHermitian(flag)

   # ==========================================================================
   # function space handling for the equation as well as the solution
   # ==========================================================================
   def setReducedOrderOn(self):
     """
     Switches reduced order on for solution and equation representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     Switches reduced order off for solution and equation representation

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     Sets order reduction state for both solution and equation representation
     according to flag.

     :param flag: if True, the order reduction is switched on for both solution
                  and equation representation, otherwise or if flag is not
                  present order reduction is switched off
     :type flag: ``bool``
     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)

   def setReducedOrderForSolutionOn(self):
     """
     Switches reduced order on for solution representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if not self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Reduced order is used for solution representation.")
         self.__reduce_solution_order=True
         self.initializeSystem()

   def setReducedOrderForSolutionOff(self):
     """
     Switches reduced order off for solution representation

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set.
     """
     if self.__reduce_solution_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Full order is used to interpolate solution.")
         self.__reduce_solution_order=False
         self.initializeSystem()

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     Sets order reduction state for solution representation according to flag.

     :param flag: if flag is True, the order reduction is switched on for
                  solution representation, otherwise or if flag is not present
                  order reduction is switched off
     :type flag: ``bool``
     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()

   def setReducedOrderForEquationOn(self):
     """
     Switches reduced order on for equation representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if not self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Reduced order is used for test functions.")
         self.__reduce_equation_order=True
         self.initializeSystem()

   def setReducedOrderForEquationOff(self):
     """
     Switches reduced order off for equation representation.

     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if self.__reduce_equation_order:
         if self.__altered_coefficients:
              raise RuntimeError("order cannot be altered after coefficients have been defined.")
         self.trace("Full order is used for test functions.")
         self.__reduce_equation_order=False
         self.initializeSystem()

   def setReducedOrderForEquationTo(self,flag=False):
     """
     Sets order reduction state for equation representation according to flag.

     :param flag: if flag is True, the order reduction is switched on for
                  equation representation, otherwise or if flag is not present
                  order reduction is switched off
     :type flag: ``bool``
     :raise RuntimeError: if order reduction is altered after a coefficient has
                          been set
     """
     if flag:
        self.setReducedOrderForEquationOn()
     else:
        self.setReducedOrderForEquationOff()

   def getOperatorType(self):
      """
      Returns the current system type.
      """
      return self.__operator_type

   def checkSymmetricTensor(self,name,verbose=True):
      """
      Tests a coefficient for symmetry.

      :param name: name of the coefficient
      :type name: ``str``
      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      :type verbose: ``bool``
      :return: True if coefficient ``name`` is symmetric
      :rtype: ``bool``
      """
      SMALL_TOLERANCE=util.EPSILON*10.
      A=self.getCoefficient(name)
      verbose=verbose or self.__debug
      out=True
      if not A.isEmpty():
         tol=util.Lsup(A)*SMALL_TOLERANCE
         s=A.getShape()
         if A.getRank() == 4:
            if s[0]==s[2] and s[1] == s[3]:
               for i in range(s[0]):
                  for j in range(s[1]):
                     for k in range(s[2]):
                        for l in range(s[3]):
                            if util.Lsup(A[i,j,k,l]-A[k,l,i,j])>tol:
                               if verbose: print(("non-symmetric problem as %s[%d,%d,%d,%d]!=%s[%d,%d,%d,%d]"%(name,i,j,k,l,name,k,l,i,j)))
                               out=False
            else:
                 if verbose: print(("non-symmetric problem because of inappropriate shape %s of coefficient %s."%(s,name)))
                 out=False
         elif A.getRank() == 2:
            if s[0]==s[1]:
               for j in range(s[0]):
                  for l in range(s[1]):
                     if util.Lsup(A[j,l]-A[l,j])>tol:
                        if verbose: print(("non-symmetric problem because %s[%d,%d]!=%s[%d,%d]"%(name,j,l,name,l,j)))
                        out=False
            else:
                 if verbose: print(("non-symmetric problem because of inappropriate shape %s of coefficient %s."%(s,name)))
                 out=False
         elif A.getRank() == 0:
            pass
         else:
             raise ValueError("Cannot check rank %s of %s."%(A.getRank(),name))
      return out

   def checkReciprocalSymmetry(self,name0,name1,verbose=True):
      """
      Tests two coefficients for reciprocal symmetry.

      :param name0: name of the first coefficient
      :type name0: ``str``
      :param name1: name of the second coefficient
      :type name1: ``str``
      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed
      :type verbose: ``bool``
      :return: True if coefficients ``name0`` and ``name1`` are reciprocally
               symmetric.
      :rtype: ``bool``
      """
      SMALL_TOLERANCE=util.EPSILON*10.
      B=self.getCoefficient(name0)
      C=self.getCoefficient(name1)
      verbose=verbose or self.__debug
      out=True
      if B.isEmpty() and not C.isEmpty():
         if verbose: print(("non-symmetric problem because %s is not present but %s is"%(name0,name1)))
         out=False
      elif not B.isEmpty() and C.isEmpty():
         if verbose: print(("non-symmetric problem because %s is not present but %s is"%(name0,name1)))
         out=False
      elif not B.isEmpty() and not C.isEmpty():
         sB=B.getShape()
         sC=C.getShape()
         tol=(util.Lsup(B)+util.Lsup(C))*SMALL_TOLERANCE/2.
         if len(sB) != len(sC):
             if verbose: print(("non-symmetric problem because ranks of %s (=%s) and %s (=%s) are different."%(name0,len(sB),name1,len(sC))))
             out=False
         else:
             if len(sB)==0:
               if util.Lsup(B-C)>tol:
                  if verbose: print(("non-symmetric problem because %s!=%s"%(name0,name1)))
                  out=False
             elif len(sB)==1:
               if sB[0]==sC[0]:
                  for j in range(sB[0]):
                     if util.Lsup(B[j]-C[j])>tol:
                        if verbose: print(("non-symmetric PDE because %s[%d]!=%s[%d]"%(name0,j,name1,j)))
                        out=False
               else:
                 if verbose: print(("non-symmetric problem because of inappropriate shapes %s and %s of coefficients %s and %s, respectively."%(sB,sC,name0,name1)))
             elif len(sB)==3:
               if sB[0]==sC[1] and sB[1]==sC[2] and sB[2]==sC[0]:
                   for i in range(sB[0]):
                      for j in range(sB[1]):
                         for k in range(sB[2]):
                            if util.Lsup(B[i,j,k]-C[k,i,j])>tol:
                                 if verbose: print(("non-symmetric problem because %s[%d,%d,%d]!=%s[%d,%d,%d]"%(name0,i,j,k,name1,k,i,j)))
                                 out=False
               else:
                 if verbose: print(("non-symmetric problem because of inappropriate shapes %s and %s of coefficients %s and %s, respectively."%(sB,sC,name0,name1)))
             else:
                 raise ValueError("Cannot check rank %s of %s and %s."%(len(sB),name0,name1))
      return out

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name``.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient
     :rtype: `Data`
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
         return self.__COEFFICIENTS[name].getValue()
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def hasCoefficient(self,name):
     """
     Returns True if ``name`` is the name of a coefficient.

     :param name: name of the coefficient enquired
     :type name: ``string``
     :return: True if ``name`` is the name of a coefficient of the general PDE,
              False otherwise
     :rtype: ``bool``
     """
     return name in self.__COEFFICIENTS

   def createCoefficient(self, name):
     """
     Creates a `Data` object corresponding to coefficient
     ``name``.

     :return: the coefficient ``name`` initialized to 0
     :rtype: `Data`
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        zero = 0.j if self.__COEFFICIENTS[name].isComplex() else 0.
        return escore.Data(zero,self.getShapeOfCoefficient(name),self.getFunctionSpaceForCoefficient(name))
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def getFunctionSpaceForCoefficient(self,name):
     """
     Returns the `FunctionSpace` to be used for
     coefficient ``name``.

     :param name: name of the coefficient enquired
     :type name: ``string``
     :return: the function space to be used for coefficient ``name``
     :rtype: `FunctionSpace`
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return self.__COEFFICIENTS[name].getFunctionSpace(self.getDomain())
     else:
        raise ValueError("unknown coefficient %s requested"%name)

   def getShapeOfCoefficient(self,name):
     """
     Returns the shape of the coefficient ``name``.

     :param name: name of the coefficient enquired
     :type name: ``string``
     :return: the shape of the coefficient ``name``
     :rtype: ``tuple`` of ``int``
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     """
     if self.hasCoefficient(name):
        return self.__COEFFICIENTS[name].getShape(self.getDomain(),self.getNumEquations(),self.getNumSolutions())
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

   def resetAllCoefficients(self):
     """
     Resets all coefficients to their default values.
     """
     for i in sorted(self.__COEFFICIENTS.keys()):
         self.__COEFFICIENTS[i].resetValue()

   def alteredCoefficient(self,name):
     """
     Announces that coefficient ``name`` has been changed.

     :param name: name of the coefficient affected
     :type name: ``string``
     :raise IllegalCoefficient: if ``name`` is not a coefficient of the PDE
     :note: if ``name`` is q or r, the method will not trigger a rebuild of the
            system as constraints are applied to the solved system.
     """
     if self.hasCoefficient(name):
        self.trace("Coefficient %s has been altered."%name)
        if not ((name=="q" or name=="r") and self.isUsingLumping()):
           if self.__COEFFICIENTS[name].isAlteringOperator(): self.invalidateOperator()
           if self.__COEFFICIENTS[name].isAlteringRightHandSide(): self.invalidateRightHandSide()
     else:
        raise IllegalCoefficient("illegal coefficient %s requested for general PDE."%name)

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
       if not self.getDomainStatus()==self.getSystemStatus(): self.invalidateSolution()
       if self.__solution_rtol>self.getSolverOptions().getTolerance() or \
          self.__solution_atol>self.getSolverOptions().getAbsoluteTolerance():
            self.invalidateSolution()
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
       if not self.getDomainStatus()==self.getSystemStatus(): self.invalidateOperator()
       if not self.getRequiredOperatorType()==self.getOperatorType(): self.invalidateOperator()
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
       self.trace("Right hand side has to be rebuilt.")
       self.invalidateSolution()
       self.__is_RHS_valid=False

   def isRightHandSideValid(self):
       """
       Returns True if the operator is still valid.
       """
       if not self.getDomainStatus()==self.getSystemStatus(): self.invalidateRightHandSide()
       return self.__is_RHS_valid

   def invalidateSystem(self):
       """
       Announces that everything has to be rebuilt.
       """
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
       self.__operator_type=None
       self.setSystemStatus()
       self.__operator=escore.Operator()
       self.__righthandside=escore.Data()
       self.__solution=escore.Data()
       self.invalidateSystem()

   def getOperator(self):
     """
     Returns the operator of the linear problem.

     :return: the operator of the problem
     """
     return self.getSystem()[0]

   def getRightHandSide(self):
     """
     Returns the right hand side of the linear problem.

     :return: the right hand side of the problem
     :rtype: `Data`
     """
     return self.getSystem()[1]

   def createRightHandSide(self):
       """
       Returns an instance of a new right hand side.
       """
       self.trace("New right hand side is allocated.")
       zero = 0.j if self.isComplex() else 0.
       if self.getNumEquations()>1:
           return escore.Data(zero,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),True)
       else:
           return escore.Data(zero,(),self.getFunctionSpaceForEquation(),True)

   def createSolution(self):
       """
       Returns an instance of a new solution.
       """
       self.trace("New solution is allocated.")
       zero = 0.j if self.isComplex() else 0.
       if self.getNumSolutions() > 1:
           return escore.Data(zero,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)
       else:
           return escore.Data(zero,(),self.getFunctionSpaceForSolution(),True)

   def resetSolution(self):
       """
       Sets the solution to zero.
       """
       if self.__solution.isEmpty():
           self.__solution=self.createSolution()
       else:
           self.__solution.setToZero()
           self.trace("Solution is reset to zero.")

   def setSolution(self,u, validate=True):
       """
       Sets the solution assuming that makes the system valid with the tolrance
       defined by the solver options
       """
       if validate:
          self.__solution_rtol=self.getSolverOptions().getTolerance()
          self.__solution_atol=self.getSolverOptions().getAbsoluteTolerance()
          self.validSolution()
       self.__solution=u

   def getCurrentSolution(self):
       """
       Returns the solution in its current state.
       """
       if self.__solution.isEmpty(): 
          self.__solution=self.createSolution()
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
       if self.getOperatorType() is None:
           if self.isUsingLumping():
               self.__operator=self.createSolution()
           else:
               self.__operator=self.createOperator()
           self.__operator_type=self.getRequiredOperatorType()
       else:
           if self.isUsingLumping():
               self.__operator.setToZero()
           else:
               if self.getOperatorType() == self.getRequiredOperatorType():
                   self.__operator.resetValues(self.shouldPreservePreconditioner())
               else:
                   self.__operator=self.createOperator()
                   self.__operator_type=self.getRequiredOperatorType()
           self.trace("Operator reset to zero")

   def getCurrentOperator(self):
       """
       Returns the operator in its current state.
       """
       return self.__operator

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      :raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      # check if the coefficients are  legal:
      for i in sorted(coefficients.keys()):
         if not self.hasCoefficient(i):
            raise IllegalCoefficient("Attempt to set unknown coefficient %s"%i)
      # if the number of unknowns or equations is still unknown we try to estimate them:
      if self.__numEquations is None or self.__numSolutions is None:
         for i,d in sorted(coefficients.items(), key=lambda x: x[0]):
            if hasattr(d,"shape"):
                s=d.shape
            elif isinstance(d, escore.Data) and not d.isEmpty():
                s=d.getShape()
            else:
                s=numpy.array(d).shape
            if s!=None:
                # get number of equations and number of unknowns:
                res=self.__COEFFICIENTS[i].estimateNumEquationsAndNumSolutions(self.getDomain(),s)
                if res is None:
                    raise IllegalCoefficientValue("Illegal shape %s of coefficient %s"%(s,i))
                else:
                    if self.__numEquations is None: self.__numEquations=res[0]
                    if self.__numSolutions is None: self.__numSolutions=res[1]
      if self.__numEquations is None: raise UndefinedPDEError("unidentified number of equations")
      if self.__numSolutions is None: raise UndefinedPDEError("unidentified number of solutions")
      # now we check the shape of the coefficient if numEquations and numSolutions are set:
      for i,d in sorted(coefficients.items(), key=lambda x: x[0]):
        try:
           self.__COEFFICIENTS[i].setValue(self.getDomain(),
                     self.getNumEquations(),self.getNumSolutions(),
                     self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
           self.alteredCoefficient(i)
        except IllegalCoefficientFunctionSpace as m:
            # if the function space is wrong then we try the reduced version:
            i_red=i+"_reduced"
            if (not i_red in list(coefficients.keys())) and i_red in list(self.__COEFFICIENTS.keys()):
                try:
                    self.__COEFFICIENTS[i_red].setValue(self.getDomain(),
                                                      self.getNumEquations(),self.getNumSolutions(),
                                                      self.reduceEquationOrder(),self.reduceSolutionOrder(),d)
                    self.alteredCoefficient(i_red)
                except IllegalCoefficientValue as m:
                    raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
                except IllegalCoefficientFunctionSpace as m:
                    raise IllegalCoefficientFunctionSpace("Coefficient %s:%s"%(i,m))
            else:
                raise IllegalCoefficientFunctionSpace("Coefficient %s:%s"%(i,m))
        except IllegalCoefficientValue as m:
           raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
      self.__altered_coefficients=True

   # ==========================================================================
   # methods that are typically overwritten when implementing a particular
   # linear problem
   # ==========================================================================
   def getRequiredOperatorType(self):
      """
      Returns the system type which needs to be used by the current set up.

      :note: Typically this method is overwritten when implementing a
             particular linear problem.
      """
      return None

   def createOperator(self):
       """
       Returns an instance of a new operator.

       :note: This method is overwritten when implementing a particular
              linear problem.
       """
       return escore.Operator()

   def checkSymmetry(self,verbose=True):
      """
      Tests the PDE for symmetry.

      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed
      :type verbose: ``bool``
      :return: True if the problem is symmetric
      :rtype: ``bool``
      :note: Typically this method is overwritten when implementing a
             particular linear problem.
      """
      out=True
      return out

   def getSolution(self,**options):
       """
       Returns the solution of the problem.

       :return: the solution
       :rtype: `Data`

       :note: This method is overwritten when implementing a particular
              linear problem.
       """
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       :return: the discrete version of the PDE
       :rtype: ``tuple`` of `Operator` and `Data`.

       :note: This method is overwritten when implementing a particular
              linear problem.
       """
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())


   def addPDEToSystem(self, operator,righthandside, A, B, C, D, X, Y,
            d, y, d_contact, y_contact, d_dirac, y_dirac):
        """
        adds a PDE to the system, results depend on domain

        :param mat:
        :type mat: `OperatorAdapter`
        :param rhs:
        :type rhs: `Data`
        :param A:
        :type A: `Data`
        :param B:
        :type B: `Data`
        :param C:
        :type C: `Data`
        :param D:
        :type D: `Data`
        :param X:
        :type X: `Data`
        :param Y:
        :type Y: `Data`
        :param d:
        :type d: `Data`
        :param y:
        :type y: `Data`
        :param d_contact:
        :type d_contact: `Data`
        :param y_contact:
        :type y_contact: `Data`
        :param d_dirac:
        :type d_dirac: `Data`
        :param y_dirac:
        :type y_dirac: `Data`
        """
        if self.domainSupportsAssemblers:
            data = [("A", A), ("B", B), ("C", C), ("D", D), ("X", X), ("Y", Y),
                    ("d", d), ("y", y), ("d_contact", d_contact),
                    ("y_contact", y_contact), ("d_dirac", d_dirac),
                    ("y_dirac", y_dirac)]
            self.addToSystem(operator,righthandside, data)
        else:
            self.getDomain().addPDEToSystem(operator,righthandside, A, B, C, D,
                    X, Y, d, y, d_contact, y_contact, d_dirac, y_dirac)

   def addToSystem(self, op, rhs, data):
        """
        adds a PDE to the system, results depend on domain

        :param mat:
        :type mat: `OperatorAdapter`
        :param rhs:
        :type rhs: `Data`
        :param data:
        :type data: `list`
        """
        self.getDomain().addToSystem(op, rhs, data, self.assembler)
        if self.hasOxley():
            self.getDomain().makeZ(self.__complex)
            self.getDomain().makeIZ(self.__complex)
            self.getDomain().finaliseA(op,self.__complex)
            rhs=self.getDomain().finaliseRhs(rhs)

   def addPDEToLumpedSystem(self, operator, a, b, c, hrz_lumping):
        """
        adds a PDE to the lumped system, results depend on domain

        :param mat:
        :type mat: `OperatorAdapter`
        :param rhs:
        :type rhs: `Data`
        :param a:
        :type a: `Data`
        :param b:
        :type b: `Data`
        :param c:
        :type c: `Data`
        :param hrz_lumping:
        :type hrz_lumping: `bool`
        """
        if self.domainSupportsAssemblers:
            self.getDomain().addPDEToLumpedSystem(operator, a, b, c, hrz_lumping, self.assembler)
        else:
            self.getDomain().addPDEToLumpedSystem(operator, a, b, c, hrz_lumping)

   def addPDEToRHS(self, righthandside, X, Y, y, y_contact, y_dirac):
        """
        adds a PDE to the right hand side, results depend on domain

        :param mat:
        :type mat: `OperatorAdapter`
        :param righthandside:
        :type righthandside: `Data`
        :param X:
        :type X: `Data`
        :param Y:
        :type Y: `Data`
        :param y:
        :type y: `Data`
        :param y_contact:
        :type y_contact: `Data`
        :param y_dirac:
        :type y_dirac: `Data`
        """
        if self.domainSupportsAssemblers:
            data = [("X", X), ("Y", Y), ("y", y), ("y_contact", y_contact),
                    ("y_dirac", y_dirac)]
            self.addToRHS(righthandside, data)
        else:
            self.getDomain().addPDEToRHS(righthandside, X, Y, y, y_contact,
                    y_dirac)

   def addToRHS(self, rhs, data):
        """
        adds a PDE to the right hand side, results depend on domain

        :param mat:
        :type mat: `OperatorAdapter`
        :param righthandside:
        :type righthandside: `Data`
        :param data:
        :type data: `list`
        """
        self.getDomain().addToRHS(rhs, data, self.assembler)

class LinearPDE(LinearProblem):
   """
   This class is used to define a general linear, steady, second order PDE
   for an unknown function *u* on a given domain defined through a
   `Domain` object.

   For a single PDE having a solution with a single component the linear PDE
   is defined in the following form:

   *-(grad(A[j,l]+A_reduced[j,l])*grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)=-grad(X+X_reduced)[j,j]+(Y+Y_reduced)*

   where *grad(F)* denotes the spatial derivative of *F*. Einstein's
   summation convention, ie. summation over indexes appearing twice in a term
   of a sum performed, is used.
   The coefficients *A*, *B*, *C*, *D*, *X* and *Y* have to be specified
   through `Data` objects in `Function` and
   the coefficients *A_reduced*, *B_reduced*, *C_reduced*, *D_reduced*,
   *X_reduced* and *Y_reduced* have to be specified through
   `Data` objects in `ReducedFunction`.
   It is also allowed to use objects that can be converted into such
   `Data` objects. *A* and *A_reduced* are rank two, *B*,
   *C*, *X*, *B_reduced*, *C_reduced* and *X_reduced* are rank one and
   *D*, *D_reduced*, *Y* and *Y_reduced* are scalar.

   The following natural boundary conditions are considered:

   *n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u)+(d+d_reduced)*u=n[j]*(X[j]+X_reduced[j])+y*

   where *n* is the outer normal field. Notice that the coefficients *A*,
   *A_reduced*, *B*, *B_reduced*, *X* and *X_reduced* are defined in the
   PDE. The coefficients *d* and *y* are each a scalar in
   `FunctionOnBoundary` and the coefficients
   *d_reduced* and *y_reduced* are each a scalar in
   `ReducedFunctionOnBoundary`.

   Constraints for the solution prescribe the value of the solution at certain
   locations in the domain. They have the form

   *u=r* where *q>0*

   *r* and *q* are each scalar where *q* is the characteristic function
   defining where the constraint is applied. The constraints override any
   other condition set by the PDE or the boundary condition.

   The PDE is symmetrical if

   *A[i,j]=A[j,i]*  and *B[j]=C[j]* and *A_reduced[i,j]=A_reduced[j,i]*
   and *B_reduced[j]=C_reduced[j]*

   For a system of PDEs and a solution with several components the PDE has the
   form

   *-grad((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])[j]+(C[i,k,l]+C_reduced[i,k,l])*grad(u[k])[l]+(D[i,k]+D_reduced[i,k]*u[k] =-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i]*

   *A* and *A_reduced* are of rank four, *B*, *B_reduced*, *C* and
   *C_reduced* are each of rank three, *D*, *D_reduced*, *X_reduced* and
   *X* are each of rank two and *Y* and *Y_reduced* are of rank one.
   The natural boundary conditions take the form:

   *n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k])+(d[i,k]+d_reduced[i,k])*u[k]=n[j]*(X[i,j]+X_reduced[i,j])+y[i]+y_reduced[i]*

   The coefficient *d* is of rank two and *y* is of rank one both in
   `FunctionOnBoundary`. The coefficients
   *d_reduced* is of rank two and *y_reduced* is of rank one both in
   `ReducedFunctionOnBoundary`.

   Constraints take the form

   *u[i]=r[i]*  where  *q[i]>0*

   *r* and *q* are each rank one. Notice that at some locations not
   necessarily all components must have a constraint.

   The system of PDEs is symmetrical if

      - *A[i,j,k,l]=A[k,l,i,j]*
      - *A_reduced[i,j,k,l]=A_reduced[k,l,i,j]*
      - *B[i,j,k]=C[k,i,j]*
      - *B_reduced[i,j,k]=C_reduced[k,i,j]*
      - *D[i,k]=D[i,k]*
      - *D_reduced[i,k]=D_reduced[i,k]*
      - *d[i,k]=d[k,i]*
      - *d_reduced[i,k]=d_reduced[k,i]*

   `LinearPDE` also supports solution discontinuities over a contact region
   in the domain. To specify the conditions across the discontinuity we are
   using the generalised flux *J* which, in the case of a system of PDEs
   and several components of the solution, is defined as

   *J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]-X[i,j]-X_reduced[i,j]*

   For the case of single solution component and single PDE *J* is defined as

   *J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u-X[i]-X_reduced[i]*

   In the context of discontinuities *n* denotes the normal on the
   discontinuity pointing from side 0 towards side 1 calculated from
   `FunctionSpace.getNormal` of `FunctionOnContactZero`.
   For a system of PDEs the contact condition takes the form

   *n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]*

   where *J0* and *J1* are the fluxes on side 0 and side 1 of the
   discontinuity, respectively. *jump(u)*, which is the difference of the
   solution at side 1 and at side 0, denotes the jump of *u* across
   discontinuity along the normal calculated by `jump`.
   The coefficient *d_contact* is of rank two and *y_contact* is of rank one
   both in `FunctionOnContactZero` or
   `FunctionOnContactOne`.
   The coefficient *d_contact_reduced* is of rank two and *y_contact_reduced*
   is of rank one both in `ReducedFunctionOnContactZero`
   or `ReducedFunctionOnContactOne`.
   In case of a single PDE and a single component solution the contact
   condition takes the form

   *n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)*

   In this case the coefficient *d_contact* and *y_contact* are each scalar
   both in `FunctionOnContactZero` or
   `FunctionOnContactOne` and the coefficient
   *d_contact_reduced* and *y_contact_reduced* are each scalar both in
   `ReducedFunctionOnContactZero` or
   `ReducedFunctionOnContactOne`.

   Typical usage::

       p = LinearPDE(dom)
       p.setValue(A=kronecker(dom), D=1, Y=0.5)
       u = p.getSolution()

   """

   def __init__(self,domain,numEquations=None,numSolutions=None, isComplex=False,debug=False):
     """
     Initializes a new linear PDE.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the PDE coefficients.
     :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the PDE
                          coefficients.
     :param debug: if True debug information is printed

     """


     super(LinearPDE, self).__init__(domain,numEquations,numSolutions,isComplex,debug)
     #
     #   the coefficients of the PDE:
     #
     self.introduceCoefficients(
       A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR, isComplex),
       B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR, isComplex),
       D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       X=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE, isComplex),
       Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       d_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       y_contact=PDECoef(PDECoef.CONTACT,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       A_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR, isComplex),
       B_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       C_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR, isComplex),
       D_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       X_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE, isComplex),
       Y_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       d_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       y_reduced=PDECoef(PDECoef.BOUNDARY_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       d_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       y_contact_reduced=PDECoef(PDECoef.CONTACT_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       d_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR, isComplex),
       y_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE, isComplex),
       r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE, isComplex),
       q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH, False) )

   def __str__(self):
     """
     Returns the string representation of the PDE.

     :return: a simple representation of the PDE
     :rtype: ``str``
     """
     return "<LinearPDE %d>"%id(self)

   def getRequiredOperatorType(self):
      """
      Returns the system type which needs to be used by the current set up.
      """
      if self.isUsingLumping():
         return "__ESCRIPT_DATA"
      else:
         return self.getDomain().getSystemMatrixTypeId(self.getSolverOptions())

   def checkSymmetry(self,verbose=True):
      """
      Tests the PDE for symmetry.

      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      :type verbose: ``bool``
      :return: True if the PDE is symmetric
      :rtype: `bool`
      :note: This is a very expensive operation. It should be used for
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
      out=out and self.checkSymmetricTensor("d_dirac", verbose)
      return out

   def createOperator(self):
       """
       Returns an instance of a new operator.
       """
       optype=self.getRequiredOperatorType()
       self.trace("New operator of type %s is allocated."%optype)
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           optype)

   def getSolution(self):
       """
       Returns the solution of the PDE.

       :return: the solution
       :rtype: `Data`
       """
       option_class=self.getSolverOptions()
       if not self.isSolutionValid():
          mat,f=self.getSystem()
          if self.isUsingLumping():
             if not util.inf(abs(mat)) > 0.:
                 raise ZeroDivisionError("Lumped mass matrix has zero entry (try order 1 elements or HRZ lumping).")
             self.setSolution(f*1/mat)
          else:
             self.trace("PDE is resolved.")
             self.trace("solver options: %s"%str(option_class))
             self.setSolution(mat.solve(f,option_class))
       return self.getCurrentSolution()

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       :return: the discrete version of the PDE
       :rtype: ``tuple`` of `Operator` and
               `Data`
       """
       if not self.isOperatorValid() or not self.isRightHandSideValid():
          if self.isUsingLumping():
              if not self.isOperatorValid():
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                      raise TypeError("Lumped matrix requires same order for equations and unknowns")
                 if not self.getCoefficient("A").isEmpty():
                      raise ValueError("coefficient A in lumped matrix may not be present.")
                 if not self.getCoefficient("B").isEmpty():
                      raise ValueError("coefficient B in lumped matrix may not be present.")
                 if not self.getCoefficient("C").isEmpty():
                      raise ValueError("coefficient C in lumped matrix may not be present.")
                 if not self.getCoefficient("d_contact").isEmpty():
                      raise ValueError("coefficient d_contact in lumped matrix may not be present.")
                 if not self.getCoefficient("A_reduced").isEmpty():
                      raise ValueError("coefficient A_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("B_reduced").isEmpty():
                      raise ValueError("coefficient B_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("C_reduced").isEmpty():
                      raise ValueError("coefficient C_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("d_contact_reduced").isEmpty():
                      raise ValueError("coefficient d_contact_reduced in lumped matrix may not be present.")
                 D=self.getCoefficient("D")
                 d=self.getCoefficient("d")
                 D_reduced=self.getCoefficient("D_reduced")
                 d_reduced=self.getCoefficient("d_reduced")
                 d_dirac=self.getCoefficient("d_dirac")

                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escore.Data()
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escore.Data()

                 if not D_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        D_reduced_times_e=util.matrix_mult(D_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_reduced_times_e=D_reduced
                 else:
                    D_reduced_times_e=escore.Data()

                 if not d_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        d_reduced_times_e=util.matrix_mult(d_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_times_e=d_reduced
                 else:
                    d_reduced_times_e=escore.Data()

                 if not d_dirac.isEmpty():
                     if self.getNumSolutions()>1:
                        d_dirac_times_e=util.matrix_mult(d_dirac,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_dirac_e=d_dirac
                 else:
                    d_dirac_times_e=escore.Data()

                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 if hasattr(self.getDomain(), "addPDEToLumpedSystem") :
                    hrz_lumping=( self.getSolverOptions().getSolverMethod() ==  SolverOptions.HRZ_LUMPING )
                    self.addPDEToLumpedSystem(operator, D_times_e, d_times_e, d_dirac_times_e,  hrz_lumping )
                    self.addPDEToLumpedSystem(operator, D_reduced_times_e, d_reduced_times_e, escore.Data(), hrz_lumping)
                 else:
                    self.addPDEToRHS(operator, \
                                                 escore.Data(), \
                                                 D_times_e, \
                                                 d_times_e,\
                                                 escore.Data(),\
                                                 d_dirac_times_e)
                    self.addPDEToRHS(operator, \
                                                 escore.Data(), \
                                                 D_reduced_times_e, \
                                                 d_reduced_times_e,\
                                                 escore.Data(), \
                                                 escore.Data())
                 self.trace("New lumped operator has been built.")
              if not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.addPDEToRHS(righthandside, \
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"),\
                               self.getCoefficient("y"),\
                               self.getCoefficient("y_contact"), \
                               self.getCoefficient("y_dirac"))
                 self.addPDEToRHS(righthandside, \
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"),\
                               self.getCoefficient("y_reduced"),\
                               self.getCoefficient("y_contact_reduced"), \
                               escore.Data())
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
                 self.addPDEToSystem(operator,righthandside, \
                               self.getCoefficient("A"), \
                               self.getCoefficient("B"), \
                               self.getCoefficient("C"), \
                               self.getCoefficient("D"), \
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"), \
                               self.getCoefficient("d"), \
                               self.getCoefficient("y"), \
                               self.getCoefficient("d_contact"), \
                               self.getCoefficient("y_contact"), \
                               self.getCoefficient("d_dirac"), \
                               self.getCoefficient("y_dirac"))
                 self.addPDEToSystem(operator,righthandside, \
                               self.getCoefficient("A_reduced"), \
                               self.getCoefficient("B_reduced"), \
                               self.getCoefficient("C_reduced"), \
                               self.getCoefficient("D_reduced"), \
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"), \
                               self.getCoefficient("d_reduced"), \
                               self.getCoefficient("y_reduced"), \
                               self.getCoefficient("d_contact_reduced"), \
                               self.getCoefficient("y_contact_reduced"), \
                               escore.Data(), \
                               escore.Data())
                 self.insertConstraint(rhs_only=False)
                 self.trace("New system has been built.")
                 self.validOperator()
                 self.validRightHandSide()
             elif not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.addPDEToRHS(righthandside,
                               self.getCoefficient("X"), \
                               self.getCoefficient("Y"),\
                               self.getCoefficient("y"),\
                               self.getCoefficient("y_contact"), \
                               self.getCoefficient("y_dirac") )
                 self.addPDEToRHS(righthandside,
                               self.getCoefficient("X_reduced"), \
                               self.getCoefficient("Y_reduced"),\
                               self.getCoefficient("y_reduced"),\
                               self.getCoefficient("y_contact_reduced"), \
                               escore.Data())
                 self.insertConstraint(rhs_only=True)
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
             elif not self.isOperatorValid():
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 self.addPDEToSystem(operator,escore.Data(), \
                            self.getCoefficient("A"), \
                            self.getCoefficient("B"), \
                            self.getCoefficient("C"), \
                            self.getCoefficient("D"), \
                            escore.Data(), \
                            escore.Data(), \
                            self.getCoefficient("d"), \
                            escore.Data(),\
                            self.getCoefficient("d_contact"), \
                            escore.Data(),                   \
                            self.getCoefficient("d_dirac"),   \
                            escore.Data())
                 self.addPDEToSystem(operator,escore.Data(), \
                            self.getCoefficient("A_reduced"), \
                            self.getCoefficient("B_reduced"), \
                            self.getCoefficient("C_reduced"), \
                            self.getCoefficient("D_reduced"), \
                            escore.Data(), \
                            escore.Data(), \
                            self.getCoefficient("d_reduced"), \
                            escore.Data(),\
                            self.getCoefficient("d_contact_reduced"), \
                            escore.Data(),  \
                            escore.Data(),  \
                            escore.Data())
                 self.insertConstraint(rhs_only=False)
                 self.trace("New operator has been built.")
                 self.validOperator()
       self.setSystemStatus()
       self.trace("System status is %s."%self.getSystemStatus())
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

   def insertConstraint(self, rhs_only=False):
      """
      Applies the constraints defined by q and r to the PDE.

      :param rhs_only: if True only the right hand side is altered by the
                       constraint
      :type rhs_only: ``bool``
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
                 operator.copyWithMask(escore.Data(1.,q.getShape(),q.getFunctionSpace()),q)
             else:
                 row_q=escore.Data(q,self.getFunctionSpaceForEquation())
                 col_q=escore.Data(q,self.getFunctionSpaceForSolution())
                 u=self.createSolution()
                 u.copyWithMask(r_s,col_q)
                 righthandside-=operator*u
                 operator.nullifyRowsAndCols(row_q,col_q,1.)
         righthandside.copyWithMask(r_s,q)

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      :param coefficients: new values assigned to coefficients
      :keyword A: value for coefficient ``A``
      :type A: any type that can be cast to a `Data` object on
               `Function`
      :keyword A_reduced: value for coefficient ``A_reduced``
      :type A_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword B: value for coefficient ``B``
      :type B: any type that can be cast to a `Data` object on
               `Function`
      :keyword B_reduced: value for coefficient ``B_reduced``
      :type B_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword C: value for coefficient ``C``
      :type C: any type that can be cast to a `Data` object on
               `Function`
      :keyword C_reduced: value for coefficient ``C_reduced``
      :type C_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword D: value for coefficient ``D``
      :type D: any type that can be cast to a `Data` object on
               `Function`
      :keyword D_reduced: value for coefficient ``D_reduced``
      :type D_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword X: value for coefficient ``X``
      :type X: any type that can be cast to a `Data` object on
               `Function`
      :keyword X_reduced: value for coefficient ``X_reduced``
      :type X_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword Y: value for coefficient ``Y``
      :type Y: any type that can be cast to a `Data` object on
               `Function`
      :keyword Y_reduced: value for coefficient ``Y_reduced``
      :type Y_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword d: value for coefficient ``d``
      :type d: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_reduced: value for coefficient ``d_reduced``
      :type d_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunctionOnBoundary`
      :keyword y: value for coefficient ``y``
      :type y: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_contact: value for coefficient ``d_contact``
      :type d_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne`
                       or `FunctionOnContactZero`
      :keyword d_contact_reduced: value for coefficient ``d_contact_reduced``
      :type d_contact_reduced: any type that can be cast to a `Data`
                               object on `ReducedFunctionOnContactOne`
                               or `ReducedFunctionOnContactZero`
      :keyword y_contact: value for coefficient ``y_contact``
      :type y_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne`
                       or `FunctionOnContactZero`
      :keyword y_contact_reduced: value for coefficient ``y_contact_reduced``
      :type y_contact_reduced: any type that can be cast to a `Data`
                               object on `ReducedFunctionOnContactOne`
                               or `ReducedFunctionOnContactZero`
      :keyword d_dirac: value for coefficient ``d_dirac``
      :type d_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`
      :keyword y_dirac: value for coefficient ``y_dirac``
      :type y_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`
      :keyword r: values prescribed to the solution at the locations of
                  constraints
      :type r: any type that can be cast to a `Data` object on
               `Solution` or `ReducedSolution`
               depending on whether reduced order is used for the solution
      :keyword q: mask for location of constraints
      :type q: any type that can be cast to a `Data` object on
               `Solution` or `ReducedSolution`
               depending on whether reduced order is used for the
               representation of the equation
      :raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      super(LinearPDE,self).setValue(**coefficients)
      # check if the systrem is inhomogeneous:
      if len(coefficients)>0 and not self.isUsingLumping():
         q=self.getCoefficient("q")
         #if not isinstance(q, list):
            #raise KeyError("q paramter '%s' shouldn't be a symbol."%q)
         r=self.getCoefficient("r")
         if not q.isEmpty() and not r.isEmpty():
             if util.Lsup(q*r)>0.:
               self.trace("Inhomogeneous constraint detected.")
               self.invalidateSystem()


   def getResidual(self,u=None):
     """
     Returns the residual of u or the current solution if u is not present.

     :param u: argument in the residual calculation. It must be representable
               in `self.getFunctionSpaceForSolution()`. If u is not present
               or equals ``None`` the current solution is used.
     :type u: `Data` or None
     :return: residual of u
     :rtype: `Data`
     """
     if u is None:
        return self.getOperator()*self.getSolution()-self.getRightHandSide()
     else:
        return self.getOperator()*escore.Data(u,self.getFunctionSpaceForSolution())-self.getRightHandSide()

   def getFlux(self,u=None):
     """
     Returns the flux *J* for a given *u*.

     *J[i,j]=(A[i,j,k,l]+A_reduced[A[i,j,k,l]]*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])u[k]-X[i,j]-X_reduced[i,j]*

     or

     *J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[l]+(B[j]+B_reduced[j])u-X[j]-X_reduced[j]*

     :param u: argument in the flux. If u is not present or equals `None` the
               current solution is used.
     :type u: `Data` or None
     :return: flux
     :rtype: `Data`
     """
     if u is None: u=self.getSolution()
     if self.getNumEquations()>1:
       out = escore.Data(0.,(self.getNumEquations(),self.getDim()),self.getFunctionSpaceForCoefficient("X"))
     else:
       out = escore.Data(0.,(self.getDim(), ),self.getFunctionSpaceForCoefficient("X"))

     A=self.getCoefficient("A")
     if not A.isEmpty():
           out+=util.tensormult(A,util.grad(u,self.getFunctionSpaceForCoefficient("A")))

     B=self.getCoefficient("B")
     if not B.isEmpty():
           if B.getRank() == 1:
               out+=B * u
           else:
               out+=util.generalTensorProduct(B,u,axis_offset=1)

     X=self.getCoefficient("X")
     if not X.isEmpty():
           out-=X

     A_reduced=self.getCoefficient("A_reduced")
     if not A_reduced.isEmpty():
           out+=util.tensormult(A_reduced, util.grad(u,self.getFunctionSpaceForCoefficient("A_reduced"))) \

     B_reduced=self.getCoefficient("B_reduced")
     if not B_reduced.isEmpty():
           if B_reduced.getRank() == 1:
                out+=B_reduced*u
           else:
                out+=util.generalTensorProduct(B_reduced,u,axis_offset=1)

     X_reduced=self.getCoefficient("X_reduced")
     if not X_reduced.isEmpty():
           out-=X_reduced
     return out

class Poisson(LinearPDE):
   """
   Class to define a Poisson equation problem. This is generally a
   `LinearPDE` of the form

   *-grad(grad(u)[j])[j] = f*

   with natural boundary conditions

   *n[j]*grad(u)[j] = 0*

   and constraints:

   *u=0* where *q>0*

   """

   def __init__(self,domain,debug=False):
     """
     Initializes a new Poisson equation.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param debug: if True debug information is printed

     """
     super(Poisson, self).__init__(domain,1,1,debug)
     self.introduceCoefficients(
                        f=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
                        f_reduced=PDECoef(PDECoef.INTERIOR_REDUCED,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE))
     self.setSymmetryOn()

   def setValue(self,**coefficients):
     """
     Sets new values to coefficients.

     :param coefficients: new values assigned to coefficients
     :keyword f: value for right hand side *f*
     :type f: any type that can be cast to a `Scalar` object
              on `Function`
     :keyword q: mask for location of constraints
     :type q: any type that can be cast to a rank zero `Data`
              object on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(Poisson, self).setValue(**coefficients)


   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name`` of the general PDE.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient ``name``
     :rtype: `Data`
     :raise IllegalCoefficient: invalid coefficient name
     :note: This method is called by the assembling routine to map the Poisson
            equation onto the general PDE.
     """
     if name == "A" :
         return escore.Data(util.kronecker(self.getDim()),escore.Function(self.getDomain()))
     elif name == "Y" :
         return self.getCoefficient("f")
     elif name == "Y_reduced" :
         return self.getCoefficient("f_reduced")
     else:
         return super(Poisson, self).getCoefficient(name)

class Helmholtz(LinearPDE):
   """
   Class to define a Helmholtz equation problem. This is generally a
   `LinearPDE` of the form

   *omega*u - grad(k*grad(u)[j])[j] = f*

   with natural boundary conditions

   *k*n[j]*grad(u)[j] = g- alphau*

   and constraints:

   *u=r* where *q>0*

   """

   def __init__(self,domain,debug=False):
     """
     Initializes a new Helmholtz equation.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param debug: if True debug information is printed

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

     :param coefficients: new values assigned to coefficients
     :keyword omega: value for coefficient *omega*
     :type omega: any type that can be cast to a `Scalar`
                  object on `Function`
     :keyword k: value for coefficient *k*
     :type k: any type that can be cast to a `Scalar` object
              on `Function`
     :keyword f: value for right hand side *f*
     :type f: any type that can be cast to a `Scalar` object
              on `Function`
     :keyword alpha: value for right hand side *alpha*
     :type alpha: any type that can be cast to a `Scalar`
                  object on `FunctionOnBoundary`
     :keyword g: value for right hand side *g*
     :type g: any type that can be cast to a `Scalar` object
              on `FunctionOnBoundary`
     :keyword r: prescribed values *r* for the solution in constraints
     :type r: any type that can be cast to a `Scalar` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :keyword q: mask for the location of constraints
     :type q: any type that can be cast to a `Scalar` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
     :raise IllegalCoefficient: if an unknown coefficient keyword is used
     """
     super(Helmholtz, self).setValue(**coefficients)

   def getCoefficient(self,name):
     """
     Returns the value of the coefficient ``name`` of the general PDE.

     :param name: name of the coefficient requested
     :type name: ``string``
     :return: the value of the coefficient ``name``
     :rtype: `Data`
     :raise IllegalCoefficient: invalid name
     """
     if name == "A" :
         if self.getCoefficient("k").isEmpty():
              return escore.Data(numpy.identity(self.getDim()),escore.Function(self.getDomain()))
         else:
              return escore.Data(numpy.identity(self.getDim()),escore.Function(self.getDomain()))*self.getCoefficient("k")
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

class WavePDE(LinearPDE):
    """
    A class specifically for waves, passes along values to native implementation
    to save computational time.
    """
    def __init__(self,domain,c,numEquations=None,numSolutions=None,debug=False):
        """
        Initializes a new linear PDE.

        :param domain: domain of the PDE
        :type domain: `Domain`
        :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the PDE coefficients.
        :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the PDE
                          coefficients.
        :param debug: if True debug information is printed

        """
        super(WavePDE, self).__init__(domain,numEquations,numSolutions,debug)
        #
        #   the coefficients of the PDE:
        #
        self.introduceCoefficients(
           A=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
           B=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           C=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION,PDECoef.BY_DIM),PDECoef.OPERATOR),
           D=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           du=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,PDECoef.BY_DIM),PDECoef.RIGHTHANDSIDE),
           Y=PDECoef(PDECoef.INTERIOR,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
           d=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           y=PDECoef(PDECoef.BOUNDARY,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
           d_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
           y_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
           r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),
           q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH))
        self.assembler = self.getDomain().createAssembler("WaveAssembler", c)


    def getSystem(self):
        """
        Returns the operator and right hand side of the PDE.

        :return: the discrete version of the PDE
        :rtype: ``tuple`` of `Operator` and
                `Data`
        """
        if not self.isOperatorValid() or not self.isRightHandSideValid():
          if self.isUsingLumping():
              if not self.isOperatorValid():
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                      raise TypeError("Lumped matrix requires same order for equations and unknowns")
                 if not self.getCoefficient("A").isEmpty():
                      raise ValueError("coefficient A in lumped matrix may not be present.")
                 if not self.getCoefficient("B").isEmpty():
                      raise ValueError("coefficient B in lumped matrix may not be present.")
                 if not self.getCoefficient("C").isEmpty():
                      raise ValueError("coefficient C in lumped matrix may not be present.")

                 D=self.getCoefficient("D")
                 d=self.getCoefficient("d")
                 d_dirac=self.getCoefficient("d_dirac")

                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escore.Data()
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escore.Data()

                 if not d_dirac.isEmpty():
                     if self.getNumSolutions()>1:
                        d_dirac_times_e=util.matrix_mult(d_dirac,numpy.ones((self.getNumSolutions(),)))
                 else:
                    d_dirac_times_e=escore.Data()

                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 if hasattr(self.getDomain(), "addPDEToLumpedSystem") :
                    hrz_lumping=( self.getSolverOptions().getSolverMethod() ==  SolverOptions.HRZ_LUMPING )
                    self.addPDEToLumpedSystem(operator, D_times_e, d_times_e, d_dirac_times_e,  hrz_lumping )
                 else:
                    self.addToRHS(operator,
                        [("Y", D_times_e), ("y", d_times_e),
                         ("y_dirac", d_dirac_times_e)])
                 self.trace("New lumped operator has been built.")
              if not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.addToRHS(righthandside,
                                [(i, self.getCoefficient(i)) for i in
                                    ["du", "Y", "y", "y_dirac"]
                                ])
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
                 data = [(i, self.getCoefficient(i)) for i in ["A", "B", "C",
                                "D", "Y", "d", "y", "d_contact",
                                "y_contact", "d_dirac", "y_dirac", "du"]
                            ]
                 self.addToSystem(operator, righthandside, data)
                 self.insertConstraint(rhs_only=False)
                 self.trace("New system has been built.")
                 self.validOperator()
                 self.validRightHandSide()
             elif not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 self.addToRHS(righthandside,
                                [(i, self.getCoefficient(i)) for i in
                                    ["du", "Y", "y", "y_contact", "y_dirac"]
                                ])
                 self.insertConstraint(rhs_only=True)
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
             elif not self.isOperatorValid():
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 data = [(i, self.getCoefficient(i)) for i in ["A", "B", "C",
                        "D", "d", "d_contact", "d_dirac", "du"]]
                 self.addToSystem(operator, escore.Data(), data)
                 self.insertConstraint(rhs_only=False)
                 self.trace("New operator has been built.")
                 self.validOperator()
        self.setSystemStatus()
        self.trace("System status is %s."%self.getSystemStatus())
        return (self.getCurrentOperator(), self.getCurrentRightHandSide())


class LameEquation(LinearPDE):
    """
    Class to define a Lame equation problem. This problem is defined as:

    *-grad(mu*(grad(u[i])[j]+grad(u[j])[i]))[j] - grad(lambda*grad(u[k])[k])[j] = F_i -grad(sigma[ij])[j]*

    with natural boundary conditions:

    *n[j]*(mu*(grad(u[i])[j]+grad(u[j])[i]) + lambda*grad(u[k])[k]) = f_i +n[j]*sigma[ij]*

    and constraints:

    *u[i]=r[i]* where *q[i]>0*

    """

    def __init__(self,domain,debug=False,useFast=True):
        """
        Initializes a new Lame equation.

        :param domain: domain of the PDE
        :type domain: `Domain`
        :param debug: if True debug information is printed

        """
        self.fastAssembler = False
        if useFast and hasattr(domain, "setAssembler"):
            self.fastAssembler = True
            self.assembler = domain.createAssembler("LameAssembler", [])
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

        :param coefficients: new values assigned to coefficients
        :keyword lame_mu: value for coefficient *mu*
        :type lame_mu: any type that can be cast to a `Scalar`
                    object on `Function`
        :keyword lame_lambda: value for coefficient *lambda*
        :type lame_lambda: any type that can be cast to a `Scalar`
                        object on `Function`
        :keyword F: value for internal force *F*
        :type F: any type that can be cast to a `Vector` object
              on `Function`
        :keyword sigma: value for initial stress *sigma*
        :type sigma: any type that can be cast to a `Tensor`
                  object on `Function`
        :keyword f: value for external force *f*
        :type f: any type that can be cast to a `Vector` object
              on `FunctionOnBoundary`
        :keyword r: prescribed values *r* for the solution in constraints
        :type r: any type that can be cast to a `Vector` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
        :keyword q: mask for the location of constraints
        :type q: any type that can be cast to a `Vector` object
              on `Solution` or
              `ReducedSolution` depending on whether
              reduced order is used for the representation of the equation
        :raise IllegalCoefficient: if an unknown coefficient keyword is used
        """
        super(LameEquation, self).setValues(**coefficients)

    def getCoefficient(self,name):
        """
        Returns the value of the coefficient ``name`` of the general PDE.

        :param name: name of the coefficient requested
        :type name: ``string``
        :return: the value of the coefficient ``name``
        :rtype: `Data`
        :raise IllegalCoefficient: invalid coefficient name
        """

        if name == "A" :
            out = self.createCoefficient("A")
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
        return super(LameEquation, self).getCoefficient(name)

    def getSystem(self):
        """
        Returns the operator and right hand side of the PDE.

        :return: the discrete version of the PDE
        :rtype: ``tuple`` of `Operator` and
               `Data`
        """

        if not self.fastAssembler:
            return super(LameEquation, self).getSystem()

        if not self.isOperatorValid() or not self.isRightHandSideValid():
          if self.isUsingLumping():
              if not self.isOperatorValid():
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                      raise TypeError("Lumped matrix requires same order for equations and unknowns")
                 if not self.getCoefficient("lame_mu").isEmpty() and not self.getCoefficient("lame_mu").isEmpty():
                      raise ValueError("coefficient A in lumped matrix may not be present.")
                 if not self.getCoefficient("B").isEmpty():
                      raise ValueError("coefficient B in lumped matrix may not be present.")
                 if not self.getCoefficient("C").isEmpty():
                      raise ValueError("coefficient C in lumped matrix may not be present.")
                 if not self.getCoefficient("d_contact").isEmpty():
                      raise ValueError("coefficient d_contact in lumped matrix may not be present.")
                 if not self.getCoefficient("A_reduced").isEmpty():
                      raise ValueError("coefficient A_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("B_reduced").isEmpty():
                      raise ValueError("coefficient B_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("C_reduced").isEmpty():
                      raise ValueError("coefficient C_reduced in lumped matrix may not be present.")
                 if not self.getCoefficient("d_contact_reduced").isEmpty():
                      raise ValueError("coefficient d_contact_reduced in lumped matrix may not be present.")
                 D=self.getCoefficient("D")
                 d=self.getCoefficient("d")
                 D_reduced=self.getCoefficient("D_reduced")
                 d_reduced=self.getCoefficient("d_reduced")
                 d_dirac=self.getCoefficient("d_dirac")

                 if not D.isEmpty():
                     if self.getNumSolutions()>1:
                        D_times_e=util.matrix_mult(D,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_times_e=D
                 else:
                    D_times_e=escore.Data()
                 if not d.isEmpty():
                     if self.getNumSolutions()>1:
                        d_times_e=util.matrix_mult(d,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_times_e=d
                 else:
                    d_times_e=escore.Data()

                 if not D_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        D_reduced_times_e=util.matrix_mult(D_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        D_reduced_times_e=D_reduced
                 else:
                    D_reduced_times_e=escore.Data()

                 if not d_reduced.isEmpty():
                     if self.getNumSolutions()>1:
                        d_reduced_times_e=util.matrix_mult(d_reduced,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_times_e=d_reduced
                 else:
                    d_reduced_times_e=escore.Data()

                 if not d_dirac.isEmpty():
                     if self.getNumSolutions()>1:
                        d_dirac_times_e=util.matrix_mult(d_dirac,numpy.ones((self.getNumSolutions(),)))
                     else:
                        d_reduced_dirac_e=d_dirac
                 else:
                    d_dirac_times_e=escore.Data()

                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 if hasattr(self.getDomain(), "addPDEToLumpedSystem") :
                    hrz_lumping=( self.getSolverOptions().getSolverMethod() ==  SolverOptions.HRZ_LUMPING )
                    self.addPDEToLumpedSystem(operator, D_times_e, d_times_e, d_dirac_times_e,  hrz_lumping )
                    self.addPDEToLumpedSystem(operator, D_reduced_times_e, d_reduced_times_e, escore.Data(), hrz_lumping)
                 else:
                    self.addToRHS(operator, [
                                                 ("Y", D_times_e),
                                                 ("y", d_times_e),
                                                 ("y_dirac", d_dirac_times_e)])
                    self.addToRHS(operator, [
                                                 ("Y",D_reduced_times_e),
                                                 ("y",d_reduced_times_e)])
                 self.trace("New lumped operator has been built.")
              if not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 data = [(i, self.getCoefficient(i)) for i in ["X", "Y", "y",
                        "y_contact", "y_dirac"]]
                 self.addToRHS(righthandside, data)
                 data = [(i, self.getCoefficient(i+"_reduced")) for i in ["X",
                        "Y", "y", "y_contact"]]
                 self.addToRHS(righthandside, data)
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
                 data = [(i, self.getCoefficient(i)) for i in ["lame_mu",
                        "lame_lambda", "B", "C", "D",
                        "X", "Y", "d", "y", "d_contact", "y_contact",
                        "d_dirac", "y_dirac"]]
                 self.addToSystem(operator,righthandside, data)
                 data = [(i, self.getCoefficient(i+"_reduced")) for i in [
                            "A", "B", "C", "D",
                            "X", "Y", "d", "y", "d_contact", "y_contact"]
                        ]
                 self.addToSystem(operator,righthandside, data)
                 self.insertConstraint(rhs_only=False)
                 self.trace("New system has been built.")
                 self.validOperator()
                 self.validRightHandSide()
             elif not self.isRightHandSideValid():
                 self.resetRightHandSide()
                 righthandside=self.getCurrentRightHandSide()
                 data = [(i, self.getCoefficient(i)) for i in ["X", "Y", "y",
                        "y_contact", "y_dirac"]]
                 self.addToRHS(righthandside, data)
                 data = [(i, self.getCoefficient(i+"_reduced")) for i in [
                        "X", "Y", "y","y_contact"]]
                 self.addToRHS(righthandside, data)
                 self.insertConstraint(rhs_only=True)
                 self.trace("New right hand side has been built.")
                 self.validRightHandSide()
             elif not self.isOperatorValid():
                 self.resetOperator()
                 operator=self.getCurrentOperator()
                 data = [(i, self.getCoefficient(i)) for i in ["lame_mu",
                        "lame_lambda", "B","C","D","d","d_contact","d_dirac"]]
                 self.addToSystem(operator, data)
                 data = [(i, self.getCoefficient(i)) for i in [
                            "lame_mu","lame_lambda"]
                        ] + [(i, self.getCoefficient(i+"_reduced")) for i in [
                            "B","C","D","d","d_contact"]
                        ]
                 self.addToSystem(operator,data)
                 self.insertConstraint(rhs_only=False)
                 self.trace("New operator has been built.")
                 self.validOperator()
        self.setSystemStatus()
        self.trace("System status is %s."%self.getSystemStatus())
        return (self.getCurrentOperator(), self.getCurrentRightHandSide())

def LinearSinglePDE(domain, isComplex=False, debug=False):
   """
   Defines a single linear PDE.

   :param domain: domain of the PDE
   :type domain: `Domain`
   :param isComplex: if true, this coefficient is part of a complex-valued
                PDE and values will be converted to complex.
   :type isComplex: ``boolean``
   :param debug: if True debug information is printed
   :rtype: `LinearPDE`
   """
   return LinearPDE(domain,numEquations=1,numSolutions=1,  isComplex=isComplex, debug=debug)

def LinearPDESystem(domain, isComplex=False, debug=False):
   """
   Defines a system of linear PDEs.

   :param domain: domain of the PDEs
   :type domain: `Domain`
   :param isComplex: if true, this coefficient is part of a complex-valued
                PDE and values will be converted to complex.
   :type isComplex: ``boolean``
   :param debug: if True debug information is printed
   :rtype: `LinearPDE`
   """
   return LinearPDE(domain,numEquations=domain.getDim(),numSolutions=domain.getDim(), isComplex=isComplex, debug=debug)


class TransportPDE(LinearProblem):
   """
   This class is used to define a transport problem given by a general linear,
   time dependent, second order PDE for an unknown, non-negative,
   time-dependent function *u* on a given domain defined through a
   `Domain` object.

   For a single equation with a solution with a single component the transport
   problem is defined in the following form:

   *(M+M_reduced)*u_t=-(grad(A[j,l]+A_reduced[j,l]) * grad(u)[l]+(B[j]+B_reduced[j])u)[j]+(C[l]+C_reduced[l])*grad(u)[l]+(D+D_reduced)-grad(X+X_reduced)[j,j]+(Y+Y_reduced)*

   where *u_t* denotes the time derivative of *u* and *grad(F)* denotes the
   spatial derivative of *F*. Einstein's summation convention,  ie. summation
   over indexes appearing twice in a term of a sum performed, is used.
   The coefficients *M*, *A*, *B*, *C*, *D*, *X* and *Y* have to be
   specified through `Data` objects in `Function`
   and the coefficients *M_reduced*, *A_reduced*, *B_reduced*, *C_reduced*,
   *D_reduced*, *X_reduced* and *Y_reduced* have to be specified through
   `Data` objects in `ReducedFunction`.
   It is also allowed to use objects that can be converted into such
   `Data` objects. *M* and *M_reduced* are scalar, *A* and
   *A_reduced* are rank two, *B*, *C*, *X*, *B_reduced*, *C_reduced* and
   *X_reduced* are rank one and *D*, *D_reduced*, *Y* and *Y_reduced* are
   scalar.

   The following natural boundary conditions are considered:

   *n[j]*((A[i,j]+A_reduced[i,j])*grad(u)[l]+(B+B_reduced)[j]*u+X[j]+X_reduced[j])+(d+d_reduced)*u+y+y_reduced=(m+m_reduced)*u_t*

   where *n* is the outer normal field. Notice that the coefficients *A*,
   *A_reduced*, *B*, *B_reduced*, *X* and *X_reduced* are defined in the
   transport problem. The coefficients *m*, *d* and *y* are each a scalar in
   `FunctionOnBoundary` and the coefficients
   *m_reduced*, *d_reduced* and *y_reduced* are each a scalar in
   `ReducedFunctionOnBoundary`.

   Constraints for the solution prescribing the value of the solution at
   certain locations in the domain have the form

   *u_t=r* where *q>0*

   *r* and *q* are each scalar where *q* is the characteristic function
   defining where the constraint is applied. The constraints override any other
   condition set by the transport problem or the boundary condition.

   The transport problem is symmetrical if

   *A[i,j]=A[j,i]* and *B[j]=C[j]* and *A_reduced[i,j]=A_reduced[j,i]* and
   *B_reduced[j]=C_reduced[j]*

   For a system and a solution with several components the transport problem
   has the form

   *(M[i,k]+M_reduced[i,k]) * u[k]_t=-grad((A[i,j,k,l]+A_reduced[i,j,k,l]) * grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k]) * u[k])[j]+(C[i,k,l]+C_reduced[i,k,l]) * grad(u[k])[l]+(D[i,k]+D_reduced[i,k] * u[k]-grad(X[i,j]+X_reduced[i,j])[j]+Y[i]+Y_reduced[i]*

   *A* and *A_reduced* are of rank four, *B*, *B_reduced*, *C* and
   *C_reduced* are each of rank three, *M*, *M_reduced*, *D*, *D_reduced*,
   *X_reduced* and *X* are each of rank two and *Y* and *Y_reduced* are of
   rank one. The natural boundary conditions take the form:

   *n[j]*((A[i,j,k,l]+A_reduced[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]+X[i,j]+X_reduced[i,j])+(d[i,k]+d_reduced[i,k])*u[k]+y[i]+y_reduced[i]= (m[i,k]+m_reduced[i,k])*u[k]_t*

   The coefficient *d* and *m* are of rank two and *y* is of rank one with
   all in `FunctionOnBoundary`. The coefficients
   *d_reduced* and *m_reduced* are of rank two and *y_reduced* is of rank
   one all in `ReducedFunctionOnBoundary`.

   Constraints take the form

   *u[i]_t=r[i]* where *q[i]>0*

   *r* and *q* are each rank one. Notice that at some locations not
   necessarily all components must have a constraint.

   The transport problem is symmetrical if

      - *M[i,k]=M[i,k]*
      - *M_reduced[i,k]=M_reduced[i,k]*
      - *A[i,j,k,l]=A[k,l,i,j]*
      - *A_reduced[i,j,k,l]=A_reduced[k,l,i,j]*
      - *B[i,j,k]=C[k,i,j]*
      - *B_reduced[i,j,k]=C_reduced[k,i,j]*
      - *D[i,k]=D[i,k]*
      - *D_reduced[i,k]=D_reduced[i,k]*
      - *m[i,k]=m[k,i]*
      - *m_reduced[i,k]=m_reduced[k,i]*
      - *d[i,k]=d[k,i]*
      - *d_reduced[i,k]=d_reduced[k,i]*
      - *d_dirac[i,k]=d_dirac[k,i]*

   `TransportPDE` also supports solution discontinuities over a contact region
   in the domain. To specify the conditions across the discontinuity we are
   using the generalised flux *J* which, in the case of a system of PDEs and
   several components of the solution, is defined as

   *J[i,j]=(A[i,j,k,l]+A_reduced[[i,j,k,l])*grad(u[k])[l]+(B[i,j,k]+B_reduced[i,j,k])*u[k]+X[i,j]+X_reduced[i,j]*

   For the case of single solution component and single PDE *J* is defined as

   *J[j]=(A[i,j]+A_reduced[i,j])*grad(u)[j]+(B[i]+B_reduced[i])*u+X[i]+X_reduced[i]*

   In the context of discontinuities *n* denotes the normal on the
   discontinuity pointing from side 0 towards side 1 calculated from
   `FunctionSpace.getNormal` of `FunctionOnContactZero`.
   For a system of transport problems the contact condition takes the form

   *n[j]*J0[i,j]=n[j]*J1[i,j]=(y_contact[i]+y_contact_reduced[i])- (d_contact[i,k]+d_contact_reduced[i,k])*jump(u)[k]*

   where *J0* and *J1* are the fluxes on side 0 and side 1 of the
   discontinuity, respectively. *jump(u)*, which is the difference of the
   solution at side 1 and at side 0, denotes the jump of *u* across
   discontinuity along the normal calculated by `jump`.
   The coefficient *d_contact* is of rank two and *y_contact* is of rank one
   both in `FunctionOnContactZero` or `FunctionOnContactOne`.
   The coefficient *d_contact_reduced* is of rank two and *y_contact_reduced*
   is of rank one both in `ReducedFunctionOnContactZero` or `ReducedFunctionOnContactOne`.
   In case of a single PDE and a single component solution the contact
   condition takes the form

   *n[j]*J0_{j}=n[j]*J1_{j}=(y_contact+y_contact_reduced)-(d_contact+y_contact_reduced)*jump(u)*

   In this case the coefficient *d_contact* and *y_contact* are each scalar
   both in `FunctionOnContactZero` or
   `FunctionOnContactOne` and the coefficient
   *d_contact_reduced* and *y_contact_reduced* are each scalar both in
   `ReducedFunctionOnContactZero` or
   `ReducedFunctionOnContactOne`.

   Typical usage::

       p = TransportPDE(dom)
       p.setValue(M=1., C=[-1.,0.])
       p.setInitialSolution(u=exp(-length(dom.getX()-[0.1,0.1])**2)
       t = 0
       dt = 0.1
       while (t < 1.):
           u = p.solve(dt)

   """
   def __init__(self,domain,numEquations=None,numSolutions=None, useBackwardEuler=None, debug=False):
     """
     Initializes a transport problem.

     :param domain: domain of the PDE
     :type domain: `Domain`
     :param numEquations: number of equations. If ``None`` the number of
                          equations is extracted from the coefficients.
     :param numSolutions: number of solution components. If ``None`` the number
                          of solution components is extracted from the
                          coefficients.
     :param debug: if True debug information is printed
     """
     super(TransportPDE, self).__init__(domain,numEquations,numSolutions,debug)
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
       d_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,PDECoef.BY_SOLUTION),PDECoef.OPERATOR),
       y_dirac=PDECoef(PDECoef.DIRACDELTA,(PDECoef.BY_EQUATION,),PDECoef.RIGHTHANDSIDE),
       r=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.RIGHTHANDSIDE),
       q=PDECoef(PDECoef.SOLUTION,(PDECoef.BY_SOLUTION,),PDECoef.BOTH) )
     if not useBackwardEuler is None:
        import warnings
        warnings.warn("Argument useBackwardEuler has expired and will be removed in a later release. Please use SolverOptions.setODESolver() instead.", PendingDeprecationWarning, stacklevel=2)
        if useBackwardEuler: self.getSolverOptions().setODESolver(SolverOptions.BACKWARD_EULER)

   def __str__(self):
     """
     Returns the string representation of the transport problem.

     :return: a simple representation of the transport problem
     :rtype: ``str``
     """
     return "<TransportPDE %d>"%id(self)

   def checkSymmetry(self,verbose=True):
      """
      Tests the transport problem for symmetry.

      :param verbose: if set to True or not present a report on coefficients
                      which break the symmetry is printed.
      :type verbose: ``bool``
      :return:  True if the PDE is symmetric
      :rtype: ``bool``
      :note: This is a very expensive operation. It should be used for
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
      out=out and self.checkSymmetricTensor("d_dirac", verbose)
      return out

   def setValue(self,**coefficients):
      """
      Sets new values to coefficients.

      :param coefficients: new values assigned to coefficients
      :keyword M: value for coefficient ``M``
      :type M: any type that can be cast to a `Data` object on
               `Function`
      :keyword M_reduced: value for coefficient ``M_reduced``
      :type M_reduced: any type that can be cast to a `Data`
                       object on `Function`
      :keyword A: value for coefficient ``A``
      :type A: any type that can be cast to a `Data` object on
               `Function`
      :keyword A_reduced: value for coefficient ``A_reduced``
      :type A_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword B: value for coefficient ``B``
      :type B: any type that can be cast to a `Data` object on
               `Function`
      :keyword B_reduced: value for coefficient ``B_reduced``
      :type B_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword C: value for coefficient ``C``
      :type C: any type that can be cast to a `Data` object on
               `Function`
      :keyword C_reduced: value for coefficient ``C_reduced``
      :type C_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword D: value for coefficient ``D``
      :type D: any type that can be cast to a `Data` object on
               `Function`
      :keyword D_reduced: value for coefficient ``D_reduced``
      :type D_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword X: value for coefficient ``X``
      :type X: any type that can be cast to a `Data` object on
               `Function`
      :keyword X_reduced: value for coefficient ``X_reduced``
      :type X_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword Y: value for coefficient ``Y``
      :type Y: any type that can be cast to a `Data` object on
               `Function`
      :keyword Y_reduced: value for coefficient ``Y_reduced``
      :type Y_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunction`
      :keyword m: value for coefficient ``m``
      :type m: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword m_reduced: value for coefficient ``m_reduced``
      :type m_reduced: any type that can be cast to a `Data`
                       object on `FunctionOnBoundary`
      :keyword d: value for coefficient ``d``
      :type d: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_reduced: value for coefficient ``d_reduced``
      :type d_reduced: any type that can be cast to a `Data`
                       object on `ReducedFunctionOnBoundary`
      :keyword y: value for coefficient ``y``
      :type y: any type that can be cast to a `Data` object on
               `FunctionOnBoundary`
      :keyword d_contact: value for coefficient ``d_contact``
      :type d_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne` or `FunctionOnContactZero`
      :keyword d_contact_reduced: value for coefficient ``d_contact_reduced``
      :type d_contact_reduced: any type that can be cast to a `Data` object on `ReducedFunctionOnContactOne` or `ReducedFunctionOnContactZero`
      :keyword y_contact: value for coefficient ``y_contact``
      :type y_contact: any type that can be cast to a `Data`
                       object on `FunctionOnContactOne` or `FunctionOnContactZero`
      :keyword y_contact_reduced: value for coefficient ``y_contact_reduced``
      :type y_contact_reduced: any type that can be cast to a `Data` object on `ReducedFunctionOnContactOne` or `ReducedFunctionOnContactZero`

      :keyword d_dirac: value for coefficient ``d_dirac``
      :type d_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`
      :keyword y_dirac: value for coefficient ``y_dirac``
      :type y_dirac: any type that can be cast to a `Data` object on `DiracDeltaFunctions`

      :keyword r: values prescribed to the solution at the locations of constraints
      :type r: any type that can be cast to a `Data` object on
               `Solution` or `ReducedSolution`
               depending on whether reduced order is used for the solution
      :keyword q: mask for the location of constraints
      :type q: any type that can be cast to a `Data` object on
               `Solution` or
               `ReducedSolution` depending on whether
               reduced order is used for the representation of the equation
      :raise IllegalCoefficient: if an unknown coefficient keyword is used
      """
      super(TransportPDE,self).setValue(**coefficients)

   def createOperator(self):
       """
       Returns an instance of a new transport operator.
       """
       optype=self.getRequiredOperatorType()
       self.trace("New Transport problem of type %s is allocated."%optype)
       return self.getDomain().newTransportProblem( \
                               self.getNumEquations(), \
                               self.getFunctionSpaceForSolution(), \
                               optype)


   def getRequiredOperatorType(self):
      """
      Returns the system type which needs to be used by the current set up.

      :return: a code to indicate the type of transport problem scheme used
      :rtype: ``float``
      """
      solver_options=self.getSolverOptions()
      return self.getDomain().getTransportTypeId(solver_options.getSolverMethod(), solver_options.getPreconditioner(),solver_options.getPackage(), solver_options.isSymmetric())

   def getUnlimitedTimeStepSize(self):
      """
      Returns the value returned by the ``getSafeTimeStepSize`` method to
      indicate no limit on the safe time step size.

       :return: the value used to indicate that no limit is set to the time
                step size
       :rtype: ``float``
       :note: Typically the biggest positive float is returned
      """
      return self.getOperator().getUnlimitedTimeStepSize()

   def getSafeTimeStepSize(self):
       """
       Returns a safe time step size to do the next time step.

       :return: safe time step size
       :rtype: ``float``
       :note: If not ``getSafeTimeStepSize()`` < ``getUnlimitedTimeStepSize()``
              any time step size can be used.
       """
       return self.getOperator().getSafeTimeStepSize()

   #====================================================================
   def getSolution(self, dt=None, u0=None):
       """
       Returns the solution by marching forward by time step dt.
       If ''u0'' is present, ''u0'' is used as the initial value otherwise
       the solution from the last call is used.

       :param dt: time step size. If ``None`` the last solution is returned.
       :type dt: positive ``float`` or ``None``
       :param u0: new initial solution or ``None``
       :type u0: any object that can be interpolated to a `Data`
                object on `Solution` or `ReducedSolution`
       :return: the solution
       :rtype: `Data`
       """
       if not dt is None:
          option_class=self.getSolverOptions()
          if dt<=0:
              raise ValueError("step size needs to be positive.")
          if u0 is None:
              u0=self.getCurrentSolution()
          else:
              u0=util.interpolate(u0,self.getFunctionSpaceForSolution())
              if self.getNumSolutions() == 1:
                if u0.getShape()!=():
                  raise ValueError("Illegal shape %s of initial solution."%(u0.getShape(),))
              else:
                 if u0.getShape()!=(self.getNumSolutions(),):
                   raise ValueError("Illegal shape %s of initial solution."%(u0.getShape(),))
          self.setSolution(self.getOperator().solve(u0, self.getRightHandSide(),dt,option_class))
          self.validSolution()
       return self.getCurrentSolution()

   def setInitialSolution(self,u):
       """
       Sets the initial solution.

       :param u: initial solution
       :type u: any object that can be interpolated to a `Data`
                object on `Solution` or `ReducedSolution`
       """
       u2=util.interpolate(u,self.getFunctionSpaceForSolution())
       if self.getNumSolutions() == 1:
          if u2.getShape()!=():
              raise ValueError("Illegal shape %s of initial solution."%(u2.getShape(),))
       else:
          if u2.getShape()!=(self.getNumSolutions(),):
              raise ValueError("Illegal shape %s of initial solution."%(u2.getShape(),))
       self.setSolution(u2,validate=False)

   def addPDEToTransportProblem(self, operator,righthandside, M, A, B, C, D, X, Y,
            d, y, d_contact, y_contact, d_dirac, y_dirac):
        """
        Adds the PDE in the given form to the system matrix
        :param tp:
        :type tp: `TransportProblemAdapter`
        :param source:
        :type source: `Data`
        :param data:
        :type data: `list`"
        :param M:
        :type M: `Data`
        :param A:
        :type A: `Data`
        :param B:
        :type B: `Data`
        :param C:
        :type C: `Data`
        :param D:
        :type D: `Data`
        :param X:
        :type X: `Data`
        :param Y:
        :type Y: `Data`
        :param d:
        :type d: `Data`
        :param y:
        :type y: `Data`
        :param d_contact:
        :type d_contact: `Data`
        :param y_contact:
        :type y_contact: `Data`
        :param d_contact:
        :type d_contact: `Data`
        :param y_contact:
        :type y_contact: `Data`
        """
        if self.domainSupportsAssemblers:
            data = [("M", M), ("A", A), ("B", B), ("C", C), ("D", D), ("X", X), ("Y", Y),
                    ("d", d), ("y", y), ("d_contact", d_contact),
                    ("y_contact", y_contact), ("d_dirac", d_dirac),
                    ("y_dirac", y_dirac)]
            self.getDomain().addPDEToTransportProblem(operator,righthandside,
                    data, self.assembler)
        else:
            self.getDomain().addPDEToTransportProblem(operator,righthandside, M, A, B, C, D,
                    X, Y, d, y, d_contact, y_contact, d_dirac, y_dirac)

   def getSystem(self):
       """
       Returns the operator and right hand side of the PDE.

       :return: the discrete version of the PDE
       :rtype: ``tuple`` of `Operator` and
               `Data`

       """
       if not self.isOperatorValid() or not self.isRightHandSideValid():
          self.resetRightHandSide()
          righthandside=self.getCurrentRightHandSide()
          self.resetOperator()
          operator=self.getCurrentOperator()
          self.addPDEToTransportProblem(
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
                            self.getCoefficient("y_contact"),
                            self.getCoefficient("d_dirac"),
                            self.getCoefficient("y_dirac") )
          self.addPDEToTransportProblem(
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
                            self.getCoefficient("y_contact_reduced"),
                            escore.Data(),
                            escore.Data() )
          operator.insertConstraint(righthandside,self.getCoefficient("q"),self.getCoefficient("r"))
          self.trace("New system has been built.")
          self.validOperator()
          self.validRightHandSide()
       self.setSystemStatus()
       self.trace("System status is %s."%self.getSystemStatus())
       return (self.getCurrentOperator(), self.getCurrentRightHandSide())

   def setDebug(self, flag):
     """
     Switches debug output on if ``flag`` is True,
     otherwise it is switched off.

     :param flag: desired debug status
     :type flag: ``bool``
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

def SingleTransportPDE(domain, debug=False):
   """
   Defines a single transport problem

   :param domain: domain of the PDE
   :type domain: `Domain`
   :param debug: if True debug information is printed
   :rtype: `TransportPDE`
   """
   return TransportPDE(domain,numEquations=1,numSolutions=1, debug=debug)


