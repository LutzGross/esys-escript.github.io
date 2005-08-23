# $Id$

## @file linearPDEs.py

"""
Functions and classes for linear PDEs
"""

import escript
import util
import numarray


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

def _CompTuple2(t1,t2):
      """
      Compare two tuples
   
      @param t1 The first tuple
      @param t2 The second tuple
   
      """
   
      dif=t1[0]+t1[1]-(t2[0]+t2[1])
      if dif<0: return 1
      elif dif>0: return -1
      else: return 0
   
class PDECoefficient:
    """
    A class for PDE coefficients
    """
    # identifier for location of Data objects defining COEFFICIENTS
    INTERIOR=0
    BOUNDARY=1
    CONTACT=2
    CONTINUOUS=3
    # identifier in the pattern of COEFFICIENTS:
    # the pattern is a tuple of EQUATION,SOLUTION,DIM where DIM represents the spatial dimension, EQUATION the number of equations and SOLUTION the
    # number of unknowns.
    EQUATION=3
    SOLUTION=4
    DIM=5
    # indicator for what is altered if the coefficient is altered:
    OPERATOR=5
    RIGHTHANDSIDE=6
    BOTH=7
    def __init__(self,where,pattern,altering):
       """
       Initialise a PDE Coefficient type
       """
       self.what=where
       self.pattern=pattern
       self.altering=altering
       self.resetValue()

    def resetValue(self):
       """
       resets coefficient value to default
       """
       self.value=escript.Data()

    def getFunctionSpace(self,domain):
       """
       defines the FunctionSpace of the coefficient on the domain

       @param domain:
       """
       if self.what==self.INTERIOR: return escript.Function(domain)
       elif self.what==self.BOUNDARY: return escript.FunctionOnBoundary(domain)
       elif self.what==self.CONTACT: return escript.FunctionOnContactZero(domain)
       elif self.what==self.CONTINUOUS: return escript.ContinuousFunction(domain)

    def getValue(self):
       """
       returns the value of the coefficient:
       """
       return self.value

    def setValue(self,domain,numEquations=1,numSolutions=1,newValue=None):
       """
       set the value of the coefficient to new value
       """
       if newValue==None:
           newValue=escript.Data()
       elif isinstance(newValue,escript.Data):
           if not newValue.isEmpty():
              newValue=escript.Data(newValue,self.getFunctionSpace(domain))
       else:
           newValue=escript.Data(newValue,self.getFunctionSpace(domain))
       if not newValue.isEmpty():
           if not self.getShape(domain,numEquations,numSolutions)==newValue.getShape():
               raise IllegalCoefficientValue,"Expected shape for coefficient %s is %s but actual shape is %s."%(self.getShape(domain,numEquations,numSolutions),newValue.getShape())
       self.value=newValue

    def isAlteringOperator(self):
        """
	return true if the operator of the PDE is changed when the coefficient is changed
	"""
        if self.altering==self.OPERATOR or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isAlteringRightHandSide(self):
        """
	return true if the right hand side of the PDE is changed when the coefficient is changed
	"""
        if self.altering==self.RIGHTHANDSIDE or self.altering==self.BOTH:
            return not None
        else:
            return None

    def estimateNumEquationsAndNumSolutions(self,domain,shape=()):
       """
       tries to estimate the number of equations in a given tensor shape for a given spatial dimension dim

       @param shape:
       @param dim:
       """
       dim=domain.getDim()
       if len(shape)>0:
           num=max(shape)+1
       else:
           num=1
       search=[]
       for u in range(num):
          for e in range(num):
             search.append((e,u))
       search.sort(_CompTuple2)
       for item in search:
             s=self.getShape(domain,item[0],item[1])
             if len(s)==0 and len(shape)==0:
                 return (1,1)
             else:
                 if s==shape: return item
       return None

    def getShape(self,domain,numEquations=1,numSolutions=1):
        """
	builds the required shape for a given number of equations e, number of unknowns u and spatial dimension dim

	@param e:
	@param u:
	@param dim:
	"""
        dim=domain.getDim()
        s=()
        for i in self.pattern:
             if i==self.EQUATION:
                if numEquations>1: s=s+(numEquations,)
             elif i==self.SOLUTION:
                if numSolutions>1: s=s+(numSolutions,)
             else:
                s=s+(dim,)
        return s

class LinearPDE:
   """
   Class to define a linear PDE of the form

   \f[
     -(A_{ijkl}u_{k,l})_{,j} -(B_{ijk}u_k)_{,j} + C_{ikl}u_{k,l} +D_{ik}u_k = - (X_{ij})_{,j} + Y_i
   \f]

   with boundary conditons:

   \f[
   n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d_{ik}u_k = - n_j*X_{ij} + y_i
   \f]

   and contact conditions

   \f[
   n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d_contact_{ik}[u_k] = - n_j*X_{ij} + y_contact_i
   \f]

   and constraints:

   \f[
   u_i=r_i \quad \mathrm{where} \quad q_i>0
   \f]

   """
   TOL=1.e-13
   # solver methods
   UNKNOWN=-1
   DEFAULT_METHOD=0
   DIRECT=1
   CHOLEVSKY=2
   PCG=3
   CR=4
   CGS=5
   BICGSTAB=6
   SSOR=7
   ILU0=8
   ILUT=9
   JACOBI=10
   GMRES=11
   PRES20=12
   LUMPING=13
   # matrix reordering:
   NO_REORDERING=0
   MINIMUM_FILL_IN=1
   NESTED_DISSECTION=2
   # important keys in the dictonary used to hand over options to the solver library:
   METHOD_KEY="method"
   SYMMETRY_KEY="symmetric"
   TOLERANCE_KEY="tolerance"


   def __init__(self,domain,numEquations=None,numSolutions=None,debug=False):
     """
     initializes a new linear PDE

     @param domain: domain of the PDE
     @type domain: L{Domain}
     @param numEquations: number of equations. If numEquations==None the number of equations
                          is exracted from the PDE coefficients.
     @param numSolutions: number of solution components. If  numSolutions==None the number of solution components
                          is exracted from the PDE coefficients.
     @param debug: if True debug informations are printed.


     """
     #
     #   the coefficients of the general PDE:
     #
     self.__COEFFICIENTS_OF_GENEARL_PDE={
       "A"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.DIM,PDECoefficient.SOLUTION,PDECoefficient.DIM),PDECoefficient.OPERATOR),
       "B"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.DIM,PDECoefficient.SOLUTION),PDECoefficient.OPERATOR),
       "C"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.SOLUTION,PDECoefficient.DIM),PDECoefficient.OPERATOR),
       "D"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.SOLUTION),PDECoefficient.OPERATOR),
       "X"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.DIM),PDECoefficient.RIGHTHANDSIDE),
       "Y"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "d"         : PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.EQUATION,PDECoefficient.SOLUTION),PDECoefficient.OPERATOR),
       "y"         : PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "d_contact" : PDECoefficient(PDECoefficient.CONTACT,(PDECoefficient.EQUATION,PDECoefficient.SOLUTION),PDECoefficient.OPERATOR),
       "y_contact" : PDECoefficient(PDECoefficient.CONTACT,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "r"         : PDECoefficient(PDECoefficient.CONTINUOUS,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "q"         : PDECoefficient(PDECoefficient.CONTINUOUS,(PDECoefficient.SOLUTION,),PDECoefficient.BOTH)}

     # COEFFICIENTS can be overwritten by subclasses:
     self.COEFFICIENTS=self.__COEFFICIENTS_OF_GENEARL_PDE
     # initialize attributes
     self.__debug=debug
     self.__domain=domain
     self.__numEquations=numEquations
     self.__numSolutions=numSolutions
     self.__resetSystem()

     # set some default values:
     self.__homogeneous_constraint=True
     self.__row_function_space=escript.Solution(self.__domain)
     self.__column_function_space=escript.Solution(self.__domain)
     self.__tolerance=1.e-8
     self.__solver_method=self.DEFAULT_METHOD
     self.__matrix_type=self.__domain.getSystemMatrixTypeId(self.DEFAULT_METHOD,False)
     self.__sym=False

     self.resetCoefficients()
     self.trace("PDE Coeffients are %s"%str(self.COEFFICIENTS.keys()))
   # =============================================================================
   #    general stuff:
   # =============================================================================
   def __str__(self):
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

     @param name: name of the coefficient enquired.
     @type name: C{string}
     """
     if self.__debug: print "%s: %s"%(str(self),text)

   # =============================================================================
   # some service functions:
   # =============================================================================
   def getDomain(self):
     """
     returns the domain of the PDE
     
     @return : the domain of the PDE
     @rtype : C{Domain}

     """
     return self.__domain

   def getDim(self):
     """
     returns the spatial dimension of the PDE

     @return : the spatial dimension of the PDE domain
     @rtype : C{int}
     """
     return self.getDomain().getDim()

   def getNumEquations(self):
     """
     returns the number of equations

     @return : the number of equations
     @rtype : C{int}
     @raise UndefinedPDEError: if the number of equations is not be specified yet.
     """
     if self.__numEquations==None:
         raise UndefinedPDEError,"Number of equations is undefined. Please specify argument numEquations."
     else:
         return self.__numEquations

   def getNumSolutions(self):
     """
     returns the number of unknowns

     @return : the number of unknowns
     @rtype : C{int}
     @raise UndefinedPDEError: if the number of unknowns is not be specified yet.
     """
     if self.__numSolutions==None:
        raise UndefinedPDEError,"Number of solution is undefined. Please specify argument numSolutions."
     else:
        return self.__numSolutions

   def getFunctionSpaceForEquation(self):
     """
     returns the L{escript.FunctionSpace} used to discretize the equation
     
     @return : representation space of equation
     @rtype : L{escript.FunctionSpace}

     """
     return self.__row_function_space

   def getFunctionSpaceForSolution(self):
     """
     returns the L{escript.FunctionSpace} used to represent the solution
     
     @return : representation space of solution
     @rtype : L{escript.FunctionSpace}

     """
     return self.__column_function_space


   def getOperator(self):
     """
     provides access to the operator of the PDE

     @return : the operator of the PDE
     @rtype : L{Operator}
     """
     m=self.getSystem()[0]
     if self.isUsingLumping():
         return self.copyConstraint(1./m)
     else:
         return m

   def getRightHandSide(self):
     """
     provides access to the right hand side of the PDE

     @return : the right hand side of the PDE
     @rtype : L{escript.Data}
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
     @type u: L{escript.Data} or None
     @return : image of u
     @rtype : L{escript.Data}
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
     @type u: L{escript.Data} or None
     @return : residual of u
     @rtype : L{escript.Data}
     """
     return self.applyOperator(u)-self.getRightHandSide()

   def checkSymmetry(self,verbose=True):
      """
      test the PDE for symmetry.


     @param verbose: if equal to True or not present a report on coefficients which are breaking the symmetry is printed.
     @type verbose: C{bool} 
     @return:  True if the PDE is symmetric. 
     @rtype : C{escript.Data}

      @note: This is a very expensive operation. It should be used for degugging only! The symmetry flag is not altered.
      """
      verbose=verbose or self.debug()
      out=True
      if self.getNumSolutions()!=self.getNumEquations():
         if verbose: print "non-symmetric PDE because of different number of equations and solutions"
         out=False
      else:
         A=self.getCoefficientOfGeneralPDE("A")
         if not A.isEmpty():
            tol=util.Lsup(A)*self.TOL
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
            tol=(util.Lsup(B)+util.Lsup(C))*self.TOL/2.
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
             tol=util.Lsup(D)*self.TOL
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(D[i,k]-D[k,i])>tol:
                      if verbose: print "non-symmetric PDE because D[%d,%d]!=D[%d,%d]"%(i,k,k,i)
                      out=False

      return out

   def getSolution(self,**options):
       """
       returns the solution of the PDE. If the solution is not valid the PDE is solved. 

       @return: the solution
       @rtype: L{escript.Data}
       @param options: solver options
       @keyword verbose: 
       @keyword reordering: reordering scheme to be used during elimination
       @keyword preconditioner: preconditioner method to be used
       @keyword iter_max: maximum number of iteration steps allowed.
       @keyword drop_tolerance:
       @keyword drop_storage:
       @keyword truncation:
       @keyword restart:
       """
       if not self.__solution_isValid:
          mat,f=self.getSystem()
          if self.isUsingLumping():
             self.__solution=self.copyConstraint(f*mat)
          else:
             options[self.TOLERANCE_KEY]=self.getTolerance()
             options[self.METHOD_KEY]=self.getSolverMethod()
             options[self.SYMMETRY_KEY]=self.isSymmetric()
             self.trace("PDE is resolved.")
             self.trace("solver options: %s"%str(options))
             self.__solution=mat.solve(f,options)
          self.__solution_isValid=True
       return self.__solution

   def getFlux(self,u=None):
     """
     returns the flux J_ij for a given u

       \f[
       J_ij=A_{ijkl}u_{k,l}+B_{ijk}u_k-X_{ij}
       \f]

     @param u: argument in the flux. If u is not present or equals L{None} the current solution is used.
     @type u: L{escript.Data} or None
     @return : flux
     @rtype : L{escript.Data}

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
       @type solver: C{int}

       """
       if solver==None: solve=self.DEFAULT_METHOD
       if not solver==self.getSolverMethod():
           self.__solver_method=solver
           self.__checkMatrixType()
           self.trace("New solver is %s"%self.getSolverMethodName())

   def getSolverMethodName(self):
       """
       returns the name of the solver currently used

       @return : the name of the solver currently used.
       @rtype: C{string}
       """

       m=self.getSolverMethod()
       if m==self.DEFAULT_METHOD: return "DEFAULT_METHOD"
       elif m==self.DIRECT: return "DIRECT"
       elif m==self.CHOLEVSKY: return "CHOLEVSKY"
       elif m==self.PCG: return "PCG"
       elif m==self.CR: return "CR"
       elif m==self.CGS: return "CGS"
       elif m==self.BICGSTAB: return "BICGSTAB"
       elif m==self.SSOR: return "SSOR"
       elif m==self.GMRES: return "GMRES"
       elif m==self.PRES20: return "PRES20"
       elif m==self.LUMPING: return "LUMPING"
       return None
       

   def getSolverMethod(self):
       """
       returns the solver method
  
       @return : the solver method currently be used.
       @rtype : C{int}
       """
       return self.__solver_method

   def isUsingLumping(self):
      """
      checks if matrix lumping is used a solver method

      @return : True is lumping is currently used a solver method.
      @rtype: C{bool}
      """
      return self.getSolverMethod()==self.LUMPING

   def setTolerance(self,tol=1.e-8):
       """
       resets the tolerance for the solver method to tol where for an appropriate norm |.|

               |self.getResidual()|<tol*|self.getRightHandSide()|

       defines the stopping criterion. 

       @param tol: new tolerance for the solver. If the tol is lower then the current tolerence
                   the system will be resolved.
       @type solver: C{float}
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
    
      @return : True is a symmetric PDE is indicated, otherwise False is returned
      @rtype : C{bool}
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
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     switches off reduced order for solution and equation representation
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     sets order reduction for both solution and equation representation according to flag. 

     @param flag: if flag is True, the order reduction is switched on for both  solution and equation representation, otherwise or 
                  if flag is not present order reduction is switched off
     @type flag: C{bool}
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)


   def setReducedOrderForSolutionOn(self):
     """
     switches on reduced order for solution representation
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         self.trace("Reduced order is used to solution representation.")
         self.__column_function_space=new_fs
         self.__resetSystem()

   def setReducedOrderForSolutionOff(self):
     """
     switches off reduced order for solution representation
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         self.trace("Full order is used to interpolate solution.")
         self.__column_function_space=new_fs
         self.__resetSystem()

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     sets order for test functions according to flag

     @param flag: if flag is True, the order reduction is switched on for solution representation, otherwise or 
                  if flag is not present order reduction is switched off
     @type flag: C{bool}
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()

   def setReducedOrderForEquationOn(self):
     """
     switches on reduced order for equation representation
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         self.trace("Reduced order is used for test functions.")
         self.__row_function_space=new_fs
         self.__resetSystem()

   def setReducedOrderForEquationOff(self):
     """
     switches off reduced order for equation representation
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         self.trace("Full order is used for test functions.")
         self.__row_function_space=new_fs
         self.__resetSystem()

   def setReducedOrderForEquationTo(self,flag=False):
     """
     sets order for test functions according to flag

     @param flag: if flag is True, the order reduction is switched on for equation representation, otherwise or 
                  if flag is not present order reduction is switched off
     @type flag: C{bool}
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
     new_matrix_type=self.getDomain().getSystemMatrixTypeId(self.getSolverMethod(),self.isSymmetric())
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
       if self.__operator_isValid: self.trace("Operator has to be rebuilt.")
       self.__invalidateSolution()
       self.__operator_isValid=False

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
       self.__operator_isValid=False
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
           self.__operator.setValue(0.)
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

     @note This method is called by the assembling routine it can be overwritten 
           to map coefficients of a particular PDE to the general PDE.

     @param name: name of the coefficient requested. 
     @type name: C{string}
     @return : the value of the coefficient  name
     @rtype : L{escript.Data}
     @raise IllegalCoefficient: if name is not one of coefficients
                  "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "r" or "q".
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
     @return : True if name is the name of a coefficient of the general PDE. Otherwise False.
     @rtype : C{bool}
     
     """
     return self.__COEFFICIENTS_OF_GENEARL_PDE.has_key(name)

   def createCoefficientOfGeneralPDE(self,name):
     """
     returns a new instance of a coefficient for coefficient name of the general PDE

     @param name: name of the coefficient requested.
     @type name: C{string}
     @return : a coefficient name initialized to 0.
     @rtype : L{escript.Data}
     @raise IllegalCoefficient: if name is not one of coefficients
                  "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "r" or "q".
     """
     if self.hasCoefficientOfGeneralPDE(name):
        return escript.Data(0,self.getShapeOfCoefficientOfGeneralPDE(name),self.getFunctionSpaceForCoefficientOfGeneralPDE(name))
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def getFunctionSpaceForCoefficientOfGeneralPDE(self,name):
     """
     return the L{escript.FunctionSpace} to be used for coefficient name of the general PDE

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return : the function space to be used for coefficient name
     @rtype : L{escript.FunctionSpace}
     @raise IllegalCoefficient: if name is not one of coefficients
                  "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "r" or "q".
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
     @return : the shape of the coefficient name
     @rtype : C{tuple} of C{int}
     @raise IllegalCoefficient: if name is not one of coefficients
                  "A", "B", "C", "D", "X", "Y", "d", "y", "d_contact", "y_contact", "r" or "q".

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
     @return : the value of the coefficient name
     @rtype : L{escript.Data}
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
     @return : True if name is the name of a coefficient of the general PDE. Otherwise False.
     @rtype : C{bool}
     """
     return self.COEFFICIENTS.has_key(name)

   def createCoefficient(self, name):
     """
     create a L{escript.Data} object corresponding to coefficient name

     @return : a coefficient name initialized to 0.
     @rtype : L{escript.Data}
     @raise IllegalCoefficient: if name is not a coefficient of the PDE.
     """
     if self.hasCoefficient(name):
        return escript.Data(0.,self.getShapeOfCoefficient(name),self.getFunctionSpaceForCoefficient(name))
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def getFunctionSpaceForCoefficient(self,name):
     """
     return the L{escript.FunctionSpace} to be used for coefficient name

     @param name: name of the coefficient enquired.
     @type name: C{string}
     @return : the function space to be used for coefficient name
     @rtype : L{escript.FunctionSpace}
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
     @return : the shape of the coefficient name
     @rtype : C{tuple} of C{int}
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
     """
     if self.hasCoefficient(name):
        self.trace("Coefficient %s has been altered."%name)
        if self.COEFFICIENTS[name].isAlteringOperator(): self.__invalidateOperator()
        if self.COEFFICIENTS[name].isAlteringRightHandSide(): self.__invalidateRightHandSide()
     else:
        raise IllegalCoefficient,"illegal coefficient %s requested for general PDE."%name

   def copyConstraint(self,u):
      """
      copies the constraint into u and returns u. 
 
      @param u: a function of rank 0 is a single PDE is solved and of shape (numSolution,) for a system of PDEs
      @type u: L{escript.Data}
      @return : the input u modified by the constraints.
      @rtype : L{escript.Data}
      @warning: u is altered if it has the appropriate L{escript.FunctionSpace}

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

      @note This method is called by the assembling routine it can be overwritten 
           to map coefficients of a particular PDE to the general PDE.

      @param name: name of the coefficient requested. 
      @type name: C{string}
      @keyword A: value for coefficient A. 
      @type A: any type that can be interpreted as L{escript.Data} object on L{escript.Function}.
      @keyword B: value for coefficient B
      @type B: any type that can be interpreted as L{escript.Data} object on L{escript.Function}.
      @keyword C: value for coefficient C
      @type C: any type that can be interpreted as L{escript.Data} object on L{escript.Function}.
      @keyword D: value for coefficient D
      @type D: any type that can be interpreted as L{escript.Data} object on L{escript.Function}.
      @keyword X: value for coefficient X
      @type X: any type that can be interpreted as L{escript.Data} object on L{escript.Function}.
      @keyword Y: value for coefficient Y
      @type Y: any type that can be interpreted as L{escript.Data} object on L{escript.Function}.
      @keyword d: value for coefficient d
      @type d: any type that can be interpreted as L{escript.Data} object on L{escript.FunctionOnBoundary}.
      @keyword y: value for coefficient y
      @type y: any type that can be interpreted as L{escript.Data} object on L{escript.FunctionOnBoundary}.
      @keyword d_contact: value for coefficient d_contact
      @type d_contact: any type that can be interpreted as L{escript.Data} object on L{escript.FunctionOnContactOne}.
                       or  L{escript.FunctionOnContactZero}.
      @keyword y_contact: value for coefficient y_contact
      @type y_contact: any type that can be interpreted as L{escript.Data} object on L{escript.FunctionOnContactOne}.
                       or  L{escript.FunctionOnContactZero}.
      @keyword r: values prescribed to the solution at the locations of constraints
      @type r: any type that can be interpreted as L{escript.Data} object on L{escript.Solution} or L{escript.ReducedSolution}
               depending of reduced order is used for the solution.
      @keyword q: mask for location of constraints
      @type q: any type that can be interpreted as L{escript.Data} object on L{escript.Solution} or L{escript.ReducedSolution}
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
           self.COEFFICIENTS[i].setValue(self.getDomain(),self.getNumEquations(),self.getNumSolutions(),d)
        except IllegalCoefficientValue,m:
           raise IllegalCoefficientValue("Coefficient %s:%s"%(i,m))
        self.alteredCoefficient(i)

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
       """
       if not self.__operator_isValid or not self.__righthandside_isValid:
          if self.isUsingLumping():
              if not self.__operator_isValid:
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution(): raise TypeError,"Lumped matrix requires same order for equations and unknowns"
                 if not self.getCoefficientOfGeneralPDE("A").isEmpty(): raise Warning,"Lumped matrix does not allow coefficient A"
                 if not self.getCoefficientOfGeneralPDE("B").isEmpty(): raise Warning,"Lumped matrix does not allow coefficient B"
                 if not self.getCoefficientOfGeneralPDE("C").isEmpty(): raise Warning,"Lumped matrix does not allow coefficient C"
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
                 self.__operator_isValid=True
              if not self.__righthandside_isValid:
                 self.getDomain().addPDEToRHS(self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfGeneralPDE("X"), \
                               self.getCoefficientOfGeneralPDE("Y"),\
                               self.getCoefficientOfGeneralPDE("y"),\
                               self.getCoefficientOfGeneralPDE("y_contact"))
                 self.trace("New right hand side as been built.")
                 self.__righthandside_isValid=True
          else:
             if not self.__operator_isValid and not self.__righthandside_isValid:
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
                 self.__operator_isValid=True
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
             elif not self.__operator_isValid:
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
                 self.__operator_isValid=True
       return (self.__operator,self.__righthandside)




class AdvectivePDE(LinearPDE):
   """
   Class to handle a linear PDE dominated by advective terms:

   class to define a linear PDE of the form

   \f[
   -(A_{ijkl}u_{k,l})_{,j} -(B_{ijk}u_k)_{,j} + C_{ikl}u_{k,l} +D_{ik}u_k = - (X_{ij})_{,j} + Y_i
   \f]

   with boundary conditons:

   \f[
   n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d_{ik}u_k = - n_j*X_{ij} + y_i
   \f]

   and contact conditions

   \f[
   n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d^{contact}_{ik}[u_k] = - n_j*X_{ij} + y^{contact}_{i}
   \f]

   and constraints:

   \f[
   u_i=r_i \quad \mathrm{where} \quad q_i>0
   \f]
   """
   def __init__(self,domain,numEquations=0,numSolutions=0,xi=None,debug=False):
      LinearPDE.__init__(self,domain,numEquations,numSolutions,debug)
      if xi==None:
         self.__xi=AdvectivePDE.ELMAN_RAMAGE
      else:
         self.__xi=xi
      self.__Xi=escript.Data()

   def __calculateXi(self,peclet_factor,Z,h):
       Z_max=util.Lsup(Z)
       if Z_max>0.:
          return h*self.__xi(Z*peclet_factor)/(Z+Z_max*self.TOL)
       else:
          return 0.

   def setValue(self,**args):
       if "A" in args.keys()   or "B" in args.keys() or "C" in args.keys(): self.__Xi=escript.Data()
       LinearPDE.setValue(**args)

   def ELMAN_RAMAGE(P):
     """   """
     return (P-1.).wherePositive()*0.5*(1.-1./(P+1.e-15))
   def SIMPLIFIED_BROOK_HUGHES(P):
     """   """
     c=(P-3.).whereNegative()
     return P/6.*c+1./2.*(1.-c)

   def HALF(P):
    """ """
    return escript.Scalar(0.5,P.getFunctionSpace())

   def getXi(self):
      if self.__Xi.isEmpty():
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         A=self.getCoefficient("A")
         h=self.getDomain().getSize()
         self.__Xi=escript.Scalar(0.,self.getFunctionSpaceForCoefficient("A"))
         if not C.isEmpty() or not B.isEmpty():
            if not C.isEmpty() and not B.isEmpty():
                Z2=escript.Scalar(0,self.getFunctionSpaceForCoefficient("A"))
                if self.getNumEquations()>1:
                   if self.getNumSolutions()>1:
                      for i in range(self.getNumEquations()):
                         for k in range(self.getNumSolutions()):
                            for l in range(self.getDim()): Z2+=(C[i,k,l]-B[i,l,k])**2
                   else:
                      for i in range(self.getNumEquations()):
                         for l in range(self.getDim()): Z2+=(C[i,l]-B[i,l])**2
                else:
                   if self.getNumSolutions()>1:
                      for k in range(self.getNumSolutions()):
                         for l in range(self.getDim()): Z2+=(C[k,l]-B[l,k])**2
                   else:
                      for l in range(self.getDim()): Z2+=(C[l]-B[l])**2
                length_of_Z=util.sqrt(Z2)
            elif C.isEmpty():
              length_of_Z=util.length(B)
            else:
              length_of_Z=util.length(C)

            Z_max=util.Lsup(length_of_Z)
            if Z_max>0.:
               length_of_A=util.length(A)
               A_max=util.Lsup(length_of_A)
               if A_max>0:
                    inv_A=1./(length_of_A+A_max*self.TOL)
               else:
                    inv_A=1./self.TOL
               peclet_number=length_of_Z*h/2*inv_A
               xi=self.__xi(peclet_number)
               self.__Xi=h*xi/(length_of_Z+Z_max*self.TOL)
               print "@ preclet number = %e"%util.Lsup(peclet_number),util.Lsup(xi),util.Lsup(length_of_Z)
      return self.__Xi


   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name:
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
            Xi=self.getXi()
            if self.getNumEquations()>1:
                for i in range(self.getNumEquations()):
                   for j in range(self.getDim()):
                      for k in range(self.getNumSolutions()):
                         for l in range(self.getDim()):
                            if not C.isEmpty() and not B.isEmpty():
                               for p in range(self.getNumEquations()): Aout[i,j,k,l]+=Xi*(C[p,i,j]-B[p,j,i])*(C[p,k,l]-B[p,l,k])
                            elif C.isEmpty():
                               for p in range(self.getNumEquations()): Aout[i,j,k,l]+=Xi*B[p,j,i]*B[p,l,k]
                            else:
                               for p in range(self.getNumEquations()): Aout[i,j,k,l]+=Xi*C[p,i,j]*C[p,k,l]
            else:
                for j in range(self.getDim()):
                   for l in range(self.getDim()):
                      if not C.isEmpty() and not B.isEmpty():
                          Aout[j,l]+=Xi*(C[j]-B[j])*(C[l]-B[l])
                      elif C.isEmpty():
                          Aout[j,l]+=Xi*B[j]*B[l]
                      else:
                          Aout[j,l]+=Xi*C[j]*C[l]
         return Aout
     elif name == "B" :
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         D=self.getCoefficient("D")
         if C.isEmpty() or D.isEmpty():
            Bout=B
         else:
            Xi=self.getXi()
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
            else:
               tmp=Xi*D
               for j in range(self.getDim()): Bout[j]+=tmp*C[j]
         return Bout
     elif name == "C" :
         B=self.getCoefficient("B")
         C=self.getCoefficient("C")
         D=self.getCoefficient("D")
         if B.isEmpty() or D.isEmpty():
            Cout=C
         else:
            Xi=self.getXi()
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
            else:
               tmp=Xi*D
               for j in range(self.getDim()): Cout[j]+=tmp*B[j]
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
            Xi=self.getXi()
            if self.getNumEquations()>1:
                 for p in range(self.getNumEquations()):
                    tmp=Xi*Y[p]
                    for i in range(self.getNumEquations()):
                       for j in range(self.getDim()):
                          if not C.isEmpty() and not B.isEmpty():
                             Xout[i,j]+=tmp*(C[p,i,j]-B[p,j,i])
                          elif C.isEmpty():
                             Xout[i,j]-=tmp*B[p,j,i]
                          else:
                             Xout[i,j]+=tmp*C[p,i,j]
            else:
                 tmp=Xi*Y
                 for j in range(self.getDim()):
                    if not C.isEmpty() and not B.isEmpty():
                       Xout[j]+=tmp*(C[j]-B[j])
                    elif C.isEmpty():
                       Xout[j]-=tmp*B[j]
                    else:
                       Xout[j]+=tmp*C[j]
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
         raise SystemError,"unknown PDE coefficient %s",name


class Poisson(LinearPDE):
   """
   Class to define a Poisson equation problem:

   class to define a linear PDE of the form
   \f[
   -u_{,jj} = f
   \f]

   with boundary conditons:

   \f[
   n_j*u_{,j} = 0
   \f]

   and constraints:

   \f[
   u=0 \quad \mathrm{where} \quad q>0
   \f]
   """

   def __init__(self,domain,f=escript.Data(),q=escript.Data(),debug=False):
       LinearPDE.__init__(self,domain,1,1,debug)
       self.COEFFICIENTS={"f": PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                          "q": PDECoefficient(PDECoefficient.CONTINUOUS,(PDECoefficient.EQUATION,),PDECoefficient.BOTH)}
       self.setSymmetryOn()
       self.setValue(f,q)

   def setValue(self,f=escript.Data(),q=escript.Data()):
       """set value of PDE parameters f and q"""
       self._LinearPDE__setValue(f=f,q=q)

   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name:
     """
     if name == "A" :
         return escript.Data(numarray.identity(self.getDim()),escript.Function(self.getDomain()))
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
         raise SystemError,"unknown PDE coefficient %s",name

class LameEquation(LinearPDE):
   """
   Class to define a Lame equation problem:

   class to define a linear PDE of the form
   \f[
   -(\mu (u_{i,j}+u_{j,i}))_{,j} - \lambda u_{j,ji}} = F_i -\sigma_{ij,j}
   \f]

   with boundary conditons:

   \f[
   n_j(\mu (u_{i,j}+u_{j,i})-sigma_{ij}) + n_i\lambda u_{j,j} = f_i
   \f]

   and constraints:

   \f[
   u_i=r_i \quad \mathrm{where} \quad q_i>0
   \f]
   """

   def __init__(self,domain,f=escript.Data(),q=escript.Data(),debug=False):
       LinearPDE.__init__(self,domain,domain.getDim(),domain.getDim(),debug)
       self.COEFFICIENTS={ "lame_lambda"  : PDECoefficient(PDECoefficient.INTERIOR,(),PDECoefficient.OPERATOR),
                          "lame_mu"      : PDECoefficient(PDECoefficient.INTERIOR,(),PDECoefficient.OPERATOR),
                          "F"            : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                          "sigma"        : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.DIM),PDECoefficient.RIGHTHANDSIDE),
                          "f"            : PDECoefficient(PDECoefficient.BOUNDARY,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
                          "r"            : PDECoefficient(PDECoefficient.CONTINUOUS,(PDECoefficient.EQUATION,),PDECoefficient.BOTH),
                          "q"            : PDECoefficient(PDECoefficient.CONTINUOUS,(PDECoefficient.EQUATION,),PDECoefficient.BOTH)}
       self.setSymmetryOn()

   def setValue(self,lame_lambda=escript.Data(),lame_mu=escript.Data(),F=escript.Data(),sigma=escript.Data(),f=escript.Data(),r=escript.Data(),q=escript.Data()):
       """set value of PDE parameters"""
       self._LinearPDE__setValue(lame_lambda=lame_lambda, \
                                 lame_mu=lame_mu, \
                                 F=F, \
                                 sigma=sigma, \
                                 f=f, \
                                 r=r, \
                                 q=q)
   def getCoefficientOfGeneralPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name:
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
         raise SystemError,"unknown PDE coefficient %s",name

# $Log$
# Revision 1.11  2005/08/23 01:24:28  jgs
# Merge of development branch dev-02 back to main trunk on 2005-08-23
#
# Revision 1.10  2005/08/12 01:45:36  jgs
# erge of development branch dev-02 back to main trunk on 2005-08-12
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
#
