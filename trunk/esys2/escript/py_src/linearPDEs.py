# $Id$

## @file linearPDEs.py

"""
Functions and classes for linear PDEs
"""

import escript
import util
import numarray


def _CompTuple2(t1,t2):
   """
   Compare two tuples

   \param t1 The first tuple
   \param t2 The second tuple
   """

   dif=t1[0]+t1[1]-(t2[0]+t2[1])
   if dif<0: return 1
   elif dif>0: return -1
   else: return 0

def ELMAN_RAMAGE(P): 
    return (P-1.).wherePositive()*0.5*(1.-1./(P+1.e-15))

def SIMPLIFIED_BROOK_HUGHES(P): 
    c=(P-3.).whereNegative()
    return P/6.*c+1./2.*(1.-c)

def HALF(P): 
    return escript.Scalar(0.5,P.getFunctionSpace())

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
     
    def setValue(self,newValue):
       """
       set the value of the coefficient to new value
       """
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

    def estimateNumEquationsAndNumSolutions(self,shape=(),dim=3):
       """
       tries to estimate the number of equations in a given tensor shape for a given spatial dimension dim

       @param shape:
       @param dim:
       """
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
             s=self.buildShape(item[0],item[1],dim)
             if len(s)==0 and len(shape)==0:
                 return (1,1)
             else:
                 if s==shape: return item
       return None

    def buildShape(self,e=1,u=1,dim=3):
        """
	builds the required shape for a given number of equations e, number of unknowns u and spatial dimension dim

	@param e:
	@param u:
	@param dim:
	"""
        s=()
        for i in self.pattern:
             if i==self.EQUATION:
                if e>1: s=s+(e,)
             elif i==self.SOLUTION:
                if u>1: s=s+(u,)
             else:
                s=s+(dim,)
        return s

class LinearPDE:
   """
   Class to handle a linear PDE
   
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
   n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d_contact_{ik}[u_k] = - n_j*X_{ij} + y_contact_i
   \f]

   and constraints:

   \f[
   u_i=r_i \quad \mathrm{where} \quad q_i>0
   \f]

   """
   TOL=1.e-13
   DEFAULT_METHOD=util.DEFAULT_METHOD
   DIRECT=util.DIRECT
   CHOLEVSKY=util.CHOLEVSKY
   PCG=util.PCG
   CR=util.CR
   CGS=util.CGS
   BICGSTAB=util.BICGSTAB
   SSOR=util.SSOR
   GMRES=util.GMRES
   PRES20=util.PRES20
   ILU0=util.ILU0
   JACOBI=util.JACOBI

   def __init__(self,domain,numEquations=0,numSolutions=0):
     """
     initializes a new linear PDE.

     @param args:
     """
     # COEFFICIENTS can be overwritten by subclasses:
     self.COEFFICIENTS={
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

     # initialize attributes
     self.__debug=None
     self.__domain=domain
     self.__numEquations=numEquations
     self.__numSolutions=numSolutions
     self.cleanCoefficients()

     self.__operator=escript.Operator()
     self.__operator_isValid=False
     self.__righthandside=escript.Data()
     self.__righthandside_isValid=False
     self.__solution=escript.Data()
     self.__solution_isValid=False

     # set some default values:
     self.__homogeneous_constraint=True
     self.__row_function_space=escript.Solution(self.__domain)
     self.__column_function_space=escript.Solution(self.__domain)
     self.__tolerance=1.e-8
     self.__solver_method=util.DEFAULT_METHOD
     self.__matrix_type=self.__domain.getSystemMatrixTypeId(util.DEFAULT_METHOD,False)
     self.__sym=False
     self.__lumping=False

   def createCoefficient(self, name):
     """
     create a data object corresponding to coefficient name
     @param name:
     """
     return escript.Data(shape = getShapeOfCoefficient(name), \
                         what = getFunctionSpaceForCoefficient(name))

   def __del__(self):
     pass

   def getCoefficient(self,name):
     """
     return the value of the parameter name

     @param name:
     """
     return self.COEFFICIENTS[name].getValue()

   def getCoefficientOfPDE(self,name):
     """
     return the value of the coefficient name of the general PDE. 
     This method is called by the assembling routine it can be 
     overwritten to map coefficients of a particualr PDE to the general PDE.

     @param name:
     """
     return self.getCoefficient(name)

   def hasCoefficient(self,name):
      """
      return true if name is the name of a coefficient

      @param name:
      """
      return self.COEFFICIENTS.has_key(name)

   def getFunctionSpaceForEquation(self):
     """
     return true if the test functions should use reduced order
     """
     return self.__row_function_space

   def getFunctionSpaceForSolution(self):
     """
     return true if the interpolation of the solution should use reduced order
     """
     return self.__column_function_space

   def setValue(self,**coefficients):
      """
      sets new values to coefficients

      @param coefficients:
      """
      self.__setValue(**coefficients)
      

   def cleanCoefficients(self):
     """
     resets all coefficients to default values. 
     """
     for i in self.COEFFICIENTS.iterkeys():
         self.COEFFICIENTS[i].resetValue()

   def createNewCoefficient(self,name):
     """
     returns a new coefficient appropriate for coefficient name:
     """
     return escript.Data(0,self.getShapeOfCoefficient(name),self.getFunctionSpaceForCoefficient(name))
      

   def getShapeOfCoefficient(self,name):
     """
     return the shape of the coefficient name

     @param name:
     """
     if self.hasCoefficient(name):
        return self.COEFFICIENTS[name].buildShape(self.getNumEquations(),self.getNumSolutions(),self.getDomain().getDim())
     else:
        raise ValueError,"Solution coefficient %s requested"%name

   def getFunctionSpaceForCoefficient(self,name):
     """
     return the atoms of the coefficient name

     @param name:
     """
     if self.hasCoefficient(name):
        return self.COEFFICIENTS[name].getFunctionSpace(self.getDomain())
     else:
        raise ValueError,"Solution coefficient %s requested"%name

   def alteredCoefficient(self,name):
     """
     announce that coefficient name has been changed

     @param name:
     """
     if self.hasCoefficient(name):
        if self.COEFFICIENTS[name].isAlteringOperator(): self.__rebuildOperator()
        if self.COEFFICIENTS[name].isAlteringRightHandSide(): self.__rebuildRightHandSide()
     else:
        raise ValueError,"unknown coefficient %s requested"%name

   # ===== debug ==============================================================
   def setDebugOn(self):
       """
       """
       self.__debug=not None

   def setDebugOff(self):
       """
       """
       self.__debug=None

   def debug(self):
       """
       returns true if the PDE is in the debug mode
       """
       return self.__debug

   #===== Lumping ===========================
   def setLumpingOn(self):
      """
      indicates to use matrix lumping
      """
      if not self.isUsingLumping():
         if self.debug() : print "PDE Debug: lumping is set on"
         self.__rebuildOperator()
         self.__lumping=True

   def setLumpingOff(self):
      """
      switches off matrix lumping
      """
      if self.isUsingLumping():
         if self.debug() : print "PDE Debug: lumping is set off"
         self.__rebuildOperator()
         self.__lumping=False

   def setLumping(self,flag=False):
      """
      set the matrix lumping flag to flag
      """
      if flag:
         self.setLumpingOn()
      else:
         self.setLumpingOff()

   def isUsingLumping(self):
      """
      
      """
      return self.__lumping

   #============ method business =========================================================
   def setSolverMethod(self,solver=util.DEFAULT_METHOD):
       """
       sets a new solver
       """
       if not solver==self.getSolverMethod():
           self.__solver_method=solver
           if self.debug() : print "PDE Debug: New solver is %s"%solver
           self.__checkMatrixType()

   def getSolverMethod(self):
       """
       returns the solver method
       """
       return self.__solver_method

   #============ tolerance business =========================================================
   def setTolerance(self,tol=1.e-8):
       """
       resets the tolerance to tol.
       """
       if not tol>0:
           raise ValueException,"Tolerance as to be positive"
       if tol<self.getTolerance(): self.__rebuildSolution()
       if self.debug() : print "PDE Debug: New tolerance %e",tol
       self.__tolerance=tol
       return
   def getTolerance(self):
       """
       returns the tolerance set for the solution
       """
       return self.__tolerance

   #===== symmetry  flag ==========================
   def isSymmetric(self):
      """
      returns true is the operator is considered to be symmetric
      """
      return self.__sym

   def setSymmetryOn(self):
      """
      sets the symmetry flag to true
      """
      if not self.isSymmetric():
         if self.debug() : print "PDE Debug: Operator is set to be symmetric"
         self.__sym=True
         self.__checkMatrixType()

   def setSymmetryOff(self):
      """
      sets the symmetry flag to false
      """
      if self.isSymmetric():
         if self.debug() : print "PDE Debug: Operator is set to be unsymmetric"
         self.__sym=False
         self.__checkMatrixType()

   def setSymmetryTo(self,flag=False):
     """
     sets the symmetry flag to flag

     @param flag:
     """
     if flag:
        self.setSymmetryOn()
     else:
        self.setSymmetryOff()

   #===== order reduction ==========================
   def setReducedOrderOn(self):
     """
     switches to on reduced order
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     switches to full order 
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     sets order according to flag

     @param flag:
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)
                                                                                                                                                           

   #===== order reduction solution ==========================
   def setReducedOrderForSolutionOn(self):
     """
     switches to reduced order to interpolate solution
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         if self.debug() : print "PDE Debug: Reduced order is used to interpolate solution."
         self.__column_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForSolutionOff(self):
     """
     switches to full order to interpolate solution
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         if self.debug() : print "PDE Debug: Full order is used to interpolate solution."
         self.__column_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     sets order for test functions according to flag

     @param flag:
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()
                                                                                                                                                           
   #===== order reduction equation ==========================
   def setReducedOrderForEquationOn(self):
     """
     switches to reduced order for test functions
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         if self.debug() : print "PDE Debug: Reduced order is used for test functions."
         self.__row_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForEquationOff(self):
     """
     switches to full order for test functions
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         if self.debug() : print "PDE Debug: Full order is used for test functions."
         self.__row_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForEquationTo(self,flag=False):
     """
     sets order for test functions according to flag

     @param flag:
     """
     if flag:
        self.setReducedOrderForEquationOn()
     else:
        self.setReducedOrderForEquationOff()
                                                                                                                                                           
   # ==== initialization =====================================================================
   def __getNewOperator(self):
       """
       """
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           self.__matrix_type)

   def __makeFreshRightHandSide(self):
       """
       """
       if self.debug() : print "PDE Debug: New right hand side allocated"
       if self.getNumEquations()>1:
           self.__righthandside=escript.Data(0.,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),True)
       else:
           self.__righthandside=escript.Data(0.,(),self.getFunctionSpaceForEquation(),True)
       return self.__righthandside

   def __getNewSolution(self):
       """
       """
       if self.debug() : print "PDE Debug: New right hand side allocated"
       if self.getNumSolutions()>1:
           return escript.Data(0.,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)
       else:
           return escript.Data(0.,(),self.getFunctionSpaceForSolution(),True)

   def __makeFreshOperator(self):
       """
       """
       if self.__operator.isEmpty():
           self.__operator=self.__getNewOperator()
           if self.debug() : print "PDE Debug: New operator allocated"
       else:
           self.__operator.setValue(0.)
           self.__operator.resetSolver()
           if self.debug() : print "PDE Debug: Operator reset to zero"
       return self.__operator

   #============ some serivice functions  =====================================================
   def getDomain(self):
     """
     returns the domain of the PDE
     """
     return self.__domain

   def getDim(self):
     """
     returns the spatial dimension of the PDE
     """
     return self.getDomain().getDim()

   def getNumEquations(self):
     """
     returns the number of equations
     """
     if self.__numEquations>0:
         return self.__numEquations
     else:
         raise ValueError,"Number of equations is undefined. Please specify argument numEquations."

   def getNumSolutions(self):
     """
     returns the number of unknowns
     """
     if self.__numSolutions>0:
        return self.__numSolutions
     else:
        raise ValueError,"Number of solution is undefined. Please specify argument numSolutions."


   def checkSymmetry(self,verbose=True):
      """
      returns if the Operator is symmetric. This is a very expensive operation!!! The symmetry flag is not altered.
      """
      verbose=verbose or self.debug()
      out=True
      if self.getNumSolutions()!=self.getNumEquations():
         if verbose: print "non-symmetric PDE because of different number of equations and solutions"
         out=False
      else:
         A=self.getCoefficientOfPDE("A")
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
         B=self.getCoefficientOfPDE("B")
         C=self.getCoefficientOfPDE("C")
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
           D=self.getCoefficientOfPDE("D")
           if not D.isEmpty():
             tol=util.Lsup(D)*self.TOL
             for i in range(self.getNumEquations()):
                for k in range(self.getNumSolutions()):
                  if util.Lsup(D[i,k]-D[k,i])>tol:
                      if verbose: print "non-symmetric PDE because D[%d,%d]!=D[%d,%d]"%(i,k,k,i)
                      out=False
            
      return out

   def getFlux(self,u):
       """
       returns the flux J_ij for a given u

       \f[
       J_ij=A_{ijkl}u_{k,l}+B_{ijk}u_k-X_{ij}
       \f]

       @param u: argument of the operator 
       """
       raise SystemError,"getFlux is not implemented yet"
       return None

   def applyOperator(self,u):
       """
       applies the operator of the PDE to a given solution u in weak from

       @param u: argument of the operator 
       """
       return self.getOperator()*escript.Data(u,self.getFunctionSpaceForSolution())
                                                                                                                                                           
   def getResidual(self,u):
       """
       return the residual of u in the weak from

       @param u:
       """
       return self.applyOperator(u)-self.getRightHandSide()

   def __setValue(self,**coefficients):
      """
      sets new values to coefficient

      @param coefficients:
      """
      # check if the coefficients are  legal:
      for i in coefficients.iterkeys():
         if not self.hasCoefficient(i):
            raise ValueError,"Attempt to set unknown coefficient %s"%i
      # if the number of unknowns or equations is still unknown we try to estimate them:
      if self.__numEquations<1 or self.__numSolutions<1:
         for i,d in coefficients.iteritems():
            if hasattr(d,"shape"):
                s=d.shape
            elif hasattr(d,"getShape"):
                s=d.getShape()
            else:
                s=numarray.array(d).shape
            if s!=None:
                # get number of equations and number of unknowns:
                res=self.COEFFICIENTS[i].estimateNumEquationsAndNumSolutions(s,self.getDim())
                if res==None:
                    raise ValueError,"Illegal shape %s of coefficient %s"%(s,i)
                else:
                    if self.__numEquations<1: self.__numEquations=res[0]
                    if self.__numSolutions<1: self.__numSolutions=res[1]
      if self.__numEquations<1: raise ValueError,"unidententified number of equations"
      if self.__numSolutions<1: raise ValueError,"unidententified number of solutions"
      # now we check the shape of the coefficient if numEquations and numSolutions are set:
      for i,d in coefficients.iteritems():
        if d==None:
             d2=escript.Data()
        elif isinstance(d,escript.Data):
             if d.isEmpty():
                d2=d
             else:
                d2=escript.Data(d,self.getFunctionSpaceForCoefficient(i))
        else:
              d2=escript.Data(d,self.getFunctionSpaceForCoefficient(i))
        if not d2.isEmpty():
           if not self.getShapeOfCoefficient(i)==d2.getShape():
               raise ValueError,"Expected shape for coefficient %s is %s but actual shape is %s."%(i,self.getShapeOfCoefficient(i),d2.getShape())
        # overwrite new values:
        if self.debug(): print "PDE Debug: Coefficient %s has been altered."%i
        self.COEFFICIENTS[i].setValue(d2)
        self.alteredCoefficient(i)
      
      # reset the HomogeneousConstraintFlag:
      self.__setHomogeneousConstraintFlag()
      if len(coefficients)>0 and not self.isUsingLumping() and not self.__homogeneous_constraint: self.__rebuildSystem()

   def __setHomogeneousConstraintFlag(self): 
      """
      checks if the constraints are homogeneous and sets self.__homogeneous_constraint accordingly.
      """
      self.__homogeneous_constraint=True
      q=self.getCoefficientOfPDE("q")
      r=self.getCoefficientOfPDE("r")
      if not q.isEmpty() and not r.isEmpty():
         if (q*r).Lsup()>=1.e-13*r.Lsup(): self.__homogeneous_constraint=False
      if self.debug():
           if self.__homogeneous_constraint:
               print "PDE Debug: Constraints are homogeneous."
           else:
               print "PDE Debug: Constraints are inhomogeneous."
 

   # ==== rebuild switches =====================================================================
   def __rebuildSolution(self,deep=False):
       """
       indicates the PDE has to be reolved if the solution is requested
       """
       if self.__solution_isValid and self.debug() : print "PDE Debug: PDE has to be resolved."
       self.__solution_isValid=False
       if deep: self.__solution=escript.Data()


   def __rebuildOperator(self,deep=False):
       """
       indicates the operator has to be rebuilt next time it is used
       """
       if self.__operator_isValid and self.debug() : print "PDE Debug: Operator has to be rebuilt."
       self.__rebuildSolution(deep)
       self.__operator_isValid=False
       if deep: self.__operator=escript.Operator()

   def __rebuildRightHandSide(self,deep=False):
       """
       indicates the right hand side has to be rebuild next time it is used
       """
       if self.__righthandside_isValid and self.debug() : print "PDE Debug: Right hand side has to be rebuilt."
       self.__rebuildSolution(deep)
       self.__righthandside_isValid=False
       if deep: self.__righthandside=escript.Data()

   def __rebuildSystem(self,deep=False):
       """
       annonced that all coefficient name has been changed
       """
       self.__rebuildSolution(deep)
       self.__rebuildOperator(deep)
       self.__rebuildRightHandSide(deep)
   
   def __checkMatrixType(self):
     """
     reassess the matrix type and, if needed, initiates an operator rebuild
     """
     new_matrix_type=self.getDomain().getSystemMatrixTypeId(self.getSolverMethod(),self.isSymmetric())
     if not new_matrix_type==self.__matrix_type:
         if self.debug() : print "PDE Debug: Matrix type is now %d."%new_matrix_type
         self.__matrix_type=new_matrix_type
         self.__rebuildOperator(deep=True)

   #============ assembling =======================================================
   def __copyConstraint(self):
      """
      copies the constrint condition into u
      """
      if not self.__righthandside.isEmpty(): 
         q=self.getCoefficientOfPDE("q")
         r=self.getCoefficientOfPDE("r")
         if not q.isEmpty():
             if r.isEmpty():
                r2=escript.Data(0,self.__righthandside.getShape(),self.__righthandside.getFunctionSpace())
             else:
                r2=escript.Data(r,self.__righthandside.getFunctionSpace())
             self.__righthandside.copyWithMask(r2,escript.Data(q,self.__righthandside.getFunctionSpace()))

   def __applyConstraint(self):
       """
       applies the constraints defined by q and r to the system
       """
       q=self.getCoefficientOfPDE("q")
       r=self.getCoefficientOfPDE("r")
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
          if self.isUsingLumping():
             self.__operator.copyWithMask(escript.Data(1,q.getShape(),self.getFunctionSpaceForEquation()),row_q)
          else:
             if not self.__righthandside.isEmpty(): self.__righthandside-=self.__operator*u
             self.__operator.nullifyRowsAndCols(row_q,col_q,1.)

   def getSystem(self):
       """
       return the operator and right hand side of the PDE
       """
       if not self.__operator_isValid or not self.__righthandside_isValid:
          if self.isUsingLumping():
              if not self.__operator_isValid:
                 if not self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                       raise TypeError,"Lumped matrix requires same order for equations and unknowns"
                 if not self.getCoefficientOfPDE("A").isEmpty():
                          raise Warning,"Lumped matrix does not allow coefficient A"
                 if not self.getCoefficientOfPDE("B").isEmpty():
                          raise Warning,"Lumped matrix does not allow coefficient B"
                 if not self.getCoefficientOfPDE("C").isEmpty():
                          raise Warning,"Lumped matrix does not allow coefficient C"
                 if self.debug() : print "PDE Debug: New lumped operator is built."
                 mat=self.__getNewOperator()
                 self.getDomain().addPDEToSystem(mat,escript.Data(), \
                           self.getCoefficientOfPDE("A"), \
                           self.getCoefficientOfPDE("B"), \
                           self.getCoefficientOfPDE("C"), \
                           self.getCoefficientOfPDE("D"), \
                           escript.Data(), \
                           escript.Data(), \
                           self.getCoefficientOfPDE("d"), \
                           escript.Data(),\
                           self.getCoefficientOfPDE("d_contact"), \
                           escript.Data())
                 self.__operator=mat*escript.Data(1,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)
                 self.__applyConstraint()
                 self.__operator_isValid=True
              if not self.__righthandside_isValid:
                 if self.debug() : print "PDE Debug: New right hand side is built."
                 self.getDomain().addPDEToRHS(self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfPDE("X"), \
                               self.getCoefficientOfPDE("Y"),\
                               self.getCoefficientOfPDE("y"),\
                               self.getCoefficientOfPDE("y_contact"))
                 self.__copyConstraint()
                 self.__righthandside_isValid=True
          else:
             if not self.__operator_isValid and not self.__righthandside_isValid:
                 if self.debug() : print "PDE Debug: New system is built."
                 self.getDomain().addPDEToSystem(self.__makeFreshOperator(),self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfPDE("A"), \
                               self.getCoefficientOfPDE("B"), \
                               self.getCoefficientOfPDE("C"), \
                               self.getCoefficientOfPDE("D"), \
                               self.getCoefficientOfPDE("X"), \
                               self.getCoefficientOfPDE("Y"), \
                               self.getCoefficientOfPDE("d"), \
                               self.getCoefficientOfPDE("y"), \
                               self.getCoefficientOfPDE("d_contact"), \
                               self.getCoefficientOfPDE("y_contact"))
                 self.__applyConstraint()
                 self.__copyConstraint()
                 self.__operator_isValid=True
                 self.__righthandside_isValid=True
             elif not self.__righthandside_isValid:
                 if self.debug() : print "PDE Debug: New right hand side is built."
                 self.getDomain().addPDEToRHS(self.__makeFreshRightHandSide(), \
                               self.getCoefficientOfPDE("X"), \
                               self.getCoefficientOfPDE("Y"),\
                               self.getCoefficientOfPDE("y"),\
                               self.getCoefficientOfPDE("y_contact"))
                 self.__copyConstraint()
                 self.__righthandside_isValid=True
             elif not self.__operator_isValid:
                 if self.debug() : print "PDE Debug: New operator is built."
                 self.getDomain().addPDEToSystem(self.__makeFreshOperator(),escript.Data(), \
                            self.getCoefficientOfPDE("A"), \
                            self.getCoefficientOfPDE("B"), \
                            self.getCoefficientOfPDE("C"), \
                            self.getCoefficientOfPDE("D"), \
                            escript.Data(), \
                            escript.Data(), \
                            self.getCoefficientOfPDE("d"), \
                            escript.Data(),\
                            self.getCoefficientOfPDE("d_contact"), \
                            escript.Data())
                 self.__applyConstraint()
                 self.__operator_isValid=True
       return (self.__operator,self.__righthandside)
   def getOperator(self):
       """
       returns the operator of the PDE
       """
       return self.getSystem()[0]

   def getRightHandSide(self):
       """
       returns the right hand side of the PDE
       """
       return self.getSystem()[1]

   def solve(self,**options):
      """
      solve the PDE

      @param options:
      """
      mat,f=self.getSystem()
      if self.isUsingLumping():
         out=f/mat
      else:
         options[util.TOLERANCE_KEY]=self.getTolerance()
         options[util.METHOD_KEY]=self.getSolverMethod()
         options[util.SYMMETRY_KEY]=self.isSymmetric()
         if self.debug() : print "PDE Debug: solver options: ",options
         out=mat.solve(f,options)
      return out

   def getSolution(self,**options):
       """
       returns the solution of the PDE

       @param options:
       """
       if not self.__solution_isValid:
           if self.debug() : print "PDE Debug: PDE is resolved."
           self.__solution=self.solve(**options)
           self.__solution_isValid=True
       return self.__solution



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
   def __init__(self,domain,numEquations=0,numSolutions=0,xi=ELMAN_RAMAGE):
      LinearPDE.__init__(self,domain,numEquations,numSolutions)
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
       self._LinearPDE__setValue(**args)
           
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
      

   def getCoefficientOfPDE(self,name):
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

   def __init__(self,domain,f=escript.Data(),q=escript.Data()):
       LinearPDE.__init__(self,domain,1,1)
       self.COEFFICIENTS={
       "f"         : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,),PDECoefficient.RIGHTHANDSIDE),
       "q"         : PDECoefficient(PDECoefficient.CONTINUOUS,(PDECoefficient.EQUATION,),PDECoefficient.BOTH)}
       self.setSymmetryOn()
       self.setValue(f,q)

   def setValue(self,f=escript.Data(),q=escript.Data()):
       """set value of PDE parameters f and q"""
       self._LinearPDE__setValue(f=f,q=q)

   def getCoefficientOfPDE(self,name):
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
   -(\lambda (u_{i,j}+u_{j,i}))_{,j} - \mu u_{j,ji}} = F_i -\sigma_{ij,j}
   \f]

   with boundary conditons:

   \f[
   n_j(\lambda(u_{i,j}+u_{j,i})-sigma_{ij}) + n_i\mu u_{j,j} = f_i
   \f]

   and constraints:

   \f[
   u_i=r_i \quad \mathrm{where} \quad q_i>0
   \f]
   """

   def __init__(self,domain,f=escript.Data(),q=escript.Data()):
       LinearPDE.__init__(self,domain,domain.getDim(),domain.getDim())
       self.COEFFICIENTS={
       "lame_lambda"  : PDECoefficient(PDECoefficient.INTERIOR,(),PDECoefficient.OPERATOR),
       "lame_mu"      : PDECoefficient(PDECoefficient.INTERIOR,(),PDECoefficient.OPERATOR),
       "F"            : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.DIM),PDECoefficient.RIGHTHANDSIDE),
       "sigma"        : PDECoefficient(PDECoefficient.INTERIOR,(PDECoefficient.EQUATION,PDECoefficient.EQUATION),PDECoefficient.RIGHTHANDSIDE),
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
   def getCoefficientOfPDE(self,name):
     """
     return the value of the coefficient name of the general PDE

     @param name:
     """
     if name == "A" : 
         A =self.createNewCoefficient("A")
         for i in range(self.getDim()):
           for j in range(self.getDim()):
             out[i,i,j,j] += self.getCoefficient("lame_mu")
             out[i,j,j,i] += self.getCoefficient("lame_lambda")
             out[i,j,i,j] += self.getCoefficient("lame_lambda")
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
# Revision 1.9  2005/07/25 05:28:13  jgs
# Merge of development branch back to main trunk on 2005-07-25
#
# Revision 1.8  2005/06/09 05:37:59  jgs
# Merge of development branch back to main trunk on 2005-06-09
#
# Revision 1.7  2005/05/06 04:26:10  jgs
# Merge of development branch back to main trunk on 2005-05-06
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
