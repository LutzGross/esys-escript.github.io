# $Id$

## @file linearPDEs.py

"""
@brief Functions and classes for linear PDEs
"""

import escript
import util
import numarray

def identifyDomain(domain=None,data={}):
     """
     @brief Return the Domain which is equal to the input domain (if not None)
     and is the domain of all Data objects in the dictionary data.
     An exception is raised if this is not possible

     @param domain
     @param data
     """
     # get the domain used by any Data object in the list data:
     data_domain=None
     for d in data.itervalues():
          if isinstance(d,escript.Data):
             if not d.isEmpty(): data_domain=d.getDomain()
     # check if domain and data_domain are identical?
     if domain == None:
         if data_domain == None:
              raise ValueError,"Undefined PDE domain. Specify a domain or use a Data class object as coefficient"
     else:
         if data_domain == None:
              data_domain=domain
         else:
           if not data_domain == domain:
                 raise ValueError,"Domain of coefficients doesnot match specified domain"
     # now we check if all Data class object coefficients are defined on data_domain:
     for i,d in data.iteritems():
         if isinstance(d,escript.Data):
            if not d.isEmpty(): 
               if not data_domain==d.getDomain():
                 raise ValueError,"Illegal domain for coefficient %s."%i
     # done:
     return data_domain

def identifyNumEquationsAndSolutions(dim,coef={}):
     # get number of equations and number of unknowns:
     numEquations=0
     numSolutions=0
     for i in coef.iterkeys():
        if not coef[i].isEmpty():
           res=_PDECoefficientTypes[i].estimateNumEquationsAndNumSolutions(coef[i].getShape(),dim)
           if res==None:
               raise ValueError,"Illegal shape %s of coefficient %s"%(coef[i].getShape().__str__(),i)
           else:
               numEquations=max(numEquations,res[0])
               numSolutions=max(numSolutions,res[1])
     return numEquations,numSolutions


def _CompTuple2(t1,t2):
   """
   @brief

   @param t1
   @param t2
   """
   dif=t1[0]+t1[1]-(t2[0]+t2[1])
   if dif<0: return 1
   elif dif>0: return -1
   else: return 0

class PDECoefficientType:
    """
    @brief
    """
    # identifier for location of Data objects defining coefficients
    INTERIOR=0
    BOUNDARY=1
    CONTACT=2
    CONTINUOUS=3
    # identifier in the pattern of coefficients:
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
       @brief Initialise a PDE Coefficient type
       """
       self.what=where
       self.pattern=pattern
       self.altering=altering

    def getFunctionSpace(self,domain):
       """
       @brief defines the FunctionSpace of the coefficient on the domain

       @param domain
       """
       if self.what==self.INTERIOR: return escript.Function(domain)
       elif self.what==self.BOUNDARY: return escript.FunctionOnBoundary(domain)
       elif self.what==self.CONTACT: return escript.FunctionOnContactZero(domain)
       elif self.what==self.CONTINUOUS: return escript.ContinuousFunction(domain)

    def isAlteringOperator(self):
        """
	@brief return true if the operator of the PDE is changed when the coefficient is changed
	"""
        if self.altering==self.OPERATOR or self.altering==self.BOTH:
            return not None
        else:
            return None

    def isAlteringRightHandSide(self):
        """
	@brief return true if the right hand side of the PDE is changed when the coefficient is changed
	"""
        if self.altering==self.RIGHTHANDSIDE or self.altering==self.BOTH:
            return not None
        else:
            return None

    def estimateNumEquationsAndNumSolutions(self,shape=(),dim=3):
       """
       @brief tries to estimate the number of equations in a given tensor shape for a given spatial dimension dim

       @param shape
       @param dim
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
	@brief builds the required shape for a given number of equations e, number of unknowns u and spatial dimension dim

	@param e
	@param u
	@param dim
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

_PDECoefficientTypes={
"A"         : PDECoefficientType(PDECoefficientType.INTERIOR,(PDECoefficientType.EQUATION,PDECoefficientType.DIM,PDECoefficientType.SOLUTION,PDECoefficientType.DIM),PDECoefficientType.OPERATOR),
"B"         : PDECoefficientType(PDECoefficientType.INTERIOR,(PDECoefficientType.EQUATION,PDECoefficientType.DIM,PDECoefficientType.SOLUTION),PDECoefficientType.OPERATOR),
"C"         : PDECoefficientType(PDECoefficientType.INTERIOR,(PDECoefficientType.EQUATION,PDECoefficientType.SOLUTION,PDECoefficientType.DIM),PDECoefficientType.OPERATOR),
"D"         : PDECoefficientType(PDECoefficientType.INTERIOR,(PDECoefficientType.EQUATION,PDECoefficientType.SOLUTION),PDECoefficientType.OPERATOR),
"X"         : PDECoefficientType(PDECoefficientType.INTERIOR,(PDECoefficientType.EQUATION,PDECoefficientType.DIM),PDECoefficientType.RIGHTHANDSIDE),
"Y"         : PDECoefficientType(PDECoefficientType.INTERIOR,(PDECoefficientType.EQUATION,),PDECoefficientType.RIGHTHANDSIDE),
"d"         : PDECoefficientType(PDECoefficientType.BOUNDARY,(PDECoefficientType.EQUATION,PDECoefficientType.SOLUTION),PDECoefficientType.OPERATOR),
"y"         : PDECoefficientType(PDECoefficientType.BOUNDARY,(PDECoefficientType.EQUATION,),PDECoefficientType.RIGHTHANDSIDE),
"d_contact" : PDECoefficientType(PDECoefficientType.CONTACT,(PDECoefficientType.EQUATION,PDECoefficientType.SOLUTION),PDECoefficientType.OPERATOR),
"y_contact" : PDECoefficientType(PDECoefficientType.CONTACT,(PDECoefficientType.EQUATION,),PDECoefficientType.RIGHTHANDSIDE),
"r"         : PDECoefficientType(PDECoefficientType.CONTINUOUS,(PDECoefficientType.EQUATION,),PDECoefficientType.RIGHTHANDSIDE),
"q"         : PDECoefficientType(PDECoefficientType.CONTINUOUS,(PDECoefficientType.SOLUTION,),PDECoefficientType.BOTH),
}

class LinearPDE:
   """
   @brief Class to define a linear PDE
   
   class to define a linear PDE of the form

     -(A_{ijkl}u_{k,l})_{,j} -(B_{ijk}u_k)_{,j} + C_{ikl}u_{k,l} +D_{ik}u_k = - (X_{ij})_{,j} + Y_i

     with boundary conditons:

        n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d_{ik}u_k = - n_j*X_{ij} + y_i

    and contact conditions

        n_j*(A_{ijkl}u_{k,l}+B_{ijk}u_k)_{,j} + d_contact_{ik}[u_k] = - n_j*X_{ij} + y_contact_i

    and constraints:

         u_i=r_i where q_i>0

   """
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

   def __init__(self,**args):
     """
     @brief initializes a new linear PDE.

     @param args
     """

     # initialize attributes
     self.__debug=None
     self.__domain=None
     self.__numEquations=0
     self.__numSolutions=0
     self.cleanCoefficients()

     self.__operator=escript.Operator()
     self.__operator_isValid=False
     self.__righthandside=escript.Data()
     self.__righthandside_isValid=False
     self.__solution=escript.Data()
     self.__solution_isValid=False

     # check the arguments
     coef={}
     for arg in args.iterkeys():
          if arg=="domain":
              self.__domain=args[arg]
          elif arg=="numEquations":
              self.__numEquations=args[arg]
          elif arg=="numSolutions":
              self.__numSolutions=args[arg]
          elif _PDECoefficientTypes.has_key(arg):
              coef[arg]=args[arg]
          else:
              raise ValueError,"Illegal argument %s"%arg

     # get the domain of the PDE
     self.__domain=identifyDomain(self.__domain,coef)

     # set some default values:
     self.__homogeneous_constraint=True
     self.__row_function_space=escript.Solution(self.__domain)
     self.__column_function_space=escript.Solution(self.__domain)
     self.__tolerance=1.e-8
     self.__solver_method=util.DEFAULT_METHOD
     self.__matrix_type=self.__domain.getSystemMatrixTypeId(util.DEFAULT_METHOD,False)
     self.__sym=False
     self.__lumping=False
     self.__numEquations=0
     self.__numSolutions=0
     # now we can set the ceofficients:
     self._setCoefficient(**coef)

   def getCoefficient(self,name):
     """
     @brief return the value of the coefficient name

     @param name
     """
     return self.__coefficient[name]

   def setValue(self,**coefficients):
      """
      @brief sets new values to coefficients

      @param coefficients
      """
      self._setCoefficient(**coefficients)
      

   def _setCoefficient(self,**coefficients):
      """
      @brief sets new values to coefficients

      @param coefficients
      """
      
      # get the dictionary of the coefficinets been altered:
      alteredCoefficients={}
      for i,d in coefficients.iteritems():
         if self.hasCoefficient(i):
            if d == None:
                alteredCoefficients[i]=escript.Data()
            elif isinstance(d,escript.Data):
                if d.isEmpty():
                  alteredCoefficients[i]=escript.Data()
                else:
                  alteredCoefficients[i]=escript.Data(d,self.getFunctionSpaceOfCoefficient(i))
            else:
                if self.__numEquations>0 and  self.__numSolutions>0:
                   alteredCoefficients[i]=escript.Data(d,self.getShapeOfCoefficient(i),self.getFunctionSpaceOfCoefficient(i))
                else:
                   alteredCoefficients[i]=escript.Data(d,self.getFunctionSpaceOfCoefficient(i))
         else:
            raise ValueError,"Attempt to set undefined coefficient %s"%i
      # if numEquations and numSolutions is undefined we try identify their values based on the coefficients:
      if self.__numEquations<1 or self.__numSolutions<1:
            numEquations,numSolutions=identifyNumEquationsAndSolutions(self.getDomain().getDim(),alteredCoefficients)
            if self.__numEquations<1 and numEquations>0: self.__numEquations=numEquations
            if self.__numSolutions<1 and numSolutions>0: self.__numSolutions=numSolutions
            if self.debug() and self.__numEquations>0: print "PDE Debug: identified number of equations is ",self.__numEquations
            if self.debug() and self.__numSolutions>0: print "PDE Debug: identified number of solutions is ",self.__numSolutions

      # now we check the shape of the coefficient if numEquations and numSolutions are set:
      if  self.__numEquations>0 and  self.__numSolutions>0:
         for i in self.__coefficient.iterkeys():
             if alteredCoefficients.has_key(i) and not alteredCoefficients[i].isEmpty():
                 if not self.getShapeOfCoefficient(i)==alteredCoefficients[i].getShape():
                    raise ValueError,"Expected shape for coefficient %s is %s but actual shape is %s."%(i,self.getShapeOfCoefficient(i),alteredCoefficients[i].getShape())
             else:
                 if not self.__coefficient[i].isEmpty():
                    if not self.getShapeOfCoefficient(i)==self.__coefficient[i].getShape():
                       raise ValueError,"Expected shape for coefficient %s is %s but actual shape is %s."%(i,self.getShapeOfCoefficient(i),self.__coefficient[i].getShape())
      # overwrite new values:
      for i,d in alteredCoefficients.iteritems():
         if self.debug(): print "PDE Debug: Coefficient %s has been altered."%i
         self.__coefficient[i]=d
         self.alteredCoefficient(i)

      # reset the HomogeneousConstraintFlag:
      self.__setHomogeneousConstraintFlag()
      if not "q" in alteredCoefficients and not self.__homogeneous_constraint: self.__rebuildSystem()

   def cleanCoefficients(self):
     """
     @brief resets all coefficients to default values. 
     """
     self.__coefficient={}
     for i in _PDECoefficientTypes.iterkeys():
         self.__coefficient[i]=escript.Data()

   def getShapeOfCoefficient(self,name):
     """
     @brief return the shape of the coefficient name

     @param name
     """
     if self.hasCoefficient(name):
        return _PDECoefficientTypes[name].buildShape(self.getNumEquations(),self.getNumSolutions(),self.getDomain().getDim())
     else:
        raise ValueError,"Solution coefficient %s requested"%name

   def getFunctionSpaceOfCoefficient(self,name):
     """
     @brief return the atoms of the coefficient name

     @param name
     """
     if self.hasCoefficient(name):
        return _PDECoefficientTypes[name].getFunctionSpace(self.getDomain())
     else:
        raise ValueError,"Solution coefficient %s requested"%name

   def alteredCoefficient(self,name):
     """
     @brief annonced that coefficient name has been changed

     @param name
     """
     if self.hasCoefficient(name):
        if _PDECoefficientTypes[name].isAlteringOperator(): self.__rebuildOperator()
        if _PDECoefficientTypes[name].isAlteringRightHandSide(): self.__rebuildRightHandSide()
     else:
        raise ValueError,"Solution coefficient %s requested"%name

   def __setHomogeneousConstraintFlag(self): 
      """
      @brief checks if the constraints are homogeneous and sets self.__homogeneous_constraint accordingly.
      """
      self.__homogeneous_constraint=True
      q=self.getCoefficient("q")
      r=self.getCoefficient("r")
      if not q.isEmpty() and not r.isEmpty():
         print (q*r).Lsup(), 1.e-13*r.Lsup()
         if (q*r).Lsup()>=1.e-13*r.Lsup(): self.__homogeneous_constraint=False
      if self.debug():
           if self.__homogeneous_constraint:
               print "PDE Debug: Constraints are homogeneous."
           else:
               print "PDE Debug: Constraints are inhomogeneous."
 

   def hasCoefficient(self,name):
      """
      @brief return true if name is the name of a coefficient

      @param name
      """
      return self.__coefficient.has_key(name)

   def getFunctionSpaceForEquation(self):
     """
     @brief return true if the test functions should use reduced order
     """
     return self.__row_function_space

   def getFunctionSpaceForSolution(self):
     """
     @brief return true if the interpolation of the solution should use reduced order
     """
     return self.__column_function_space

   # ===== debug ==============================================================
   def setDebugOn(self):
       """
       @brief
       """
       self.__debug=not None

   def setDebugOff(self):
       """
       @brief
       """
       self.__debug=None

   def debug(self):
       """
       @brief returns true if the PDE is in the debug mode
       """
       return self.__debug

   #===== Lumping ===========================
   def setLumpingOn(self):
      """
      @brief indicates to use matrix lumping
      """
      if not self.isUsingLumping():
         raise SystemError,"Lumping is not working yet! Talk to the experts"
         if self.debug() : print "PDE Debug: lumping is set on"
         self.__rebuildOperator()
         self.__lumping=True

   def setLumpingOff(self):
      """
      @brief switches off matrix lumping
      """
      if self.isUsingLumping():
         if self.debug() : print "PDE Debug: lumping is set off"
         self.__rebuildOperator()
         self.__lumping=False

   def setLumping(self,flag=False):
      """
      @brief set the matrix lumping flag to flag
      """
      if flag:
         self.setLumpingOn()
      else:
         self.setLumpingOff()

   def isUsingLumping(self):
      """
      @brief 
      """
      return self.__lumping

   #============ method business =========================================================
   def setSolverMethod(self,solver=util.DEFAULT_METHOD):
       """
       @brief sets a new solver
       """
       if not solver==self.getSolverMethod():
           self.__solver_method=solver
           if self.debug() : print "PDE Debug: New solver is %s"%solver
           self.__checkMatrixType()

   def getSolverMethod(self):
       """
       @brief returns the solver method
       """
       return self.__solver_method

   #============ tolerance business =========================================================
   def setTolerance(self,tol=1.e-8):
       """
       @brief resets the tolerance to tol.
       """
       if not tol>0:
           raise ValueException,"Tolerance as to be positive"
       if tol<self.getTolerance(): self.__rebuildSolution()
       if self.debug() : print "PDE Debug: New tolerance %e",tol
       self.__tolerance=tol
       return
   def getTolerance(self):
       """
       @brief returns the tolerance set for the solution
       """
       return self.__tolerance

   #===== symmetry  flag ==========================
   def isSymmetric(self):
      """
      @brief returns true is the operator is considered to be symmetric
      """
      return self.__sym

   def setSymmetryOn(self):
      """
      @brief sets the symmetry flag to true
      """
      if not self.isSymmetric():
         if self.debug() : print "PDE Debug: Operator is set to be symmetric"
         self.__sym=True
         self.__checkMatrixType()

   def setSymmetryOff(self):
      """
      @brief sets the symmetry flag to false
      """
      if self.isSymmetric():
         if self.debug() : print "PDE Debug: Operator is set to be unsymmetric"
         self.__sym=False
         self.__checkMatrixType()

   def setSymmetryTo(self,flag=False):
     """
     @brief sets the symmetry flag to flag

     @param flag
     """
     if flag:
        self.setSymmetryOn()
     else:
        self.setSymmetryOff()

   #===== order reduction ==========================
   def setReducedOrderOn(self):
     """
     @brief switches to on reduced order
     """
     self.setReducedOrderForSolutionOn()
     self.setReducedOrderForEquationOn()

   def setReducedOrderOff(self):
     """
     @brief switches to full order 
     """
     self.setReducedOrderForSolutionOff()
     self.setReducedOrderForEquationOff()

   def setReducedOrderTo(self,flag=False):
     """
     @brief sets order according to flag

     @param flag
     """
     self.setReducedOrderForSolutionTo(flag)
     self.setReducedOrderForEquationTo(flag)
                                                                                                                                                           

   #===== order reduction solution ==========================
   def setReducedOrderForSolutionOn(self):
     """
     @brief switches to reduced order to interpolate solution
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         if self.debug() : print "PDE Debug: Reduced order is used to interpolate solution."
         self.__column_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForSolutionOff(self):
     """
     @brief switches to full order to interpolate solution
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         if self.debug() : print "PDE Debug: Full order is used to interpolate solution."
         self.__column_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForSolutionTo(self,flag=False):
     """
     @brief sets order for test functions according to flag

     @param flag
     """
     if flag:
        self.setReducedOrderForSolutionOn()
     else:
        self.setReducedOrderForSolutionOff()
                                                                                                                                                           
   #===== order reduction equation ==========================
   def setReducedOrderForEquationOn(self):
     """
     @brief switches to reduced order for test functions
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         if self.debug() : print "PDE Debug: Reduced order is used for test functions."
         self.__row_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForEquationOff(self):
     """
     @brief switches to full order for test functions
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         if self.debug() : print "PDE Debug: Full order is used for test functions."
         self.__row_function_space=new_fs
         self.__rebuildSystem(deep=True)

   def setReducedOrderForEquationTo(self,flag=False):
     """
     @brief sets order for test functions according to flag

     @param flag
     """
     if flag:
        self.setReducedOrderForEquationOn()
     else:
        self.setReducedOrderForEquationOff()
                                                                                                                                                           

   # ==== initialization =====================================================================
   def __makeNewOperator(self):
       """
       @brief
       """
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           self.__matrix_type)

   def __makeNewRightHandSide(self):
       """
       @brief
       """
       return escript.Data(0.,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),True)

   def __makeNewSolution(self):
       """
       @brief
       """
       return escript.Data(0.,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),True)

   def __getFreshOperator(self):
       """
       @brief
       """
       if self.__operator.isEmpty():
           self.__operator=self.__makeNewOperator()
           if self.debug() : print "PDE Debug: New operator allocated"
       else:
           self.__operator.setValue(0.)
           if self.debug() : print "PDE Debug: Operator reset to zero"
       return self.__operator

   def __getFreshRightHandSide(self):
       """
       @brief
       """
       if self.__righthandside.isEmpty():
           self.__righthandside=self.__makeNewRightHandSide()
           if self.debug() : print "PDE Debug: New right hand side allocated"
       else:
           print "fix self.__righthandside*=0"
           self.__righthandside*=0.
           if self.debug() : print "PDE Debug: Right hand side reset to zero"
       return  self.__righthandside

   # ==== rebuild switches =====================================================================
   def __rebuildSolution(self,deep=False):
       """
       @brief indicates the PDE has to be reolved if the solution is requested
       """
       if self.__solution_isValid and self.debug() : print "PDE Debug: PDE has to be resolved."
       self.__solution_isValid=False
       if deep: self.__solution=escript.Data(deep)


   def __rebuildOperator(self,deep=False):
       """
       @brief indicates the operator has to be rebuilt next time it is used
       """
       if self.__operator_isValid and self.debug() : print "PDE Debug: Operator has to be rebuilt."
       self.__rebuildSolution(deep)
       self.__operator_isValid=False
       if deep: self.__operator=escript.Operator()

   def __rebuildRightHandSide(self,deep=False):
       """
       @brief indicates the right hand side has to be rebuild next time it is used
       """
       if self.__righthandside_isValid and self.debug() : print "PDE Debug: Right hand side has to be rebuilt."
       self.__rebuildSolution(deep)
       self.__righthandside_isValid=False
       if not self.__homogeneous_constraint: self.__rebuildOperator()
       if deep: self.__righthandside=escript.Data()

   def __rebuildSystem(self,deep=False):
       """
       @brief annonced that all coefficient name has been changed
       """
       self.__rebuildSolution(deep)
       self.__rebuildOperator(deep)
       self.__rebuildRightHandSide(deep)
   
   def __checkMatrixType(self):
     """
     @brief reassess the matrix type and, if needed, initiates an operator rebuild
     """
     new_matrix_type=self.getDomain().getSystemMatrixTypeId(self.getSolverMethod(),self.isSymmetric())
     if not new_matrix_type==self.__matrix_type:
         if self.debug() : print "PDE Debug: Matrix type is now %d."%new_matrix_type
         self.__matrix_type=new_matrix_type
         self.__rebuildOperator(deep=True)

   #============ assembling =======================================================
   def __copyConstraint(self,u):
      """
      @brief copies the constrint condition into u
      """
      q=self.getCoefficient("q")
      r=self.getCoefficient("r")
      if not q.isEmpty():
          if r.isEmpty():
             r2=escript.Data(0,u.getShape(),u.getFunctionSpace())
          else:
             r2=escript.Data(r,u.getFunctionSpace())
          u.copyWithMask(r2,escript.Data(q,u.getFunctionSpace()))

   def __applyConstraint(self,rhs_update=True):
       """
       @brief applies the constraints  defined by q and r to the system
       """
       q=self.getCoefficient("q")
       r=self.getCoefficient("r")
       if not q.isEmpty() and not self.__operator.isEmpty():
          # q is the row and column mask to indicate where constraints are set:
          row_q=escript.Data(q,self.getFunctionSpaceForEquation())
          col_q=escript.Data(q,self.getFunctionSpaceForSolution())
          u=self.__makeNewSolution()
          if r.isEmpty():
             r_s=self.__makeNewSolution()
          else:
             r_s=escript.Data(r,self.getFunctionSpaceForSolution())
          u.copyWithMask(r_s,col_q)
          if not self.__righthandside.isEmpty() and rhs_update: self.__righthandside-=self.__operator*u
          self.__operator.nullifyRowsAndCols(row_q,col_q,1.)

   def getOperator(self):
       """
       @brief returns the operator of the PDE
       """
       if not self.__operator_isValid:
           # some Constraints are applying for a lumpled stifness matrix:
           if self.isUsingLumping():
              if self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                       raise TypeError,"Lumped matrix requires same order for equations and unknowns"
              if not self.getCoefficient("A").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient A"
              if not self.getCoefficient("B").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient B"
              if not self.getCoefficient("C").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient C"
              if not self.getCoefficient("D").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient D"
              if self.debug() : print "PDE Debug: New lumped operator is built."
              mat=self.__makeNewOperator(self)
           else:
              if self.debug() : print "PDE Debug: New operator is built."
              mat=self.__getFreshOperator()

           self.getDomain().addPDEToSystem(mat,escript.Data(), \
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
           if self.isUsingLumping():
              self.__operator=mat*escript.Data(1,(self.getNumSolutions(),),self.getFunctionSpaceOfSolution(),True)
           else:
              self.__applyConstraint(rhs_update=False)
           self.__operator_isValid=True
       return self.__operator

   def getRightHandSide(self,ignoreConstraint=False):
       """
       @brief returns the right hand side of the PDE

       @param ignoreConstraint
       """
       if not self.__righthandside_isValid:
           if self.debug() : print "PDE Debug: New right hand side is built."
           self.getDomain().addPDEToRHS(self.__getFreshRightHandSide(), \
                         self.getCoefficient("X"), \
                         self.getCoefficient("Y"),\
                         self.getCoefficient("y"),\
                         self.getCoefficient("y_contact"))
           self.__righthandside_isValid=True
           if ignoreConstraint: self.__copyConstraint(self.__righthandside)
       return self.__righthandside

   def getSystem(self):
       """
       @brief
       """
       if not self.__operator_isValid and not self.__righthandside_isValid:
          if self.debug() : print "PDE Debug: New PDE is built."
          if self.isUsingLumping():
              self.getRightHandSide(ignoreConstraint=True)
              self.getOperator()
          else:
              self.getDomain().addPDEToSystem(self.__getFreshOperator(),self.__getFreshRightHandSide(), \
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
          self.__operator_isValid=True
          self.__righthandside_isValid=True
          self.__applyConstraint()
          self.__copyConstraint(self.__righthandside)
       elif not self.__operator_isValid:
          self.getOperator()
       elif not self.__righthandside_isValid:
          self.getRightHandSide()
       return (self.__operator,self.__righthandside)

   def solve(self,**options):
      """
      @brief solve the PDE

      @param options
      """
      mat,f=self.getSystem()
      if self.isUsingLumping():
         out=f/mat
         self.__copyConstraint(out)
      else:
         options[util.TOLERANCE_KEY]=self.getTolerance()
         options[util.METHOD_KEY]=self.getSolverMethod()
         options[util.SYMMETRY_KEY]=self.isSymmetric()
         if self.debug() : print "PDE Debug: solver options: ",options
         out=mat.solve(f,options)
      return out

   def getSolution(self,**options):
       """
       @brief returns the solution of the PDE

       @param options
       """
       if not self.__solution_isValid:
           if self.debug() : print "PDE Debug: PDE is resolved."
           self.__solution=self.solve(**options)
           self.__solution_isValid=True
       return self.__solution
   #============ some serivice functions  =====================================================
   def getDomain(self):
     """
     @brief returns the domain of the PDE
     """
     return self.__domain

   def getNumEquations(self):
     """
     @brief returns the number of equations
     """
     if self.__numEquations>0:
         return self.__numEquations
     else:
         raise ValueError,"Number of equations is undefined. Please specify argument numEquations."

   def getNumSolutions(self):
     """
     @brief returns the number of unknowns
     """
     if self.__numSolutions>0:
        return self.__numSolutions
     else:
        raise ValueError,"Number of solution is undefined. Please specify argument numSolutions."


   def checkSymmetry(self):
      """
      @brief returns if the Operator is symmetric. This is a very expensive operation!!! The symmetry flag is not altered.
      """
      raise SystemError,"checkSymmetry is not implemented yet"

      return None

   def getFlux(self,u):
       """
       @brief returns the flux J_ij for a given u

            J_ij=A_{ijkl}u_{k,l}+B_{ijk}u_k-X_{ij}

       @param u argument of the operator 

       """
       raise SystemError,"getFlux is not implemented yet"
       return None

   def applyOperator(self,u):
       """
       @brief applies the operator of the PDE to a given solution u in weak from

       @param u argument of the operator 

       """
       return self.getOperator()*escript.Data(u,self.getFunctionSpaceForSolution())
                                                                                                                                                           
   def getResidual(self,u):
       """
       @brief return the residual of u in the weak from

       @param u 
       """
       return self.applyOperator(u)-self.getRightHandSide()

class Poisson(LinearPDE):
   """
   @brief Class to define a Poisson equstion problem:
                                                                                                                                                             
   class to define a linear PDE of the form
                                                                                                                                                             
        -u_{,jj} = f
                                                                                                                                                             
     with boundary conditons:
                                                                                                                                                             
        n_j*u_{,j} = 0
                                                                                                                                                             
    and constraints:
                                                                                                                                                             
         u=0 where q>0
                                                                                                                                                             
   """

   def __init__(self,domain=None,f=escript.Data(),q=escript.Data()):
       LinearPDE.__init__(self,domain=identifyDomain(domain,{"f":f, "q":q}))
       self._setCoefficient(A=numarray.identity(self.getDomain().getDim()))
       self.setSymmetryOn()
       self.setValue(f,q)

   def setValue(self,f=escript.Data(),q=escript.Data()):
       self._setCoefficient(Y=f,q=q)
 
                                                                                                                                                           
# $Log$
# Revision 1.2  2004/12/15 07:08:27  jgs
# *** empty log message ***
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
