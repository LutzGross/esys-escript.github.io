# $Id$

## @file linearPDE.py

"""
@brief Functions and classes for linear PDEs
"""

import escript
import util

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
             data_domain=d.getDomain()
     print "identifyDomain: data_domain: ", data_domain
     print "identifyDomain: domain: ", domain
     # check if domain and data_domain are identical?
     if domain == None:
         if data_domain == None:
              raise ValueError,"Undefined PDE domain. Specify a domain or use a Data class object as coefficient"
     else:
         if data_domain == None:
              data_domain=domain
         else:
              if not data_domain==domain:
                 raise ValueError,"Domain of coefficients doesnot match specified domain"
     # now we check if all Data class object coefficients are defined on data_domain:
     for i,d in data.iteritems():
         if isinstance(d,escript.Data):
            if not data_domain==d.getDomain():
                 raise ValueError,"Illegal domain for coefficient %s."%i
     # done:
     return data_domain


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

class linearPDE:
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

   def __init__(self,**args):
     """
     @brief initializes a new linear PDE.

     @param args
     """

     print "Creating new linearPDE"

     # initialize attributes
     self.__debug=None
     self.__domain=None
     self.__numEquations=0
     self.__numSolutions=0
     self.__coefficient={}
     for i in _PDECoefficientTypes.iterkeys():
         self.__coefficient[i]=escript.Data()
     self.__operator=escript.Operator()
     self.__righthandside=escript.Data()
     self.__solution=escript.Data()
     self.__solveroptions=None
     self.__matrixtype=util.UNKNOWN
     self.__sym=False
     coef={}

     # check the arguments
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

     # now each coefficient is changed into a Data object if it is not already one or is None
     for key in coef.iterkeys():
           self.__coefficient[key]=escript.Data(coef[key],_PDECoefficientTypes[key].getFunctionSpace(self.__domain))

     # get number of equations and number of unknowns:
     numEquations=0
     numSolutions=0
     numDim=self.__domain.getDim()
     for i in self.__coefficient.iterkeys():
        if not self.__coefficient[i].isEmpty():
           res=_PDECoefficientTypes[i].estimateNumEquationsAndNumSolutions(self.__coefficient[i].getShape(),numDim)
           if res==None:
               raise ValueError,"Illegal shape %s of coefficient %s"%(self.__coefficient[i].getShape().__str__(),i)
           else:
               numEquations=max(numEquations,res[0])
               numSolutions=max(numSolutions,res[1])

     # check the number of equations:
     if numEquations>0:
        if self.__numEquations>0:
           if self.__numEquations!=numEquations:
              raise ValueError,"Shape of coefficients is inconsistent with the specified number of equations (=%d)"%self.__numEquations
        else:
           self.__numEquations=numEquations
     else:
        if self.__numEquations<1:
            raise ValueError,"Number of equations could not be identified. Please specify argument numEquations."

     # check the number of unknowns:
     if numSolutions>0:
        if self.__numSolutions>0:
           if self.__numSolutions!=numSolutions:
              raise ValueError,"Shape of coefficients is inconsistent with the specified number of unknowns (=%d)"%self.__numSolutions
        else:
           self.__numSolutions=numSolutions
     else:
        if self.__numSolutions<1: self.__numSolutions=self.__numEquations

     print "Identified number of equations and unknowns: (",self.__numEquations,self.__numSolutions,")"

     # check the shape of all coefficients:

     for i,d in self.__coefficient.iteritems():
          if not d.isEmpty():
            if not self.getShapeOfCoefficient(i)==d.getShape():
               raise ValueError,"Expected shape for coefficient %s is %s but actual shape is %s."%(i,self.getShape(i).__str__(),d.getShape().__str__())
     #
     self.__row_function_space=escript.Solution(self.__domain)
     self.__column_function_space=escript.Solution(self.__domain)

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

   def getCoefficient(self,name):
     """
     @brief return the value of the coefficient name

     @param name
     """
     return self.__coefficient[name]

   def setCoefficient(self,**coefficients):
      """
      @brief sets new values to coefficients

      @param coefficients
      """
      alteredCoefficients=[]
      for i,d in coefficients.iteritems():
         if self.hasCoefficient(i):
            if d == None:
                 if not self.__coefficient[i].isEmpty():
                     self.__coefficient[i]=escript.Data()
                     alteredCoefficients.append(i)
            elif isinstance(d,escript.Data):
                if d.isEmpty():
                     if not self.__coefficient[i].isEmpty():
                         self.__coefficient[i]=escript.Data()
                         alteredCoefficients.append(i)
                else:
                     if not self.getShapeOfCoefficient(i)==d.getShape():
                             raise ValueError,"Illegal shape for coefficient %s"%i
                     if not self.getDomain()==d.getDomain():
                             raise ValueError,"Illegal domain for coefficient %s"%i
                     self.__coefficient[i]=escript.Data(d,self.getFunctionSpaceOfCoefficient(i))
                     alteredCoefficients.append(i)
            else:
                self.__coefficient[i]=escript.Data(d,self.getShapeOfCoefficient(i),self.getFunctionSpaceOfCoefficient(i))
                alteredCoefficients.append(i)
            if self.debug(): print "PDE Debug: Coefficient %s has been altered."%i
         else:
            raise ValueError,"Attempt to set unknown coefficient %s"%i

      if len(alteredCoefficients) > 0:
         inhomogeneous=None
         # even if q has not been altered, the system has to be rebuilt if the constraint is not homogeneous:
         if not "q" in alteredCoefficients:
            q=self.getCoefficient("q")
            r=self.getCoefficient("r")
            if not q.isEmpty() and not r.isEmpty():
                       if (q*r).Lsup()>1.e-13*r.Lsup(): inhomogeneous=not None

         if inhomogeneous:
            if self.debug() : print "PDE Debug: System has to be rebuilt because of inhomogeneous constraints."
            self.rebuildSystem()
         else:
            for i in alteredCoefficients:
               self.alteredCoefficient(i)

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

   def getMatrixType(self):
     """
     @brief returns the matrix type to be used to store the operator
     """
     return self.__matrixtype

   def setMatrixType(self,type=util.UNKNOWN):
     """
     @brief sets the new matrix type

     @param type
     """
     if not type==self.__matrixtype:
         if self.debug() : print "PDE Debug: Matrix type is now %d."%type
         self.__matrixtype=type
         self.rebuildOperator()

   def isSymmetric(self):
      """
      @brief returns true is the operator is considered to be symmetric
      """
      return self.__sym

   def setReducedOrderForSolutionsOn(self):
     """
     @brief switches to reduced order to interpolate solution
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         if self.debug() : print "PDE Debug: Reduced order is used to interpolate solution."
         self.__column_function_space=new_fs
         self.rebuildOperator()

   def setReducedOrderForSolutionsOff(self):
     """
     @brief switches to full order to interpolate solution
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForSolution()!=new_fs:
         if self.debug() : print "PDE Debug: Full order is used to interpolate solution."
         self.__column_function_space=new_fs
         self.rebuildOperator()

   def setReducedOrderForEquationsOn(self):
     """
     @brief switches to reduced order for test functions
     """
     new_fs=escript.ReducedSolution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         if self.debug() : print "PDE Debug: Reduced order is used for test functions."
         self.__row_function_space=new_fs
         self.rebuildOperator()
         self.rebuildRightHandSide()

   def setReducedOrderForSolutionsOff(self):
     """
     @brief switches to full order for test functions
     """
     new_fs=escript.Solution(self.getDomain())
     if self.getFunctionSpaceForEquation()!=new_fs:
         if self.debug() : print "PDE Debug: Full order is used for test functions."
         self.__row_function_space=new_fs
         self.rebuildOperator()
         self.rebuildRightHandSide()

   def getDomain(self):
     """
     @brief returns the domain of the PDE
     """
     return self.__domain

   def getNumEquations(self):
     """
     @brief returns the number of equations
     """
     return self.__numEquations

   def getNumSolutions(self):
     """
     @brief returns the number of unknowns
     """
     return self.__numSolutions

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
        if _PDECoefficientTypes[name].isAlteringOperator(): self.rebuildOperator()
        if _PDECoefficientTypes[name].isAlteringRightHandSide(): self.rebuildRightHandSide()
     else:
        raise ValueError,"Solution coefficient %s requested"%name

   def initialiseNewOperator(self):
       """
       @brief
       """
       return self.getDomain().newOperator( \
                           self.getNumEquations(), \
                           self.getFunctionSpaceForEquation(), \
                           self.getNumSolutions(), \
                           self.getFunctionSpaceForSolution(), \
                           self.getMatrixType(), \
                           self.isSymmetric())

   def initialiseNewRightHandSide(self,expanded=True):
       """
       @brief

       @param expanded
       """
       return escript.Data(0.,(self.getNumEquations(),),self.getFunctionSpaceForEquation(),expanded)

   def initialiseNewSolution(self,expanded=True):
       """
       @brief

       @param expanded
       """
       return escript.Data(0.,(self.getNumSolutions(),),self.getFunctionSpaceForSolution(),expanded)

   def rebuildSolution(self):
       """
       @brief indicates the PDE has to be reolved if the solution is requested
       """
       if not self.__solution.isEmpty() and self.debug() : print "PDE Debug: PDE has to be resolved."
       self.__solution=escript.Data()

   def rebuildOperator(self):
       """
       @brief indicates the operator has to be rebuilt next time it is used
       """
       if not self.__operator.isEmpty() and self.debug() : print "PDE Debug: Operator has to be rebuilt."
       self.rebuildSolution()
       self.__operator=escript.Operator()

   def rebuildRightHandSide(self):
       """
       @brief indicates the right hand side has to be rebuild next time it is used
       """
       if not self.__righthandside.isEmpty() and self.debug() : print "PDE Debug: Right hand side has to be rebuilt."
       self.rebuildSolution()
       self.__righthandside=escript.Data()

   def rebuildSystem(self):
       """
       @brief annonced that all coefficient name has been changed
       """
       self.rebuildSolution()
       self.rebuildOperator()
       self.rebuildRightHandSide()

   def applyConstraint(self):
       """
       @brief applies the constraints  defined by q and r to the system
       """
       q=self.getCoefficient("q")
       r=self.getCoefficient("r")
       rhs=self.getRightHandSide()
       mat=self.getOperator()
       if not q.isEmpty():
          # q is the row and column mask to indicate where constraints are set:
          row_q=escript.Data(q,self.getFunctionSpaceForEquation())
          col_q=escript.Data(q,self.getFunctionSpaceForSolution())
          if r.isEmpty():
              r2=self.initialiseNewRightHandSide()
          else:
              u=self.initialiseNewSolution()
              src=escript.Data(r,self.getFunctionSpaceForEquation())
              u.copyWithMask(src,row_q)
              if not rhs.isEmpty(): rhs-=mat*u
              r2=escript.Data(r,self.getFunctionSpaceForEquation())
          if not rhs.isEmpty(): rhs.copyWithMask(r2,row_q)
          if not mat.isEmpty(): mat.nullifyRowsAndCols(row_q,col_q,1.)

   def getOperator(self):
       """
       @brief returns the operator of the PDE
       """
       if self.__operator.isEmpty():
           if self.debug() : print "PDE Debug: New operator is built."
           # some constrainst are applying for a lumpled stifness matrix:
           if self.getMatrixType()==util.LUMPED:
              if self.getFunctionSpaceForEquation()==self.getFunctionSpaceForSolution():
                       raise Warning,"Lumped matrix requires same order for equations and unknowns"
              if not self.getCoefficient("A").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient A"
              if not self.getCoefficient("B").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient B"
              if not self.getCoefficient("C").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient C"
              if not self.getCoefficient("D").isEmpty():
                       raise Warning,"Lumped matrix does not allow coefficient D"
           # get a new matrix:
           mat=self.initialiseNewOperator()
           #
           #   assemble stiffness matrix:
           #
           self.getDomain().addPDEToSystem(mat,escript.Data(), \
                         self.getCoefficient("A"), \
                         self.getCoefficient("B"), \
                         self.getCoefficient("C"), \
                         self.getCoefficient("D"), \
                         escript.Data(), \
                         escript.Data())
           self.getDomain().addRobinConditionsToSystem(mat,escript.Data(), \
                         self.getCoefficient("d"), \
                         escript.Data())
           self.getDomain().addContactToSystem(mat,escript.Data(), \
                         self.getCoefficient("d_contact"), \
                         escript.Data())
           self.applyConstraint()
           if self.getMatrixType()==util.LUMPED:
                self.__operator=mat*escript.Data(1,(self.getNumSolutions(),),self.getFunctionSpaceOfSolution(),True)
           else:
                self.__operator=mat
       return self.__operator

   def getRightHandSide(self,ignoreConstraint=None):
       """
       @brief returns the right hand side of the PDE

       @param ignoreConstraint
       """
       if self.__righthandside.isEmpty():
           if self.debug() : print "PDE Debug: New right hand side is built."
           self.__righthandside=self.initialiseNewRightHandSide()
           self.getDomain().addPDEToSystem(escript.Operator(),self.__righthandside, \
                         escript.Data(), \
                         escript.Data(), \
                         escript.Data(), \
                         escript.Data(), \
                         self.getCoefficient("X"), \
                         self.getCoefficient("Y"))
           self.getDomain().addRobinConditionsToSystem(escript.Operator(),self.__righthandside, \
                         escript.Data(), \
                         self.getCoefficient("y"))
           self.getDomain().addContactToSystem(escript.Operator(),self.__righthandside, \
                         escript.Data(), \
                         self.getCoefficient("y_contact"))
           if not ignoreConstraint: self.applyConstraint()
       return self.__righthandside

   def getSystem(self):
       """
       @brief
       """
       if self.__operator.isEmpty() and self.__righthandside.isEmpty():
          if self.debug() : print "PDE Debug: New PDE is built."
          if self.getMatrixType()==util.LUMPED:
              self.getRightHandSide(ignoreConstraint=not None)
              self.getOperator()
          else:
              self.__righthandside=self.initialiseNewRightHandSide()
              self.__operator=self.initialiseNewOperator()
              self.getDomain().addPDEToSystem(self.__operator,self.__righthandside, \
                            self.getCoefficient("A"), \
                            self.getCoefficient("B"), \
                            self.getCoefficient("C"), \
                            self.getCoefficient("D"), \
                            self.getCoefficient("X"), \
                            self.getCoefficient("Y"))
              self.getDomain().addRobinConditionsToSystem(self.__operator,self.__righthandside, \
                            self.getCoefficient("d"), \
                            self.getCoefficient("y"))
              self.getDomain().addContactToSystem(self.__operator,self.__righthandside, \
                            self.getCoefficient("d_contact"), \
                            self.getCoefficient("y_contact"))
              self.applyConstraint()
       elif self.__operator.isEmpty():
            self.getOperator()
       elif self.__righthandside.isEmpty():
            self.getRightHandSide()
       return (self.__operator,self.__righthandside)

   def solve(self,**options):
      """
      @brief solve the PDE

      @param options
      """
      if not options.has_key("iterative"): options["iterative"]=True
      mat,f=self.getSystem()
      if isinstance(mat,escript.Data):
         return f/mat
      else:
         return mat.solve(f,options)

   def getSolution(self,**options):
       """
       @brief returns the solution of the PDE

       @param options
       """
       if self.__solution.isEmpty() or not self.__solveroptions==options:
           self.__solveroptions=options
           self.__solution=self.solve(**options)
           if self.debug() : print "PDE Debug: PDE is resolved."
       return self.__solution

# $Log$
# Revision 1.1  2004/10/26 06:53:56  jgs
# Initial revision
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
