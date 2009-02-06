
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
Provides some tools related to PDEs.

Currently includes:
    - Projector - to project a discontinuous function onto a continuous function
    - Locator - to trace values in data objects at a certain location
    - TimeIntegrationManager - to handle extrapolation in time
    - SaddlePointProblem - solver for Saddle point problems using the inexact uszawa scheme

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"


import escript
import linearPDEs
import numarray
import util
import math

##### Added by Artak
# from Numeric import zeros,Int,Float64
###################################


class TimeIntegrationManager:
  """
  A simple mechanism to manage time dependend values.

  Typical usage is::

     dt=0.1 # time increment
     tm=TimeIntegrationManager(inital_value,p=1)
     while t<1.
         v_guess=tm.extrapolate(dt) # extrapolate to t+dt
         v=...
         tm.checkin(dt,v)
         t+=dt

  @note: currently only p=1 is supported.
  """
  def __init__(self,*inital_values,**kwargs):
     """
     Sets up the value manager where C{inital_values} are the initial values
     and p is the order used for extrapolation.
     """
     if kwargs.has_key("p"):
            self.__p=kwargs["p"]
     else:
            self.__p=1
     if kwargs.has_key("time"):
            self.__t=kwargs["time"]
     else:
            self.__t=0.
     self.__v_mem=[inital_values]
     self.__order=0
     self.__dt_mem=[]
     self.__num_val=len(inital_values)

  def getTime(self):
      return self.__t

  def getValue(self):
      out=self.__v_mem[0]
      if len(out)==1:
          return out[0]
      else:
          return out

  def checkin(self,dt,*values):
      """
      Adds new values to the manager. The p+1 last values are lost.
      """
      o=min(self.__order+1,self.__p)
      self.__order=min(self.__order+1,self.__p)
      v_mem_new=[values]
      dt_mem_new=[dt]
      for i in range(o-1):
         v_mem_new.append(self.__v_mem[i])
         dt_mem_new.append(self.__dt_mem[i])
      v_mem_new.append(self.__v_mem[o-1])
      self.__order=o
      self.__v_mem=v_mem_new
      self.__dt_mem=dt_mem_new
      self.__t+=dt

  def extrapolate(self,dt):
      """
      Extrapolates to C{dt} forward in time.
      """
      if self.__order==0:
         out=self.__v_mem[0]
      else:
        out=[]
        for i in range(self.__num_val):
           out.append((1.+dt/self.__dt_mem[0])*self.__v_mem[0][i]-dt/self.__dt_mem[0]*self.__v_mem[1][i])

      if len(out)==0:
         return None
      elif len(out)==1:
         return out[0]
      else:
         return out


class Projector:
  """
  The Projector is a factory which projects a discontinuous function onto a
  continuous function on a given domain.
  """
  def __init__(self, domain, reduce=True, fast=True):
    """
    Creates a continuous function space projector for a domain.

    @param domain: Domain of the projection.
    @param reduce: Flag to reduce projection order
    @param fast: Flag to use a fast method based on matrix lumping
    """
    self.__pde = linearPDEs.LinearPDE(domain)
    if fast:
        self.__pde.setSolverMethod(linearPDEs.LinearPDE.LUMPING)
    self.__pde.setSymmetryOn()
    self.__pde.setReducedOrderTo(reduce)
    self.__pde.setValue(D = 1.)
    return

  def __call__(self, input_data):
    """
    Projects C{input_data} onto a continuous function.

    @param input_data: the data to be projected
    """
    out=escript.Data(0.,input_data.getShape(),self.__pde.getFunctionSpaceForSolution())
    self.__pde.setValue(Y = escript.Data(), Y_reduced = escript.Data())
    if input_data.getRank()==0:
        self.__pde.setValue(Y = input_data)
        out=self.__pde.getSolution()
    elif input_data.getRank()==1:
        for i0 in range(input_data.getShape()[0]):
           self.__pde.setValue(Y = input_data[i0])
           out[i0]=self.__pde.getSolution()
    elif input_data.getRank()==2:
        for i0 in range(input_data.getShape()[0]):
           for i1 in range(input_data.getShape()[1]):
              self.__pde.setValue(Y = input_data[i0,i1])
              out[i0,i1]=self.__pde.getSolution()
    elif input_data.getRank()==3:
        for i0 in range(input_data.getShape()[0]):
           for i1 in range(input_data.getShape()[1]):
              for i2 in range(input_data.getShape()[2]):
                 self.__pde.setValue(Y = input_data[i0,i1,i2])
                 out[i0,i1,i2]=self.__pde.getSolution()
    else:
        for i0 in range(input_data.getShape()[0]):
           for i1 in range(input_data.getShape()[1]):
              for i2 in range(input_data.getShape()[2]):
                 for i3 in range(input_data.getShape()[3]):
                    self.__pde.setValue(Y = input_data[i0,i1,i2,i3])
                    out[i0,i1,i2,i3]=self.__pde.getSolution()
    return out

class NoPDE:
     """
     Solves the following problem for u:

     M{kronecker[i,j]*D[j]*u[j]=Y[i]}

     with constraint

     M{u[j]=r[j]}  where M{q[j]>0}

     where M{D}, M{Y}, M{r} and M{q} are given functions of rank 1.

     In the case of scalars this takes the form

     M{D*u=Y}

     with constraint

     M{u=r} where M{q>0}

     where M{D}, M{Y}, M{r} and M{q} are given scalar functions.

     The constraint overwrites any other condition.

     @note: This class is similar to the L{linearPDEs.LinearPDE} class with
            A=B=C=X=0 but has the intention that all input parameters are given
            in L{Solution} or L{ReducedSolution}.
     """
     # The whole thing is a bit strange and I blame Rob Woodcock (CSIRO) for
     # this.
     def __init__(self,domain,D=None,Y=None,q=None,r=None):
         """
         Initializes the problem.

         @param domain: domain of the PDE
         @type domain: L{Domain}
         @param D: coefficient of the solution
         @type D: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         @param Y: right hand side
         @type Y: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         @param q: location of constraints
         @type q: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         @param r: value of solution at locations of constraints
         @type r: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         """
         self.__domain=domain
         self.__D=D
         self.__Y=Y
         self.__q=q
         self.__r=r
         self.__u=None
         self.__function_space=escript.Solution(self.__domain)

     def setReducedOn(self):
         """
         Sets the L{FunctionSpace} of the solution to L{ReducedSolution}.
         """
         self.__function_space=escript.ReducedSolution(self.__domain)
         self.__u=None

     def setReducedOff(self):
         """
         Sets the L{FunctionSpace} of the solution to L{Solution}.
         """
         self.__function_space=escript.Solution(self.__domain)
         self.__u=None

     def setValue(self,D=None,Y=None,q=None,r=None):
         """
         Assigns values to the parameters.

         @param D: coefficient of the solution
         @type D: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         @param Y: right hand side
         @type Y: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         @param q: location of constraints
         @type q: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         @param r: value of solution at locations of constraints
         @type r: C{float}, C{int}, C{numarray.NumArray}, L{Data}
         """
         if not D==None:
            self.__D=D
            self.__u=None
         if not Y==None:
            self.__Y=Y
            self.__u=None
         if not q==None:
            self.__q=q
            self.__u=None
         if not r==None:
            self.__r=r
            self.__u=None

     def getSolution(self):
         """
         Returns the solution.

         @return: the solution of the problem
         @rtype: L{Data} object in the L{FunctionSpace} L{Solution} or
                 L{ReducedSolution}
         """
         if self.__u==None:
            if self.__D==None:
               raise ValueError,"coefficient D is undefined"
            D=escript.Data(self.__D,self.__function_space)
            if D.getRank()>1:
               raise ValueError,"coefficient D must have rank 0 or 1"
            if self.__Y==None:
               self.__u=escript.Data(0.,D.getShape(),self.__function_space)
            else:
               self.__u=util.quotient(self.__Y,D)
            if not self.__q==None:
                q=util.wherePositive(escript.Data(self.__q,self.__function_space))
                self.__u*=(1.-q)
                if not self.__r==None: self.__u+=q*self.__r
         return self.__u

class Locator:
     """
     Locator provides access to the values of data objects at a given spatial
     coordinate x.

     In fact, a Locator object finds the sample in the set of samples of a
     given function space or domain which is closest to the given point x.
     """

     def __init__(self,where,x=numarray.zeros((3,))):
       """
       Initializes a Locator to access values in Data objects on the Doamin
       or FunctionSpace for the sample point which is closest to the given
       point x.

       @param where: function space
       @type where: L{escript.FunctionSpace}
       @param x: coefficient of the solution
       @type x: C{numarray.NumArray} or C{list} of C{numarray.NumArray}
       """
       if isinstance(where,escript.FunctionSpace):
          self.__function_space=where
       else:
          self.__function_space=escript.ContinuousFunction(where)
       if isinstance(x, list):
           self.__id=[]
           for p in x:
              self.__id.append(util.length(self.__function_space.getX()-p[:self.__function_space.getDim()]).minGlobalDataPoint())
       else:
           self.__id=util.length(self.__function_space.getX()-x[:self.__function_space.getDim()]).minGlobalDataPoint()

     def __str__(self):
       """
       Returns the coordinates of the Locator as a string.
       """
       x=self.getX()
       if instance(x,list):
          out="["
          first=True
          for xx in x:
            if not first:
                out+=","
            else:
                first=False
            out+=str(xx)
          out+="]>"
       else:
          out=str(x)
       return out

     def getX(self):
        """
        Returns the exact coordinates of the Locator.
        """
        return self(self.getFunctionSpace().getX())

     def getFunctionSpace(self):
        """
        Returns the function space of the Locator.
        """
        return self.__function_space

     def getId(self,item=None):
        """
        Returns the identifier of the location.
        """
        if item == None:
           return self.__id
        else:
           if isinstance(self.__id,list):
              return self.__id[item]
           else:
              return self.__id


     def __call__(self,data):
        """
        Returns the value of data at the Locator of a Data object.
        """
        return self.getValue(data)

     def getValue(self,data):
        """
        Returns the value of C{data} at the Locator if C{data} is a L{Data}
        object otherwise the object is returned.
        """
        if isinstance(data,escript.Data):
           if data.getFunctionSpace()==self.getFunctionSpace():
             dat=data
           else:
             dat=data.interpolate(self.getFunctionSpace())
           id=self.getId()
           r=data.getRank()
           if isinstance(id,list):
               out=[]
               for i in id:
                  o=data.getValueOfGlobalDataPoint(*i)
                  if data.getRank()==0:
                     out.append(o[0])
                  else:
                     out.append(o)
               return out
           else:
             out=data.getValueOfGlobalDataPoint(*id)
             if data.getRank()==0:
                return out[0]
             else:
                return out
        else:
           return data

class SolverSchemeException(Exception):
   """
   This is a generic exception thrown by solvers.
   """
   pass

class IndefinitePreconditioner(SolverSchemeException):
   """
   Exception thrown if the preconditioner is not positive definite.
   """
   pass

class MaxIterReached(SolverSchemeException):
   """
   Exception thrown if the maximum number of iteration steps is reached.
   """
   pass

class CorrectionFailed(SolverSchemeException):
   """
   Exception thrown if no convergence has been achieved in the solution
   correction scheme.
   """
   pass

class IterationBreakDown(SolverSchemeException):
   """
   Exception thrown if the iteration scheme encountered an incurable breakdown.
   """
   pass

class NegativeNorm(SolverSchemeException):
   """
   Exception thrown if a norm calculation returns a negative norm.
   """
   pass

def PCG(r, Aprod, x, Msolve, bilinearform, atol=0, rtol=1.e-8, iter_max=100, initial_guess=True, verbose=False):
   """
   Solver for

   M{Ax=b}

   with a symmetric and positive definite operator A (more details required!).
   It uses the conjugate gradient method with preconditioner M providing an
   approximation of A.

   The iteration is terminated if

   M{|r| <= atol+rtol*|r0|}

   where M{r0} is the initial residual and M{|.|} is the energy norm. In fact

   M{|r| = sqrt( bilinearform(Msolve(r),r))}

   For details on the preconditioned conjugate gradient method see the book:

   I{Templates for the Solution of Linear Systems by R. Barrett, M. Berry,
   T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo,
   C. Romine, and H. van der Vorst}.

   @param r: initial residual M{r=b-Ax}. C{r} is altered.
   @type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   @param x: an initial guess for the solution
   @type x: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   @param Aprod: returns the value Ax
   @type Aprod: function C{Aprod(x)} where C{x} is of the same object like
                argument C{x}. The returned object needs to be of the same type
                like argument C{r}.
   @param Msolve: solves Mx=r
   @type Msolve: function C{Msolve(r)} where C{r} is of the same type like
                 argument C{r}. The returned object needs to be of the same
                 type like argument C{x}.
   @param bilinearform: inner product C{<x,r>}
   @type bilinearform: function C{bilinearform(x,r)} where C{x} is of the same
                       type like argument C{x} and C{r} is. The returned value
                       is a C{float}.
   @param atol: absolute tolerance
   @type atol: non-negative C{float}
   @param rtol: relative tolerance
   @type rtol: non-negative C{float}
   @param iter_max: maximum number of iteration steps
   @type iter_max: C{int}
   @return: the solution approximation and the corresponding residual
   @rtype: C{tuple}
   @warning: C{r} and C{x} are altered.
   """
   iter=0
   rhat=Msolve(r)
   d = rhat
   rhat_dot_r = bilinearform(rhat, r)
   if rhat_dot_r<0: raise NegativeNorm,"negative norm."
   norm_r0=math.sqrt(rhat_dot_r)
   atol2=atol+rtol*norm_r0
   if atol2<=0:
      raise ValueError,"Non-positive tolarance."
   atol2=max(atol2, 100. * util.EPSILON * norm_r0)

   if verbose: print "PCG: initial residual norm = %e (absolute tolerance = %e)"%(norm_r0, atol2)


   while not math.sqrt(rhat_dot_r) <= atol2:
       iter+=1
       if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

       q=Aprod(d)
       alpha = rhat_dot_r / bilinearform(d, q)
       x += alpha * d
       r += (-alpha) * q

       rhat=Msolve(r)
       rhat_dot_r_new = bilinearform(rhat, r)
       beta = rhat_dot_r_new / rhat_dot_r
       rhat+=beta * d
       d=rhat

       rhat_dot_r = rhat_dot_r_new
       if rhat_dot_r<0: raise NegativeNorm,"negative norm."
       if verbose: print "PCG: iteration step %s: residual norm = %e"%(iter, math.sqrt(rhat_dot_r))
   if verbose: print "PCG: tolerance reached after %s steps."%iter
   return x,r

class Defect(object):
    """
    Defines a non-linear defect F(x) of a variable x.
    """
    def __init__(self):
        """
        Initializes defect.
        """
        self.setDerivativeIncrementLength()

    def bilinearform(self, x0, x1):
        """
        Returns the inner product of x0 and x1

        @param x0: value for x0
        @param x1: value for x1
        @return: the inner product of x0 and x1
        @rtype: C{float}
        """
        return 0

    def norm(self,x):
        """
        Returns the norm of argument C{x}.

        @param x: a value
        @return: norm of argument x
        @rtype: C{float}
        @note: by default C{sqrt(self.bilinearform(x,x)} is returned.
        """
        s=self.bilinearform(x,x)
        if s<0: raise NegativeNorm,"negative norm."
        return math.sqrt(s)

    def eval(self,x):
        """
        Returns the value F of a given C{x}.

        @param x: value for which the defect C{F} is evaluated
        @return: value of the defect at C{x}
        """
        return 0

    def __call__(self,x):
        return self.eval(x)

    def setDerivativeIncrementLength(self,inc=math.sqrt(util.EPSILON)):
        """
        Sets the relative length of the increment used to approximate the
        derivative of the defect. The increment is inc*norm(x)/norm(v)*v in the
        direction of v with x as a starting point.

        @param inc: relative increment length
        @type inc: positive C{float}
        """
        if inc<=0: raise ValueError,"positive increment required."
        self.__inc=inc

    def getDerivativeIncrementLength(self):
        """
        Returns the relative increment length used to approximate the
        derivative of the defect.
        @return: value of the defect at C{x}
        @rtype: positive C{float}
        """
        return self.__inc

    def derivative(self, F0, x0, v, v_is_normalised=True):
        """
        Returns the directional derivative at C{x0} in the direction of C{v}.

        @param F0: value of this defect at x0
        @param x0: value at which derivative is calculated
        @param v: direction
        @param v_is_normalised: True to indicate that C{v} is nomalized
                                (self.norm(v)=0)
        @return: derivative of this defect at x0 in the direction of C{v}
        @note: by default numerical evaluation (self.eval(x0+eps*v)-F0)/eps is
               used but this method maybe overwritten to use exact evaluation.
        """
        normx=self.norm(x0)
        if normx>0:
             epsnew = self.getDerivativeIncrementLength() * normx
        else:
             epsnew = self.getDerivativeIncrementLength()
        if not v_is_normalised:
            normv=self.norm(v)
            if normv<=0:
               return F0*0
            else:
               epsnew /= normv
        F1=self.eval(x0 + epsnew * v)
        return (F1-F0)/epsnew

######################################
def NewtonGMRES(defect, x, iter_max=100, sub_iter_max=20, atol=0,rtol=1.e-4, sub_tol_max=0.5, gamma=0.9, verbose=False):
   """
   Solves a non-linear problem M{F(x)=0} for unknown M{x} using the stopping
   criterion:

   M{norm(F(x) <= atol + rtol * norm(F(x0)}

   where M{x0} is the initial guess.

   @param defect: object defining the function M{F}. C{defect.norm} defines the
                  M{norm} used in the stopping criterion.
   @type defect: L{Defect}
   @param x: initial guess for the solution, C{x} is altered.
   @type x: any object type allowing basic operations such as
            C{numarray.NumArray}, L{Data}
   @param iter_max: maximum number of iteration steps
   @type iter_max: positive C{int}
   @param sub_iter_max: maximum number of inner iteration steps
   @type sub_iter_max: positive C{int}
   @param atol: absolute tolerance for the solution
   @type atol: positive C{float}
   @param rtol: relative tolerance for the solution
   @type rtol: positive C{float}
   @param gamma: tolerance safety factor for inner iteration
   @type gamma: positive C{float}, less than 1
   @param sub_tol_max: upper bound for inner tolerance
   @type sub_tol_max: positive C{float}, less than 1
   @return: an approximation of the solution with the desired accuracy
   @rtype: same type as the initial guess
   """
   lmaxit=iter_max
   if atol<0: raise ValueError,"atol needs to be non-negative."
   if rtol<0: raise ValueError,"rtol needs to be non-negative."
   if rtol+atol<=0: raise ValueError,"rtol or atol needs to be non-negative."
   if gamma<=0 or gamma>=1: raise ValueError,"tolerance safety factor for inner iteration (gamma =%s) needs to be positive and less than 1."%gamma
   if sub_tol_max<=0 or sub_tol_max>=1: raise ValueError,"upper bound for inner tolerance for inner iteration (sub_tol_max =%s) needs to be positive and less than 1."%sub_tol_max

   F=defect(x)
   fnrm=defect.norm(F)
   stop_tol=atol + rtol*fnrm
   sub_tol=sub_tol_max
   if verbose: print "NewtonGMRES: initial residual = %e."%fnrm
   if verbose: print "             tolerance = %e."%sub_tol
   iter=1
   #
   # main iteration loop
   #
   while not fnrm<=stop_tol:
            if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max
            #
	    #   adjust sub_tol_
	    #
            if iter > 1:
	       rat=fnrm/fnrmo
               sub_tol_old=sub_tol
	       sub_tol=gamma*rat**2
	       if gamma*sub_tol_old**2 > .1: sub_tol=max(sub_tol,gamma*sub_tol_old**2)
	       sub_tol=max(min(sub_tol,sub_tol_max), .5*stop_tol/fnrm)
	    #
	    # calculate newton increment xc
            #     if iter_max in __FDGMRES is reached MaxIterReached is thrown
            #     if iter_restart -1 is returned as sub_iter
            #     if  atol is reached sub_iter returns the numer of steps performed to get there
            #
            #
            if verbose: print "             subiteration (GMRES) is called with relative tolerance %e."%sub_tol
            try:
               xc, sub_iter=__FDGMRES(F, defect, x, sub_tol*fnrm, iter_max=iter_max-iter, iter_restart=sub_iter_max)
            except MaxIterReached:
               raise MaxIterReached,"maximum number of %s steps reached."%iter_max
            if sub_iter<0:
               iter+=sub_iter_max
            else:
               iter+=sub_iter
            # ====
	    x+=xc
            F=defect(x)
	    iter+=1
            fnrmo, fnrm=fnrm, defect.norm(F)
            if verbose: print "             step %s: residual %e."%(iter,fnrm)
   if verbose: print "NewtonGMRES: completed after %s steps."%iter
   return x

def __givapp(c,s,vin):
    """
    Applies a sequence of Givens rotations (c,s) recursively to the vector
    C{vin}

    @warning: C{vin} is altered.
    """
    vrot=vin
    if isinstance(c,float):
        vrot=[c*vrot[0]-s*vrot[1],s*vrot[0]+c*vrot[1]]
    else:
        for i in range(len(c)):
            w1=c[i]*vrot[i]-s[i]*vrot[i+1]
	    w2=s[i]*vrot[i]+c[i]*vrot[i+1]
            vrot[i:i+2]=w1,w2
    return vrot

def __FDGMRES(F0, defect, x0, atol, iter_max=100, iter_restart=20):
   h=numarray.zeros((iter_restart,iter_restart),numarray.Float64)
   c=numarray.zeros(iter_restart,numarray.Float64)
   s=numarray.zeros(iter_restart,numarray.Float64)
   g=numarray.zeros(iter_restart,numarray.Float64)
   v=[]

   rho=defect.norm(F0)
   if rho<=0.: return x0*0

   v.append(-F0/rho)
   g[0]=rho
   iter=0
   while rho > atol and iter<iter_restart-1:
        if iter  >= iter_max:
            raise MaxIterReached,"maximum number of %s steps reached."%iter_max

        p=defect.derivative(F0,x0,v[iter], v_is_normalised=True)
        v.append(p)

        v_norm1=defect.norm(v[iter+1])

        # Modified Gram-Schmidt
        for j in range(iter+1):
            h[j,iter]=defect.bilinearform(v[j],v[iter+1])
            v[iter+1]-=h[j,iter]*v[j]

        h[iter+1,iter]=defect.norm(v[iter+1])
        v_norm2=h[iter+1,iter]

        # Reorthogonalize if needed
        if v_norm1 + 0.001*v_norm2 == v_norm1:   #Brown/Hindmarsh condition (default)
            for j in range(iter+1):
                hr=defect.bilinearform(v[j],v[iter+1])
                h[j,iter]=h[j,iter]+hr
                v[iter+1] -= hr*v[j]

            v_norm2=defect.norm(v[iter+1])
            h[iter+1,iter]=v_norm2
        #   watch out for happy breakdown
        if not v_norm2 == 0:
            v[iter+1]=v[iter+1]/h[iter+1,iter]

        #   Form and store the information for the new Givens rotation
        if iter > 0 :
            hhat=numarray.zeros(iter+1,numarray.Float64)
            for i in range(iter+1) : hhat[i]=h[i,iter]
            hhat=__givapp(c[0:iter],s[0:iter],hhat);
            for i in range(iter+1) : h[i,iter]=hhat[i]

        mu=math.sqrt(h[iter,iter]*h[iter,iter]+h[iter+1,iter]*h[iter+1,iter])

        if mu!=0 :
            c[iter]=h[iter,iter]/mu
            s[iter]=-h[iter+1,iter]/mu
            h[iter,iter]=c[iter]*h[iter,iter]-s[iter]*h[iter+1,iter]
            h[iter+1,iter]=0.0
            g[iter:iter+2]=__givapp(c[iter],s[iter],g[iter:iter+2])

        # Update the residual norm
        rho=abs(g[iter+1])
        iter+=1

   # At this point either iter > iter_max or rho < tol.
   # It's time to compute x and leave.
   if iter > 0 :
     y=numarray.zeros(iter,numarray.Float64)
     y[iter-1] = g[iter-1] / h[iter-1,iter-1]
     if iter > 1 :
        i=iter-2
        while i>=0 :
          y[i] = ( g[i] - numarray.dot(h[i,i+1:iter], y[i+1:iter])) / h[i,i]
          i=i-1
     xhat=v[iter-1]*y[iter-1]
     for i in range(iter-1):
	xhat += v[i]*y[i]
   else :
      xhat=v[0] * 0

   if iter<iter_restart-1:
      stopped=iter
   else:
      stopped=-1

   return xhat,stopped

def GMRES(r, Aprod, x, bilinearform, atol=0, rtol=1.e-8, iter_max=100, iter_restart=20, verbose=False):
   """
   Solver for

   M{Ax=b}

   with a general operator A (more details required!).
   It uses the generalized minimum residual method (GMRES).

   The iteration is terminated if

   M{|r| <= atol+rtol*|r0|}

   where M{r0} is the initial residual and M{|.|} is the energy norm. In fact

   M{|r| = sqrt( bilinearform(r,r))}

   @param r: initial residual M{r=b-Ax}. C{r} is altered.
   @type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   @param x: an initial guess for the solution
   @type x: same like C{r}
   @param Aprod: returns the value Ax
   @type Aprod: function C{Aprod(x)} where C{x} is of the same object like
                argument C{x}. The returned object needs to be of the same
                type like argument C{r}.
   @param bilinearform: inner product C{<x,r>}
   @type bilinearform: function C{bilinearform(x,r)} where C{x} is of the same
                       type like argument C{x} and C{r}. The returned value is
                       a C{float}.
   @param atol: absolute tolerance
   @type atol: non-negative C{float}
   @param rtol: relative tolerance
   @type rtol: non-negative C{float}
   @param iter_max: maximum number of iteration steps
   @type iter_max: C{int}
   @param iter_restart: in order to save memory the orthogonalization process
                        is terminated after C{iter_restart} steps and the
                        iteration is restarted.
   @type iter_restart: C{int}
   @return: the solution approximation and the corresponding residual
   @rtype: C{tuple}
   @warning: C{r} and C{x} are altered.
   """
   m=iter_restart
   restarted=False
   iter=0
   if rtol>0:
      r_dot_r = bilinearform(r, r)
      if r_dot_r<0: raise NegativeNorm,"negative norm."
      atol2=atol+rtol*math.sqrt(r_dot_r)
      if verbose: print "GMRES: norm of right hand side = %e (absolute tolerance = %e)"%(math.sqrt(r_dot_r), atol2)
   else:
      atol2=atol
      if verbose: print "GMRES: absolute tolerance = %e"%atol2
   if atol2<=0:
      raise ValueError,"Non-positive tolarance."

   while True:
      if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached"%iter_max
      if restarted:
         r2 = r-Aprod(x-x2)
      else:
         r2=1*r
      x2=x*1.
      x,stopped=_GMRESm(r2, Aprod, x, bilinearform, atol2, iter_max=iter_max-iter, iter_restart=m, verbose=verbose)
      iter+=iter_restart
      if stopped: break
      if verbose: print "GMRES: restart."
      restarted=True
   if verbose: print "GMRES: tolerance has been reached."
   return x

def _GMRESm(r, Aprod, x, bilinearform, atol, iter_max=100, iter_restart=20, verbose=False):
   iter=0

   h=numarray.zeros((iter_restart+1,iter_restart),numarray.Float64)
   c=numarray.zeros(iter_restart,numarray.Float64)
   s=numarray.zeros(iter_restart,numarray.Float64)
   g=numarray.zeros(iter_restart+1,numarray.Float64)
   v=[]

   r_dot_r = bilinearform(r, r)
   if r_dot_r<0: raise NegativeNorm,"negative norm."
   rho=math.sqrt(r_dot_r)

   v.append(r/rho)
   g[0]=rho

   if verbose: print "GMRES: initial residual %e (absolute tolerance = %e)"%(rho,atol)
   while not (rho<=atol or iter==iter_restart):

	if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

	p=Aprod(v[iter])
	v.append(p)

	v_norm1=math.sqrt(bilinearform(v[iter+1], v[iter+1]))

# Modified Gram-Schmidt
	for j in range(iter+1):
	  h[j,iter]=bilinearform(v[j],v[iter+1])
	  v[iter+1]-=h[j,iter]*v[j]

	h[iter+1,iter]=math.sqrt(bilinearform(v[iter+1],v[iter+1]))
	v_norm2=h[iter+1,iter]

# Reorthogonalize if needed
	if v_norm1 + 0.001*v_norm2 == v_norm1:   #Brown/Hindmarsh condition (default)
   	 for j in range(iter+1):
	    hr=bilinearform(v[j],v[iter+1])
      	    h[j,iter]=h[j,iter]+hr
      	    v[iter+1] -= hr*v[j]

   	 v_norm2=math.sqrt(bilinearform(v[iter+1], v[iter+1]))
	 h[iter+1,iter]=v_norm2

#   watch out for happy breakdown
        if not v_norm2 == 0:
         v[iter+1]=v[iter+1]/h[iter+1,iter]

#   Form and store the information for the new Givens rotation
	if iter > 0: h[:iter+1,iter]=__givapp(c[:iter],s[:iter],h[:iter+1,iter])
	mu=math.sqrt(h[iter,iter]*h[iter,iter]+h[iter+1,iter]*h[iter+1,iter])

	if mu!=0 :
		c[iter]=h[iter,iter]/mu
		s[iter]=-h[iter+1,iter]/mu
		h[iter,iter]=c[iter]*h[iter,iter]-s[iter]*h[iter+1,iter]
		h[iter+1,iter]=0.0
		g[iter:iter+2]=__givapp(c[iter],s[iter],g[iter:iter+2])
# Update the residual norm

        rho=abs(g[iter+1])
        if verbose: print "GMRES: iteration step %s: residual %e"%(iter,rho)
	iter+=1

# At this point either iter > iter_max or rho < tol.
# It's time to compute x and leave.

   if verbose: print "GMRES: iteration stopped after %s step."%iter
   if iter > 0 :
     y=numarray.zeros(iter,numarray.Float64)
     y[iter-1] = g[iter-1] / h[iter-1,iter-1]
     if iter > 1 :
        i=iter-2
        while i>=0 :
          y[i] = ( g[i] - numarray.dot(h[i,i+1:iter], y[i+1:iter])) / h[i,i]
          i=i-1
     xhat=v[iter-1]*y[iter-1]
     for i in range(iter-1):
	xhat += v[i]*y[i]
   else:
     xhat=v[0] * 0

   x += xhat
   if iter<iter_restart-1:
      stopped=True
   else:
      stopped=False

   return x,stopped

def MINRES(r, Aprod, x, Msolve, bilinearform, atol=0, rtol=1.e-8, iter_max=100):
    """
    Solver for

    M{Ax=b}

    with a symmetric and positive definite operator A (more details required!).
    It uses the minimum residual method (MINRES) with preconditioner M
    providing an approximation of A.

    The iteration is terminated if

    M{|r| <= atol+rtol*|r0|}

    where M{r0} is the initial residual and M{|.|} is the energy norm. In fact

    M{|r| = sqrt( bilinearform(Msolve(r),r))}

    For details on the preconditioned conjugate gradient method see the book:

    I{Templates for the Solution of Linear Systems by R. Barrett, M. Berry,
    T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo,
    C. Romine, and H. van der Vorst}.

    @param r: initial residual M{r=b-Ax}. C{r} is altered.
    @type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
    @param x: an initial guess for the solution
    @type x: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
    @param Aprod: returns the value Ax
    @type Aprod: function C{Aprod(x)} where C{x} is of the same object like
                 argument C{x}. The returned object needs to be of the same
                 type like argument C{r}.
    @param Msolve: solves Mx=r
    @type Msolve: function C{Msolve(r)} where C{r} is of the same type like
                  argument C{r}. The returned object needs to be of the same
                  type like argument C{x}.
    @param bilinearform: inner product C{<x,r>}
    @type bilinearform: function C{bilinearform(x,r)} where C{x} is of the same
                        type like argument C{x} and C{r} is. The returned value
                        is a C{float}.
    @param atol: absolute tolerance
    @type atol: non-negative C{float}
    @param rtol: relative tolerance
    @type rtol: non-negative C{float}
    @param iter_max: maximum number of iteration steps
    @type iter_max: C{int}
    @return: the solution approximation and the corresponding residual
    @rtype: C{tuple}
    @warning: C{r} and C{x} are altered.
    """
    #------------------------------------------------------------------
    # Set up y and v for the first Lanczos vector v1.
    # y  =  beta1 P' v1,  where  P = C**(-1).
    # v is really P' v1.
    #------------------------------------------------------------------
    r1    = r
    y = Msolve(r)
    beta1 = bilinearform(y,r)

    if beta1< 0: raise NegativeNorm,"negative norm."

    #  If r = 0 exactly, stop with x
    if beta1==0: return x

    if beta1> 0: beta1  = math.sqrt(beta1)

    #------------------------------------------------------------------
    # Initialize quantities.
    # ------------------------------------------------------------------
    iter   = 0
    Anorm = 0
    ynorm = 0
    oldb   = 0
    beta   = beta1
    dbar   = 0
    epsln  = 0
    phibar = beta1
    rhs1   = beta1
    rhs2   = 0
    rnorm  = phibar
    tnorm2 = 0
    ynorm2 = 0
    cs     = -1
    sn     = 0
    w      = r*0.
    w2     = r*0.
    r2     = r1
    eps    = 0.0001

    #---------------------------------------------------------------------
    # Main iteration loop.
    # --------------------------------------------------------------------
    while not rnorm<=atol+rtol*Anorm*ynorm:    #  checks ||r|| < (||A|| ||x||) * TOL

	if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max
        iter    = iter  +  1

        #-----------------------------------------------------------------
        # Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
        # The general iteration is similar to the case k = 1 with v0 = 0:
        #
        #   p1      = Operator * v1  -  beta1 * v0,
        #   alpha1  = v1'p1,
        #   q2      = p2  -  alpha1 * v1,
        #   beta2^2 = q2'q2,
        #   v2      = (1/beta2) q2.
        #
        # Again, y = betak P vk,  where  P = C**(-1).
        #-----------------------------------------------------------------
        s = 1/beta                 # Normalize previous vector (in y).
        v = s*y                    # v = vk if P = I

        y      = Aprod(v)

        if iter >= 2:
          y = y - (beta/oldb)*r1

        alfa   = bilinearform(v,y)              # alphak
        y      += (- alfa/beta)*r2
        r1     = r2
        r2     = y
        y = Msolve(r2)
        oldb   = beta                           # oldb = betak
        beta   = bilinearform(y,r2)             # beta = betak+1^2
        if beta < 0: raise NegativeNorm,"negative norm."

        beta   = math.sqrt( beta )
        tnorm2 = tnorm2 + alfa*alfa + oldb*oldb + beta*beta

        if iter==1:                 # Initialize a few things.
          gmax   = abs( alfa )      # alpha1
          gmin   = gmax             # alpha1

        # Apply previous rotation Qk-1 to get
        #   [deltak epslnk+1] = [cs  sn][dbark    0   ]
        #   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

        oldeps = epsln
        delta  = cs * dbar  +  sn * alfa  # delta1 = 0         deltak
        gbar   = sn * dbar  -  cs * alfa  # gbar 1 = alfa1     gbar k
        epsln  =               sn * beta  # epsln2 = 0         epslnk+1
        dbar   =            -  cs * beta  # dbar 2 = beta2     dbar k+1

        # Compute the next plane rotation Qk

        gamma  = math.sqrt(gbar*gbar+beta*beta)  # gammak
        gamma  = max(gamma,eps)
        cs     = gbar / gamma             # ck
        sn     = beta / gamma             # sk
        phi    = cs * phibar              # phik
        phibar = sn * phibar              # phibark+1

        # Update  x.

        denom = 1/gamma
        w1    = w2
        w2    = w
        w     = (v - oldeps*w1 - delta*w2) * denom
        x     +=  phi*w

        # Go round again.

        gmax   = max(gmax,gamma)
        gmin   = min(gmin,gamma)
        z      = rhs1 / gamma
        ynorm2 = z*z  +  ynorm2
        rhs1   = rhs2 -  delta*z
        rhs2   =      -  epsln*z

        # Estimate various norms and test for convergence.

        Anorm  = math.sqrt( tnorm2 )
        ynorm  = math.sqrt( ynorm2 )

        rnorm  = phibar

    return x

def TFQMR(r, Aprod, x, bilinearform, atol=0, rtol=1.e-8, iter_max=100):
  """
  Solver for

  M{Ax=b}

  with a general operator A (more details required!).
  It uses the Transpose-Free Quasi-Minimal Residual method (TFQMR).

  The iteration is terminated if

  M{|r| <= atol+rtol*|r0|}

  where M{r0} is the initial residual and M{|.|} is the energy norm. In fact

  M{|r| = sqrt( bilinearform(r,r))}

  @param r: initial residual M{r=b-Ax}. C{r} is altered.
  @type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
  @param x: an initial guess for the solution
  @type x: same like C{r}
  @param Aprod: returns the value Ax
  @type Aprod: function C{Aprod(x)} where C{x} is of the same object like
               argument C{x}. The returned object needs to be of the same type
               like argument C{r}.
  @param bilinearform: inner product C{<x,r>}
  @type bilinearform: function C{bilinearform(x,r)} where C{x} is of the same
                      type like argument C{x} and C{r}. The returned value is
                      a C{float}.
  @param atol: absolute tolerance
  @type atol: non-negative C{float}
  @param rtol: relative tolerance
  @type rtol: non-negative C{float}
  @param iter_max: maximum number of iteration steps
  @type iter_max: C{int}
  @rtype: C{tuple}
  @warning: C{r} and C{x} are altered.
  """
  u1=0
  u2=0
  y1=0
  y2=0

  w = r
  y1 = r
  iter = 0
  d = 0
  v = Aprod(y1)
  u1 = v

  theta = 0.0;
  eta = 0.0;
  rho=bilinearform(r,r)
  if rho < 0: raise NegativeNorm,"negative norm."
  tau = math.sqrt(rho)
  norm_r0=tau
  while tau>atol+rtol*norm_r0:
    if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

    sigma = bilinearform(r,v)
    if sigma == 0.0: raise IterationBreakDown,'TFQMR breakdown, sigma=0'

    alpha = rho / sigma

    for j in range(2):
#
#   Compute y2 and u2 only if you have to
#
      if ( j == 1 ):
        y2 = y1 - alpha * v
        u2 = Aprod(y2)

      m = 2 * (iter+1) - 2 + (j+1)
      if j==0:
         w = w - alpha * u1
         d = y1 + ( theta * theta * eta / alpha ) * d
      if j==1:
         w = w - alpha * u2
         d = y2 + ( theta * theta * eta / alpha ) * d

      theta = math.sqrt(bilinearform(w,w))/ tau
      c = 1.0 / math.sqrt ( 1.0 + theta * theta )
      tau = tau * theta * c
      eta = c * c * alpha
      x = x + eta * d
#
#  Try to terminate the iteration at each pass through the loop
#
    if rho == 0.0: raise IterationBreakDown,'TFQMR breakdown, rho=0'

    rhon = bilinearform(r,w)
    beta = rhon / rho;
    rho = rhon;
    y1 = w + beta * y2;
    u1 = Aprod(y1)
    v = u1 + beta * ( u2 + beta * v )

    iter += 1

  return x


#############################################

class ArithmeticTuple(object):
   """
   Tuple supporting inplace update x+=y and scaling x=a*y where C{x,y} is an
   ArithmeticTuple and C{a} is a float.

   Example of usage::

       from esys.escript import Data
       from numarray import array
       a=Data(...)
       b=array([1.,4.])
       x=ArithmeticTuple(a,b)
       y=5.*x

   """
   def __init__(self,*args):
       """
       Initializes object with elements C{args}.

       @param args: tuple of objects that support inplace add (x+=y) and
                    scaling (x=a*y)
       """
       self.__items=list(args)

   def __len__(self):
       """
       Returns the number of items.

       @return: number of items
       @rtype: C{int}
       """
       return len(self.__items)

   def __getitem__(self,index):
       """
       Returns item at specified position.

       @param index: index of item to be returned
       @type index: C{int}
       @return: item with index C{index}
       """
       return self.__items.__getitem__(index)

   def __mul__(self,other):
       """
       Scales by C{other} from the right.

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise self*other
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       try:
           l=len(other)
           if l!=len(self):
               raise ValueError,"length of arguments don't match."
           for i in range(l): out.append(self[i]*other[i])
       except TypeError:
           for i in range(len(self)): out.append(self[i]*other)
       return ArithmeticTuple(*tuple(out))

   def __rmul__(self,other):
       """
       Scales by C{other} from the left.

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise other*self
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       try:
           l=len(other)
           if l!=len(self):
               raise ValueError,"length of arguments don't match."
           for i in range(l): out.append(other[i]*self[i])
       except TypeError:
           for i in range(len(self)): out.append(other*self[i])
       return ArithmeticTuple(*tuple(out))

   def __div__(self,other):
       """
       Scales by (1/C{other}) from the right.

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise self/other
       @rtype: L{ArithmeticTuple}
       """
       return self*(1/other)

   def __rdiv__(self,other):
       """
       Scales by (1/C{other}) from the left.

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise other/self
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       try:
           l=len(other)
           if l!=len(self):
               raise ValueError,"length of arguments don't match."
           for i in range(l): out.append(other[i]/self[i])
       except TypeError:
           for i in range(len(self)): out.append(other/self[i])
       return ArithmeticTuple(*tuple(out))

   def __iadd__(self,other):
       """
       Inplace addition of C{other} to self.

       @param other: increment
       @type other: C{ArithmeticTuple}
       """
       if len(self) != len(other):
           raise ValueError,"tuple lengths must match."
       for i in range(len(self)):
           self.__items[i]+=other[i]
       return self

   def __add__(self,other):
       """
       Adds C{other} to self.

       @param other: increment
       @type other: C{ArithmeticTuple}
       """
       out=[]
       try:
           l=len(other)
           if l!=len(self):
               raise ValueError,"length of arguments don't match."
           for i in range(l): out.append(self[i]+other[i])
       except TypeError:
           for i in range(len(self)): out.append(self[i]+other)
       return ArithmeticTuple(*tuple(out))

   def __sub__(self,other):
       """
       Subtracts C{other} from self.

       @param other: decrement
       @type other: C{ArithmeticTuple}
       """
       out=[]
       try:
           l=len(other)
           if l!=len(self):
               raise ValueError,"length of arguments don't match."
           for i in range(l): out.append(self[i]-other[i])
       except TypeError:
           for i in range(len(self)): out.append(self[i]-other)
       return ArithmeticTuple(*tuple(out))

   def __isub__(self,other):
       """
       Inplace subtraction of C{other} from self.

       @param other: decrement
       @type other: C{ArithmeticTuple}
       """
       if len(self) != len(other):
           raise ValueError,"tuple length must match."
       for i in range(len(self)):
           self.__items[i]-=other[i]
       return self

   def __neg__(self):
       """
       Negates values.
       """
       out=[]
       for i in range(len(self)):
           out.append(-self[i])
       return ArithmeticTuple(*tuple(out))


class HomogeneousSaddlePointProblem(object):
      """
      This class provides a framework for solving linear homogeneous saddle
      point problems of the form::

          M{Av+B^*p=f}
          M{Bv     =0}

      for the unknowns M{v} and M{p} and given operators M{A} and M{B} and
      given right hand side M{f}. M{B^*} is the adjoint operator of M{B}.
      """
      def __init__(self,**kwargs):
        self.setTolerance()
        self.setAbsoluteTolerance()
        self.setSubProblemTolerance()

      #=============================================================
      def initialize(self):
        """
        Initializes the problem (overwrite).
        """
        pass

      def inner_pBv(self,p,v):
         """
         Returns inner product of element p and Bv (overwrite).

         @param p: a pressure increment
         @param v: a residual
         @return: inner product of element p and Bv
         @rtype: C{float}
         @note: used if PCG is applied.
         """
         raise NotImplementedError,"no inner product for p implemented."

      def inner_p(self,p0,p1):
         """
         Returns inner product of p0 and p1 (overwrite).

         @param p0: a pressure
         @param p1: a pressure
         @return: inner product of p0 and p1
         @rtype: C{float}
         """
         raise NotImplementedError,"no inner product for p implemented."
   
      def norm_v(self,v):
         """
         Returns the norm of v (overwrite).

         @param v: a velovity
         @return: norm of v
         @rtype: non-negative C{float}
         """
         raise NotImplementedError,"no norm of v implemented."


      def getV(self, p, v0):
         """
         return the value for v for a given p (overwrite)

         @param p: a pressure
         @param v0: a initial guess for the value v to return. 
         @return: v given as M{v= A^{-1} (f-B^*p)}
         """
         raise NotImplementedError,"no v calculation implemented."
        
      def norm_Bv(self,v):
        """
        Returns Bv (overwrite).

        @rtype: equal to the type of p
        @note: boundary conditions on p should be zero!
        """
        raise NotImplementedError, "no operator B implemented."

      def solve_AinvBt(self,p):
         """
         Solves M{Av=B^*p} with accuracy L{self.getSubProblemTolerance()}
         (overwrite).

         @param p: a pressure increment
         @return: the solution of M{Av=B^*p} 
         @note: boundary conditions on v should be zero!
         """
         raise NotImplementedError,"no operator A implemented."

      def solve_precB(self,v):
         """
         Provides a preconditioner for M{BA^{-1}B^*} with accuracy
         L{self.getSubProblemTolerance()} (overwrite).

         @rtype: equal to the type of p
         @note: boundary conditions on p should be zero!
         """
         raise NotImplementedError,"no preconditioner for Schur complement implemented."
      #=============================================================
      def __Aprod_PCG(self,p):
          return self.solve_AinvBt(p)

      def __inner_PCG(self,p,v):
         return self.inner_pBv(p,v)

      def __Msolve_PCG(self,v):
          return self.solve_precB(v)
      #=============================================================
      def __Aprod_GMRES(self,p):
          return self.solve_precB(self.solve_AinvBt(p))
      def __inner_GMRES(self,p0,p1):
         return self.inner_p(p0,p1)
      #=============================================================
      def norm_p(self,p):
          """
          calculates the norm of C{p}
          
          @param p: a pressure
          @return: the norm of C{p} using the inner product for pressure
          @rtype: C{float}
          """
          f=self.inner_p(p,p)
          if f<0: raise ValueError,"negative pressure norm."
          return math.sqrt(f)
          

      def solve(self,v,p,max_iter=20, verbose=False, show_details=False, usePCG=True, iter_restart=20, max_correction_steps=10):
         """
         Solves the saddle point problem using initial guesses v and p.

         @param v: initial guess for velocity
         @param p: initial guess for pressure
         @type v: L{Data}
         @type p: L{Data}
         @param usePCG: indicates the usage of the PCG rather than GMRES scheme.
         @param max_iter: maximum number of iteration steps per correction
                          attempt
         @param verbose: if True, shows information on the progress of the
                         saddlepoint problem solver.
         @param show_details: if True, shows details of the sub problem solver
         @param iter_restart: restart the iteration after C{iter_restart} steps
                              (only used if useUzaw=False)
         @type usePCG: C{bool}
         @type max_iter: C{int}
         @type verbose: C{bool}
         @type show_details: C{bool}
         @type iter_restart: C{int}
         @rtype: C{tuple} of L{Data} objects
         """
         self.verbose=verbose
         self.show_details=show_details and self.verbose
         rtol=self.getTolerance()
         atol=self.getAbsoluteTolerance()
         correction_step=0
         converged=False
         while not converged:
              # calculate velocity for current pressure:
              v=self.getV(p,v)
              #
              norm_v=self.norm_v(v)
              norm_Bv=self.norm_Bv(v)
              ATOL=norm_v*rtol+atol
              if self.verbose: print "saddle point solver: norm v= %e, norm_Bv= %e, tolerance = %e."%(norm_v, norm_Bv,ATOL)
              if not ATOL>0: raise ValueError,"overall absolute tolerance needs to be positive."
              if norm_Bv <= ATOL:
                 converged=True
              else:
                 correction_step+=1
                 if correction_step>max_correction_steps:
                      raise CorrectionFailed,"Given up after %d correction steps."%correction_step
                 dp=self.solve_precB(v)
                 if usePCG:
                   norm2=self.inner_pBv(dp,v)
                   if norm2<0: raise ValueError,"negative PCG norm."
                   norm2=math.sqrt(norm2)
                 else:
                   norm2=self.norm_p(dp)
                 ATOL_ITER=ATOL/norm_Bv*norm2
                 if self.verbose: print "saddle point solver: tolerance for solver: %e"%ATOL_ITER
                 if usePCG:
                       p,v0=PCG(v,self.__Aprod_PCG,p,self.__Msolve_PCG,self.__inner_PCG,atol=ATOL_ITER, rtol=0.,iter_max=max_iter, verbose=self.verbose)
                 else:
                       p=GMRES(dp,self.__Aprod_GMRES, p, self.__inner_GMRES,atol=ATOL_ITER, rtol=0.,iter_max=max_iter, iter_restart=iter_restart, verbose=self.verbose)
         if self.verbose: print "saddle point solver: tolerance reached."
 	 return v,p

      #========================================================================
      def setTolerance(self,tolerance=1.e-4):
         """
         Sets the relative tolerance for (v,p).

         @param tolerance: tolerance to be used
         @type tolerance: non-negative C{float}
         """
         if tolerance<0:
             raise ValueError,"tolerance must be positive."
         self.__rtol=tolerance
         self.setSubProblemTolerance()

      def getTolerance(self):
         """
         Returns the relative tolerance.

         @return: relative tolerance
         @rtype: C{float}
         """
         return self.__rtol

      def setAbsoluteTolerance(self,tolerance=0.):
         """
         Sets the absolute tolerance.

         @param tolerance: tolerance to be used
         @type tolerance: non-negative C{float}
         """
         if tolerance<0:
             raise ValueError,"tolerance must be non-negative."
         self.__atol=tolerance

      def getAbsoluteTolerance(self):
         """
         Returns the absolute tolerance.

         @return: absolute tolerance
         @rtype: C{float}
         """
         return self.__atol

      def setSubProblemTolerance(self,rtol=None):
         """
         Sets the relative tolerance to solve the subproblem(s).

         @param rtol: relative tolerence
         @type rtol: positive C{float}
         """
         if rtol == None:
              rtol=max(200.*util.EPSILON,self.getTolerance()**2)
         if rtol<=0:
             raise ValueError,"tolerance must be positive."
         self.__sub_tol=rtol

      def getSubProblemTolerance(self):
         """
         Returns the subproblem reduction factor.

         @return: subproblem reduction factor
         @rtype: C{float}
         """
         return self.__sub_tol

def MaskFromBoundaryTag(domain,*tags):
   """
   Creates a mask on the Solution(domain) function space where the value is
   one for samples that touch the boundary tagged by tags.

   Usage: m=MaskFromBoundaryTag(domain, "left", "right")

   @param domain: domain to be used
   @type domain: L{escript.Domain}
   @param tags: boundary tags
   @type tags: C{str}
   @return: a mask which marks samples that are touching the boundary tagged
            by any of the given tags
   @rtype: L{escript.Data} of rank 0
   """
   pde=linearPDEs.LinearPDE(domain,numEquations=1, numSolutions=1)
   d=escript.Scalar(0.,escript.FunctionOnBoundary(domain))
   for t in tags: d.setTaggedValue(t,1.)
   pde.setValue(y=d)
   return util.whereNonZero(pde.getRightHandSide())

#==============================================================================
class SaddlePointProblem(object):
   """
   This implements a solver for a saddle point problem

   M{f(u,p)=0}
   M{g(u)=0}

   for u and p. The problem is solved with an inexact Uszawa scheme for p:

   M{Q_f (u^{k+1}-u^{k}) = - f(u^{k},p^{k})}
   M{Q_g (p^{k+1}-p^{k}) =   g(u^{k+1})}

   where Q_f is an approximation of the Jacobian A_f of f with respect to u and
   Q_f is an approximation of A_g A_f^{-1} A_g with A_g is the Jacobian of g
   with respect to p. As a the construction of a 'proper' Q_g can be difficult,
   non-linear conjugate gradient method is applied to solve for p, so Q_g plays
   in fact the role of a preconditioner.
   """
   def __init__(self,verbose=False,*args):
       """
       Initializes the problem.

       @param verbose: if True, some information is printed in the process
       @type verbose: C{bool}
       @note: this method may be overwritten by a particular saddle point
              problem
       """
       print "SaddlePointProblem should not be used anymore!"
       if not isinstance(verbose,bool):
            raise TypeError("verbose needs to be of type bool.")
       self.__verbose=verbose
       self.relaxation=1.
       DeprecationWarning("SaddlePointProblem should not be used anymore.")

   def trace(self,text):
       """
       Prints C{text} if verbose has been set.

       @param text: a text message
       @type text: C{str}
       """
       if self.__verbose: print "%s: %s"%(str(self),text)

   def solve_f(self,u,p,tol=1.e-8):
       """
       Solves

       A_f du = f(u,p)

       with tolerance C{tol} and returns du. A_f is Jacobian of f with respect
       to u.

       @param u: current approximation of u
       @type u: L{escript.Data}
       @param p: current approximation of p
       @type p: L{escript.Data}
       @param tol: tolerance expected for du
       @type tol: C{float}
       @return: increment du
       @rtype: L{escript.Data}
       @note: this method has to be overwritten by a particular saddle point
              problem
       """
       pass

   def solve_g(self,u,tol=1.e-8):
       """
       Solves

       Q_g dp = g(u)

       where Q_g is a preconditioner for A_g A_f^{-1} A_g with A_g is the
       Jacobian of g with respect to p.

       @param u: current approximation of u
       @type u: L{escript.Data}
       @param tol: tolerance expected for dp
       @type tol: C{float}
       @return: increment dp
       @rtype: L{escript.Data}
       @note: this method has to be overwritten by a particular saddle point
              problem
       """
       pass

   def inner(self,p0,p1):
       """
       Inner product of p0 and p1 approximating p. Typically this returns
       C{integrate(p0*p1)}.
       @return: inner product of p0 and p1
       @rtype: C{float}
       """
       pass

   subiter_max=3
   def solve(self,u0,p0,tolerance=1.e-6,tolerance_u=None,iter_max=100,accepted_reduction=0.995,relaxation=None):
        """
        Runs the solver.

        @param u0: initial guess for C{u}
        @type u0: L{esys.escript.Data}
        @param p0: initial guess for C{p}
        @type p0: L{esys.escript.Data}
        @param tolerance: tolerance for relative error in C{u} and C{p}
        @type tolerance: positive C{float}
        @param tolerance_u: tolerance for relative error in C{u} if different
                            from C{tolerance}
        @type tolerance_u: positive C{float}
        @param iter_max: maximum number of iteration steps
        @type iter_max: C{int}
        @param accepted_reduction: if the norm g cannot be reduced by
                                   C{accepted_reduction} backtracking to adapt
                                   the relaxation factor. If
                                   C{accepted_reduction=None} no backtracking
                                   is used.
        @type accepted_reduction: positive C{float} or C{None}
        @param relaxation: initial relaxation factor. If C{relaxation==None},
                           the last relaxation factor is used.
        @type relaxation: C{float} or C{None}
        """
        tol=1.e-2
        if tolerance_u==None: tolerance_u=tolerance
        if not relaxation==None: self.relaxation=relaxation
        if accepted_reduction ==None:
              angle_limit=0.
        elif accepted_reduction>=1.:
              angle_limit=0.
        else:
              angle_limit=util.sqrt(1-accepted_reduction**2)
        self.iter=0
        u=u0
        p=p0
        #
        #   initialize things:
        #
        converged=False
        #
        #  start loop:
        #
        #  initial search direction is g
        #
        while not converged :
            if self.iter>iter_max:
                raise ArithmeticError("no convergence after %s steps."%self.iter)
            f_new=self.solve_f(u,p,tol)
            norm_f_new = util.Lsup(f_new)
            u_new=u-f_new
            g_new=self.solve_g(u_new,tol)
            self.iter+=1
            norm_g_new = util.sqrt(self.inner(g_new,g_new))
            if norm_f_new==0. and norm_g_new==0.: return u, p
            if self.iter>1 and not accepted_reduction==None:
               #
               #   did we manage to reduce the norm of G? I
               #   if not we start a backtracking procedure
               #
               # print "new/old norm = ",norm_g_new, norm_g, norm_g_new/norm_g
               if norm_g_new > accepted_reduction * norm_g:
                  sub_iter=0
                  s=self.relaxation
                  d=g
                  g_last=g
                  self.trace("    start substepping: f = %s, g = %s, relaxation = %s."%(norm_f_new, norm_g_new, s))
                  while sub_iter < self.subiter_max and  norm_g_new > accepted_reduction * norm_g:
                     dg= g_new-g_last
                     norm_dg=abs(util.sqrt(self.inner(dg,dg))/self.relaxation)
                     rad=self.inner(g_new,dg)/self.relaxation
                     # print "   ",sub_iter,": rad, norm_dg:",abs(rad), norm_dg*norm_g_new * angle_limit
                     # print "   ",sub_iter,": rad, norm_dg:",rad, norm_dg, norm_g_new, norm_g
                     if abs(rad) < norm_dg*norm_g_new * angle_limit:
                         if sub_iter>0: self.trace("    no further improvements expected from backtracking.")
                         break
                     r=self.relaxation
                     self.relaxation= - rad/norm_dg**2
                     s+=self.relaxation
                     #####
                     # a=g_new+self.relaxation*dg/r
                     # print "predicted new norm = ",util.sqrt(self.inner(a,a)),util.sqrt(self.inner(g_new,g_new)), self.relaxation
                     #####
                     g_last=g_new
                     p+=self.relaxation*d
                     f_new=self.solve_f(u,p,tol)
                     u_new=u-f_new
                     g_new=self.solve_g(u_new,tol)
                     self.iter+=1
                     norm_f_new = util.Lsup(f_new)
                     norm_g_new = util.sqrt(self.inner(g_new,g_new))
                     # print "   ",sub_iter," new g norm",norm_g_new
                     self.trace("    %s th sub-step: f = %s, g = %s, relaxation = %s."%(sub_iter, norm_f_new, norm_g_new, s))
                     #
                     #   can we expect reduction of g?
                     #
                     # u_last=u_new
                     sub_iter+=1
                  self.relaxation=s
            #
            #  check for convergence:
            #
            norm_u_new = util.Lsup(u_new)
            p_new=p+self.relaxation*g_new
            norm_p_new = util.sqrt(self.inner(p_new,p_new))
            self.trace("%s th step: f = %s, f/u = %s, g = %s, g/p = %s, relaxation = %s."%(self.iter, norm_f_new ,norm_f_new/norm_u_new, norm_g_new, norm_g_new/norm_p_new, self.relaxation))

            if self.iter>1:
               dg2=g_new-g
               df2=f_new-f
               norm_dg2=util.sqrt(self.inner(dg2,dg2))
               norm_df2=util.Lsup(df2)
               # print norm_g_new, norm_g, norm_dg, norm_p, tolerance
               tol_eq_g=tolerance*norm_dg2/(norm_g*abs(self.relaxation))*norm_p_new
               tol_eq_f=tolerance_u*norm_df2/norm_f*norm_u_new
               if norm_g_new <= tol_eq_g and norm_f_new <= tol_eq_f:
                   converged=True
            f, norm_f, u, norm_u, g, norm_g, p, norm_p = f_new, norm_f_new, u_new, norm_u_new, g_new, norm_g_new, p_new, norm_p_new
        self.trace("convergence after %s steps."%self.iter)
        return u,p

