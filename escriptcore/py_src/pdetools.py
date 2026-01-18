##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""
Provides some tools related to PDEs.

Currently includes:
    - Projector - to project a discontinuous function onto a continuous function
    - Locator - to trace values in data objects at a certain location
    - TimeIntegrationManager - to handle extrapolation in time
    - SaddlePointProblem - solver for Saddle point problems using the inexact uszawa scheme

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"
__author__="Lutz Gross, l.gross@uq.edu.au"


from . import escriptcpp as escore
from . import linearPDEs
from . import util
import math
import numpy

class TimeIntegrationManager(object):
  """
  A simple mechanism to manage time dependent values.

  Typical usage is::

     dt=0.1 # time increment
     tm=TimeIntegrationManager(inital_value,p=1)
     while t<1.
         v_guess=tm.extrapolate(dt) # extrapolate to t+dt
         v=...
         tm.checkin(dt,v)
         t+=dt

  :note: currently only p=1 is supported.
  """
  def __init__(self,*inital_values,**kwargs):
     """
     Sets up the value manager where ``inital_values`` are the initial values
     and p is the order used for extrapolation.
     """
     if "p" in kwargs:
            self.__p=kwargs["p"]
     else:
            self.__p=1
     if "time" in kwargs:
            self.__t=kwargs["time"]
     else:
            self.__t=0.
     self.__v_mem=[inital_values]
     self.__order=0
     self.__dt_mem=[]
     self.__num_val=len(inital_values)

  def getTime(self):
      """
      Returns the current time.

      :return: the current time value
      :rtype: ``float``
      """
      return self.__t

  def getValue(self):
      """
      Returns the current value(s).

      :return: the current stored value(s)
      :rtype: single value or ``tuple`` of values
      """
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
      Extrapolates to ``dt`` forward in time.
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


class Projector(object):
  """
  The Projector is a factory which projects a discontinuous function onto a
  continuous function on a given domain.
  """
  def __init__(self, domain, reduce=True, fast=True):
    """
    Creates a continuous function space projector for a domain.

    :param domain: Domain of the projection.
    :type domain: `Domain`
    :param reduce: Flag to reduce projection order
    :type reduce: ``bool``
    :param fast: Flag to use a fast method based on matrix lumping
    :type fast: ``bool``
    """
    self.__pde = linearPDEs.LinearPDE(domain)
    if fast:
        self.__pde.getSolverOptions().setSolverMethod(linearPDEs.SolverOptions.LUMPING)
    self.__pde.setSymmetryOn()
    self.__pde.setReducedOrderTo(reduce)
    self.__pde.setValue(D = 1.)
    return
  def getSolverOptions(self):
    """
    Returns the solver options of the PDE solver.
    
    :rtype: `linearPDEs.SolverOptions`
    """
    return self.__pde.getSolverOptions()

  def getValue(self, input_data):
    """
    Projects ``input_data`` onto a continuous function.

    :param input_data: the data to be projected
    """
    return self(input_data)

  def __call__(self, input_data):
    """
    Projects ``input_data`` onto a continuous function.

    :param input_data: the data to be projected
    """
    out=escore.Data(0.,input_data.getShape(),self.__pde.getFunctionSpaceForSolution())
    self.__pde.setValue(Y = escore.Data(), Y_reduced = escore.Data())
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

class NoPDE(object):
     """
     Solves the following problem for u:

     *kronecker[i,j]*D[j]*u[j]=Y[i]*

     with constraint

     *u[j]=r[j]*  where *q[j]>0*

     where *D*, *Y*, *r* and *q* are given functions of rank 1.

     In the case of scalars this takes the form

     *D*u=Y*

     with constraint

     *u=r* where *q>0*

     where *D*, *Y*, *r* and *q* are given scalar functions.

     The constraint overwrites any other condition.

     :note: This class is similar to the `linearPDEs.LinearPDE` class with
            A=B=C=X=0 but has the intention that all input parameters are given
            in `Solution` or `ReducedSolution`.
     """
     # The whole thing is a bit strange and I blame Rob Woodcock (CSIRO) for
     # this.
     def __init__(self,domain,D=None,Y=None,q=None,r=None):
         """
         Initializes the problem.

         :param domain: domain of the PDE
         :type domain: `Domain`
         :param D: coefficient of the solution
         :type D: ``float``, ``int``, ``numpy.ndarray``, `Data`
         :param Y: right hand side
         :type Y: ``float``, ``int``, ``numpy.ndarray``, `Data`
         :param q: location of constraints
         :type q: ``float``, ``int``, ``numpy.ndarray``, `Data`
         :param r: value of solution at locations of constraints
         :type r: ``float``, ``int``, ``numpy.ndarray``, `Data`
         """
         self.__domain=domain
         self.__D=D
         self.__Y=Y
         self.__q=q
         self.__r=r
         self.__u=None
         self.__function_space=escore.Solution(self.__domain)

     def setReducedOn(self):
         """
         Sets the `FunctionSpace` of the solution to `ReducedSolution`.
         """
         self.__function_space=escore.ReducedSolution(self.__domain)
         self.__u=None

     def setReducedOff(self):
         """
         Sets the `FunctionSpace` of the solution to `Solution`.
         """
         self.__function_space=escore.Solution(self.__domain)
         self.__u=None

     def setValue(self,D=None,Y=None,q=None,r=None):
         """
         Assigns values to the parameters.

         :param D: coefficient of the solution
         :type D: ``float``, ``int``, ``numpy.ndarray``, `Data`
         :param Y: right hand side
         :type Y: ``float``, ``int``, ``numpy.ndarray``, `Data`
         :param q: location of constraints
         :type q: ``float``, ``int``, ``numpy.ndarray``, `Data`
         :param r: value of solution at locations of constraints
         :type r: ``float``, ``int``, ``numpy.ndarray``, `Data`
         """
         if not D is None:
            self.__D=D
            self.__u=None
         if not Y is None:
            self.__Y=Y
            self.__u=None
         if not q is None:
            self.__q=q
            self.__u=None
         if not r is None:
            self.__r=r
            self.__u=None

     def getSolution(self):
         """
         Returns the solution.

         :return: the solution of the problem
         :rtype: `Data` object in the `FunctionSpace` `Solution` or
                 `ReducedSolution`
         """
         if self.__u is None:
            if self.__D is None:
               raise ValueError("coefficient D is undefined")
            D=escore.Data(self.__D,self.__function_space)
            if D.getRank()>1:
               raise ValueError("coefficient D must have rank 0 or 1")
            if self.__Y is None:
               self.__u=escore.Data(0.,D.getShape(),self.__function_space)
            else:
               self.__u=1./D*self.__Y
            if not self.__q is None:
                q=util.wherePositive(escore.Data(self.__q,self.__function_space))
                self.__u*=(1.-q)
                if not self.__r is None: self.__u+=q*self.__r
         return self.__u

class Locator(object):
     """
     Locator provides access to the values of data objects at a given spatial
     coordinate x.

     In fact, a Locator object finds the sample in the set of samples of a
     given function space or domain which is closest to the given point x.
     """

     def __init__(self,where,x=numpy.zeros((3,))):
       """
       Initializes a Locator to access values in Data objects on the Doamin
       or FunctionSpace for the sample point which is closest to the given
       point x.

       :param where: function space
       :type where: `escript.FunctionSpace`
       :param x: location(s) of the Locator
       :type x: ``numpy.ndarray`` or ``list`` of ``numpy.ndarray``
       """
       if isinstance(where,escore.FunctionSpace):
          self.__function_space=where
       else:
          self.__function_space=escore.ContinuousFunction(where)
       iterative=False
       if isinstance(x, list):
           if len(x)==0: 
              raise ValueError("At least one point must be given.")
           try:
             iter(x[0])
             iterative=True
           except TypeError:
             iterative=False
       xxx=self.__function_space.getX()
       if iterative:
           self.__id=[]
           for p in x:
              self.__id.append(util.length(xxx-p[:self.__function_space.getDim()]).internal_minGlobalDataPoint())
       else:
           self.__id=util.length(xxx-x[:self.__function_space.getDim()]).internal_minGlobalDataPoint()

     def __str__(self):
       """
       Returns the coordinates of the Locator as a string.
       """
       x=self.getX()
       if isinstance(x,list):
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
        if item is None:
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
        Returns the value of ``data`` at the Locator if ``data`` is a `Data`
        object otherwise the object is returned.
        """
        if isinstance(data,escore.Data):
           dat=util.interpolate(data,self.getFunctionSpace())
           ii=self.getId()
           r=data.getRank()
           if isinstance(ii,list):
               out=[]
               for i in ii:
                  o=numpy.array(dat.getTupleForGlobalDataPoint(*i))
                  if data.getRank()==0:
                     out.append(o[0])
                  else:
                     out.append(o)
               return out
           else:
             out=numpy.array(dat.getTupleForGlobalDataPoint(*ii))
             if data.getRank()==0:
                return out[0]
             else:
                return out
        else:
           return data
           
     def setValue(self, data, v):
      """
      Sets the value of the ``data`` at the Locator.
      """
      if isinstance(data, escore.Data):
         if data.getFunctionSpace()!=self.getFunctionSpace():
           raise TypeError("setValue: FunctionSpace of Locator and Data object must match.")
         data.expand()  
         ii=self.getId()
         if isinstance(ii, list):
          for i in ii:
           data._setTupleForGlobalDataPoint(i[1], i[0], v)
         else:
           data._setTupleForGlobalDataPoint(ii[1], ii[0], v)
      else:
           raise TypeError("setValue: Invalid argument type.")


def getInfLocator(arg):
    """
    Return a Locator for a point with the inf value over all arg.
    """
    if not isinstance(arg, escore.Data):
       raise TypeError("getInfLocator: Unknown argument type.")
    a_inf=util.inf(arg)
    loc=util.length(arg-a_inf).internal_minGlobalDataPoint()     # This gives us the location but not coords
    x=arg.getFunctionSpace().getX()
    x_min=x.getTupleForGlobalDataPoint(*loc)
    return Locator(arg.getFunctionSpace(),x_min)

def getSupLocator(arg):
    """
    Return a Locator for a point with the sup value over all arg.
    """
    if not isinstance(arg, escore.Data):
       raise TypeError("getSupLocator: Unknown argument type.")
    a_inf=util.sup(arg)
    loc=util.length(arg-a_inf).internal_minGlobalDataPoint()     # This gives us the location but not coords
    x=arg.getFunctionSpace().getX()
    x_min=x.getTupleForGlobalDataPoint(*loc)
    return Locator(arg.getFunctionSpace(),x_min)
        

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

   *Ax=b*

   with a symmetric and positive definite operator A (more details required!).
   It uses the conjugate gradient method with preconditioner M providing an
   approximation of A.

   The iteration is terminated if

   *|r| <= atol+rtol*|r0|*

   where *r0* is the initial residual and *|.|* is the energy norm. In fact

   *|r| = sqrt( bilinearform(Msolve(r),r))*

   For details on the preconditioned conjugate gradient method see the book:

   "Templates for the Solution of Linear Systems by R. Barrett, M. Berry,
   T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo,
   C. Romine, and H. van der Vorst".

   :param r: initial residual *r=b-Ax*. ``r`` is altered.
   :type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   :param x: an initial guess for the solution
   :type x: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   :param Aprod: returns the value Ax
   :type Aprod: function ``Aprod(x)`` where ``x`` is of the same object like
                argument ``x``. The returned object needs to be of the same type
                like argument ``r``.
   :param Msolve: solves Mx=r
   :type Msolve: function ``Msolve(r)`` where ``r`` is of the same type like
                 argument ``r``. The returned object needs to be of the same
                 type like argument ``x``.
   :param bilinearform: inner product ``<x,r>``
   :type bilinearform: function ``bilinearform(x,r)`` where ``x`` is of the same
                       type like argument ``x`` and ``r`` is. The returned value
                       is a ``float``.
   :param atol: absolute tolerance
   :type atol: non-negative ``float``
   :param rtol: relative tolerance
   :type rtol: non-negative ``float``
   :param iter_max: maximum number of iteration steps
   :type iter_max: ``int``
   :return: the solution approximation and the corresponding residual
   :rtype: ``tuple``
   :warning: ``r`` and ``x`` are altered.
   """
   iter=0
   rhat=Msolve(r)
   d = rhat
   rhat_dot_r = bilinearform(rhat, r)
   if rhat_dot_r<0: raise NegativeNorm("negative norm.")
   norm_r0=math.sqrt(rhat_dot_r)
   atol2=atol+rtol*norm_r0
   if atol2<=0:
      raise ValueError("Non-positive tolerance.")
   atol2=max(atol2, 100. * util.EPSILON * norm_r0)

   if verbose: print(("PCG: initial residual norm = %e (absolute tolerance = %e)"%(norm_r0, atol2)))


   while not math.sqrt(rhat_dot_r) <= atol2:
       iter+=1
       if iter  >= iter_max: raise MaxIterReached("maximum number of %s steps reached."%iter_max)

       q=Aprod(d)
       alpha = rhat_dot_r / bilinearform(d, q)
       x += alpha * d
       if isinstance(q,ArithmeticTuple):
          r += q * (-alpha)      # Doing it the other way calls the float64.__mul__ not AT.__rmul__ 
       else:
           r += (-alpha) * q
       rhat=Msolve(r)
       rhat_dot_r_new = bilinearform(rhat, r)
       beta = rhat_dot_r_new / rhat_dot_r
       rhat+=beta * d
       d=rhat

       rhat_dot_r = rhat_dot_r_new
       if rhat_dot_r<0: raise NegativeNorm("negative norm.")
       if verbose: print(("PCG: iteration step %s: residual norm = %e"%(iter, math.sqrt(rhat_dot_r))))
   if verbose: print(("PCG: tolerance reached after %s steps."%iter))
   return x,r,math.sqrt(rhat_dot_r)

class Defect(object):
    """
    Defines a non-linear defect F(x) of a variable x. This class includes
    two functions (bilinearform and eval) that must be overridden by subclassing
    before use.
    """
    def __init__(self):
        """
        Initializes defect.
        """
        self.setDerivativeIncrementLength()

    def bilinearform(self, x0, x1):
        """
        Returns the inner product of x0 and x1
        
        NOTE: MUST BE OVERRIDDEN BY A SUBCLASS

        :param x0: value for x0
        :param x1: value for x1
        :return: the inner product of x0 and x1
        :rtype: ``float``
        """
        raise NotImplementedError("Defect bilinearform method not overridden")

    def norm(self,x):
        """
        Returns the norm of argument ``x``.

        :param x: a value
        :return: norm of argument x
        :rtype: ``float``
        :note: by default ``sqrt(self.bilinearform(x,x)`` is returned.
        """
        s=self.bilinearform(x,x)
        if s<0: raise NegativeNorm("negative norm.")
        return math.sqrt(s)

    def eval(self,x):
        """
        Returns the value F of a given ``x``.

        NOTE: MUST BE OVERRIDDEN BY A SUBCLASS

        :param x: value for which the defect ``F`` is evaluated
        :return: value of the defect at ``x``
        """
        raise NotImplementedError("Defect eval() method not overridden")

    def __call__(self,x):
        return self.eval(x)

    def setDerivativeIncrementLength(self,inc=1000.*math.sqrt(util.EPSILON)):
        """
        Sets the relative length of the increment used to approximate the
        derivative of the defect. The increment is inc*norm(x)/norm(v)*v in the
        direction of v with x as a starting point.

        :param inc: relative increment length
        :type inc: positive ``float``
        """
        if inc<=0: raise ValueError("positive increment required.")
        self.__inc=inc

    def getDerivativeIncrementLength(self):
        """
        Returns the relative increment length used to approximate the
        derivative of the defect.
        :return: value of the defect at ``x``
        :rtype: positive ``float``
        """
        return self.__inc

    def derivative(self, F0, x0, v, v_is_normalised=True):
        """
        Returns the directional derivative at ``x0`` in the direction of ``v``.

        :param F0: value of this defect at x0
        :param x0: value at which derivative is calculated
        :param v: direction
        :param v_is_normalised: True to indicate that ``v`` is nomalized
                                (self.norm(v)=0)
        :return: derivative of this defect at x0 in the direction of ``v``
        :note: by default numerical evaluation (self.eval(x0+eps*v)-F0)/eps is
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
def NewtonGMRES(defect, x, iter_max=100, sub_iter_max=20, atol=0,rtol=1.e-4, subtol_max=0.5, gamma=0.9, verbose=False):
   """
   Solves a non-linear problem *F(x)=0* for unknown *x* using the stopping
   criterion:

   *norm(F(x) <= atol + rtol * norm(F(x0)*

   where *x0* is the initial guess.

   :param defect: object defining the function *F*. ``defect.norm`` defines the
                  *norm* used in the stopping criterion.
   :type defect: `Defect`
   :param x: initial guess for the solution, ``x`` is altered.
   :type x: any object type allowing basic operations such as
            ``numpy.ndarray``, `Data`
   :param iter_max: maximum number of iteration steps
   :type iter_max: positive ``int``
   :param sub_iter_max: maximum number of inner iteration steps
   :type sub_iter_max: positive ``int``
   :param atol: absolute tolerance for the solution
   :type atol: positive ``float``
   :param rtol: relative tolerance for the solution
   :type rtol: positive ``float``
   :param gamma: tolerance safety factor for inner iteration
   :type gamma: positive ``float``, less than 1
   :param subtol_max: upper bound for inner tolerance
   :type subtol_max: positive ``float``, less than 1
   :return: an approximation of the solution with the desired accuracy
   :rtype: same type as the initial guess
   """
   lmaxit=iter_max
   if atol<0: raise ValueError("atol needs to be non-negative.")
   if rtol<0: raise ValueError("rtol needs to be non-negative.")
   if rtol+atol<=0: raise ValueError("rtol or atol needs to be non-negative.")
   if gamma<=0 or gamma>=1: raise ValueError("tolerance safety factor for inner iteration (gamma =%s) needs to be positive and less than 1."%gamma)
   if subtol_max<=0 or subtol_max>=1: raise ValueError("upper bound for inner tolerance for inner iteration (subtol_max =%s) needs to be positive and less than 1."%subtol_max)

   F=defect(x)
   fnrm=defect.norm(F)
   stop_tol=atol + rtol*fnrm
   subtol=subtol_max
   if verbose: print(("NewtonGMRES: initial residual = %e."%fnrm))
   if verbose: print(("             tolerance = %e."%subtol))
   iter=1
   #
   # main iteration loop
   #
   while not fnrm<=stop_tol:
            if iter  >= iter_max: raise MaxIterReached("maximum number of %s steps reached."%iter_max)
            #
            #   adjust subtol_
            #
            if iter > 1:
               rat=fnrm/fnrmo
               subtol_old=subtol
               subtol=gamma*rat**2
               if gamma*subtol_old**2 > .1: subtol=max(subtol,gamma*subtol_old**2)
               subtol=max(min(subtol,subtol_max), .5*stop_tol/fnrm)
            #
            # calculate newton increment xc
            #     if iter_max in __FDGMRES is reached MaxIterReached is thrown
            #     if iter_restart -1 is returned as sub_iter
            #     if  atol is reached sub_iter returns the numer of steps performed to get there
            #
            #
            if verbose: print(("             subiteration (GMRES) is called with relative tolerance %e."%subtol))
            try:
               xc, sub_iter=__FDGMRES(F, defect, x, subtol*fnrm, iter_max=iter_max-iter, iter_restart=sub_iter_max)
            except MaxIterReached:
               raise MaxIterReached("maximum number of %s steps reached."%iter_max)
            if sub_iter<0:
               iter+=sub_iter_max
            else:
               iter+=sub_iter
            # ====
            x+=xc
            F=defect(x)
            iter+=1
            fnrmo, fnrm=fnrm, defect.norm(F)
            if verbose: print(("             step %s: residual %e."%(iter,fnrm)))
   if verbose: print(("NewtonGMRES: completed after %s steps."%iter))
   return x

def __givapp(c,s,vin):
    """
    Applies a sequence of Givens rotations (c,s) recursively to the vector
    ``vin``

    :warning: ``vin`` is altered.
    """
    vrot=vin
    if isinstance(c,float):
        vrot=[c*vrot[0]-s*vrot[1],s*vrot[0]+c*vrot[1]]
    else:
        for i in range(len(c)):
            w1=c[i]*vrot[i]-s[i]*vrot[i+1]
            w2=s[i]*vrot[i]+c[i]*vrot[i+1]
            vrot[i]=w1
            vrot[i+1]=w2
    return vrot

def __FDGMRES(F0, defect, x0, atol, iter_max=100, iter_restart=20):
   h=numpy.zeros((iter_restart,iter_restart),numpy.float64)
   c=numpy.zeros(iter_restart,numpy.float64)
   s=numpy.zeros(iter_restart,numpy.float64)
   g=numpy.zeros(iter_restart,numpy.float64)
   v=[]

   rho=defect.norm(F0)
   if rho<=0.: return x0*0

   v.append(-F0/rho)
   g[0]=rho
   iter=0
   while rho > atol and iter<iter_restart-1:
        if iter  >= iter_max:
            raise MaxIterReached("maximum number of %s steps reached."%iter_max)

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
            hhat=numpy.zeros(iter+1,numpy.float64)
            for i in range(iter+1) : hhat[i]=h[i,iter]
            hhat=__givapp(c[0:iter],s[0:iter],hhat)
            for i in range(iter+1) : h[i,iter]=hhat[i]

        mu=math.sqrt(h[iter,iter]*h[iter,iter]+h[iter+1,iter]*h[iter+1,iter])

        if mu!=0 :
            c[iter]=h[iter,iter]/mu
            s[iter]=-h[iter+1,iter]/mu
            h[iter,iter]=c[iter]*h[iter,iter]-s[iter]*h[iter+1,iter]
            h[iter+1,iter]=0.0
            gg=__givapp(c[iter],s[iter],[g[iter],g[iter+1]])
            g[iter]=gg[0] 
            g[iter+1]=gg[1]

        # Update the residual norm
        rho=abs(g[iter+1])
        iter+=1

   # At this point either iter > iter_max or rho < tol.
   # It's time to compute x and leave.
   if iter > 0 :
     y=numpy.zeros(iter,numpy.float64)
     y[iter-1] = g[iter-1] / h[iter-1,iter-1]
     if iter > 1 :
        i=iter-2
        while i>=0 :
          y[i] = ( g[i] - numpy.dot(h[i,i+1:iter], y[i+1:iter])) / h[i,i]
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

def GMRES(r, Aprod, x, bilinearform, atol=0, rtol=1.e-8, iter_max=100, iter_restart=20, verbose=False,P_R=None):
   """
   Solver for

   *Ax=b*

   with a general operator A (more details required!).
   It uses the generalized minimum residual method (GMRES).

   The iteration is terminated if

   *|r| <= atol+rtol*|r0|*

   where *r0* is the initial residual and *|.|* is the energy norm. In fact

   *|r| = sqrt( bilinearform(r,r))*

   :param r: initial residual *r=b-Ax*. ``r`` is altered.
   :type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   :param x: an initial guess for the solution
   :type x: same like ``r``
   :param Aprod: returns the value Ax
   :type Aprod: function ``Aprod(x)`` where ``x`` is of the same object like
                argument ``x``. The returned object needs to be of the same
                type like argument ``r``.
   :param bilinearform: inner product ``<x,r>``
   :type bilinearform: function ``bilinearform(x,r)`` where ``x`` is of the same
                       type like argument ``x`` and ``r``. The returned value is
                       a ``float``.
   :param atol: absolute tolerance
   :type atol: non-negative ``float``
   :param rtol: relative tolerance
   :type rtol: non-negative ``float``
   :param iter_max: maximum number of iteration steps
   :type iter_max: ``int``
   :param iter_restart: in order to save memory the orthogonalization process
                        is terminated after ``iter_restart`` steps and the
                        iteration is restarted.
   :type iter_restart: ``int``
   :return: the solution approximation and the corresponding residual
   :rtype: ``tuple``
   :warning: ``r`` and ``x`` are altered.
   """
   m=iter_restart
   restarted=False
   iter=0
   if rtol>0:
      r_dot_r = bilinearform(r, r)
      if r_dot_r<0: raise NegativeNorm("negative norm.")
      atol2=atol+rtol*math.sqrt(r_dot_r)
      if verbose: print(("GMRES: norm of right hand side = %e (absolute tolerance = %e)"%(math.sqrt(r_dot_r), atol2)))
   else:
      atol2=atol
      if verbose: print(("GMRES: absolute tolerance = %e"%atol2))
   if atol2<=0:
      raise ValueError("Non-positive tolarance.")

   while True:
      if iter  >= iter_max: raise MaxIterReached("maximum number of %s steps reached"%iter_max)
      if restarted:
         r2 = r-Aprod(x-x2)
      else:
         r2=1*r
      x2=x*1.
      x,stopped=_GMRESm(r2, Aprod, x, bilinearform, atol2, iter_max=iter_max-iter, iter_restart=m, verbose=verbose,P_R=P_R)
      iter+=iter_restart
      if stopped: break
      if verbose: print("GMRES: restart.")
      restarted=True
   if verbose: print("GMRES: tolerance has been reached.")
   return x

def _GMRESm(r, Aprod, x, bilinearform, atol, iter_max=100, iter_restart=20, verbose=False, P_R=None):
   iter=0

   h=numpy.zeros((iter_restart+1,iter_restart),numpy.float64)
   c=numpy.zeros(iter_restart,numpy.float64)
   s=numpy.zeros(iter_restart,numpy.float64)
   g=numpy.zeros(iter_restart+1,numpy.float64)
   v=[]

   r_dot_r = bilinearform(r, r)
   if r_dot_r<0: raise NegativeNorm("negative norm.")
   rho=math.sqrt(r_dot_r)

   v.append(r/rho)
   g[0]=rho

   if verbose: print(("GMRES: initial residual %e (absolute tolerance = %e)"%(rho,atol)))
   while not (rho<=atol or iter==iter_restart):

        if iter  >= iter_max: raise MaxIterReached("maximum number of %s steps reached."%iter_max)

        if P_R!=None:
            p=Aprod(P_R(v[iter]))
        else:
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
                gg=__givapp(c[iter],s[iter],[g[iter],g[iter+1]])
                g[iter]=gg[0] 
                g[iter+1]=gg[1]
# Update the residual norm

        rho=abs(g[iter+1])
        if verbose: print("GMRES: iteration step %s: residual %e"%(iter,numpy.linalg.norm(rho)))
        iter+=1

# At this point either iter > iter_max or rho < tol.
# It's time to compute x and leave.

   if verbose: print(("GMRES: iteration stopped after %s step."%iter))
   if iter > 0 :
     y=numpy.zeros(iter,numpy.float64)
     y[iter-1] = g[iter-1] / h[iter-1,iter-1]
     if iter > 1 :
        i=iter-2
        while i>=0 :
          y[i] = ( g[i] - numpy.dot(h[i,i+1:iter], y[i+1:iter])) / h[i,i]
          i=i-1
     xhat=v[iter-1]*y[iter-1]
     for i in range(iter-1):
       xhat += v[i]*y[i]
   else:
     xhat=v[0] * 0
   if P_R!=None:
      x += P_R(xhat)
   else:
      x += xhat
   if iter<iter_restart-1:
      stopped=True
   else:
      stopped=False

   return x,stopped

def MINRES(r, Aprod, x, Msolve, bilinearform, atol=0, rtol=1.e-8, iter_max=100):
    """
    Solver for

    *Ax=b*

    with a symmetric and positive definite operator A (more details required!).
    It uses the minimum residual method (MINRES) with preconditioner M
    providing an approximation of A.

    The iteration is terminated if

    *|r| <= atol+rtol*|r0|*

    where *r0* is the initial residual and *|.|* is the energy norm. In fact

    *|r| = sqrt( bilinearform(Msolve(r),r))*

    For details on the preconditioned conjugate gradient method see the book:

    "Templates for the Solution of Linear Systems by R. Barrett, M. Berry,
    T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo,
    C. Romine, and H. van der Vorst".

    :param r: initial residual *r=b-Ax*. ``r`` is altered.
    :type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
    :param x: an initial guess for the solution
    :type x: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
    :param Aprod: returns the value Ax
    :type Aprod: function ``Aprod(x)`` where ``x`` is of the same object like
                 argument ``x``. The returned object needs to be of the same
                 type like argument ``r``.
    :param Msolve: solves Mx=r
    :type Msolve: function ``Msolve(r)`` where ``r`` is of the same type like
                  argument ``r``. The returned object needs to be of the same
                  type like argument ``x``.
    :param bilinearform: inner product ``<x,r>``
    :type bilinearform: function ``bilinearform(x,r)`` where ``x`` is of the same
                        type like argument ``x`` and ``r`` is. The returned value
                        is a ``float``.
    :param atol: absolute tolerance
    :type atol: non-negative ``float``
    :param rtol: relative tolerance
    :type rtol: non-negative ``float``
    :param iter_max: maximum number of iteration steps
    :type iter_max: ``int``
    :return: the solution approximation and the corresponding residual
    :rtype: ``tuple``
    :warning: ``r`` and ``x`` are altered.
    """
    #------------------------------------------------------------------
    # Set up y and v for the first Lanczos vector v1.
    # y  =  beta1 P' v1,  where  P = C**(-1).
    # v is really P' v1.
    #------------------------------------------------------------------
    r1    = r
    y = Msolve(r)
    beta1 = bilinearform(y,r)

    if beta1< 0: raise NegativeNorm("negative norm.")

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

        if iter  >= iter_max: raise MaxIterReached("maximum number of %s steps reached."%iter_max)
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
        if beta < 0: raise NegativeNorm("negative norm.")

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

  *Ax=b*

  with a general operator A (more details required!).
  It uses the Transpose-Free Quasi-Minimal Residual method (TFQMR).

  The iteration is terminated if

  *|r| <= atol+rtol*|r0|*

  where *r0* is the initial residual and *|.|* is the energy norm. In fact

  *|r| = sqrt( bilinearform(r,r))*

  :param r: initial residual *r=b-Ax*. ``r`` is altered.
  :type r: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
  :param x: an initial guess for the solution
  :type x: same like ``r``
  :param Aprod: returns the value Ax
  :type Aprod: function ``Aprod(x)`` where ``x`` is of the same object like
               argument ``x``. The returned object needs to be of the same type
               like argument ``r``.
  :param bilinearform: inner product ``<x,r>``
  :type bilinearform: function ``bilinearform(x,r)`` where ``x`` is of the same
                      type like argument ``x`` and ``r``. The returned value is
                      a ``float``.
  :param atol: absolute tolerance
  :type atol: non-negative ``float``
  :param rtol: relative tolerance
  :type rtol: non-negative ``float``
  :param iter_max: maximum number of iteration steps
  :type iter_max: ``int``
  :rtype: ``tuple``
  :warning: ``r`` and ``x`` are altered.
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
  if rho < 0: raise NegativeNorm("negative norm.")
  tau = math.sqrt(rho)
  norm_r0=tau
  while tau>atol+rtol*norm_r0:
    if iter  >= iter_max: raise MaxIterReached("maximum number of %s steps reached."%iter_max)

    sigma = bilinearform(r,v)
    if sigma == 0.0: raise IterationBreakDown('TFQMR breakdown, sigma=0')

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
    if rho == 0.0: raise IterationBreakDown('TFQMR breakdown, rho=0')

    rhon = bilinearform(r,w)
    beta = rhon / rho
    rho = rhon
    y1 = w + beta * y2;
    u1 = Aprod(y1)
    v = u1 + beta * ( u2 + beta * v )

    iter += 1

  return x


#############################################

class ArithmeticTuple(object):
   """
   Tuple supporting inplace update x+=y and scaling x=a*y where ``x,y`` is an
   ArithmeticTuple and ``a`` is a float.

   Example of usage::

       from esys.escript import Data
       from numpy import array
       a=eData(...)
       b=array([1.,4.])
       x=ArithmeticTuple(a,b)
       y=5.*x

   """
   def __init__(self,*args):
       """
       Initializes object with elements ``args``.

       :param args: tuple of objects that support inplace add (x+=y) and
                    scaling (x=a*y)
       """
       self.__items=list(args)

   def __len__(self):
       """
       Returns the number of items.

       :return: number of items
       :rtype: ``int``
       """
       return len(self.__items)

   def __getitem__(self,index):
       """
       Returns item at specified position.

       :param index: index of item to be returned
       :type index: ``int``
       :return: item with index ``index``
       """
       return self.__items.__getitem__(index)

   def __mul__(self,other):
       """
       Scales by ``other`` from the right.

       :param other: scaling factor
       :type other: ``float``
       :return: itemwise self*other
       :rtype: `ArithmeticTuple`
       """
       out=[]
       try:
           l=len(other)
           if l!=len(self):
               raise ValueError("length of arguments don't match.")
           for i in range(l): 
                if self.__isEmpty(self[i]) or self.__isEmpty(other[i]):
                    out.append(escore.Data())
                else:
                    out.append(self[i]*other[i])
       except TypeError:
           for i in range(len(self)):  
                if self.__isEmpty(self[i]) or self.__isEmpty(other):
                    out.append(escore.Data())
                else:
                    out.append(self[i]*other)
       return ArithmeticTuple(*tuple(out))

   def __rmul__(self,other):
      """
      Scales by ``other`` from the left.

      :param other: scaling factor
      :type other: ``float``
      :return: itemwise other*self
      :rtype: `ArithmeticTuple`
      """
      out=[]
      try:
          l=len(other)
          if l!=len(self):
              raise ValueError("length of arguments don't match.")
          for i in range(l): 
                if self.__isEmpty(self[i]) or self.__isEmpty(other[i]):
                    out.append(escore.Data())
                else:
                    out.append(other[i]*self[i])
      except TypeError:
          for i in range(len(self)):  
                if self.__isEmpty(self[i]) or self.__isEmpty(other):
                    out.append(escore.Data())
                else:
                    out.append(other*self[i])
      return ArithmeticTuple(*tuple(out))

   def __div__(self,other):
       """
       Scales by (1/``other``) from the right.

       :param other: scaling factor
       :type other: ``float``
       :return: itemwise self/other
       :rtype: `ArithmeticTuple`
       """
       return self*(1/other)

   def __rdiv__(self,other):
      """
      Scales by (1/``other``) from the left.

      :param other: scaling factor
      :type other: ``float``
      :return: itemwise other/self
      :rtype: `ArithmeticTuple`
      """
      out=[]
      try:
          l=len(other)
          if l!=len(self):
              raise ValueError("length of arguments don't match.")
          
          for i in range(l): 
                if self.__isEmpty(self[i]):
                    raise ZeroDivisionError("in component %s"%i)
                else:
                    if self.__isEmpty(other[i]):
                        out.append(escore.Data())
                    else:
                        out.append(other[i]/self[i])
      except TypeError:
          for i in range(len(self)):
                if self.__isEmpty(self[i]):
                    raise ZeroDivisionError("in component %s"%i)
                else:
                    if self.__isEmpty(other):
                        out.append(escore.Data())
                    else:
                        out.append(other/self[i])
      return ArithmeticTuple(*tuple(out))

   def __iadd__(self,other):
      """
      Inplace addition of ``other`` to self.

      :param other: increment
      :type other: ``ArithmeticTuple``
      """
      if len(self) != len(other):
          raise ValueError("tuple lengths must match.")
      for i in range(len(self)):
          if self.__isEmpty(self.__items[i]):
              self.__items[i]=other[i]
          else:
              self.__items[i]+=other[i]
              
      return self

   def __add__(self,other):
      """
      Adds ``other`` to self.

      :param other: increment
      :type other: ``ArithmeticTuple``
      """
      out=[]
      try:
          l=len(other)
          if l!=len(self):
              raise ValueError("length of arguments don't match.")
          for i in range(l): 
                if self.__isEmpty(self[i]):
                    out.append(other[i])
                elif self.__isEmpty(other[i]):
                    out.append(self[i])
                else:
                    out.append(self[i]+other[i])
      except TypeError:
            for i in range(len(self)):     
                if self.__isEmpty(self[i]):
                    out.append(other)
                elif self.__isEmpty(other):
                    out.append(self[i])
                else:
                    out.append(self[i]+other)
      return ArithmeticTuple(*tuple(out))

   def __sub__(self,other):
      """
      Subtracts ``other`` from self.

      :param other: decrement
      :type other: ``ArithmeticTuple``
      """
      out=[]
      try:
          l=len(other)
          if l!=len(self):
              raise ValueError("length of arguments don't match.")
          for i in range(l): 
                if self.__isEmpty(other[i]):
                    out.append(self[i])
                elif self.__isEmpty(self[i]):
                    out.append(-other[i])
                else:
                    out.append(self[i]-other[i])
      except TypeError:
            for i in range(len(self)):     
                if  self.__isEmpty(other):
                    out.append(self[i])
                elif self.__isEmpty(self[i]):
                    out.append(-other)
                else:
                    out.append(self[i]-other)
                    
      return ArithmeticTuple(*tuple(out))

   def __isub__(self,other):
      """
      Inplace subtraction of ``other`` from self.

      :param other: decrement
      :type other: ``ArithmeticTuple``
      """
      if len(self) != len(other):
          raise ValueError("tuple length must match.")
      for i in range(len(self)):
          if not self.__isEmpty(other[i]):
              if self.__isEmpty(self.__items[i]):
                  self.__items[i]=-other[i]
              else:
                  self.__items[i]=other[i]
      return self

   def __neg__(self):
      """
      Negates values.
      """
      out=[]
      for i in range(len(self)):
          if self.__isEmpty(self[i]):
              out.append(escore.Data())
          else:
              out.append(-self[i])
          
      return ArithmeticTuple(*tuple(out))
   def __isEmpty(self, d):
    if isinstance(d, escore.Data):
        return d.isEmpty()
    else:
        return False
                
   def __str__(self):
    s="("
    for i in self:
      s=s+str(i)+", "
    s=s+")"
    return s

class HomogeneousSaddlePointProblem(object):
      """
      This class provides a framework for solving linear homogeneous saddle
      point problems of the form::

          *Av+B^*p=f*
          *Bv     =0*

      for the unknowns *v* and *p* and given operators *A* and *B* and
      given right hand side *f*. *B^** is the adjoint operator of *B*.
      *A* may depend weakly on *v* and *p*.
      """
      def __init__(self, **kwargs):
        """
        initializes the saddle point problem
        """
        self.resetControlParameters()
        self.setTolerance()
        self.setAbsoluteTolerance()
      def resetControlParameters(self, K_p=1., K_v=1., rtol_max=0.01, rtol_min = 1.e-7, chi_max=0.5, reduction_factor=0.3, theta = 0.1):
         """
         sets a control parameter

         :param K_p: initial value for constant to adjust pressure tolerance
         :type K_p: ``float``
         :param K_v: initial value for constant to adjust velocity tolerance
         :type K_v: ``float``
         :param rtol_max: maximuim relative tolerance used to calculate presssure and velocity increment.
         :type rtol_max: ``float``
         :param chi_max: maximum tolerable converegence rate.
         :type chi_max: ``float``
         :param reduction_factor: reduction factor for adjustment factors.
         :type reduction_factor: ``float``
         """
         self.setControlParameter(K_p, K_v, rtol_max, rtol_min, chi_max, reduction_factor, theta)

      def setControlParameter(self,K_p=None, K_v=None, rtol_max=None, rtol_min=None, chi_max=None, reduction_factor=None, theta=None):
         """
         sets a control parameter


         :param K_p: initial value for constant to adjust pressure tolerance
         :type K_p: ``float``
         :param K_v: initial value for constant to adjust velocity tolerance
         :type K_v: ``float``
         :param rtol_max: maximuim relative tolerance used to calculate presssure and velocity increment.
         :type rtol_max: ``float``
         :param chi_max: maximum tolerable converegence rate.
         :type chi_max: ``float``
         :type reduction_factor: ``float``
         """
         if not K_p is None:
            if K_p<1:
               raise ValueError("K_p need to be greater or equal to 1.")
         else:
            K_p=self.__K_p

         if not K_v is None:
            if K_v<1:
               raise ValueError("K_v need to be greater or equal to 1.")
         else:
            K_v=self.__K_v

         if not rtol_max is None:
            if rtol_max<=0 or rtol_max>=1: 
               raise ValueError("rtol_max needs to be positive and less than 1.")
         else:
            rtol_max=self.__rtol_max

         if not rtol_min is None:
            if rtol_min<=0 or rtol_min>=1: 
               raise ValueError("rtol_min needs to be positive and less than 1.")
         else:
            rtol_min=self.__rtol_min

         if not chi_max is None:
            if chi_max<=0 or chi_max>=1: 
               raise ValueError("chi_max needs to be positive and less than 1.")
         else:
            chi_max = self.__chi_max 

         if not reduction_factor is None:
            if reduction_factor<=0 or reduction_factor>1:
               raise ValueError("reduction_factor need to be between zero and one.")
         else:
            reduction_factor=self.__reduction_factor

         if not theta is None:
            if theta<=0 or theta>1:
               raise ValueError("theta need to be between zero and one.")
         else:
            theta=self.__theta

         if rtol_min>=rtol_max:
             raise ValueError("rtol_max = %e needs to be greater than rtol_min = %e"%(rtol_max,rtol_min))
         self.__chi_max = chi_max
         self.__rtol_max = rtol_max
         self.__K_p = K_p
         self.__K_v = K_v
         self.__reduction_factor = reduction_factor
         self.__theta = theta
         self.__rtol_min=rtol_min

      #=============================================================
      def inner_pBv(self,p,Bv):
         """
         Returns inner product of element p and Bv (overwrite).

         :param p: a pressure increment
         :param Bv: a residual
         :return: inner product of element p and Bv
         :rtype: ``float``
         :note: used if PCG is applied.
         """
         raise NotImplementedError("no inner product for p and Bv implemented.")

      def inner_p(self,p0,p1):
         """
         Returns inner product of p0 and p1 (overwrite).

         :param p0: a pressure
         :param p1: a pressure
         :return: inner product of p0 and p1
         :rtype: ``float``
         """
         raise NotImplementedError("no inner product for p implemented.")
   
      def norm_v(self,v):
         """
         Returns the norm of v (overwrite).

         :param v: a velovity
         :return: norm of v
         :rtype: non-negative ``float``
         """
         raise NotImplementedError("no norm of v implemented.")
      def getDV(self, p, v, tol):
         """
         return a correction to the value for a given v and a given p with accuracy `tol` (overwrite) 

         :param p: pressure
         :param v: pressure
         :return: dv given as *dv= A^{-1} (f-A v-B^*p)*
         :note: Only *A* may depend on *v* and *p*
         """
         raise NotImplementedError("no dv calculation implemented.")

        
      def Bv(self,v, tol):
        """
        Returns Bv with accuracy `tol` (overwrite)

        :rtype: equal to the type of p
        :note: boundary conditions on p should be zero!
        """
        raise NotImplementedError("no operator B implemented.")

      def norm_Bv(self,Bv):
        """
        Returns the norm of Bv (overwrite).

        :rtype: equal to the type of p
        :note: boundary conditions on p should be zero!
        """
        raise NotImplementedError("no norm of Bv implemented.")

      def solve_AinvBt(self,dp, tol):
         """
         Solves *A dv=B^*dp* with accuracy `tol`

         :param dp: a pressure increment
         :return: the solution of *A dv=B^*dp* 
         :note: boundary conditions on dv should be zero! *A* is the operator used in ``getDV`` and must not be altered.
         """
         raise NotImplementedError("no operator A implemented.")

      def solve_prec(self,Bv, tol):
         """
         Provides a preconditioner for *(BA^{-1}B^ * )* applied to Bv with accuracy `tol`

         :rtype: equal to the type of p
         :note: boundary conditions on p should be zero!
         """
         raise NotImplementedError("no preconditioner for Schur complement implemented.")
      #=============================================================
      def __Aprod_PCG(self,dp):
          dv=self.solve_AinvBt(dp, self.__subtol)
          return ArithmeticTuple(dv,self.Bv(dv, self.__subtol))

      def __inner_PCG(self,p,r):
         return self.inner_pBv(p,r[1])

      def __Msolve_PCG(self,r):
          return self.solve_prec(r[1], self.__subtol)
      #=============================================================
      def __Aprod_GMRES(self,p):
          return self.solve_prec(self.Bv(self.solve_AinvBt(p, self.__subtol), self.__subtol), self.__subtol)
      def __inner_GMRES(self,p0,p1):
         return self.inner_p(p0,p1)

      #=============================================================
      def norm_p(self,p):
          """
          calculates the norm of ``p``
          
          :param p: a pressure
          :return: the norm of ``p`` using the inner product for pressure
          :rtype: ``float``
          """
          f=self.inner_p(p,p)
          if f<0: raise ValueError("negative pressure norm.")
          return math.sqrt(f)
          
      def solve(self,v,p,max_iter=20, verbose=False, usePCG=True, iter_restart=20, max_correction_steps=10):
         """
         Solves the saddle point problem using initial guesses v and p.

         :param v: initial guess for velocity
         :param p: initial guess for pressure
         :type v: `Data`
         :type p: `Data`
         :param usePCG: indicates the usage of the PCG rather than GMRES scheme.
         :param max_iter: maximum number of iteration steps per correction
                          attempt
         :param verbose: if True, shows information on the progress of the
                         saddlepoint problem solver.
         :param iter_restart: restart the iteration after ``iter_restart`` steps
                              (only used if useUzaw=False)
         :type usePCG: ``bool``
         :type max_iter: ``int``
         :type verbose: ``bool``
         :type iter_restart: ``int``
         :rtype: ``tuple`` of `Data` objects
         :note: typically this method is overwritten by a subclass. It provides a wrapper for the ``_solve`` method.
         """
         return self._solve(v=v,p=p,max_iter=max_iter,verbose=verbose, usePCG=usePCG, iter_restart=iter_restart, max_correction_steps=max_correction_steps)

      def _solve(self,v,p,max_iter=20, verbose=False, usePCG=True, iter_restart=20, max_correction_steps=10):
         """
         see `_solve` method.
         """
         self.verbose=verbose
         rtol=self.getTolerance()
         atol=self.getAbsoluteTolerance()

         K_p=self.__K_p
         K_v=self.__K_v
         correction_step=0
         converged=False
         chi=None
         eps=None

         if self.verbose: print(("HomogeneousSaddlePointProblem: start iteration: rtol= %e, atol=%e"%(rtol, atol)))
         while not converged:

             # get tolerance for velecity increment:
             if chi is None:
                rtol_v=self.__rtol_max 
             else:
                rtol_v=min(chi/K_v,self.__rtol_max)
             rtol_v=max(rtol_v, self.__rtol_min)
             if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: rtol_v= %e"%(correction_step,rtol_v)))
             # get velocity increment:
             dv1=self.getDV(p,v,rtol_v)
             v1=v+dv1
             Bv1=self.Bv(v1, rtol_v)
             norm_Bv1=self.norm_Bv(Bv1)
             norm_dv1=self.norm_v(dv1)
             if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: norm_Bv1 = %e, norm_dv1 = %e"%(correction_step, norm_Bv1, norm_dv1)))
             if norm_dv1*self.__theta < norm_Bv1:
                # get tolerance for pressure increment:
                large_Bv1=True
                if chi is None or eps is None:
                   rtol_p=self.__rtol_max 
                else:
                   rtol_p=min(chi**2*eps/K_p/norm_Bv1, self.__rtol_max)
                self.__subtol=max(rtol_p**2, self.__rtol_min)
                if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: rtol_p= %e"%(correction_step,rtol_p)))
                # now we solve for the pressure increment dp from B*A^{-1}B^* dp = Bv1
                if usePCG:
                    dp,r,a_norm=PCG(ArithmeticTuple(v1,Bv1),self.__Aprod_PCG,0*p,self.__Msolve_PCG,self.__inner_PCG,atol=0, rtol=rtol_p,iter_max=max_iter, verbose=self.verbose)
                    v2=r[0]
                    Bv2=r[1]
                else:
                    # don't use!!!!
                    dp=GMRES(self.solve_prec(Bv1,self.__subtol),self.__Aprod_GMRES, 0*p, self.__inner_GMRES,atol=0, rtol=rtol_p,iter_max=max_iter, iter_restart=iter_restart, verbose=self.verbose)
                    dv2=self.solve_AinvBt(dp, self.__subtol)
                    v2=v1-dv2
                    Bv2=self.Bv(v2, self.__subtol)
                p2=p+dp
             else:
                large_Bv1=False
                v2=v1
                p2=p
             # update business:
             norm_dv2=self.norm_v(v2-v)
             norm_v2=self.norm_v(v2)
             if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: v2 = %e, norm_dv2 = %e"%(correction_step, norm_v2, self.norm_v(v2-v))))
             eps, eps_old = max(norm_Bv1, norm_dv2), eps
             if eps_old is None:
                  chi, chi_old = None, chi
             else:
                  chi, chi_old = min(eps/ eps_old, self.__chi_max), chi
             if eps != None:
                 if chi !=None:
                    if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: convergence rate = %e, correction = %e"%(correction_step,chi, eps)))
                 else:
                    if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: correction = %e"%(correction_step, eps)))
             if eps <= rtol*norm_v2+atol :
                 converged = True
             else:
                 if correction_step>=max_correction_steps:
                      raise CorrectionFailed("Given up after %d correction steps."%correction_step)
                 if chi_old!=None:
                    K_p=max(1,self.__reduction_factor*K_p,(chi-chi_old)/chi_old**2*K_p)
                    K_v=max(1,self.__reduction_factor*K_v,(chi-chi_old)/chi_old**2*K_p)
                    if self.verbose: print(("HomogeneousSaddlePointProblem: step %s: new adjustment factor K = %e"%(correction_step,K_p)))
             correction_step+=1
             v,p =v2, p2
         if self.verbose: print(("HomogeneousSaddlePointProblem: tolerance reached after %s steps."%correction_step))
         return v,p
      #========================================================================
      def setTolerance(self,tolerance=1.e-4):
         """
         Sets the relative tolerance for (v,p).

         :param tolerance: tolerance to be used
         :type tolerance: non-negative ``float``
         """
         if tolerance<0:
             raise ValueError("tolerance must be positive.")
         self.__rtol=tolerance

      def getTolerance(self):
         """
         Returns the relative tolerance.

         :return: relative tolerance
         :rtype: ``float``
         """
         return self.__rtol

      def setAbsoluteTolerance(self,tolerance=0.):
         """
         Sets the absolute tolerance.

         :param tolerance: tolerance to be used
         :type tolerance: non-negative ``float``
         """
         if tolerance<0:
             raise ValueError("tolerance must be non-negative.")
         self.__atol=tolerance

      def getAbsoluteTolerance(self):
         """
         Returns the absolute tolerance.

         :return: absolute tolerance
         :rtype: ``float``
         """
         return self.__atol


def getMaskFromBoundaryTag(domain,*tags):
   """
   Creates a mask on the Solution(domain) function space where the value is
   one for samples that touch the **boundary** tagged by tags.

   Usage: m=getMaskFromBoundaryTag(domain, "left", "right")

   :param domain: domain to be used
   :type domain: `escript.Domain`
   :param tags: boundary tags
   :type tags: ``str``
   :return: a mask which marks samples that are touching the boundary tagged
            by any of the given tags
   :rtype: `escript.Data` of rank 0
   """
   pde=linearPDEs.LinearPDE(domain,numEquations=1, numSolutions=1)
   d=escore.Scalar(0.,escore.FunctionOnBoundary(domain))
   for t in tags: d.setTaggedValue(t,1.)
   pde.setValue(y=d)
   return util.whereNonZero(pde.getRightHandSide())

def MaskFromBoundaryTag(domain,*tags):
   """
    Deprecates: use getMaskFromBoundaryTag instead,

   Creates a mask on the Solution(domain) function space where the value is
   one for samples that touch the **boundary** tagged by tags.

   Usage: m=MaskFromBoundaryTag(domain, "left", "right")

   :param domain: domain to be used
   :type domain: `escript.Domain`
   :param tags: boundary tags
   :type tags: ``str``
   :return: a mask which marks samples that are touching the boundary tagged
            by any of the given tags
   :rtype: `escript.Data` of rank 0
   """
   import warnings
   warnings.warn("MaskFromTag is depreciated, use getMaskFromBoundaryTag")
   return getMaskFromBoundaryTag(domain)

def MaskFromTag(domain,*tags):
   """
     Deprecates: use getMaskFromBoundaryTag instead,

   Creates a mask on the Solution(domain) function space where the value is
   one for samples that touch regions tagged by tags.

   Usage: m=MaskFromTag(domain, "ham")

   :param domain: domain to be used
   :type domain: `escript.Domain`
   :param tags: boundary tags
   :type tags: ``str``
   :return: a mask which marks samples that are touching the boundary tagged
            by any of the given tags
   :rtype: `escript.Data` of rank 0


   """
   import warnings
   warnings.warn("MaskFromTag is depreciated, use getMaskFromBoundaryTag")
   return getMaskFromBoundaryTag(domain)

def getBoundaryValuesFromVolumeTag(domain, **values):
       """
       Creates values as the Solution(domain) function space where the value is
       one for samples that touch regions tagged by tags.

       Usage: m=getBoundaryValuesFromVolumeTag(domain, ham=1, f=6)

       :param domain: domain to be used
       :type domain: `escript.Domain`
       :return: a mask which marks samples that are touching the boundary tagged
                by any of the given tags
       :rtype: `escript.Data` of rank 0
       """
       pde = linearPDEs.LinearPDE(domain, numEquations=1, numSolutions=1)
       out = escore.Scalar(0., escore.FunctionOnBoundary(domain))
       for t, v in values.items():
           d = escore.Scalar(0., escore.Function(domain))
           d.setTaggedValue(t, 1.)
           pde.setValue(Y=d)
           out += v * util.whereZero(
               util.interpolate(util.whereNonZero(pde.getRightHandSide()), escore.FunctionOnBoundary(domain)) - 1.)
       return out



def BoundaryValuesFromVolumeTag(domain,**values):
   """
   depreciated: use getBoundaryMaskFromVolumeTag instead.
   Creates a mask on the Solution(domain) function space where the value is
   one for samples that touch regions tagged by tags.

   Usage: m=getBoundaryValuesFromVolumeTag(domain, ham=1, f=6)

   :param domain: domain to be used
   :type domain: `escript.Domain`
   :return: a mask which marks samples that are touching the boundary tagged
            by any of the given tags
   :rtype: `escript.Data` of rank 0
   """
   import warnings
   warnings.warn("BoundaryValuesFromVolumeTag is depreciated, use getBoundaryValuesFromVolumeTag")
   return getBoundaryValuesFromVolumeTag(domain, **values)

class Wavelet(object):
       """
       place holder for source wavelet
       """
       pass


class Ricker(Wavelet):
    """
    The Ricker Wavelet w=f(t)
    """

    def __init__(self, f_dom=40, t_dom=None):
        """
        Sets up a Ricker wavelet wih dominant frequence `f_dom` and
        center at time `t_dom`. If `t_dom` is not given an estimate
        for suitable `t_dom` is calculated so f(0)~0.

        :note: maximum frequence is about 2 x the dominant frequence.
        """
        drop = 18
        self.__f = f_dom
        self.__f_max = numpy.sqrt(7) * f_dom
        self.__s = math.pi * self.__f
        if t_dom == None:
            t_dom = numpy.sqrt(drop) / self.__s
        self.__t0 = t_dom

    def getCenter(self):
        """
        Return value of wavelet center
        """
        return self.__t0

    def getTimeScale(self):
        """
        Returns the time scale which is the inverse of the largest
        frequence with a significant spectral component.
        """
        return 1 / self.__f_max

    def getValue(self, t):
        """
        get value of wavelet at time `t`
        """
        e2 = (self.__s * (t - self.__t0)) ** 2
        return (1 - 2 * e2) * numpy.exp(-e2)

    def getVelocity(self, t):
        """
        get the velocity f'(t) at time `t`
        """
        e2 = (self.__s * (t - self.__t0)) ** 2
        return (-3 + 2 * e2) * numpy.exp(-e2) * 2 * self.__s ** 2 * (t - self.__t0)

    def getAcceleration(self, t):
        """
        get the acceleration f''(t) at time `t`
        """
        e2 = (self.__s * (t - self.__t0)) ** 2
        return 2 * self.__s ** 2 * (-4 * e2 ** 2 + 12 * e2 - 3) * numpy.exp(-e2)
