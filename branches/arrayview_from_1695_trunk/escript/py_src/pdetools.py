#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
Provides some tools related to PDEs. 

Currently includes:
    - Projector - to project a discontinuous
    - Locator - to trace values in data objects at a certain location
    - TimeIntegrationManager - to handel extraplotion in time
    - SaddlePointProblem - solver for Saddle point problems using the inexact uszawa scheme

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


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
  a simple mechanism to manage time dependend values. 

  typical usage is::

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
     sets up the value manager where inital_value is the initial value and p is order used for extrapolation
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
      adds new values to the manager. the p+1 last value get lost
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
      extrapolates to dt forward in time.
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
  The Projector is a factory which projects a discontiuous function onto a
  continuous function on the a given domain.
  """
  def __init__(self, domain, reduce = True, fast=True):
    """
    Create a continuous function space projector for a domain.

    @param domain: Domain of the projection.
    @param reduce: Flag to reduce projection order (default is True)
    @param fast: Flag to use a fast method based on matrix lumping (default is true)
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
    Projects input_data onto a continuous function

    @param input_data: The input_data to be projected.
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
     solves the following problem for u:

     M{kronecker[i,j]*D[j]*u[j]=Y[i]} 

     with constraint

     M{u[j]=r[j]}  where M{q[j]>0}

     where D, Y, r and q are given functions of rank 1.

     In the case of scalars this takes the form

     M{D*u=Y} 

     with constraint

     M{u=r}  where M{q>0}

     where D, Y, r and q are given scalar functions.

     The constraint is overwriting any other condition.

     @note: This class is similar to the L{linearPDEs.LinearPDE} class with A=B=C=X=0 but has the intention
            that all input parameter are given in L{Solution} or L{ReducedSolution}. The whole
            thing is a bit strange and I blame Robert.Woodcock@csiro.au for this.
     """
     def __init__(self,domain,D=None,Y=None,q=None,r=None):
         """
         initialize the problem

         @param domain: domain of the PDE.
         @type domain: L{Domain}
         @param D: coefficient of the solution. 
         @type D: C{float}, C{int}, L{numarray.NumArray}, L{Data}
         @param Y: right hand side
         @type Y: C{float}, C{int}, L{numarray.NumArray}, L{Data}
         @param q: location of constraints
         @type q: C{float}, C{int}, L{numarray.NumArray}, L{Data}
         @param r: value of solution at locations of constraints
         @type r: C{float}, C{int}, L{numarray.NumArray}, L{Data}
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
         sets the L{FunctionSpace} of the solution to L{ReducedSolution}
         """
         self.__function_space=escript.ReducedSolution(self.__domain)
         self.__u=None

     def setReducedOff(self):
         """
         sets the L{FunctionSpace} of the solution to L{Solution}
         """
         self.__function_space=escript.Solution(self.__domain)
         self.__u=None
         
     def setValue(self,D=None,Y=None,q=None,r=None):
         """
         assigns values to the parameters.

         @param D: coefficient of the solution. 
         @type D: C{float}, C{int}, L{numarray.NumArray}, L{Data}
         @param Y: right hand side
         @type Y: C{float}, C{int}, L{numarray.NumArray}, L{Data}
         @param q: location of constraints
         @type q: C{float}, C{int}, L{numarray.NumArray}, L{Data}
         @param r: value of solution at locations of constraints
         @type r: C{float}, C{int}, L{numarray.NumArray}, L{Data}
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
         returns the solution
        
         @return: the solution of the problem
         @rtype: L{Data} object in the L{FunctionSpace} L{Solution} or L{ReducedSolution}.
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
     Locator provides access to the values of data objects at a given
     spatial coordinate x.  
     
     In fact, a Locator object finds the sample in the set of samples of a
     given function space or domain where which is closest to the given
     point x. 
     """

     def __init__(self,where,x=numarray.zeros((3,))):
       """
       Initializes a Locator to access values in Data objects on the Doamin 
       or FunctionSpace where for the sample point which
       closest to the given point x.

       @param where: function space
       @type where: L{escript.FunctionSpace} 
       @param x: coefficient of the solution. 
       @type x: L{numarray.NumArray} or C{list} of L{numarray.NumArray}
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
	Returns the value of data at the Locator of a Data object otherwise 
	the object is returned.
	"""
        return self.getValue(data)

     def getValue(self,data):
        """
	Returns the value of data at the Locator if data is a Data object 
	otherwise the object is returned.
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
   exceptions thrown by solvers 
   """
   pass

class IndefinitePreconditioner(SolverSchemeException):
   """
   the preconditioner is not positive definite.
   """
   pass
class MaxIterReached(SolverSchemeException):
   """
   maxium number of iteration steps is reached.
   """
   pass
class IterationBreakDown(SolverSchemeException):
   """
   iteration scheme econouters an incurable breakdown.
   """
   pass
class NegativeNorm(SolverSchemeException):
   """
   a norm calculation returns a negative norm.
   """
   pass

class IterationHistory(object):
   """
   The IterationHistory class is used to define a stopping criterium. It keeps track of the 
   residual norms. The stoppingcriterium indicates termination if the residual norm has been reduced by 
   a given tolerance.
   """
   def __init__(self,tolerance=math.sqrt(util.EPSILON),verbose=False):
      """
      Initialization

      @param tolerance: tolerance
      @type tolerance: positive C{float}
      @param verbose: switches on the printing out some information
      @type verbose: C{bool}
      """
      if not tolerance>0.:
          raise ValueError,"tolerance needs to be positive."
      self.tolerance=tolerance
      self.verbose=verbose
      self.history=[]
   def stoppingcriterium(self,norm_r,r,x):
       """
       returns True if the C{norm_r} is C{tolerance}*C{norm_r[0]} where C{norm_r[0]}  is the residual norm at the first call.

       
       @param norm_r: current residual norm
       @type norm_r: non-negative C{float}
       @param r: current residual (not used)
       @param x: current solution approximation (not used)
       @return: C{True} is the stopping criterium is fullfilled. Otherwise C{False} is returned.
       @rtype: C{bool}

       """
       self.history.append(norm_r)
       if self.verbose: print "iter: %s:  inner(rhat,r) = %e"#(len(self.history)-1, self.history[-1])
       return self.history[-1]<=self.tolerance * self.history[0]

   def stoppingcriterium2(self,norm_r,norm_b,solver="GMRES",TOL=None):
       """
       returns True if the C{norm_r} is C{tolerance}*C{norm_b} 

       
       @param norm_r: current residual norm
       @type norm_r: non-negative C{float}
       @param norm_b: norm of right hand side
       @type norm_b: non-negative C{float}
       @return: C{True} is the stopping criterium is fullfilled. Otherwise C{False} is returned.
       @rtype: C{bool}

       """
       if TOL==None:
          TOL=self.tolerance
       self.history.append(norm_r)
       if self.verbose: print "iter: %s:  norm(r) = %e"#(len(self.history)-1, self.history[-1])
       return self.history[-1]<=TOL * norm_b

def PCG(b, Aprod, Msolve, bilinearform, stoppingcriterium, x=None, iter_max=100):
   """
   Solver for 

   M{Ax=b}

   with a symmetric and positive definite operator A (more details required!). 
   It uses the conjugate gradient method with preconditioner M providing an approximation of A. 

   The iteration is terminated if the C{stoppingcriterium} function return C{True}.

   For details on the preconditioned conjugate gradient method see the book:

   Templates for the Solution of Linear Systems by R. Barrett, M. Berry, 
   T.F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo, 
   C. Romine, and H. van der Vorst.

   @param b: the right hand side of the liner system. C{b} is altered.
   @type b: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   @param Aprod: returns the value Ax
   @type Aprod: function C{Aprod(x)} where C{x} is of the same object like argument C{x}. The returned object needs to be of the same type like argument C{b}.
   @param Msolve: solves Mx=r 
   @type Msolve: function C{Msolve(r)} where C{r} is of the same type like argument C{b}. The returned object needs to be of the same 
type like argument C{x}.
   @param bilinearform: inner product C{<x,r>} 
   @type bilinearform: function C{bilinearform(x,r)} where C{x} is of the same type like argument C{x} and C{r} is . The returned value is a C{float}.
   @param stoppingcriterium: function which returns True if a stopping criterium is meet. C{stoppingcriterium} has the arguments C{norm_r}, C{r} and C{x} giving the current norm of the residual (=C{sqrt(bilinearform(Msolve(r),r)}), the current residual and the current solution approximation. C{stoppingcriterium} is called in each iteration step.
   @type stoppingcriterium: function that returns C{True} or C{False}
   @param x: an initial guess for the solution. If no C{x} is given 0*b is used.
   @type x: any object supporting inplace add (x+=y) and scaling (x=scalar*y)
   @param iter_max: maximum number of iteration steps.
   @type iter_max: C{int}
   @return: the solution approximation and the corresponding residual
   @rtype: C{tuple} 
   @warning: C{b} and C{x} are altered.
   """
   iter=0
   if x==None:
      x=0*b
   else:
      b += (-1)*Aprod(x) 
   r=b
   rhat=Msolve(r)
   d = rhat 
   rhat_dot_r = bilinearform(rhat, r)
   if rhat_dot_r<0: raise NegativeNorm,"negative norm."

   while not stoppingcriterium(math.sqrt(rhat_dot_r),r,x):
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

   return x,r


############################
# Added by Artak
#################################3

#Apply a sequence of k Givens rotations, used within gmres codes
# vrot=givapp(c, s, vin, k)
def givapp(c,s,vin):
    vrot=vin # warning: vin is altered!!!!
    if isinstance(c,float):
        vrot=[c*vrot[0]-s*vrot[1],s*vrot[0]+c*vrot[1]]
    else:
        for i in range(len(c)):
            w1=c[i]*vrot[i]-s[i]*vrot[i+1]
	    w2=s[i]*vrot[i]+c[i]*vrot[i+1]
            vrot[i:i+2]=w1,w2
    return vrot

##############################################
def GMRES(b, Aprod, Msolve, bilinearform, stoppingcriterium, x=None, iter_max=100, iter_restart=20):
################################################
   m=iter_restart
   iter=0
   xc=x
   while True:
      if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached"%iter_max
      xc,stopped=GMRESm(b*1, Aprod, Msolve, bilinearform, stoppingcriterium, x=xc*1, iter_max=iter_max-iter, iter_restart=m)
      if stopped: break
      iter+=iter_restart	
   return xc

def GMRESm(b, Aprod, Msolve, bilinearform, stoppingcriterium, x=None, iter_max=100, iter_restart=20):
   iter=0
   r=Msolve(b)
   r_dot_r = bilinearform(r, r)
   if r_dot_r<0: raise NegativeNorm,"negative norm."
   norm_b=math.sqrt(r_dot_r)

   if x==None:
      x=0*b
   else: 
      r=Msolve(b-Aprod(x))
      r_dot_r = bilinearform(r, r)
      if r_dot_r<0: raise NegativeNorm,"negative norm."
   
   h=numarray.zeros((iter_restart,iter_restart),numarray.Float64)
   c=numarray.zeros(iter_restart,numarray.Float64)
   s=numarray.zeros(iter_restart,numarray.Float64)
   g=numarray.zeros(iter_restart,numarray.Float64)
   v=[]

   rho=math.sqrt(r_dot_r)
   
   v.append(r/rho)
   g[0]=rho

   while not (stoppingcriterium(rho,norm_b) or iter==iter_restart-1):

	if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

	p=Msolve(Aprod(v[iter]))

	v.append(p)

	v_norm1=math.sqrt(bilinearform(v[iter+1], v[iter+1]))  

# Modified Gram-Schmidt	
	for j in range(iter+1): 
	  h[j][iter]=bilinearform(v[j],v[iter+1])   
	  v[iter+1]-=h[j][iter]*v[j]
       
	h[iter+1][iter]=math.sqrt(bilinearform(v[iter+1],v[iter+1])) 
	v_norm2=h[iter+1][iter]

# Reorthogonalize if needed
	if v_norm1 + 0.001*v_norm2 == v_norm1:   #Brown/Hindmarsh condition (default)
   	 for j in range(iter+1):  
	    hr=bilinearform(v[j],v[iter+1])
      	    h[j][iter]=h[j][iter]+hr 
      	    v[iter+1] -= hr*v[j]

   	 v_norm2=math.sqrt(bilinearform(v[iter+1], v[iter+1]))  
	 h[iter+1][iter]=v_norm2

#   watch out for happy breakdown 
        if not v_norm2 == 0:
         v[iter+1]=v[iter+1]/h[iter+1][iter]

#   Form and store the information for the new Givens rotation
	if iter > 0 :
		hhat=numarray.zeros(iter+1,numarray.Float64)
		for i in range(iter+1) : hhat[i]=h[i][iter]
		hhat=givapp(c[0:iter],s[0:iter],hhat);
	        for i in range(iter+1) : h[i][iter]=hhat[i]

	mu=math.sqrt(h[iter][iter]*h[iter][iter]+h[iter+1][iter]*h[iter+1][iter])

	if mu!=0 :
		c[iter]=h[iter][iter]/mu
		s[iter]=-h[iter+1][iter]/mu
		h[iter][iter]=c[iter]*h[iter][iter]-s[iter]*h[iter+1][iter]
		h[iter+1][iter]=0.0
		g[iter:iter+2]=givapp(c[iter],s[iter],g[iter:iter+2])

# Update the residual norm
               
        rho=abs(g[iter+1])
	iter+=1

# At this point either iter > iter_max or rho < tol.
# It's time to compute x and leave.        

   if iter > 0 : 
     y=numarray.zeros(iter,numarray.Float64)	
     y[iter-1] = g[iter-1] / h[iter-1][iter-1]
     if iter > 1 :	
        i=iter-2   
        while i>=0 :
          y[i] = ( g[i] - numarray.dot(h[i][i+1:iter], y[i+1:iter])) / h[i][i]
          i=i-1
     xhat=v[iter-1]*y[iter-1]
     for i in range(iter-1):
	xhat += v[i]*y[i]
   else : xhat=v[0] 

   x += xhat
   if iter<iter_restart-1: 
      stopped=True 
   else: 
      stopped=False

   return x,stopped


######################################################
def dirder(x, w, bilinearform, Aprod, Msolve, f0, b ):
######################################################

# DIRDER estimates the directional derivative of a function F.


# Hardwired difference increment.
#
  epsnew = 1.0e-07
#
#  Scale the step.
#
  norm_w=math.sqrt(bilinearform(w,w))
  if norm_w== 0.0:
    return x/x

  epsnew = epsnew / norm_w

  if norm_w > 0.0: 
    epsnew = epsnew * math.sqrt(bilinearform(x,x))
#
#  DEL and F1 could share the same space if storage
#  is more important than clarity.
#

  DEL = x + epsnew * w
  f1 = -Msolve(Aprod(DEL))
  z = ( f1 - f0 ) / epsnew
  return z

######################################################
def FDGMRES(f0, Aprod, Msolve, bilinearform, stoppingcriterium, xc=None, x=None, iter_max=100, iter_restart=20,TOL=None):
######################################################
   b=-f0
   b_dot_b = bilinearform(b, b)
   if b_dot_b<0: raise NegativeNorm,"negative norm."
   norm_b=math.sqrt(b_dot_b)

   r=b

   if x==None:
      x=0*f0
   else:
      r=-dirder(xc,x,bilinearform,Aprod,Msolve,f0,b)-f0   
      
   r_dot_r = bilinearform(r, r)
   if r_dot_r<0: raise NegativeNorm,"negative norm."
   
   h=numarray.zeros((iter_restart,iter_restart),numarray.Float64)
   c=numarray.zeros(iter_restart,numarray.Float64)
   s=numarray.zeros(iter_restart,numarray.Float64)
   g=numarray.zeros(iter_restart,numarray.Float64)
   v=[]

   rho=math.sqrt(r_dot_r)
   
   v.append(r/rho)
   g[0]=rho
   iter=0

   while not (stoppingcriterium(rho,norm_b,solver="FDGMRES",TOL=TOL) or iter==iter_restart-1):

	if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

	
        p=dirder(xc, v[iter], bilinearform,Aprod,Msolve,f0,b)

	v.append(p)

	v_norm1=math.sqrt(bilinearform(v[iter+1], v[iter+1]))  

# Modified Gram-Schmidt	
	for j in range(iter+1):
	  h[j][iter]=bilinearform(v[j],v[iter+1])   
	  v[iter+1]+=(-1.)*h[j][iter]*v[j]
       
	h[iter+1][iter]=math.sqrt(bilinearform(v[iter+1],v[iter+1])) 
	v_norm2=h[iter+1][iter]


# Reorthogonalize if needed
	if v_norm1 + 0.001*v_norm2 == v_norm1:   #Brown/Hindmarsh condition (default)
   	 for j in range(iter+1):
	    hr=bilinearform(v[j],v[iter+1])
      	    h[j][iter]=h[j][iter]+hr #vhat
      	    v[iter+1] +=(-1.)*hr*v[j]

   	 v_norm2=math.sqrt(bilinearform(v[iter+1], v[iter+1]))  
	 h[iter+1][iter]=v_norm2

#   watch out for happy breakdown 
        if v_norm2 != 0:
         v[iter+1]=v[iter+1]/h[iter+1][iter]

#   Form and store the information for the new Givens rotation
	if iter > 0 :
		hhat=[]
		for i in range(iter+1) : hhat.append(h[i][iter])
		hhat=givapp(c[0:iter],s[0:iter],hhat);
	        for i in range(iter+1) : h[i][iter]=hhat[i]

	mu=math.sqrt(h[iter][iter]*h[iter][iter]+h[iter+1][iter]*h[iter+1][iter])
	if mu!=0 :
		c[iter]=h[iter][iter]/mu
		s[iter]=-h[iter+1][iter]/mu
		h[iter][iter]=c[iter]*h[iter][iter]-s[iter]*h[iter+1][iter]
		h[iter+1][iter]=0.0
		g[iter:iter+2]=givapp(c[iter],s[iter],g[iter:iter+2])

# Update the residual norm
        rho=abs(g[iter+1])
	iter+=1

   if iter > 0 : 
     y=numarray.zeros(iter,numarray.Float64)	
     y[iter-1] = g[iter-1] / h[iter-1][iter-1]
     if iter > 1 :	
        i=iter-2   
        while i>=0 :
          y[i] = ( g[i] - numarray.dot(h[i][i+1:iter], y[i+1:iter])) / h[i][i]
          i=i-1
     xhat=v[iter-1]*y[iter-1]
     for i in range(iter-1):
	xhat += v[i]*y[i]
   else : xhat=v[0] 
    
   x += xhat
   if iter<iter_restart-1: 
      stopped=True 
   else: 
      stopped=False

   return x,stopped

#################################################
def MINRES(b, Aprod, Msolve, bilinearform, stoppingcriterium, x=None, iter_max=100):
#################################################
    #
    #  minres solves the system of linear equations Ax = b
    #  where A is a symmetric matrix (possibly indefinite or singular)
    #  and b is a given vector.
    #  
    #  "A" may be a dense or sparse matrix (preferably sparse!)
    #  or the name of a function such that
    #               y = A(x)
    #  returns the product y = Ax for any given vector x.
    #
    #  "M" defines a positive-definite preconditioner M = C C'.
    #  "M" may be a dense or sparse matrix (preferably sparse!)
    #  or the name of a function such that
    #  solves the system My = x for any given vector x.
    #
    #
    
    #------------------------------------------------------------------
    # Set up y and v for the first Lanczos vector v1.
    # y  =  beta1 P' v1,  where  P = C**(-1).
    # v is really P' v1.
    #------------------------------------------------------------------
    if x==None:
      x=0*b
    else:
      b += (-1)*Aprod(x) 

    r1    = b
    y = Msolve(b)
    beta1 = bilinearform(y,b)
 
    if beta1< 0: raise NegativeNorm,"negative norm."

    #  If b = 0 exactly, stop with x = 0.
    if beta1==0: return x*0.

    if beta1> 0:
      beta1  = math.sqrt(beta1)       

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
    w      = b*0.
    w2     = b*0.
    r2     = r1
    eps    = 0.0001

    #---------------------------------------------------------------------
    # Main iteration loop.
    # --------------------------------------------------------------------
    while not stoppingcriterium(rnorm,Anorm*ynorm,'MINRES'):    #  checks ||r|| < (||A|| ||x||) * TOL

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

######################################    
def NewtonGMRES(b, Aprod, Msolve, bilinearform, stoppingcriterium,x=None, iter_max=100,iter_restart=20,atol=1.e-2,rtol=1.e-4):
#####################################
    gamma=.9
    lmaxit=100
    etamax=.5

    n = 1 #len(x)
    iter=0
    
    # evaluate f at the initial iterate
    # compute the stop tolerance
    #
    r=b
    if x==None:
      x=0*b
    else:
      r += (-1)*Aprod(x) 

    f0=-Msolve(r)
    fnrm=math.sqrt(bilinearform(f0,f0))/math.sqrt(n)
    fnrmo=1
    stop_tol=atol + rtol*fnrm
    #
    # main iteration loop
    #
    while not stoppingcriterium(fnrm*1,stop_tol,'NewtonGMRES',TOL=1.):

            if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max 
	    #
	    # keep track of the ratio (rat = fnrm/frnmo)
	    # of successive residual norms and 
	    # the iteration counter (iter)
	    #
	    #rat=fnrm/fnrmo
	    fnrmo=fnrm 
	    iter+=1
	    #
    	    # compute the step using a GMRES(m) routine especially designed for this purpose
	    #
            initer=0
            while True:
               xc,stopped=FDGMRES(f0*1, Aprod, Msolve, bilinearform, stoppingcriterium, xc=x*1, iter_max=lmaxit-initer, iter_restart=iter_restart, TOL=etamax)
               if stopped: break
               initer+=iter_restart
	    x+=xc
	    f0=-Msolve(Aprod(x))
	    fnrm=math.sqrt(bilinearform(f0,f0))/math.sqrt(n)
	    rat=fnrm/fnrmo


	#   adjust eta 
	#
	    if etamax > 0:
	        etaold=etamax
	        etanew=gamma*rat*rat
	        if gamma*etaold*etaold > .1 :
	            etanew=max(etanew,gamma*etaold*etaold)
	        etamax=min(etanew,etamax)
	        etamax=max(etamax,.5*stop_tol/fnrm)
    return x

def TFQMR(b, Aprod, Msolve, bilinearform, stoppingcriterium, x=None, iter_max=100):

# TFQMR solver for linear systems
#
#
# initialization
#
  errtol = math.sqrt(bilinearform(b,b))
  norm_b=errtol
  kmax  = iter_max
  error = []

  if math.sqrt(bilinearform(x,x)) != 0.0:
    r = b - Aprod(x)
  else:
    r = b

  r=Msolve(r)

  u1=0
  u2=0
  y1=0
  y2=0

  w = r
  y1 = r 
  iter = 0 
  d = 0
  
  v = Msolve(Aprod(y1))
  u1 = v
  
  theta = 0.0;
  eta = 0.0;
  tau = math.sqrt(bilinearform(r,r))
  error = [ error, tau ] 
  rho = tau * tau
  m=1
#
#  TFQMR iteration
#
#  while ( iter < kmax-1 ):
   
  while not stoppingcriterium(tau*math.sqrt ( m + 1 ),norm_b,'TFQMR'):
    if iter  >= iter_max: raise MaxIterReached,"maximum number of %s steps reached."%iter_max

    sigma = bilinearform(r,v)

    if ( sigma == 0.0 ):
      raise 'TFQMR breakdown, sigma=0'
    

    alpha = rho / sigma

    for j in range(2):
#
#   Compute y2 and u2 only if you have to
#
      if ( j == 1 ):
        y2 = y1 - alpha * v
        u2 = Msolve(Aprod(y2))
 
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
     # if ( tau * math.sqrt ( m + 1 ) <= errtol ):
     #   error = [ error, tau ]
     #   total_iters = iter
     #   break
      

    if ( rho == 0.0 ):
      raise 'TFQMR breakdown, rho=0'
    

    rhon = bilinearform(r,w)
    beta = rhon / rho;
    rho = rhon;
    y1 = w + beta * y2;
    u1 = Msolve(Aprod(y1))
    v = u1 + beta * ( u2 + beta * v )
    error = [ error, tau ]
    total_iters = iter
    
    iter = iter + 1

  return x


#############################################

class ArithmeticTuple(object):
   """
   tuple supporting inplace update x+=y and scaling x=a*y where x,y is an ArithmeticTuple and a is a float.

   example of usage:

   from esys.escript import Data
   from numarray import array
   a=Data(...)
   b=array([1.,4.])
   x=ArithmeticTuple(a,b)
   y=5.*x

   """
   def __init__(self,*args):
       """
       initialize object with elements args.

       @param args: tuple of object that support implace add (x+=y) and scaling (x=a*y)
       """
       self.__items=list(args)

   def __len__(self):
       """
       number of items

       @return: number of items
       @rtype: C{int}
       """
       return len(self.__items)

   def __getitem__(self,index):
       """
       get an item

       @param index: item to be returned
       @type index: C{int}
       @return: item with index C{index}
       """
       return self.__items.__getitem__(index)

   def __mul__(self,other):
       """
       scaling from the right 

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise self*other
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       other=1.*other
       if isinstance(other,float):
	for i in range(len(self)):
           out.append(self[i]*other)
       else:
        for i in range(len(self)):
           out.append(self[i]*other[i])
       return ArithmeticTuple(*tuple(out))

   def __rmul__(self,other):
       """
       scaling from the left

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise other*self
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       other=1.*other
       if isinstance(other,float):
	for i in range(len(self)):
           out.append(other*self[i])
       else:
        for i in range(len(self)):
           out.append(other[i]*self[i])
       return ArithmeticTuple(*tuple(out))

#########################
# Added by Artak
#########################
   def __div__(self,other):
       """
       dividing from the right 

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise self/other
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       other=1.*other
       if isinstance(other,float):
	for i in range(len(self)):
           out.append(self[i]*(1./other))
       else:
        for i in range(len(self)):
           out.append(self[i]*(1./other[i]))
       return ArithmeticTuple(*tuple(out))

   def __rdiv__(self,other):
       """
       dividing from the left

       @param other: scaling factor
       @type other: C{float}
       @return: itemwise other/self
       @rtype: L{ArithmeticTuple}
       """
       out=[]
       other=1.*other
       if isinstance(other,float):
        for i in range(len(self)):
           out.append(other*(1./self[i]))
       else:
        for i in range(len(self)):
           out.append(other[i]*(1./self[i]))
       return ArithmeticTuple(*tuple(out))
  
##########################################33

   def __iadd__(self,other):
       """
       in-place add of other to self

       @param other: increment
       @type other: C{ArithmeticTuple}
       """
       if len(self) != len(other):
           raise ValueError,"tuple length must match."
       for i in range(len(self)):
           self.__items[i]+=other[i]
       return self

   def __add__(self,other):
       """
       add other to self

       @param other: increment
       @type other: C{ArithmeticTuple}
       """
#      if not isinstance(other,float):
       if len(self) != len(other):
          raise ValueError,"tuple length must match."
       for i in range(len(self)):
          self.__items[i]+=other[i]
#       else:
#        for i in range(len(self)):
#           self.__items[i]+=other

       return self

   def __sub__(self,other):
       """
       subtract other from self

       @param other: increment
       @type other: C{ArithmeticTuple}
       """
       if len(self) != len(other):
           raise ValueError,"tuple length must match."
       for i in range(len(self)):
           self.__items[i]-=other[i]
       return self
   
   def __isub__(self,other):
       """
       subtract other from self

       @param other: increment
       @type other: C{ArithmeticTuple}
       """
       if len(self) != len(other):
           raise ValueError,"tuple length must match."
       for i in range(len(self)):
           self.__items[i]-=other[i]
       return self

   def __neg__(self):
       """
       negate 

       """
       for i in range(len(self)):
           self.__items[i]=-self.__items[i]
       return self


class HomogeneousSaddlePointProblem(object):
      """
      This provides a framwork for solving homogeneous saddle point problem of the form

             Av+B^*p=f
             Bv    =0

      for the unknowns v and p and given operators A and B and given right hand side f.
      B^* is the adjoint operator of B is the given inner product.

      """
      def __init__(self,**kwargs):
        self.setTolerance()
        self.setToleranceReductionFactor()

      def initialize(self):
        """
        initialize the problem (overwrite)
        """
        pass
      def B(self,v):
         """
         returns Bv (overwrite)
         @rtype: equal to the type of p

         @note: boundary conditions on p should be zero!
         """
         pass

      def inner(self,p0,p1):
         """
         returns inner product of two element p0 and p1  (overwrite)
         
         @type p0: equal to the type of p
         @type p1: equal to the type of p
         @rtype: C{float}

         @rtype: equal to the type of p
         """
         pass

      def solve_A(self,u,p):
         """
         solves Av=f-Au-B^*p with accuracy self.getReducedTolerance() (overwrite) 

         @rtype: equal to the type of v
         @note: boundary conditions on v should be zero!
         """
         pass

      def solve_prec(self,p):
         """
         provides a preconditioner for BA^{-1}B^* with accuracy self.getReducedTolerance() (overwrite)

         @rtype: equal to the type of p
         """
         pass

      def stoppingcriterium(self,Bv,v,p):
         """
         returns a True if iteration is terminated. (overwrite)

         @rtype: C{bool}
         """
         pass
             
      def __inner(self,p,r):
         return self.inner(p,r[1])

      def __inner_p(self,p1,p2):
         return self.inner(p1,p2)
      
      def __inner_a(self,a1,a2):
         return self.inner_a(a1,a2)

      def __inner_a1(self,a1,a2):
         return self.inner(a1[1],a2[1])

      def __stoppingcriterium(self,norm_r,r,p):
          return self.stoppingcriterium(r[1],r[0],p)

      def __stoppingcriterium2(self,norm_r,norm_b,solver='GMRES',TOL=None):
          return self.stoppingcriterium2(norm_r,norm_b,solver,TOL)

      def setTolerance(self,tolerance=1.e-8):
              self.__tol=tolerance
      def getTolerance(self):
              return self.__tol
      def setToleranceReductionFactor(self,reduction=0.01):
              self.__reduction=reduction
      def getSubProblemTolerance(self):
              return self.__reduction*self.getTolerance()

      def solve(self,v,p,max_iter=20, verbose=False, show_details=False, solver='PCG',iter_restart=20):
              """
              solves the saddle point problem using initial guesses v and p.

              @param max_iter: maximum number of iteration steps.
              """
              self.verbose=verbose
              self.show_details=show_details and self.verbose

              # assume p is known: then v=A^-1(f-B^*p) 
              # which leads to BA^-1B^*p = BA^-1f  

	      # Az=f is solved as A(z-v)=f-Av (z-v = 0 on fixed_u_mask)	     
	      self.__z=v+self.solve_A(v,p*0)
              Bz=self.B(self.__z) 
              #
	      #   solve BA^-1B^*p = Bz 
              #
              #
              #
              self.iter=0
	      if solver=='GMRES':   	
                if self.verbose: print "enter GMRES method (iter_max=%s)"%max_iter
                p=GMRES(Bz,self.__Aprod2,self.__Msolve2,self.__inner_p,self.__stoppingcriterium2,iter_max=max_iter, x=p*1.,iter_restart=iter_restart)
                # solve Au=f-B^*p 
                #       A(u-v)=f-B^*p-Av
                #       u=v+(u-v)
		u=v+self.solve_A(v,p)

	      if solver=='NewtonGMRES':   	
                if self.verbose: print "enter NewtonGMRES method (iter_max=%s)"%max_iter
                p=NewtonGMRES(Bz,self.__Aprod_Newton,self.__Msolve2,self.__inner_p,self.__stoppingcriterium2,iter_max=max_iter, x=p*1.,atol=0,rtol=self.getTolerance())
                # solve Au=f-B^*p 
                #       A(u-v)=f-B^*p-Av
                #       u=v+(u-v)
		u=v+self.solve_A(v,p)
                

	      if solver=='TFQMR':   	
                if self.verbose: print "enter TFQMR method (iter_max=%s)"%max_iter
                p=TFQMR(Bz,self.__Aprod2,self.__Msolve2,self.__inner_p,self.__stoppingcriterium2,iter_max=max_iter, x=p*1.)
                # solve Au=f-B^*p 
                #       A(u-v)=f-B^*p-Av
                #       u=v+(u-v)
		u=v+self.solve_A(v,p)

	      if solver=='MINRES':   	
                if self.verbose: print "enter MINRES method (iter_max=%s)"%max_iter
                p=MINRES(Bz,self.__Aprod2,self.__Msolve2,self.__inner_p,self.__stoppingcriterium2,iter_max=max_iter, x=p*1.)
                # solve Au=f-B^*p 
                #       A(u-v)=f-B^*p-Av
                #       u=v+(u-v)
		u=v+self.solve_A(v,p)
              
	      if solver=='GMRESC':   	
                if self.verbose: print "enter GMRES coupled method (iter_max=%s)"%max_iter
                p0=self.solve_prec1(Bz)
	        #alfa=(1/self.vol)*util.integrate(util.interpolate(p,escript.Function(self.domain)))
                #p-=alfa
                x=GMRES(ArithmeticTuple(self.__z*1.,p0*1),self.__Anext,self.__Mempty,self.__inner_a,self.__stoppingcriterium2,iter_max=max_iter, x=ArithmeticTuple(v*1,p*1),iter_restart=20)
                #x=NewtonGMRES(ArithmeticTuple(self.__z*1.,p0*1),self.__Aprod_Newton2,self.__Mempty,self.__inner_a,self.__stoppingcriterium2,iter_max=max_iter, x=ArithmeticTuple(v*1,p*1),atol=0,rtol=self.getTolerance())

                # solve Au=f-B^*p 
                #       A(u-v)=f-B^*p-Av
                #       u=v+(u-v)
         	p=x[1]
		u=v+self.solve_A(v,p)		
		#p=x[1]
		#u=x[0]

              if solver=='PCG':
                #   note that the residual r=Bz-BA^-1B^*p = B(z-A^-1B^*p) = Bv
                #
                #   with                    Av=Az-B^*p = f - B^*p (v=z on fixed_u_mask)
                #                           A(v-z)= f -Az - B^*p (v-z=0 on fixed_u_mask)
                if self.verbose: print "enter PCG method (iter_max=%s)"%max_iter
                p,r=PCG(ArithmeticTuple(self.__z*1.,Bz),self.__Aprod,self.__Msolve,self.__inner,self.__stoppingcriterium,iter_max=max_iter, x=p)
	        u=r[0]  
                # print "DIFF=",util.integrate(p)

              # print "RESULT div(u)=",util.Lsup(self.B(u)),util.Lsup(u)

 	      return u,p

      def __Msolve(self,r):
          return self.solve_prec(r[1])

      def __Msolve2(self,r):
          return self.solve_prec(r*1)

      def __Mempty(self,r):
          return r


      def __Aprod(self,p):
          # return BA^-1B*p 
          #solve Av =B^*p as Av =f-Az-B^*(-p)
          v=self.solve_A(self.__z,-p)
          return ArithmeticTuple(v, self.B(v))

      def __Anext(self,x):
          # return next v,p 
          #solve Av =-B^*p as Av =f-Az-B^*p

	  pc=x[1]
          v=self.solve_A(self.__z,-pc)
	  p=self.solve_prec1(self.B(v))

          return ArithmeticTuple(v,p)


      def __Aprod2(self,p):
          # return BA^-1B*p 
          #solve Av =B^*p as Av =f-Az-B^*(-p)
	  v=self.solve_A(self.__z,-p)
          return self.B(v)

      def __Aprod_Newton(self,p):
          # return BA^-1B*p - Bz 
          #solve Av =-B^*p as Av =f-Az-B^*p
	  v=self.solve_A(self.__z,-p)
          return self.B(v-self.__z)

      def __Aprod_Newton2(self,x):
          # return BA^-1B*p - Bz 
          #solve Av =-B^*p as Av =f-Az-B^*p
          pc=x[1]
	  v=self.solve_A(self.__z,-pc)
          p=self.solve_prec1(self.B(v-self.__z))
          return ArithmeticTuple(v,p)

class SaddlePointProblem(object):
   """
   This implements a solver for a saddlepoint problem

   M{f(u,p)=0}
   M{g(u)=0}

   for u and p. The problem is solved with an inexact Uszawa scheme for p:

   M{Q_f (u^{k+1}-u^{k}) = - f(u^{k},p^{k})}
   M{Q_g (p^{k+1}-p^{k}) =   g(u^{k+1})}

   where Q_f is an approximation of the Jacobiean A_f of f with respect to u  and Q_f is an approximation of
   A_g A_f^{-1} A_g with A_g is the jacobiean of g with respect to p. As a the construction of a 'proper'
   Q_g can be difficult, non-linear conjugate gradient method is applied to solve for p, so Q_g plays
   in fact the role of a preconditioner.
   """
   def __init__(self,verbose=False,*args):
       """
       initializes the problem

       @param verbose: switches on the printing out some information
       @type verbose: C{bool}
       @note: this method may be overwritten by a particular saddle point problem
       """
       print "SaddlePointProblem should not be used anymore!"
       if not isinstance(verbose,bool):
            raise TypeError("verbose needs to be of type bool.")
       self.__verbose=verbose
       self.relaxation=1.

   def trace(self,text):
       """
       prints text if verbose has been set

       @param text: a text message
       @type text: C{str}
       """
       if self.__verbose: print "%s: %s"%(str(self),text)

   def solve_f(self,u,p,tol=1.e-8):
       """
       solves 

       A_f du = f(u,p) 

       with tolerance C{tol} and return du. A_f is Jacobiean of f with respect to u.

       @param u: current approximation of u
       @type u: L{escript.Data}
       @param p: current approximation of p
       @type p: L{escript.Data}
       @param tol: tolerance expected for du
       @type tol: C{float}
       @return: increment du
       @rtype: L{escript.Data}
       @note: this method has to be overwritten by a particular saddle point problem
       """
       pass

   def solve_g(self,u,tol=1.e-8):
       """
       solves 

       Q_g dp = g(u) 

       with Q_g is a preconditioner for A_g A_f^{-1} A_g with  A_g is the jacobiean of g with respect to p.

       @param u: current approximation of u
       @type u: L{escript.Data}
       @param tol: tolerance expected for dp
       @type tol: C{float}
       @return: increment dp
       @rtype: L{escript.Data}
       @note: this method has to be overwritten by a particular saddle point problem
       """
       pass

   def inner(self,p0,p1):
       """
       inner product of p0 and p1 approximating p. Typically this returns integrate(p0*p1)
       @return: inner product of p0 and p1
       @rtype: C{float}
       """
       pass

   subiter_max=3
   def solve(self,u0,p0,tolerance=1.e-6,tolerance_u=None,iter_max=100,accepted_reduction=0.995,relaxation=None):
        """
        runs the solver

        @param u0: initial guess for C{u}
        @type u0: L{esys.escript.Data}
        @param p0: initial guess for C{p}
        @type p0: L{esys.escript.Data}
        @param tolerance: tolerance for relative error in C{u} and C{p}
        @type tolerance: positive C{float}
        @param tolerance_u: tolerance for relative error in C{u} if different from C{tolerance}
        @type tolerance_u: positive C{float}
        @param iter_max: maximum number of iteration steps.
        @type iter_max: C{int}
        @param accepted_reduction: if the norm  g cannot be reduced by C{accepted_reduction} backtracking to adapt the 
                                   relaxation factor. If C{accepted_reduction=None} no backtracking is used.
        @type accepted_reduction: positive C{float} or C{None}
        @param relaxation: initial relaxation factor. If C{relaxation==None}, the last relaxation factor is used.
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
#   def solve(self,u0,p0,tolerance=1.e-6,iter_max=10,self.relaxation=1.):
#      tol=1.e-2
#      iter=0
#      converged=False
#      u=u0*1.
#      p=p0*1.
#      while not converged and iter<iter_max:
#          du=self.solve_f(u,p,tol)
#          u-=du
#          norm_du=util.Lsup(du)
#          norm_u=util.Lsup(u)
#        
#          dp=self.relaxation*self.solve_g(u,tol)
#          p+=dp
#          norm_dp=util.sqrt(self.inner(dp,dp))
#          norm_p=util.sqrt(self.inner(p,p))
#          print iter,"-th step rel. errror u,p= (%s,%s) (%s,%s)(%s,%s)"%(norm_du,norm_dp,norm_du/norm_u,norm_dp/norm_p,norm_u,norm_p)
#          iter+=1
#
#          converged = (norm_du <= tolerance*norm_u) and  (norm_dp <= tolerance*norm_p)
#      if converged:
#          print "convergence after %s steps."%iter
#      else:
#          raise ArithmeticError("no convergence after %s steps."%iter)
#
#      return u,p
          
def MaskFromBoundaryTag(function_space,*tags):
   """
   create a mask on the given function space which one for samples 
   that touch the boundary tagged by tags.

   usage: m=MaskFromBoundaryTag(Solution(domain),"left", "right")

   @param function_space: a given function space 
   @type function_space: L{escript.FunctionSpace}
   @param tags: boundray tags
   @type tags: C{str}
   @return: a mask which marks samples used by C{function_space} that are touching the
            boundary tagged by any of the given tags.
   @rtype: L{escript.Data} of rank 0
   """
   pde=linearPDEs.LinearPDE(function_space.getDomain(),numEquations=1, numSolutions=1)
   d=escript.Scalar(0.,escript.FunctionOnBoundary(function_space.getDomain()))
   for t in tags: d.setTaggedValue(t,1.)
   pde.setValue(y=d)
   out=util.whereNonZero(pde.getRightHandSide())
   if out.getFunctionSpace() == function_space:
      return out
   else:
      return util.whereNonZero(util.interpolate(out,function_space))



