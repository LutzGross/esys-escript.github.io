# $Id$

"""
Provides some tools related to PDEs. 

Currently includes:
    - Projector - to project a discontinuous
    - Locator - to trace values in data objects at a certain location
    - TimeIntegrationManager - to handel extraplotion in time
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__licence__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licences/osl-3.0.php"""

import escript
import linearPDEs
import numarray
import util

class TimeIntegrationManager:
  """
  a simple mechanism to manage time dependend values. 

  typical usage is:

  dt=0.1 # time increment
  tm=TimeIntegrationManager(inital_value,p=1)
  while t<1.
      v_guess=tm.extrapolate(dt) # extrapolate to t+dt
      v=...
      tm.checkin(dt,v)
      t+=dt

  @remark: currently only p=1 is supported.
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

  def __del__(self):
    return

  def __call__(self, input_data):
    """
    Projects input_data onto a continuous function

    @param input_data: The input_data to be projected.
    """
    out=escript.Data(0.,input_data.getShape(),self.__pde.getFunctionSpaceForSolution())
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

     @remark: This class is similar to the L{LinearPDE} class with A=B=C=X=0 but has the intention
              that all input parameter are given in L{Solution} or L{ReducedSolution}. The whole
              thing is a bit strange and I blame Robert.Woodcock@csiro.au for this.
     """
     def __init__(self,domain,D=None,Y=None,q=None,r=None):
         """
         initialize the problem

         @param domain: domain of the PDE.
         @type domain: L{Domain}
         @param D: coefficient of the solution. 
         @type D: C{float}, C{int}, L{NumArray}, L{Data}
         @param Y: right hand side
         @type Y: C{float}, C{int}, L{NumArray}, L{Data}
         @param q: location of constraints
         @type q: C{float}, C{int}, L{NumArray}, L{Data}
         @param r: value of solution at locations of constraints
         @type r: C{float}, C{int}, L{NumArray}, L{Data}
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
         @type D: C{float}, C{int}, L{NumArray}, L{Data}
         @param Y: right hand side
         @type Y: C{float}, C{int}, L{NumArray}, L{Data}
         @param q: location of constraints
         @type q: C{float}, C{int}, L{NumArray}, L{Data}
         @param r: value of solution at locations of constraints
         @type r: C{float}, C{int}, L{NumArray}, L{Data}
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
       """
       if isinstance(where,escript.FunctionSpace):
          self.__function_space=where
       else:
          self.__function_space=escript.ContinuousFunction(where)
       self.__id=util.length(self.__function_space.getX()-x[:self.__function_space.getDim()]).mindp()

     def __str__(self):
       """
       Returns the coordinates of the Locator as a string.
       """
       return "<Locator %s>"%str(self.getX())

     def getFunctionSpace(self):
        """
	Returns the function space of the Locator.
	"""
        return self.__function_space

     def getId(self):
        """
	Returns the identifier of the location.
	"""
        return self.__id

     def getX(self):
        """
	Returns the exact coordinates of the Locator.
	"""
        return self(self.getFunctionSpace().getX())

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
             out=data.convertToNumArrayFromDPNo(self.getId()[0],self.getId()[1])
           else:
             out=data.interpolate(self.getFunctionSpace()).convertToNumArrayFromDPNo(self.getId()[0],self.getId()[1])
           if data.getRank()==0:
              return out[0]
           else:
              return out
        else:
           return data

# vim: expandtab shiftwidth=4:
