# $Id$

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
              self.__id.append(util.length(self.__function_space.getX()-p[:self.__function_space.getDim()]).mindp())
       else:
           self.__id=util.length(self.__function_space.getX()-x[:self.__function_space.getDim()]).mindp()

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
                  o=data.convertToNumArrayFromDPNo(*i)
                  if data.getRank()==0:
                     out.append(o[0])
                  else:
                     out.append(o)
               return out
           else:
             out=data.convertToNumArrayFromDPNo(*id)
             if data.getRank()==0:
                return out[0]
             else:
                return out
        else:
           return data

class SaddlePointProblem(object):
   """
   This implements a solver for a saddlepoint problem

   M{f(u,p)=0}
   M{g(u)=0}

   for u and p. The problem is solved with an inexact Uszawa scheme for p:

   M{Q_f (u^{k+1}-u^{k}) = - f(u^{k},p^{k})
   M{Q_g (p^{k+1}-p^{k}) =   g(u^{k+1})}

   where Q_f is an approximation of the Jacobiean A_f of f with respect to u  and Q_f is an approximation of
   A_g A_f^{-1} A_g with A_g is the jacobiean of g with respect to p. As a the construction of a 'proper'
   Q_g can be difficult, non-linear conjugate gradient method is applied to solve for p, so Q_g plays
   in fact the role of a preconditioner.
   """
   def __init__(self,verbose=False,*args):
       """
       initializes the problem

       @parm verbose: switches on the printing out some information
       @type verbose: C{bool}
       @note: this method may be overwritten by a particular saddle point problem
       """
       self.__verbose=verbose
       self.relaxation=1.

   def trace(self,text):
       """
       prints text if verbose has been set

       @parm text: a text message
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
            self.trace("%s th step: f/u = %s, g/p = %s, relaxation = %s."%(self.iter,norm_f_new/norm_u_new, norm_g_new/norm_p_new, self.relaxation))

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
                   break
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
          
# vim: expandtab shiftwidth=4:
