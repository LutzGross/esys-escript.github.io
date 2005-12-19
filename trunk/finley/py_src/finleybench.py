# $Id:$

#
#      COPYRIGHT ACcESS 2004 -  All Rights Reserved
#
#   This software is the property of ACcESS.  No part of this code
#   may be copied in any form or by any means without the expressed written
#   consent of ACcESS.  Copying, use or modification of this software
#   by any unauthorised person is illegal unless that
#   person has a software license agreement with ACcESS.
#

"""
some benchmarks for tetsing the finley solver. The idea is to develop a set of standart benchmarks

  * Laplace2Dorder1_?k
  * Laplace3Dorder2_?k

where ? is approximatively the number of unknowns in 1000.

@var __author__: name of author
@var __licence__: licence agreement
var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"

from esys.escript.benchmark import BenchmarkProblem, Options, BenchmarkFilter
from esys.escript import Lsup
import esys.finley 
import os

class FinleyFilter(BenchmarkFilter):
   """
   defines a filter for L{FinleyProblem} characteristics
   """
   TIME="t [sec]"
   ERROR="rel. error"


   def __init__(self,args=None):
      """
      sets up the filter

      @param args: list of value names to be filtered
      @type args: C{list} of L{TIME}, L{ERROR}
      """
      if args==None: args=[FinleyFilter.TIME,FinleyFilter.ERROR]
      super(FinleyFilter,self).__init__()
      self.__args=args

   def getResultNames(self):
       """
       return the names of the results produced when run() is called.
       
       @return: names the list of the names to be used when the results of the run() call are printed
       @rtype: C{list} of C{str}
       """
       return self.__args

   def __call__(self,result):
       """
       filters out the characteristic values
       
       @param result: characteristics rturned by a L{FinleyProblem} run
       @type result: C{dict}
       @return: filtered values
       @rtype: C{list} of C{str}
       """
       out=[]
       for a in self.__args:
         out.append(result[a])
       return out
       


class FinleyProblem(BenchmarkProblem):
   """
   The general benchmark problem for Finley
   """
   def run(self,options):
       """
       creates a domain and a PDE on this domain, solves it (with the given options) and returns the 
       elapsed time and the error.
       
       @param options: paso solver options
       @type options:  L{PasoOptions}
       @return:  elapsed time and the error of calculated solution
       @rtype: pair of C{float}
       """
       domain=self.getDomain()
       pde,u=self.getTestProblem(domain)
       pde.setTolerance(options.getTolerance())
       pde.setSolverMethod(options.getSolverMethod())
       pde.setSolverPackage(options.getSolverPackage())
       a=os.times()[4]
       u_h=pde.getSolution(options.getPasoOptions())
       a=os.times()[4]-a
       if u==None:
          return {FinleyFilter.TIME : a ,FinleyFilter.TIME : None }
       else:
          error=Lsup(u-uh)/Lsup(u)
          return {FinleyFilter.TIME : a ,FinleyFilter.TIME : error }

   def getTestProblem(self,domain):
       """
       returns a PDEto be solved and an exact solution

       @param domain: the PDE domain
       @type domain: L{escript.Domain}
       @return: a linear PDE to be solved an a reference solution
       @rtype: L{LinearPDE},L{escript.Data}
       @remark: must be overwritten by a particular problem
       """
       raise NotImplementedError

   def getDomain(self):
       """
       returns the domain of the problem

       @return: a domain
       @rtype: L{escript.Domain}
       @remark: must be overwritten by a particular problem
       """
       raise NotImplementedError

class RegularFinleyProblem(FinleyProblem):
    """
    base class for finley problem on a rectangular mesh
    """
    def __init__(self,n=1,order=1,dim=2):
       """
       sets up a recangular mesh in finley on a unit cube/square
 
       @param n: number of elements in each spactial direction
       @type n: C{int}
       @param order: element order
       @type order: 1 or 2
       @param dim: spatial dimension
       @type n: 2 or 3
       """
       super(RegularFinleyProblem,self).__init__(name=str((order*n+1)**dim))
       self.__n=n
       self.__order=order
       self.__dim=dim

    def getDomain(self):
       """
       returns the unit square/cube with a rectangular mesh

       @return: a domain
       @rtype: L{escript.Domain}
       """
       if dim==2:
          domain=esys.finley.Rectangle(n0=self.__n,n1=self.__n,order=order)
       else:
          domain=esys.finley.Brick(n0=self.__n,n1=self.__n,n2=self.__n,order=order)
       return domain

class LaplaceProblem(RegularFinleyProblem):
    """
    base class for the Lapalce eqaution on a rectangular mesh
    """
    def getTestProblem(self,domain):
         """
         returns a PDE and a test solution on the given domain
     
         @param doamin: a domain
         @type domain: L{escript.Domain}
         @return: the Laplace equation and a test solution
         @rtype: C{tuple} of C{LinearPDE} and C{escript.Data}
         """
         x=domain.getX()
         msk=whereZero(x[0])+whereZero(x[0]-1.)
         u=x[0]
         for i in range(1,domain.getDim()):
            msk+=whereZero(x[i])+whereZero(x[i]-1.)
            u*=(x[i]-i)
         pde=LinearPDE(mydomain)
         pde.setSymmetryOn() 
         pde.setValue(A=kronnecker,q=msk,r=u)
         return pde,u

class Laplace2DOrder1_30k(LaplaceProblem): 
    def __init__(self):
         super(Laplace2DOrder1_30k,self).__init__(n=176,order=1,dim=2)
class Laplace2DOrder2_30k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder2_30k,self).__init__(n=88,order=2,dim=2)
class Laplace2DOrder1_60k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder1_60k,self).__init__(n=248,order=1,dim=2)
class Laplace2DOrder2_60k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder2_60k,self).__init__(n=124,order=2,dim=2)
class Laplace2DOrder1_120k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder1_120k,self).__init__(n=349,order=1,dim=2)
class Laplace2DOrder2_120k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder2_120k,self).__init__(n=175,order=2,dim=2)
class Laplace2DOrder1_240k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder1_240k,self).__init__(n=492,order=1,dim=2)
class Laplace2DOrder2_240k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder2_240k,self).__init__(n=246,order=2,dim=2)
class Laplace2DOrder1_480k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder1_480k,self).__init__(n=694,order=1,dim=2)
class Laplace2DOrder2_480k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder2_480k,self).__init__(n=347,order=2,dim=2)
class Laplace2DOrder1_960k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder1_960k,self).__init__(n=978,order=1,dim=2)
class Laplace2DOrder2_960k(LaplaceProblem):
    def __init__(self):
         super(Laplace2DOrder2_960k,self).__init__(n=489,order=2,dim=2)
