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
  * Laplace2Dorder2_?k
  * Laplace3Dorder1_?k
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

from esys.escript import Lsup,whereZero,kronecker
from esys.escript.benchmark import BenchmarkProblem, Options, BenchmarkFilter
import esys.finley 
from esys.escript.linearPDEs import LinearPDE
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

class FinleyOptions(Options):
   """
   finley solver options to be handed over to paso

   """
   def __init__(self,solver_method=None,
                     preconditioner=None,
                     package=None,
                     tolerance=None,
                     verbose=False):
       self.strmap={
                      LinearPDE.DIRECT : "DIRECT",
                      LinearPDE.PCG:  "PCG",
                      LinearPDE.CR:  "CR",
                      LinearPDE.CGS: "CGS",
                      LinearPDE.BICGSTAB: "BICGSTAB",
                      LinearPDE.SSOR: "SSOR",
                      LinearPDE.ILU0: "ILU0",
                      LinearPDE.ILUT: "ILUT",
                      LinearPDE.JACOBI: "JACOBI",
                      LinearPDE.GMRES:  "GMRES",
                      LinearPDE.PRES20:  "PRES20",
                      LinearPDE.LUMPING:  "LUMPIMG",
                      LinearPDE.NO_REORDERING:  "NO_REORDERING",
                      LinearPDE.MINIMUM_FILL_IN:  "MINIMUM_FILL_IN",
                      LinearPDE.NESTED_DISSECTION: "NESTED_DISSECTION",
                      LinearPDE.SCSL:  "SCSL",
                      LinearPDE.MKL:  "MKL",
                      LinearPDE.UMFPACK: "UMFPACK",
                      LinearPDE.PASO:  "PASO"
                  }
       name=""
       if solver_method==None: 
             solver_method==LinearPDE.PRES20
       else:
             name+=self.strmap[solver_method]
       if preconditioner==None: 
             preconditioner==LinearPDE.JACOBI
       else:
             if not name=="": name+="+"
             name+=self.strmap[preconditioner]
       if package==None: 
             package==LinearPDE.PASO
       else:
             if not name=="": name+=" with "
             name+=self.strmap[package]
       if tolerance==None: 
             tolerance=1.e-8
       else:
             if not name=="": name+=", "
             name+="tol = %s"%tolerance
       self.solver_method=solver_method
       self.preconditioner=preconditioner
       self.tolerance=tolerance
       self.package=package
       self.verbose=verbose
       super(FinleyOptions,self).__init__(name=name)



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
       pde.setTolerance(options.tolerance)
       pde.setSolverMethod(options.solver_method,options.preconditioner)
       pde.setSolverPackage(options.package)
       a=os.times()[4]
       uh=pde.getSolution(verbose=options.verbose)
       a=os.times()[4]-a
       if u==None:
          return {FinleyFilter.TIME : a , FinleyFilter.ERROR : None }
       else:
          error=Lsup(u-uh)/Lsup(u)
          return {FinleyFilter.TIME : a , FinleyFilter.ERROR : error }

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
       if self.__dim==2:
          domain=esys.finley.Rectangle(n0=self.__n,n1=self.__n,order=self.__order)
       else:
          domain=esys.finley.Brick(n0=self.__n,n1=self.__n,n2=self.__n,order=self.__order)
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
         pde=LinearPDE(domain)
         pde.setSymmetryOn() 
         pde.setValue(A=kronecker(domain),q=msk,r=u)
         return pde,u

class Laplace2DOrder1_30k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_30k,self).__init__(n=172,order=1,dim=2)
class Laplace2DOrder1_60k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_60k,self).__init__(n=244,order=1,dim=2)
class Laplace2DOrder1_120k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_120k,self).__init__(n=345,order=1,dim=2)
class Laplace2DOrder1_240k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_240k,self).__init__(n=489,order=1,dim=2)
class Laplace2DOrder1_480k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_480k,self).__init__(n=692,order=1,dim=2)
class Laplace2DOrder1_960k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_960k,self).__init__(n=979,order=1,dim=2)
class Laplace2DOrder1_1920k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_1920k,self).__init__(n=1385,order=1,dim=2)
class Laplace2DOrder1_3840k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_3840k,self).__init__(n=1959,order=1,dim=2)
class Laplace2DOrder1_7680k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_7680k,self).__init__(n=2770,order=1,dim=2)
class Laplace2DOrder1_15360k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder1_15360k,self).__init__(n=3918,order=1,dim=2)
class Laplace2DOrder2_30k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_30k,self).__init__(n=86,order=2,dim=2)
class Laplace2DOrder2_60k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_60k,self).__init__(n=122,order=2,dim=2)
class Laplace2DOrder2_120k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_120k,self).__init__(n=173,order=2,dim=2)
class Laplace2DOrder2_240k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_240k,self).__init__(n=244,order=2,dim=2)
class Laplace2DOrder2_480k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_480k,self).__init__(n=346,order=2,dim=2)
class Laplace2DOrder2_960k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_960k,self).__init__(n=489,order=2,dim=2)
class Laplace2DOrder2_1920k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_1920k,self).__init__(n=692,order=2,dim=2)
class Laplace2DOrder2_3840k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_3840k,self).__init__(n=979,order=2,dim=2)
class Laplace2DOrder2_7680k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_7680k,self).__init__(n=1385,order=2,dim=2)
class Laplace2DOrder2_15360k(LaplaceProblem):
   def __init__(self):
      super(Laplace2DOrder2_15360k,self).__init__(n=1959,order=2,dim=2)
class Laplace3DOrder1_30k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_30k,self).__init__(n=30,order=1,dim=3)
class Laplace3DOrder1_60k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_60k,self).__init__(n=38,order=1,dim=3)
class Laplace3DOrder1_120k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_120k,self).__init__(n=48,order=1,dim=3)
class Laplace3DOrder1_240k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_240k,self).__init__(n=61,order=1,dim=3)
class Laplace3DOrder1_480k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_480k,self).__init__(n=77,order=1,dim=3)
class Laplace3DOrder1_960k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_960k,self).__init__(n=98,order=1,dim=3)
class Laplace3DOrder1_1920k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_1920k,self).__init__(n=123,order=1,dim=3)
class Laplace3DOrder1_3840k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_3840k,self).__init__(n=156,order=1,dim=3)
class Laplace3DOrder1_7680k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_7680k,self).__init__(n=196,order=1,dim=3)
class Laplace3DOrder1_15360k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder1_15360k,self).__init__(n=248,order=1,dim=3)
class Laplace3DOrder2_30k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_30k,self).__init__(n=15,order=2,dim=3)
class Laplace3DOrder2_60k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_60k,self).__init__(n=19,order=2,dim=3)
class Laplace3DOrder2_120k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_120k,self).__init__(n=24,order=2,dim=3)
class Laplace3DOrder2_240k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_240k,self).__init__(n=31,order=2,dim=3)
class Laplace3DOrder2_480k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_480k,self).__init__(n=39,order=2,dim=3)
class Laplace3DOrder2_960k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_960k,self).__init__(n=49,order=2,dim=3)
class Laplace3DOrder2_1920k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_1920k,self).__init__(n=62,order=2,dim=3)
class Laplace3DOrder2_3840k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_3840k,self).__init__(n=78,order=2,dim=3)
class Laplace3DOrder2_7680k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_7680k,self).__init__(n=98,order=2,dim=3)
class Laplace3DOrder2_15360k(LaplaceProblem):
   def __init__(self):
      super(Laplace3DOrder2_15360k,self).__init__(n=124,order=2,dim=3)

if __name__=="__main__":
   test=""
   n0=30000
   for d in [2,3]:
     for o in [1,2]:
        for i in range(10):
             dofs=n0*2**i
             n=int((float(dofs)**(1./float(d))-1)/o+0.5)
             name="Laplace%sDOrder%s_%sk"%(d,o,dofs/1000)
             print "class %s(LaplaceProblem):"%name
             print "   def __init__(self):"
             print "      super(%s,self).__init__(n=%s,order=%s,dim=%s)"%(name,n,o,d)
             test+="addProblem(%s())\n"%name
   print test

