
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file in the repository root for contributors and development history
# https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS
#
##############################################################################


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
some benchmarks for tetsing the finley solver. The idea is to develop a set of standart benchmarks.

:var __author__: name of author
:var __licence__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from esys.escript import *
from esys.escript.benchmark import BenchmarkProblem, Options, BenchmarkFilter
import esys.finley 
from esys.escript.linearPDEs import LinearPDE, SolverOptions
import os
import math
import numpy

class FinleyFilter(BenchmarkFilter):
   """
   defines a filter for `FinleyProblem` characteristics
   """
   TIME="t [sec]"
   ERROR="rel. error"


   def __init__(self,args=None):
      """
      sets up the filter

      :param args: list of value names to be filtered
      :type args: ``list`` of `TIME`, `ERROR`
      """
      if args==None: args=[FinleyFilter.TIME,FinleyFilter.ERROR]
      super(FinleyFilter,self).__init__()
      self.__args=args

   def getResultNames(self):
       """
       return the names of the results produced when run() is called.
       
       :return: names the list of the names to be used when the results of the run() call are printed
       :rtype: ``list`` of ``str``
       """
       return self.__args

   def __call__(self,result):
       """
       filters out the characteristic values
       
       :param result: characteristics rturned by a `FinleyProblem` run
       :type result: ``dict``
       :return: filtered values
       :rtype: ``list`` of ``str``
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
                      SolverOptions.DIRECT : "DIRECT",
                      SolverOptions.PCG:  "PCG",
                      SolverOptions.CR:  "CR",
                      SolverOptions.CGS: "CGS",
                      SolverOptions.BICGSTAB: "BICGSTAB",
                      SolverOptions.SSOR: "SSOR",
                      SolverOptions.ILU0: "ILU0",
                      SolverOptions.ILUT: "ILUT",
                      SolverOptions.JACOBI: "JACOBI",
                      SolverOptions.GMRES:  "GMRES",
                      SolverOptions.PRES20:  "PRES20",
                      SolverOptions.LUMPING:  "LUMPINremG",
                      SolverOptions.NO_REORDERING:  "NO_REORDERING",
                      SolverOptions.MINIMUM_FILL_IN:  "MINIMUM_FILL_IN",
                      SolverOptions.NESTED_DISSECTION: "NESTED_DISSECTION",
                      SolverOptions.MKL:  "MKL",
                      SolverOptions.UMFPACK: "UMFPACK",
                      SolverOptions.MUMPS: "MUMPS",
                      SolverOptions.TRILINOS: "TRILINOS",
                      SolverOptions.PASO:  "PASO",
                      SolverOptions.RILU: "RILU",
                      SolverOptions.AMG:  "AMG"
                  }
       name=""
       if solver_method==None: 
             solver_method==SolverOptions.PRES20
       else:
             name+=self.strmap[solver_method]
       if preconditioner==None: 
             preconditioner==SolverOptions.JACOBI
       else:
             if not name=="": name+="+"
             name+=self.strmap[preconditioner]
       if package==None: 
             package==SolverOptions.PASO
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
       
       :param options: solver options
       :type options:  `FinleyOptions`
       :return:  elapsed time and the error of calculated solution
       :rtype: pair of ``float``
       """
       domain=self.getDomain()
       pde,u=self.getTestProblem(domain)
       pde.getSolverOptions().setTolerance(options.tolerance)
       pde.getSolverOptions().setPreconditioner(options.preconditioner)
       pde.getSolverOptions().setSolverMethod(options.solver_method)
       pde.getSolverOptions().setPackage(options.package)
       pde.getSolverOptions().setVerbosity(options.verbose)
       pde.getSolverOptions().setIterMax(6000)
       a=os.times()[4]
       uh=pde.getSolution()
       a=os.times()[4]-a
       if u==None:
          return {FinleyFilter.TIME : a , FinleyFilter.ERROR : None }
       else:
          error=Lsup(u-uh)/Lsup(u)
          return {FinleyFilter.TIME : a , FinleyFilter.ERROR : error }

   def getTestProblem(self,domain):
       """
       returns a PDEto be solved and an exact solution

       :param domain: the PDE domain
       :type domain: `escript.Domain`
       :return: a linear PDE to be solved an a reference solution
       :rtype: `LinearPDE`,`escript.Data`
       :note: must be overwritten by a particular problem
       """
       raise NotImplementedError

   def getDomain(self):
       """
       returns the domain of the problem

       :return: a domain
       :rtype: `escript.Domain`
       :note: must be overwritten by a particular problem
       """
       raise NotImplementedError

class RegularFinleyProblem(FinleyProblem):
    """
    base class for finley problem on a rectangular mesh
    """
    def __init__(self,n=1,order=1,dim=2,num_equations=1):
       """
       sets up a recangular mesh in finley on a unit cube/square
 
       :param n: number of elements in each spactial direction
       :type n: ``int``
       :param order: element order
       :type order: 1 or 2
       :param dim: spatial dimension
       :type dim: 2 or 3
       :param num_equations: number of equations
       :type num_equations: ``int``
       """
       super(RegularFinleyProblem,self).__init__(name=str(num_equations*(order*n+1)**dim))
       self.__n=n
       self.__order=order
       self.__dim=dim
       self.__num_equations=num_equations

    def getDomain(self):
       """
       returns the unit square/cube with a rectangular mesh

       :return: a domain
       :rtype: `escript.Domain`
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
     
         :param domain: a domain
         :type domain: `escript.Domain`
         :return: the Laplace equation and a test solution
         :rtype: ``tuple`` of ``LinearPDE`` and ``escript.Data``
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

class AnisotropicProblem(RegularFinleyProblem):
    """
    base class for the Anisotropic scalar problem on a rectangular mesh
    """
    def __init__(self,n,order,dim,gamma,c):
        self.c=c
        self.gamma=gamma
        super(AnisotropicProblem,self).__init__(n,order,dim)


    def getTestProblem(self,domain):
         """
         returns a PDE and a test solution on the given domain
     
         :param domain: a domain
         :type domain: `escript.Domain`
         :return: the Laplace equation and a test solution
         :rtype: ``tuple`` of ``LinearPDE`` and ``escript.Data``
         """
         x=domain.getX()
         msk=whereZero(x[0])+whereZero(x[0]-1.)
         u=x[0]
         for i in range(1,domain.getDim()):
            msk+=whereZero(x[i])+whereZero(x[i]-1.)
            u*=(x[i]-i)

         gamma_rad=self.gamma/360.*8*math.atan(1.)
         cg=math.cos(gamma_rad)
         sg=math.sin(gamma_rad)
         C=kronecker(domain)
         C[0,0]=cg**2+self.c*sg**2
         C[1,0]=(self.c-1.)*cg*sg
         C[0,1]=C[1,0]
         C[1,1]=sg**2+self.c*cg**2
         F=2*(1.-self.c)*cg*sg
         if domain.getDim()==3: F*=x[2]-2.
         pde=LinearPDE(domain)
         pde.setSymmetryOn() 
         pde.setValue(A=C,Y=F,q=msk,r=u)
         return pde,u

class AnisotropicSystem(RegularFinleyProblem):
    """
    base class for the Anisotropic system problem on a rectangular mesh
    with an anisotropic 

    *- (mu*(u_{i,j}+u_{j,i}))_j+lam*u_{k,k})_j=X_{ij,j}*
       
    where 
           - *u_i=x_i+1/d*\prod{i!=j} x_j*
           - *mu(x) = 1* for inner(x,normal)<inner(1.,normal)/2.
             and *mu(x) =  mu0* for inner(x,normal)>inner(1.,normal)/2.
           - lam(x)=max(1,mu0)*alpha (constant)
    plus constraints on the boundary
    """
    def __init__(self,n,order,dim,mu0,normal,alpha):
        self.mu0=mu0
        self.normal=numpy.array(normal)
        self.alpha=alpha
        super(AnisotropicSystem,self).__init__(n,order,dim,dim)


    def getTestProblem(self,domain):
         """
         returns a PDE and a test solution on the given domain
     
         :param domain: a domain
         :type domain: `escript.Domain`
         :return: the Laplace equation and a test solution
         :rtype: ``tuple`` of ``LinearPDE`` and ``escript.Data``
         """
         x=domain.getX()
         d=domain.getDim()

         msk=whereZero(x[0])+whereZero(x[0]-1.)
         for i in range(1,d):
            msk+=whereZero(x[i])+whereZero(x[i]-1.)
         msk=msk*numpy.ones((d,),numpy.float64)

         u=x[:]
         for i in range(d):
            s=1.
            for k in range(d):
               if not i==k: s=s*x[k]
            u[i]+=1./d*s

         s=whereNegative(inner(x-numpy.ones((d,),numpy.float64)/2,self.normal))
         mu=s+self.mu0*(1.-s)
         lam=max(1.,self.mu0)*self.alpha

         F=Tensor(0.,Function(domain))
         for i in range(d):
            for j in range(d):
                 if i==j:
                    F[i,i]+=mu*d*lam
                 else:
                    s=1.
                    for k in range(d): 
                       if not i==k and not j==k:
                          s*=x[k]
                    F[i,j]+=mu/d*s
         C=Tensor4(0.,Function(domain))
         for i in range(domain.getDim()):
           for j in range(domain.getDim()):
                C[i,i,j,j]+=lam
                C[j,i,j,i]+=mu
                C[j,i,i,j]+=mu
         pde=LinearPDE(domain)
         pde.setSymmetryOn() 
         pde.setValue(A=C,X=F,q=msk,r=u)
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

class Anisotropic2DOrder1Gamma30_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_30k,self).__init__(n=172,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_60k,self).__init__(n=244,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_120k,self).__init__(n=345,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_240k,self).__init__(n=489,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_480k,self).__init__(n=692,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_960k,self).__init__(n=979,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_1920k,self).__init__(n=1385,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_3840k,self).__init__(n=1959,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_7680k,self).__init__(n=2770,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma30_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma30_15360k,self).__init__(n=3918,order=1,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder1Gamma45_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_30k,self).__init__(n=172,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_60k,self).__init__(n=244,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_120k,self).__init__(n=345,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_240k,self).__init__(n=489,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_480k,self).__init__(n=692,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_960k,self).__init__(n=979,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_1920k,self).__init__(n=1385,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_3840k,self).__init__(n=1959,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_7680k,self).__init__(n=2770,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder1Gamma45_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder1Gamma45_15360k,self).__init__(n=3918,order=1,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma30_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_30k,self).__init__(n=86,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_60k,self).__init__(n=122,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_120k,self).__init__(n=173,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_240k,self).__init__(n=244,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_480k,self).__init__(n=346,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_960k,self).__init__(n=489,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_1920k,self).__init__(n=692,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_3840k,self).__init__(n=979,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_7680k,self).__init__(n=1385,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma30_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma30_15360k,self).__init__(n=1959,order=2,dim=2,gamma=30,c=0.001)
class Anisotropic2DOrder2Gamma45_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_30k,self).__init__(n=86,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_60k,self).__init__(n=122,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_120k,self).__init__(n=173,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_240k,self).__init__(n=244,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_480k,self).__init__(n=346,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_960k,self).__init__(n=489,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_1920k,self).__init__(n=692,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_3840k,self).__init__(n=979,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_7680k,self).__init__(n=1385,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic2DOrder2Gamma45_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic2DOrder2Gamma45_15360k,self).__init__(n=1959,order=2,dim=2,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma30_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_30k,self).__init__(n=30,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_60k,self).__init__(n=38,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_120k,self).__init__(n=48,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_240k,self).__init__(n=61,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_480k,self).__init__(n=77,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_960k,self).__init__(n=98,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_1920k,self).__init__(n=123,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_3840k,self).__init__(n=156,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_7680k,self).__init__(n=196,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma30_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma30_15360k,self).__init__(n=248,order=1,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder1Gamma45_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_30k,self).__init__(n=30,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_60k,self).__init__(n=38,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_120k,self).__init__(n=48,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_240k,self).__init__(n=61,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_480k,self).__init__(n=77,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_960k,self).__init__(n=98,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_1920k,self).__init__(n=123,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_3840k,self).__init__(n=156,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_7680k,self).__init__(n=196,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder1Gamma45_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder1Gamma45_15360k,self).__init__(n=248,order=1,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma30_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_30k,self).__init__(n=15,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_60k,self).__init__(n=19,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_120k,self).__init__(n=24,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_240k,self).__init__(n=31,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_480k,self).__init__(n=39,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_960k,self).__init__(n=49,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_1920k,self).__init__(n=62,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_3840k,self).__init__(n=78,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_7680k,self).__init__(n=98,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma30_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma30_15360k,self).__init__(n=124,order=2,dim=3,gamma=30,c=0.001)
class Anisotropic3DOrder2Gamma45_30k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_30k,self).__init__(n=15,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_60k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_60k,self).__init__(n=19,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_120k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_120k,self).__init__(n=24,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_240k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_240k,self).__init__(n=31,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_480k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_480k,self).__init__(n=39,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_960k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_960k,self).__init__(n=49,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_1920k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_1920k,self).__init__(n=62,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_3840k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_3840k,self).__init__(n=78,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_7680k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_7680k,self).__init__(n=98,order=2,dim=3,gamma=45,c=0.001)
class Anisotropic3DOrder2Gamma45_15360k(AnisotropicProblem):
   def __init__(self):
      super(Anisotropic3DOrder2Gamma45_15360k,self).__init__(n=124,order=2,dim=3,gamma=45,c=0.001)
class Lame2DOrder1_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder1Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder1Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder1Alpha100_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE2Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE2Alpha100_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder1JumpE6Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder1JumpE6Alpha100_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE2Normal45Alpha100_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder1JumpE6Normal45Alpha100_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class Lame2DOrder2_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=1.0)
class Lame2DOrder2Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class Lame2DOrder2Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame2DOrder2Alpha100_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE2Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE2Alpha100_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=1.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomgeneousLame2DOrder2JumpE6Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame2DOrder2JumpE6Alpha100_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE2Normal45Alpha100_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=1.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame2DOrder2JumpE6Normal45Alpha100_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=100.0)
class Lame3DOrder1_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder1Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder1Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder1Alpha100_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE2Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE2Alpha100_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder1JumpE6Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder1JumpE6Alpha100_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE2Normal45Alpha100_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder1JumpE6Normal45Alpha100_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class Lame3DOrder2_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=1.0)
class Lame3DOrder2Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class Lame3DOrder2Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(Lame3DOrder2Alpha100_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE2Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE2Alpha100_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=1.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomgeneousLame3DOrder2JumpE6Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomgeneousLame3DOrder2JumpE6Alpha100_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE2Normal45Alpha100_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=1.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_30k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_60k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_120k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_240k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_480k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_960k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_1920k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_3840k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_7680k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)
class InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_15360k(AnisotropicSystem):
   def __init__(self):
      super(InhomogeneousLame3DOrder2JumpE6Normal45Alpha100_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=100.0)

class CompressibleLame2DOrder1_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder1_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder1_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE2_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder1JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder1JumpE6_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE2Normal45_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_30k,self).__init__(n=121,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_60k,self).__init__(n=172,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_120k,self).__init__(n=244,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_240k,self).__init__(n=345,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_480k,self).__init__(n=489,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_960k,self).__init__(n=692,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_1920k,self).__init__(n=979,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_3840k,self).__init__(n=1385,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_7680k,self).__init__(n=1959,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder1JumpE6Normal45_15360k,self).__init__(n=2770,order=1,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleLame2DOrder2_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleLame2DOrder2_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame2DOrder2_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+00,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE2_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+02,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomgeneousLame2DOrder2JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame2DOrder2JumpE6_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+06,normal=[1.,0.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE2Normal45_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+02,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_30k,self).__init__(n=61,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_60k,self).__init__(n=86,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_120k,self).__init__(n=122,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_240k,self).__init__(n=173,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_480k,self).__init__(n=244,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_960k,self).__init__(n=346,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_1920k,self).__init__(n=489,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_3840k,self).__init__(n=692,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_7680k,self).__init__(n=979,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame2DOrder2JumpE6Normal45_15360k,self).__init__(n=1385,order=2,dim=2,mu0=1.000000e+06,normal=[1.,1.],alpha=0.0)
class CompressibleLame3DOrder1_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder1_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder1_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE2_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder1JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder1JumpE6_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE2Normal45_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_30k,self).__init__(n=21,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_60k,self).__init__(n=26,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_120k,self).__init__(n=33,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_240k,self).__init__(n=42,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_480k,self).__init__(n=53,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_960k,self).__init__(n=67,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_1920k,self).__init__(n=85,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_3840k,self).__init__(n=108,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_7680k,self).__init__(n=136,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder1JumpE6Normal45_15360k,self).__init__(n=171,order=1,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleLame3DOrder2_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleLame3DOrder2_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleLame3DOrder2_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+00,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE2_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE2_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+02,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomgeneousLame3DOrder2JumpE6_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomgeneousLame3DOrder2JumpE6_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+06,normal=[1.,0.,0.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE2Normal45_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+02,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_30k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_30k,self).__init__(n=10,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_60k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_60k,self).__init__(n=13,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_120k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_120k,self).__init__(n=17,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_240k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_240k,self).__init__(n=21,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_480k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_480k,self).__init__(n=27,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_960k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_960k,self).__init__(n=34,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_1920k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_1920k,self).__init__(n=43,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_3840k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_3840k,self).__init__(n=54,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_7680k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_7680k,self).__init__(n=68,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)
class CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_15360k(AnisotropicSystem):
   def __init__(self):
      super(CompressibleInhomogeneousLame3DOrder2JumpE6Normal45_15360k,self).__init__(n=86,order=2,dim=3,mu0=1.000000e+06,normal=[1.,1.,1.],alpha=0.0)

if __name__=="__main__":
   test=""
   n0=30000
   for d in [2,3]:
    for o in [1,2]:
      for g in [0,45]:
       for jump in [0,2,6]:
         for alpha in [0]:   # [1,100]:
           for i in range(10):
            if not jump==0 or g==0:
             dofs=(n0*2**i)/float(d)
             n=int((float(dofs)**(1./float(d))-1)/o+0.5)
             if jump==0:
                if g==0:
                   if alpha==1:
                      name="Lame%sDOrder%s_%sk"%(d,o,int(d*dofs/1000))
                   elif alpha==0:
                      name="CompressibleLame%sDOrder%s_%sk"%(d,o,int(d*dofs/1000))
                   else:
                      name="Lame%sDOrder%sAlpha%s_%sk"%(d,o,alpha,int(d*dofs/1000))
                else:
                   if alpha==1:
                      name="InhomogeneousLame%sDOrder%sNormal%s_%sk"%(d,o,g,int(d*dofs/1000))
                   elif alpha==0:
                      name="ComressibleInhomogeneousLame%sDOrder%sNormal%s_%sk"%(d,o,g,int(d*dofs/1000))
                   else:
                      name="InhomogeneousLame%sDOrder%sNormal%sAlpha%s_%sk"%(d,o,g,alpha,int(d*dofs/1000))
             else:
                if g==0:
                   if alpha==1:
                      name="InhomgeneousLame%sDOrder%sJumpE%s_%sk"%(d,o,jump,int(d*dofs/1000))
                   elif alpha==0:
                      name="CompressibleInhomgeneousLame%sDOrder%sJumpE%s_%sk"%(d,o,jump,int(d*dofs/1000))
                   else:
                      name="InhomgeneousLame%sDOrder%sJumpE%sAlpha%s_%sk"%(d,o,jump,alpha,int(d*dofs/1000))
                else:
                   if alpha==1:
                     name="InhomogeneousLame%sDOrder%sJumpE%sNormal%s_%sk"%(d,o,jump,g,int(d*dofs/1000))
                   elif alpha==0:
                     name="CompressibleInhomogeneousLame%sDOrder%sJumpE%sNormal%s_%sk"%(d,o,jump,g,int(d*dofs/1000))
                   else:
                     name="InhomogeneousLame%sDOrder%sJumpE%sNormal%sAlpha%s_%sk"%(d,o,jump,g,alpha,int(d*dofs/1000))

             if g==45:
                if d==2:
                    normal="[1.,1.]"
                else:
                    normal="[1.,1.,1.]"
             else:
                if d==2:
                    normal="[1.,0.]"
                else:
                    normal="[1.,0.,0.]"

             print("class %s(AnisotropicSystem):"%name)
             print("   def __init__(self):")
             print("      super(%s,self).__init__(n=%s,order=%s,dim=%s,mu0=%e,normal=%s,alpha=%s)"%(name,n,o,d,10.**jump,normal,float(alpha)))
             test+="addProblem(%s())\n"%name
   print(test)

