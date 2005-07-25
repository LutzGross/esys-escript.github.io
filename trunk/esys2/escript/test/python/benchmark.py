# $Id$

from os import getenv
from time import time
from escript.escript import *
from finley.finley import Brick
from numarray import ones

cache_size=3 # Mbytes
nterm_max=256
ne=20000
np=int(getenv("OMP_NUM_THREADS"))
print "np =",np

print "nthreads, num. terms, num. comp, volume[Mbytes], cache usage, Gflops/thread"

for var in [True,False]:

     n=int(ne**(1./3.))
     dom=Brick(n,n,n)
     l=(n+1)**3
     sx=dom.getX()[0]+dom.getX()[1]+dom.getX()[2]

     for ncomp in [1, 2, 4, 8, 16, 32, 128, 256]:

        fac=ones([ncomp])
        dataset=[]
        for i in range(nterm_max+1): 
             dataset.append(sx*fac*(i+1))

        for nterm in [2, 4, 8, 16, 32, 64, 128, 256]:

           total_length=(ncomp*float(l)*(nterm+1)*8.)*1.e-6

           if var:

              opcount=(nterm-1)*float(l)*ncomp*1.e-9
              a=sx*0.*fac
              rv=DataVariable(a)
              dv=[]
              for d in dataset[:nterm]:
                dv.append(DataVariable(d))
              t=time()
              for i in range(1,nterm):
                dv[i].sum(dv[i-1])
              r=dv[nterm-1].evaluate()
              t=time()-t

           else:

              opcount=nterm*float(l)*ncomp*1.e-9
              r=sx*0.*fac
              t=time()
              for d in dataset[:nterm]:
                r+=d
              t=time()-t

           err=Lsup(sx*fac*(nterm+1.)*nterm/2.-r)
           if err>1.e-8:
               raise SystemError,"result wrong (error %e)!"%err
           print "%d,%d,%d,%f,%f,%f"%(np,nterm,ncomp,total_length,total_length/(np*cache_size),opcount/t/np)
