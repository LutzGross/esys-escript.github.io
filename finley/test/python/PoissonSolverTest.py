
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript import *
from esys.escript.linearPDEs import Poisson
from esys import finley

ne_list=[10,15,22,33,50,75]
height_list=[0.25,0.5,1.]


def getDomain(dim,ne,height):

    if dim==2:
     ne1=int(ne*height+0.5)
     mydomain=finley.Rectangle(n0=ne,n1=ne1,l1=height,order=1)
     totne=ne1*ne
    else:
     ne2=int(ne*height+0.5)
     mydomain=finley.Brick(n0=ne,n1=ne,n2=ne2,l2=height,order=2)
     totne=ne2*ne*ne
    print("%d -dimensional domain generated."%dim)
    print("height of the domain is ",height)
    print("total number of elements is ",totne)
    return mydomain


def Solve1(mydomain,height):
    print("Fully constraint solution")
    l=[1.,1.,1.]
    l[mydomain.getDim()-1]=height
    cf=ContinuousFunction(mydomain)
    x=cf.getX()
    #construct exact solution:
    u_ex=Scalar(1.,cf)
    for i in range(mydomain.getDim()):
      u_ex*=x[i]*(x[i]-l[i])
    #construct mask:
    msk=Scalar(0.,cf)
    for i in range(mydomain.getDim()):
      msk+=whereZero(x[i])+whereZero(x[i]-l[i])
    #construct right hand side 
    f=Scalar(0,cf)
    for i in range(mydomain.getDim()):
       f_p=Scalar(1,cf)
       for j in range(mydomain.getDim()):
          if i==j:
             f_p*=-2.
          else:
             f_p*=x[j]*(x[j]-l[j])
       f+=f_p

    mypde=Poisson(mydomain)
    mypde.setTolerance(1.e-10)
    mypde.setValue(f=f,q=msk)
    u=mypde.getSolution()
    error=Lsup(u-u_ex)/Lsup(u_ex)
    print("error = ",error)
    return error

def Solve2(mydomain,height):
    print("Partially constraint solution")
    l=[1.,1.,1.]
    l[mydomain.getDim()-1]=height
    print(l)
    cf=ContinuousFunction(mydomain)
    x=cf.getX()
    #construct exact solution:
    u_ex=Scalar(1.,cf)
    for i in range(mydomain.getDim()):
      u_ex*=x[i]*(2*l[i]-x[i])
    #construct mask:
    msk=Scalar(0.,cf)
    for i in range(mydomain.getDim()):
      msk+=whereZero(x[i])
    #construct right hand side 
    f=Scalar(0,cf)
    for i in range(mydomain.getDim()):
       f_p=Scalar(1,cf)
       for j in range(mydomain.getDim()):
          if i==j:
             f_p*=2.
          else:
             f_p*=x[j]*(2*l[j]-x[j])
       f+=f_p
    mypde=Poisson(mydomain)
    mypde.setTolerance(1.e-10)
    mypde.setValue(f=f,q=msk)
    u=mypde.getSolution()
    error=Lsup(u-u_ex)/Lsup(u_ex)
    print("error = ",error)
    return error


def main() :
    error=0
    for ne in ne_list:
       for dim in [2,3]:
       # for dim in [2]:
          for height in height_list:
             print("***************************************************************")
             mydomain= getDomain(dim,ne,height)
             print("---------------------------------------------------------------")
             error=max(error,Solve1(mydomain,height))
             print("---------------------------------------------------------------")
             error=max(error,Solve2(mydomain,height))
             print("***************************************************************")

    print("***************************************************************")
    print("maximum error: ",error)
    print("***************************************************************")



import profile as Pr, pstats as Ps


if __name__ == "__main__":
    pr = Pr.Profile()
    pr.calibrate(10000)
    Pr.run('main()','eos_stats')
    stats = Ps.Stats('eos_stats')
    stats.strip_dirs()
    stats.sort_stats('time')
    stats.print_stats()
