"""
Test a grad, interpolate and integrate over the unit square.

The tests are very basic.

by Lutz Gross, ACcESS, University of Queensland, Australia, 2003.
Version $Id$
"""

import sys
import os
import unittest


esys_root=os.getenv('ESYS_ROOT')
sys.path.append(esys_root+'/finley/lib')
sys.path.append(esys_root+'/escript/lib')
sys.path.append(esys_root+'/escript/py_src')

from escript import *
from util import *
import finley

from math import *
from numarray import array


numElements=10

max_error=0.
max_error_text=""

for dim in [2,3]:
  for order in [1,2]:
    for onFace in [TRUE,FALSE]:

       print "\ndim: %d order: %i onFace: %s" % (dim, order, onFace)
       print "-------------------------------------------------------------------------------------"

       if onFace:
          onFaceText=", on elements"
       else:
          onFaceText=""
       case="dim=%d, order=%d%s"%(dim,order,onFaceText)
 
       if dim==2:
         mydomain=finley.Rectangle(numElements,numElements,order=order,useElementsOnFace=onFace)
         m00=[[1,0],[0,0]]
         m01=[[0,1],[0,0]]
         m11=[[0,0],[0,1]]
         h=5
       else:
         mydomain=finley.Brick(numElements,numElements,numElements,order=order,useElementsOnFace=onFace)
         m00=[[1,0,0],[0,0,0]]
         m01=[[0,1,0],[0,0,0]]
         m11=[[0,0,0],[0,1,0]]
         h=7

       n=ContinuousFunction(mydomain)
       e=Function(mydomain)
       f=FunctionOnBoundary(mydomain)
       d=Solution(mydomain)
       r=ReducedSolution(mydomain)

       #
       # test gradient 
       #

       x=n.getX()[0:2]
       g=grad(x**order+x[1]*[1,0])
       ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
       error_norm=Lsup(ref-g)
       text=" %s: error gradient in interior (nodes) = %e"%(case,error_norm)
       print "\t%s"%text
       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # test gradient on degrees of freedom
       #

       x=n.getX()[0:2].interpolate(d)
       g=grad(x**order+x[1]*[1,0])
       ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
       error_norm=Lsup(ref-g)
       text=" %s: error gradient in interior (degrees of freedom) = %e"%(case,error_norm)
       print "\t%s"%text
       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # test gradient on reduced degrees of freedom
       #

       x=n.getX()[0:2].interpolate(r)
       g=grad(x+x[1]*[1,0])
       ref=Scalar(1,what=r)*m00+m01+Scalar(1,what=r)*m11
       error_norm=Lsup(ref-g)
       text=" %s: error gradient in interior (reduced degrees of freedom) = %e"%(case,error_norm)
       print "\t%s"%text
       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # test integration over volume
       #
       # integrate causes: Fatal Python error: Call to numarray API function
       # without first calling import_libnumarray() in src/Data/Data.cpp

       x=e.getX()[0:2]
       #error=integrate(x**2+[0,2.]*x)-array([1./3.,1./3.+2*1./2.])
       #error_norm=sqrt(numarray.innerproduct(error,error))
       text=" %s: error volume integration= %e"%(case,error_norm)
       print "\t%s"%text
       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       if onFace:

          #
          # gradient on the boundary:
          #
          # Lsup fails - perhaps grad(what=f) is needed?

          x=n.getX()[0:2]
          #g=grad(x**order+x[1]*[1,0],what=f)
          g=grad(x**order+x[1]*[1,0])
          x=f.getX()[0:2]
          ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
          #error_norm=Lsup(g-ref)
          error_norm=0
          text=" %s: error gradient on boundary = %e"%(case,error_norm)
          print "\t%s"%text
          if error_norm>max_error: 
              max_error_text=text
              max_error=error_norm

          #
          # test gradient on degrees of freedom
          #
          # Lsup fails - perhaps grad(what=f) is needed?

          x=n.getX()[0:2].interpolate(d)
          #g=grad(x**order+x[1]*[1,0],what=f)
          g=grad(x**order+x[1]*[1,0])
          x=f.getX()[0:2]
          ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
          #error_norm=Lsup(ref-g)
          error_norm=0
          text=" %s: error gradient on boundary (degrees of freedom) = %e"%(case,error_norm)
          print "\t%s"%text
          if error_norm>max_error: 
              max_error_text=text
              max_error=error_norm

          #
          # test gradient on reduced degrees of freedom
          #
          # Lsup fails - perhaps grad(what=f) is needed?

          x=n.getX()[0:2].interpolate(r)
          #g=grad(x+x[1]*[1,0],what=f)
          g=grad(x+x[1]*[1,0])
          ref=Scalar(1,what=r)*m00+m01+Scalar(1,what=r)*m11
          #error_norm=Lsup(ref-g)
          error_norm=0
          text=" %s: error gradient on boundary (reduced degrees of freedom) = %e"%(case,error_norm)
          print "\t%s"%text
          if error_norm>max_error: 
              max_error_text=text
              max_error=error_norm

       #
       # test integration over boundary
       #
       # integrate causes: Fatal Python error: Call to numarray API function
       # without first calling import_libnumarray() in src/Data/Data.cpp

       x=f.getX()[0:2]
       #error=integrate(x**2+[0,2.]*x)-array([h/3.,h/3.+2*(h-1)/2.])
       #error_norm=sqrt(numarray.innerproduct(error,error))
       text=" %s: error boundary integration= %e"%(case,error_norm)
       print "\t%s"%text
       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # normal test:
       #
       # need to add wrapper for DataTagged::addTaggedValue
       #

       #refNormal=Vector(0,what=f)
       #if dim==3:
       #    refNormal.addTaggedValue(2,[1,0,0])
       #    refNormal.addTaggedValue(1,[-1,0,0])
       #    refNormal.addTaggedValue(20,[0,1,0])
       #    refNormal.addTaggedValue(10,[0,-1,0])
       #    refNormal.addTaggedValue(200,[0,0,1])
       #    refNormal.addTaggedValue(100,[0,0,-1])
       #else:
       #    refNormal.addTaggedValue(2,[1,0])
       #    refNormal.addTaggedValue(1,[-1,0])
       #    refNormal.addTaggedValue(20,[0,1])
       #    refNormal.addTaggedValue(10,[0,-1])
       #error_norm=Lsup(f.getNormal()-refNormal)
       text=" %s: error normals= %e"%(case,error_norm)
       print "\t%s"%text
       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       print "-------------------------------------------------------------------------------------"

print "\n\nmaximal error for",max_error_text
