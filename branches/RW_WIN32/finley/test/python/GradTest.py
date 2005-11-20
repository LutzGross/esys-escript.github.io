# $Id$
"""
Test a grad, interpolate and integrate over the unit square.

The tests are very basic.

by Lutz Gross, ACcESS, University of Queensland, Australia, 2003.
"""

import sys
import os
import unittest
import time

from esys.escript import *
from esys.escript.linearPDEs import *
from esys import finley

from math import *
from numarray import array

starttime = time.clock()

numElements=10
max_error=0.
max_error_text=""
error_tol=pow(10,-10)

for dim in [2,3]:
  for order in [1,2]:
    for onFace in [True,False]:

       print "\n"
       print "-----------------------------------------------------------------------------------------------------"
       print "dim: %d, order: %i, onFace: %s." % (dim, order, onFace)
       print "-----------------------------------------------------------------------------------------------------"

       if onFace:
          onFaceText=", on elements"
       else:
          onFaceText=""

       case="dim=%d, order=%d%s" % (dim,order,onFaceText)
 
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

       test="error gradient in interior (nodes)"

       x=n.getX()[0:2]
       g=grad(x**order+x[1]*[1,0])
       ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
       error_norm=Lsup(ref-g)

       text="%s: %55s = %e" % (case, test, error_norm)
       print "%s" % text

       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # test gradient on degrees of freedom
       #

       test="error gradient in interior (degrees of freedom)"

       x=n.getX()[0:2].interpolate(d)
       g=grad(x**order+x[1]*[1,0])
       ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
       error_norm=Lsup(ref-g)/Lsup(ref)

       text="%s: %55s = %e" % (case, test, error_norm)
       print "%s" % text

       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # test gradient on reduced degrees of freedom
       #

       test="error gradient in interior (reduced degrees of freedom)"

       x=n.getX()[0:2].interpolate(r)
       g=grad(x+x[1]*[1,0])
       ref=Scalar(1,what=r)*m00+m01+Scalar(1,what=r)*m11
       error_norm=Lsup(ref-g)/Lsup(ref)

       text="%s: %55s = %e" % (case, test, error_norm)
       print "%s" % text

       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # test integration over volume
       #

       test="error volume integration"

       x=e.getX()[0:2]
       error=integrate(x**2+[0,2.]*x)-array([1./3.,1./3.+2*1./2.])
       error_norm=sqrt(numarray.innerproduct(error,error))

       text="%s: %55s = %e" % (case, test, error_norm)
       print "%s" % text

       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       if onFace:

          #
          # gradient on the boundary:
          #

          test="error gradient on boundary"

          x=n.getX()[0:2]
          g=grad(x**order+x[1]*[1,0],where=f)
          x=f.getX()[0:2]
          ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
          error_norm=Lsup(g-ref)

          text="%s: %55s = %e" % (case, test, error_norm)
          print "%s" % text

          if error_norm>max_error: 
              max_error_text=text
              max_error=error_norm

          #
          # test gradient on degrees of freedom
          #

          test="error gradient on boundary (degrees of freedom)"

          x=n.getX()[0:2].interpolate(d)
          g=grad(x**order+x[1]*[1,0],where=f)
          x=f.getX()[0:2]
          ref=order*x[0]**(order-1)*m00+m01+order*x[1]**(order-1)*m11
          error_norm=Lsup(ref-g)/Lsup(ref)

          text="%s: %55s = %e" % (case, test, error_norm)
          print "%s" % text

          if error_norm>max_error: 
              max_error_text=text
              max_error=error_norm

          #
          # test gradient on reduced degrees of freedom
          #

          test="error gradient on boundary (reduced degrees of freedom)"

          x=n.getX()[0:2].interpolate(r)
          g=grad(x+x[1]*[1,0],where=f)
          ref=Scalar(1,what=r)*m00+m01+Scalar(1,what=r)*m11
          error_norm=Lsup(ref-g)/Lsup(ref)

          text="%s: %55s = %e" % (case, test, error_norm)
          print "%s" % text

          if error_norm>max_error: 
              max_error_text=text
              max_error=error_norm

       #
       # test integration over boundary
       #

       test="error boundary integration"

       x=f.getX()[0:2]
       error=integrate(x**2+[0,2.]*x)-array([h/3.,h/3.+2*(h-1)/2.])
       error_norm=sqrt(numarray.innerproduct(error,error))

       text="%s: %55s = %e" % (case, test, error_norm)
       print "%s" % text

       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       #
       # normal test:
       #

       test="error normals"

       refNormal=Vector(0,what=f)
       if dim==3:
           refNormal.setTaggedValue(2,[1,0,0])
           refNormal.setTaggedValue(1,[-1,0,0])
           refNormal.setTaggedValue(20,[0,1,0])
           refNormal.setTaggedValue(10,[0,-1,0])
           refNormal.setTaggedValue(200,[0,0,1])
           refNormal.setTaggedValue(100,[0,0,-1])
       else:
           refNormal.setTaggedValue(2,[1,0])
           refNormal.setTaggedValue(1,[-1,0])
           refNormal.setTaggedValue(20,[0,1])
           refNormal.setTaggedValue(10,[0,-1])
       error_norm=Lsup(f.getNormal()-refNormal)

       text="%s: %55s = %e" % (case, test, error_norm)
       print "%s" % text

       if error_norm>max_error: 
           max_error_text=text
           max_error=error_norm

       print "-----------------------------------------------------------------------------------------------------"

print "\n\n"
print "******************************************************************************************************************"
print "maximal error:", max_error_text
print "******************************************************************************************************************"

stoptime = time.clock()
elapsed = stoptime - starttime
print "\nElapsed time: ", elapsed, "\n"

if max_error > error_tol:
  print "max error exceeds tolerance"
  sys.exit(1)

sys.exit(0)
