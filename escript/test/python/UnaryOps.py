"""

Test unary ops on Data objects.

Version $Id$

"""

import sys
import unittest
import os

from esys.escript import *
from esys import bruce

import numarray

from numarray import array,Float64,ones,greater

Tag1=10
Tag2=15

tol=1.E-15

#
#  list of arguments: a list item has the form [a0,a1,a2]
#  what a0 is the default value and a1 is used for tag Tag1
#  and a2 for Tag2.
#  

arglist = [ \
[ [3,4], [-5,6.], [2,3] ], \
[ [[1,2],[3,4]], [[5,6],[7,8]], [[-5,-6],[7,8]] ], \
[ [[15,8],[12,8]], [[-9,9],[13,8]], [[7,34],[19,7]] ], \
[ [[[15,8],[12,8]],[[-9,9],[13,8]]], [[[3,4],[-9,4]],[[1,-9],[7,4]]], [[[5,2],[6,2]],[[-6,4],[7,5]]] ], \
[ [3.0], [6.0], [3] ] \
]

def prepareArg(val,ex,wh):
    if ex=="Array":
        out=val[0]
    else:
        if ex=="Expanded":
            exx=True
        else:
            exx=False
        out=Data(val[0],what=wh,expand=exx)
        if ex=="Tagged1":
            out.setTaggedValue(Tag1,val[1])
        elif ex=="Tagged2":
            out.setTaggedValue(Tag1,val[1])
            out.setTaggedValue(Tag2,val[2])
    return out

def turnToArray(val,tagged):
    if tagged=="Tagged1":
        out=[array(val[0],Float64),array(val[1],Float64),array(val[0],Float64)]
    elif tagged=="Tagged2":
        out=[array(val[0],Float64),array(val[1],Float64),array(val[2],Float64)]
    else: 
        out=[array(val[0],Float64),array(val[0],Float64),array(val[0],Float64)]
    return out

def checkResult(text,res,val0,val1,val2,wh):
    ref=Data(val0,what=wh,expand=False)
    ref.setTaggedValue(Tag1,val1)
    ref.setTaggedValue(Tag2,val2)
    norm=Lsup(ref)+tol
    error=Lsup(ref-res)/norm
    print "@@ %s, shape %s: error = %e"%(text,ref.getShape(),error)
    if error>tol:
        print "**** %s: error is too large"%(text)
        raise SystemError,"@@ %s at %s: error is too large"%(text,wh)
        sys.exit(1)

def getRank(arg):
    if isinstance(arg,Data):
        return arg.getRank()
    else:
        g=array(arg)
        if g.rank==0:
            return 1
        else:
            return g.rank

def isScalar(arg):
    if isinstance(arg,Data):
        if arg.getRank()==1 and arg.getShape()[0]==1:
            return not None
        else:
            return None
    else:
        g=array(arg)
        if g.rank==0:
            return not None
        else:
            if g.rank==1 and g.shape[0]==1:
                return not None
            else:
                return None

#
# ==============================================================
#
# test unary operators:
#

msh=bruce.Rectangle(20,6)
for wh in [ContinuousFunction(msh),Function(msh)]:

  print wh

  for ex1 in ["Constant","Expanded","Tagged1","Tagged2"]:

    print "Unary Ops:", ex1

    for a1 in arglist:

      arg1=prepareArg(a1,ex1,wh)
      arrays1=turnToArray(a1,ex1)
      if isScalar(arg1):
         t1="(scalar)"
      else:
         t1=""

      # + identity:
      ref=checkResult("+"+ex1, \
                      +arg1, \
                      arrays1[0], \
                      arrays1[1], \
                      arrays1[2], \
                      wh)

      # - negation:
      ref=checkResult("-"+ex1, \
                      -arg1, \
                      -arrays1[0], \
                      -arrays1[1], \
                      -arrays1[2], \
                      wh)

      # where positive:
      ref=checkResult("where positive("+ex1+")", \
                      (arg1-3)._wherePositive(), \
                      numarray.greater(arrays1[0],3.), \
                      numarray.greater(arrays1[1],3.), \
                      numarray.greater(arrays1[2],3.), \
                      wh)

      # where negative:
      ref=checkResult("where negative("+ex1+")", \
                      (arg1-3)._whereNegative(), \
                      numarray.greater(3.,arrays1[0]), \
                      numarray.greater(3.,arrays1[1]), \
                      numarray.greater(3.,arrays1[2]), \
                      wh)

      # where non-negative:
      ref=checkResult("where nonnegative("+ex1+")", \
                      (arg1-3)._whereNonNegative(), \
                      numarray.greater_equal(arrays1[0],3.), \
                      numarray.greater_equal(arrays1[1],3.), \
                      numarray.greater_equal(arrays1[2],3.), \
                      wh)

      # where non-positive:
      ref=checkResult("where nonpositive("+ex1+")", \
                      (arg1-3)._whereNonPositive(), \
                      numarray.greater_equal(3.,arrays1[0]), \
                      numarray.greater_equal(3.,arrays1[1]), \
                      numarray.greater_equal(3.,arrays1[2]), \
                      wh)

      # where zero:
      ref=checkResult("where zero("+ex1+")", \
                      (arg1-3)._whereZero(), \
                      numarray.less_equal(numarray.abs(arrays1[0]-3.),0.0), \
                      numarray.less_equal(numarray.abs(arrays1[1]-3.),0.0), \
                      numarray.less_equal(numarray.abs(arrays1[2]-3.),0.0), \
                      wh)

      # where non-zero:
      ref=checkResult("where nonzero("+ex1+")", \
                      (arg1-3)._whereNonZero(), \
                      numarray.greater(numarray.abs(arrays1[0]-3.),0.0), \
                      numarray.greater(numarray.abs(arrays1[1]-3.),0.0), \
                      numarray.greater(numarray.abs(arrays1[2]-3.),0.0), \
                      wh)

      # exponential function:
      ref=checkResult("exp("+ex1+")", \
                      arg1._exp(), \
                      numarray.exp(arrays1[0]), \
                      numarray.exp(arrays1[1]), \
                      numarray.exp(arrays1[2]), \
                      wh)

      # sqrt
      ref=checkResult("sqrt("+ex1+")", \
                      arg1._sqrt(), \
                      numarray.sqrt(numarray.abs(arrays1[0])), \
                      numarray.sqrt(numarray.abs(arrays1[1])), \
                      numarray.sqrt(numarray.abs(arrays1[2])), \
                      wh)

      # sin:
      ref=checkResult("sin("+ex1+")", \
                      arg1._sin(), \
                      numarray.sin(arrays1[0]), \
                      numarray.sin(arrays1[1]), \
                      numarray.sin(arrays1[2]), \
                      wh)

      # cos:
      ref=checkResult("cos("+ex1+")", \
                      arg1._cos(), \
                      numarray.cos(arrays1[0]), \
                      numarray.cos(arrays1[1]), \
                      numarray.cos(arrays1[2]), \
                      wh)

      # tan:
      ref=checkResult("tan("+ex1+")", \
                      arg1._tan(), \
                      numarray.tan(arrays1[0]), \
                      numarray.tan(arrays1[1]), \
                      numarray.tan(arrays1[2]), \
                      wh)

      # numarray has no asin/acos/atan funcs

      # asin:
      #ref=checkResult("asin("+ex1+")", \
      #                arg1._asin(), \
      #                numarray.asin(arrays1[0]), \
      #                numarray.asin(arrays1[1]), \
      #                numarray.asin(arrays1[2]), \
      #                wh)

      # acos:
      #ref=checkResult("acos("+ex1+")", \
      #                arg1._acos(), \
      #                numarray.acos(arrays1[0]), \
      #                numarray.acos(arrays1[1]), \
      #                numarray.acos(arrays1[2]), \
      #                wh)

      # atan:
      #ref=checkResult("atan("+ex1+")", \
      #                arg1._atan(), \
      #                numarray.atan(arrays1[0]), \
      #                numarray.atan(arrays1[1]), \
      #                numarray.atan(arrays1[2]), \
      #                wh)

      # sinh:
      ref=checkResult("sinh("+ex1+")", \
                      arg1._sinh(), \
                      numarray.sinh(arrays1[0]), \
                      numarray.sinh(arrays1[1]), \
                      numarray.sinh(arrays1[2]), \
                      wh)

      # cosh:
      ref=checkResult("cosh("+ex1+")", \
                      arg1._cosh(), \
                      numarray.cosh(arrays1[0]), \
                      numarray.cosh(arrays1[1]), \
                      numarray.cosh(arrays1[2]), \
                      wh)

      # tanh:
      ref=checkResult("tanh("+ex1+")", \
                      arg1._tanh(), \
                      numarray.tanh(arrays1[0]), \
                      numarray.tanh(arrays1[1]), \
                      numarray.tanh(arrays1[2]), \
                      wh)

      # numarray has no asinh/acosh/atanh funcs

      # asinh:
      #ref=checkResult("asinh("+ex1+")", \
      #                arg1._asinh(), \
      #                numarray.asinh(arrays1[0]), \
      #                numarray.asinh(arrays1[1]), \
      #                numarray.asinh(arrays1[2]), \
      #                wh)

      # acosh:
      #ref=checkResult("acosh("+ex1+")", \
      #                arg1._acosh(), \
      #                numarray.acosh(arrays1[0]), \
      #                numarray.acosh(arrays1[1]), \
      #                numarray.acosh(arrays1[2]), \
      #                wh)

      # atanh:
      #ref=checkResult("atanh("+ex1+")", \
      #                arg1._atanh(), \
      #                numarray.atanh(arrays1[0]), \
      #                numarray.atanh(arrays1[1]), \
      #                numarray.atanh(arrays1[2]), \
      #                wh)

      # get the maximum value at each data point
      ref=checkResult("maxval("+ex1+")", \
                      arg1._maxval(), \
                      arrays1[0].max(), \
                      arrays1[1].max(), \
                      arrays1[2].max(), \
                      wh)

      # get the minimum value at each data point
      ref=checkResult("minval("+ex1+")", \
                      arg1._minval(), \
                      arrays1[0].min(), \
                      arrays1[1].min(), \
                      arrays1[2].min(), \
                      wh)

      # length/trace/transpose not yet implemented for Data

      # get the length at each data point = sqrt(sum_{i,j,k,l} A[i,j,k,l]^2)
      #ref=checkResult("length("+ex1+")", \
      #                arg1._length(), \
      #                numarray.sqrt((arrays1[0]**2).sum()), \
      #                numarray.sqrt((arrays1[1]**2).sum()), \
      #                numarray.sqrt((arrays1[2]**2).sum()), \
      #                wh)

      # trace:
      #ref=checkResult("trace("+ex1+")", \
      #                arg1._trace(), \
      #                numarray.trace(arrays1[0]), \
      #                numarray.trace(arrays1[1]), \
      #                numarray.trace(arrays1[2]), \
      #                wh)

      # transpose:
      #axis=arrays1[0]/2
      #ref=checkResult("transpose("+ex1+")", \
      #                arg1._transpose(), \
      #                numarray.transpose(arrays1[0],axis), \
      #                numarray.transpose(arrays1[1],axis), \
      #                numarray.transpose(arrays1[2],axis), \
      #                wh)

      # get the signs of the values:
      ref=checkResult("sign("+ex1+")", \
                      arg1._sign(), \
                      numarray.greater(arrays1[0],numarray.zeros(arrays1[0].shape)) \
                        -numarray.less(arrays1[0],numarray.zeros(arrays1[0].shape)),\
                      numarray.greater(arrays1[1],numarray.zeros(arrays1[1].shape)) \
                        -numarray.less(arrays1[1],numarray.zeros(arrays1[1].shape)),\
                      numarray.greater(arrays1[2],numarray.zeros(arrays1[2].shape)) \
                        -numarray.less(arrays1[2],numarray.zeros(arrays1[2].shape)),\
                      wh)

sys.exit(0)
# end
