import sys
import unittest
import os

from esys.escript import *
from esys import bruce

import numarray

"""

Test unary ops on Data objects.

Version $Id$

"""

from numarray import array,Float64,ones,greater

Tag1=10
Tag2=11

tol=1.E-15

#
#  list of arguments: a list item has the form [a0,a1,a2]
#  what a0 is the default value and a1 is used for tag Tag1
#  and a2 for tag2. a0,a1,a2 are converted into numarrays.
#  

arglist = [ \
[ [3,4], [-5,6.], [2,3] ], \
[ [[1,2],[3,4]], [[5,6],[7,8]], [[-5,-6],[7,8]] ], \
[ [[15,8],[12,8]], [[-9,9],[13,8]], [[7,34],[19,7]] ], \
[ [[[15,8],[12,8]],[[-9,9],[13,8]]], [[[3,4],[-9,4]],[[1,-9],[7,4]]], [[[5,2],[6,2]],[[-6,4],[7,5]]] ], \
[ 3.0, 6.0, 3 ] \
]

# these are used to test slicing:
a_r1=[ [1,2,3], [-1,-2,-3], [100,200,300] ]
a_r1_in=[ [1./1,2,3], [-1./1,-1./2,-1./3], [1./100,1./200,1./300] ]
a_r4=[ \
[ [ [[ 1,2,3],[11,12,13]], [[21,22,23],[31,32,33]], [[41,42,43],[51,52,53]]  ], [ [[101,102,103],[111,112,113]], [[121,122,123],[131,132,133]], [[141,142,143],[151,152,153]]  ], [ [[201,202,203],[211,212,213]], [[221,222,223],[231,232,233]], [[241,242,243],[251,252,253]]  ] ], \
[ [ [[ -1,-2,-3],[-11,-12,-13]], [[-21,-22,-23],[-31,-32,-33]], [[-41,-42,-43],[-51,-52,-53]]  ], [ [[-101,-102,-103],[-111,-112,-113]], [[-121,-122,-123],[-131,-132,-133]], [[-141,-142,-143],[-151,-152,-153]]  ], [ [[-201,-202,-203],[-211,-212,-213]], [[-221,-222,-223],[-231,-232,-233]], [[-241,-242,-243],[-251,-252,-253]]  ] ], \
[ [[[ 11,12,13],[111,112,113]], [[121,122,123],[131,132,133]], [[141,142,143],[151,152,153]]  ], [ [[1101,1102,1103],[1111,1112,1113]], [[1121,1122,1123],[1131,1132,1133]], [[1141,1142,1143],[1151,1152,1153]]  ], [ [[1201,1202,1203],[1211,1212,1213]], [[1221,1222,1223],[1231,1232,1233]], [[1241,1242,1243],[1251,1252,1253]]  ] ] ]
a_r4_in=[ \
[ [ [[ 1./1,1./2,1./3],[1./11,1./12,1./13]], [[1./21,1./22,1./23],[1./31,1./32,1./33]], [[1./41,1./42,1./43],[1./51,1./52,1./53]]  ], [ [[1./101,1./102,1./103],[1./111,1./112,1./113]], [[1./121,1./122,1./123],[1./131,1./132,1./133]], [[1./141,1./142,1./143],[1./151,1./152,1./153]]  ], [ [[1./201,1./202,1./203],[1./211,1./212,1./213]], [[1./221,1./222,1./223],[1./231,1./232,1./233]], [[1./241,1./242,1./243],[1./251,1./252,1./253]]  ] ], \
[ [ [[ -1./1,-1./2,-1./3],[-1./11,-1./12,-1./13]], [[-1./21,-1./22,-1./23],[-1./31,-1./32,-1./33]], [[-1./41,-1./42,-1./43],[-1./51,-1./52,-1./53]]  ], [ [[-1./101,-1./102,-1./103],[-1./111,-1./112,-1./113]], [[-1./121,-1./122,-1./123],[-1./131,-1./132,-1./133]], [[-1./141,-1./142,-1./143],[1./-151,-1./152,-1./153]]  ], [ [[-1./201,-1./202,-1./203],[-1./211,-1./212,-1./213]], [[-1./221,-1./222,-1./223],[-1./231,-1./232,-1./233]], [[-1./241,-1./242,-1./243],[-1./251,-1./252,-1./253]]  ] ], \
[ [[[ 1./11,1./12,1./13],[1./111,1./112,1./113]], [[1./121,1./122,1./123],[1./131,1./132,1./133]], [[1./141,1./142,1./143],[1./151,1./152,1./153]]  ], [ [[1./1101,1./1102,1./1103],[1./1111,1./1112,1./1113]], [[1./1121,1./1122,1./1123],[1./1131,1./1132,1./1133]], [[1./1141,1./1142,1./1143],[1./1151,1./1152,1./1153]]  ], [ [[1./1201,1./1202,1./1203],[1./1211,1./1212,1./1213]], [[1./1221,1./1222,1./1223],[1./1231,1./1232,1./1233]], [[1./1241,1./1242,1./1243],[1./1251,1./1252,1./1253]]  ] ] ]

def turnToArray(val,tagged):
     if tagged=="Tagged1":
         out=[array(val[0],Float64),array(val[1],Float64),array(val[0],Float64)]
     elif tagged=="Tagged2":
         out=[array(val[0],Float64),array(val[1],Float64),array(val[2],Float64)]
     else: 
         out=[array(val[0],Float64),array(val[0],Float64),array(val[0],Float64)]
     return out

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

  #for ex1 in ["Constant","Expanded","Tagged1","Tagged2"]:
  for ex1 in ["Constant","Expanded"]:

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
                      +arrays1[0], \
                      +arrays1[1], \
                      +arrays1[2], \
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
                      (arg1-3).wherePositive(), \
                      numarray.greater(arrays1[0],3.), \
                      numarray.greater(arrays1[1],3.), \
                      numarray.greater(arrays1[2],3.), \
                      wh)

      # where negative:
      ref=checkResult("where negative("+ex1+")", \
                      (arg1-3).whereNegative(), \
                      numarray.greater(3.,arrays1[0]), \
                      numarray.greater(3.,arrays1[1]), \
                      numarray.greater(3.,arrays1[2]), \
                      wh)

      # where non-negative:
      ref=checkResult("where nonnegative("+ex1+")", \
                      (arg1-3).whereNonNegative(), \
                      numarray.greater_equal(arrays1[0],3.), \
                      numarray.greater_equal(arrays1[1],3.), \
                      numarray.greater_equal(arrays1[2],3.), \
                      wh)

      # where non-positive:
      ref=checkResult("where nonpositive("+ex1+")", \
                      (arg1-3).whereNonPositive(), \
                      numarray.greater_equal(3.,arrays1[0]), \
                      numarray.greater_equal(3.,arrays1[1]), \
                      numarray.greater_equal(3.,arrays1[2]), \
                      wh)

      # where zero:
      ref=checkResult("where zero("+ex1+")", \
                      (arg1-3).whereZero(), \
                      numarray.less_equal(numarray.abs(arrays1[0]-3.),0.0), \
                      numarray.less_equal(numarray.abs(arrays1[1]-3.),0.0), \
                      numarray.less_equal(numarray.abs(arrays1[2]-3.),0.0), \
                      wh)

      # where non-zero:
      ref=checkResult("where nonzero("+ex1+")", \
                      (arg1-3).whereNonZero(), \
                      numarray.greater(numarray.abs(arrays1[0]-3.),0.0), \
                      numarray.greater(numarray.abs(arrays1[1]-3.),0.0), \
                      numarray.greater(numarray.abs(arrays1[2]-3.),0.0), \
                      wh)

      # exponential function:
      ref=checkResult("exp("+ex1+")", \
                      arg1.exp(), \
                      numarray.exp(arrays1[0]), \
                      numarray.exp(arrays1[1]), \
                      numarray.exp(arrays1[2]), \
                      wh)

      # sqrt
      ref=checkResult("sqrt("+ex1+")", \
                      arg1.abs().sqrt(), \
                      numarray.sqrt(numarray.abs(arrays1[0])), \
                      numarray.sqrt(numarray.abs(arrays1[1])), \
                      numarray.sqrt(numarray.abs(arrays1[2])), \
                      wh)

      # sin:
      ref=checkResult("sin("+ex1+")", \
                      arg1.sin(), \
                      numarray.sin(arrays1[0]), \
                      numarray.sin(arrays1[1]), \
                      numarray.sin(arrays1[2]), \
                      wh)

      # cos:
      ref=checkResult("cos("+ex1+")", \
                      arg1.cos(), \
                      numarray.cos(arrays1[0]), \
                      numarray.cos(arrays1[1]), \
                      numarray.cos(arrays1[2]), \
                      wh)

      # tan:
      ref=checkResult("tan("+ex1+")", \
                      arg1.tan(), \
                      numarray.tan(arrays1[0]), \
                      numarray.tan(arrays1[1]), \
                      numarray.tan(arrays1[2]), \
                      wh)

      # asin:
      #ref=checkResult("asin("+ex1+")", \
      #                arg1.asin(), \
      #                numarray.asin(arrays1[0]), \
      #                numarray.asin(arrays1[1]), \
      #                numarray.asin(arrays1[2]), \
      #                wh)

      # acos:
      #ref=checkResult("acos("+ex1+")", \
      #                arg1.acos(), \
      #                numarray.acos(arrays1[0]), \
      #                numarray.acos(arrays1[1]), \
      #                numarray.acos(arrays1[2]), \
      #                wh)

      # atan:
      #ref=checkResult("atan("+ex1+")", \
      #                arg1.atan(), \
      #                numarray.atan(arrays1[0]), \
      #                numarray.atan(arrays1[1]), \
      #                numarray.atan(arrays1[2]), \
      #                wh)

      # sinh:
      ref=checkResult("sinh("+ex1+")", \
                      arg1.sinh(), \
                      numarray.sinh(arrays1[0]), \
                      numarray.sinh(arrays1[1]), \
                      numarray.sinh(arrays1[2]), \
                      wh)

      # cosh:
      ref=checkResult("cosh("+ex1+")", \
                      arg1.cosh(), \
                      numarray.cosh(arrays1[0]), \
                      numarray.cosh(arrays1[1]), \
                      numarray.cosh(arrays1[2]), \
                      wh)

      # tanh:
      ref=checkResult("tanh("+ex1+")", \
                      arg1.tanh(), \
                      numarray.tanh(arrays1[0]), \
                      numarray.tanh(arrays1[1]), \
                      numarray.tanh(arrays1[2]), \
                      wh)

      # asinh:
      #ref=checkResult("asinh("+ex1+")", \
      #                arg1.asinh(), \
      #                numarray.asinh(arrays1[0]), \
      #                numarray.asinh(arrays1[1]), \
      #                numarray.asinh(arrays1[2]), \
      #                wh)

      # acosh:
      #ref=checkResult("acosh("+ex1+")", \
      #                arg1.acosh(), \
      #                numarray.acosh(arrays1[0]), \
      #                numarray.acosh(arrays1[1]), \
      #                numarray.acosh(arrays1[2]), \
      #                wh)

      # atanh:
      #ref=checkResult("atanh("+ex1+")", \
      #                arg1.atanh(), \
      #                numarray.atanh(arrays1[0]), \
      #                numarray.atanh(arrays1[1]), \
      #                numarray.atanh(arrays1[2]), \
      #                wh)

      # get the maximum value at each data point
      #ref=checkResult("maxval("+ex1+")", \
      #                arg1.maxval(), \
      #                arrays1[0].max(), \
      #                arrays1[1].max(), \
      #                arrays1[2].max(), \
      #                wh)

      # get the minimum value at each data point
      #ref=checkResult("minval("+ex1+")", \
      #                arg1.minval(), \
      #                arrays1[0].min(), \
      #                arrays1[1].min(), \
      #                arrays1[2].min(), \
      #                wh)

      # get the length at each data point = sqrt(sum_{i,j,k,l} A[i,j,k,l]^2)
      #ref=checkResult("length("+ex1+")", \
      #                arg1.length(), \
      #                numarray.sqrt((arrays1[0]**2).sum()), \
      #                numarray.sqrt((arrays1[1]**2).sum()), \
      #                numarray.sqrt((arrays1[2]**2).sum()), \
      #                wh)

      # trace:
      #ref=checkResult("trace("+ex1+")", \
      #                arg1.trace(), \
      #                numarray.trace(arrays1[0]), \
      #                numarray.trace(arrays1[1]), \
      #                numarray.trace(arrays1[2]), \
      #                wh)

      # transpose:
      #axis=arrays1[0]/2
      #ref=checkResult("transpose("+ex1+")", \
      #                arg1.transpose(), \
      #                numarray.transpose(arrays1[0],axis), \
      #                numarray.transpose(arrays1[1],axis), \
      #                numarray.transpose(arrays1[2],axis), \
      #                wh)

      # get the signs of the values:
      ref=checkResult("sign("+ex1+")", \
                      arg1.sign(), \
                      numarray.greater(arrays1[0],numarray.zeros(arrays1[0].shape)) \
                        -numarray.less(arrays1[0],numarray.zeros(arrays1[0].shape)),\
                      numarray.greater(arrays1[1],numarray.zeros(arrays1[1].shape)) \
                        -numarray.less(arrays1[1],numarray.zeros(arrays1[1].shape)),\
                      numarray.greater(arrays1[2],numarray.zeros(arrays1[2].shape)) \
                        -numarray.less(arrays1[2],numarray.zeros(arrays1[2].shape)),\
                      wh)

sys.exit(0)
# end
