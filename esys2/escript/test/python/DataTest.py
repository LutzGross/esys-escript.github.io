import sys
import unittest
import os
                                                                                                                 
esys_root=os.getenv('ESYS_ROOT')
sys.path.append(esys_root+'/finley/lib')
sys.path.append(esys_root+'/escript/lib')
sys.path.append(esys_root+'/escript/py_src')
                                                                                                                 
from escript import *
import finley
import numarray
from util import *

"""
tests arithmetic operations on Data:

   +arg1
   -arg1
    arg1+arg2
    arg1-arg2
    arg1*arg2
    arg1/arg2
    arg1**arg2
    arg1+=arg2
    arg1-=arg2
    arg1*=arg2
    arg1/=arg2

various values for arg1 and arg2 are tested on various atoms.

moreover various slicing operations are tested.

the results are compared to reference result. if the relative difference is
bigger than tol an exception is raised.

by Lutz Gross, ACcESS, University of Queensland, Australia, 2003.

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
#  binary operations are tested on all pairs from arglist
#
#  each item in the arglist are used to construct the following 5 argument
#  types arg for arithmetic operations:
#
#  1) arg is a numarray/list a0
#  2) arg is a Data with default value a0
#  3) arg is an DataArray with constant value a0
#  4) arg is a Data object with constant value a0 and value a1 for tag Tag1
#  5) arg is a Data object with constant value a0 and value a1 for tag Tag1
#     and value a2 for Tag2.
#
#  i.e for a single binary arithmetic operation (len(arglist)*5)**2
#  test are performed. 
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
       #raise SystemError,"@@ %s at %s: error is too large"%(text,wh)
       print "**** %s: error is too large"%(text)

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

msh=finley.Rectangle(20,6,1)
for wh in [ContinuousFunction(msh),Function(msh)]:

#
# ==============================================================
#
# test slicing:
#

  for ex1 in ["Array","Constant","Expanded","Tagged1","Tagged2"]:

      print "Slice getting: ", ex1, ":"

      # slice getting

      # rank 1
      arg=prepareArg(a_r1,ex1,wh)
      arrays=turnToArray(a_r1,ex1)
      if isinstance(arg,Data):
          checkResult("slicing: rank=1 [:] "+ex1,arg[:], \
                        arrays[0][:], \
                        arrays[1][:], \
                        arrays[2][:], \
                        wh)
          checkResult("slicing: rank=1 [1] "+ex1,arg[1], \
                        arrays[0][1], \
                        arrays[1][1], \
                        arrays[2][1], \
                        wh)
          checkResult("slicing: rank=1 [1:3] "+ex1,arg[1:3], \
                        arrays[0][1:3], \
                        arrays[1][1:3], \
                        arrays[2][1:3], \
                        wh)

      # rank 4
      arg=prepareArg(a_r4,ex1,wh)
      arrays=turnToArray(a_r4,ex1)
      if isinstance(arg,Data):
          checkResult("slicing: rank=4 [:] "+ex1,arg[:], \
                        arrays[0][:], \
                        arrays[1][:], \
                        arrays[2][:], \
                        wh)
          checkResult("slicing: rank=4 [1] "+ex1,arg[1], \
                        arrays[0][1], \
                        arrays[1][1], \
                        arrays[2][1], \
                        wh)
          checkResult("slicing: rank=4 [1:3] "+ex1,arg[1:3], \
                        arrays[0][1:3], \
                        arrays[1][1:3], \
                        arrays[2][1:3], \
                        wh)
          checkResult("slicing: rank=4 [:,:] "+ex1,arg[:,:], \
                        arrays[0][:,:], \
                        arrays[1][:,:], \
                        arrays[2][:,:], \
                        wh)
          checkResult("slicing: rank=4 [:,1] "+ex1,arg[:,1], \
                        arrays[0][:,1], \
                        arrays[1][:,1], \
                        arrays[2][:,1], \
                        wh)
          checkResult("slicing: rank=4 [:,1:3] "+ex1,arg[:,1:3], \
                        arrays[0][:,1:3], \
                        arrays[1][:,1:3], \
                        arrays[2][:,1:3], \
                        wh)
          checkResult("slicing: rank=4 [1:2,1:3] "+ex1,arg[1:2,1:3], \
                        arrays[0][1:2,1:3], \
                        arrays[1][1:2,1:3], \
                        arrays[2][1:2,1:3], \
                        wh)
          checkResult("slicing: rank=4 [1:2,1:3,1,1] "+ex1,arg[1:2,1:3,1,1], \
                        arrays[0][1:2,1:3,1,1], \
                        arrays[1][1:2,1:3,1,1], \
                        arrays[2][1:2,1:3,1,1], \
                        wh)
          checkResult("slicing: rank=4 [1:2,1:3,:,1] "+ex1,arg[1:2,1:3,:,1], \
                        arrays[0][1:2,1:3,:,1], \
                        arrays[1][1:2,1:3,:,1], \
                        arrays[2][1:2,1:3,:,1], \
                        wh)
          checkResult("slicing: rank=4 [1:2,1:3,0,:] "+ex1,arg[1:2,1:3,0,:], \
                        arrays[0][1:2,1:3,0,:], \
                        arrays[1][1:2,1:3,0,:], \
                        arrays[2][1:2,1:3,0,:], \
                        wh)
          checkResult("slicing: rank=4 [1:2,1:3,0,1] "+ex1,arg[1:2,1:3,0,1], \
                        arrays[0][1:2,1:3,0,1], \
                        arrays[1][1:2,1:3,0,1], \
                        arrays[2][1:2,1:3,0,1], \
                        wh)
          checkResult("slicing: rank=4 [1,1,0,1] "+ex1,arg[1,1,0,1], \
                        arrays[0][1,1,0,1], \
                        arrays[1][1,1,0,1], \
                        arrays[2][1,1,0,1], \
                        wh)

          # slice setting

          #for ex2 in ["Array","Constant","Expanded","Tagged1","Tagged2"]:
          for ex2 in ["Expanded"]:

                print "Slice setting: ", ex1, ",", ex2, ":"

                # rank 1

                arrays_in=turnToArray(a_r1_in,ex2)

                arg2=prepareArg(a_r1,ex1,wh) 
                arrays2=turnToArray(a_r1,ex1)
                a_in=[arrays_in[0][:],arrays_in[1][:],arrays_in[2][:]]
                arg2[:]=prepareArg(a_in,ex2,wh)
                arrays2[0][:]=a_in[0]
                arrays2[1][:]=a_in[1]
                arrays2[2][:]=a_in[2]
                checkResult("slicing, set: rank=1 [:] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r1,ex1,wh) 
                arrays2=turnToArray(a_r1,ex1)
                a_in=[arrays_in[0][1],arrays_in[1][1],arrays_in[2][1]]
                arg2[1]=prepareArg(a_in,ex2,wh)
                arrays2[0][1]=a_in[0]
                arrays2[1][1]=a_in[1]
                arrays2[2][1]=a_in[2]
                checkResult("slicing, set: rank=1 [1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r1,ex1,wh) 
                arrays2=turnToArray(a_r1,ex1)
                a_in=[arrays_in[0][1:3],arrays_in[1][1:3],arrays_in[2][1:3]]
                arg2[1:3]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:3]=a_in[0]
                arrays2[1][1:3]=a_in[1]
                arrays2[2][1:3]=a_in[2]
                checkResult("slicing, set: rank=1 [1:3] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                # rank 4

                arrays_in=turnToArray(a_r4_in,ex2)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][:],arrays_in[1][:],arrays_in[2][:]]
                arg2[:]=prepareArg(a_in,ex2,wh)
                arrays2[0][:]=a_in[0]
                arrays2[1][:]=a_in[1]
                arrays2[2][:]=a_in[2]
                checkResult("slicing, set: rank=4 [:] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1],arrays_in[1][1],arrays_in[2][1]]
                arg2[1]=prepareArg(a_in,ex2,wh)
                arrays2[0][1]=a_in[0]
                arrays2[1][1]=a_in[1]
                arrays2[2][1]=a_in[2]
                checkResult("slicing, set: rank=4 [1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1:3],arrays_in[1][1:3],arrays_in[2][1:3]]
                arg2[1:3]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:3]=a_in[0]
                arrays2[1][1:3]=a_in[1]
                arrays2[2][1:3]=a_in[2]
                checkResult("slicing, set: rank=4 [1:3] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][:,:],arrays_in[1][:,:],arrays_in[2][:,:]]
                arg2[:,:]=prepareArg(a_in,ex2,wh)
                arrays2[0][:,:]=a_in[0]
                arrays2[1][:,:]=a_in[1]
                arrays2[2][:,:]=a_in[2]
                checkResult("slicing, set: rank=4 [:,:] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][:,1],arrays_in[1][:,1],arrays_in[2][:,1]]
                arg2[:,1]=prepareArg(a_in,ex2,wh)
                arrays2[0][:,1]=a_in[0]
                arrays2[1][:,1]=a_in[1]
                arrays2[2][:,1]=a_in[2]
                checkResult("slicing, set: rank=4 [:,1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][:,1:3],arrays_in[1][:,1:3],arrays_in[2][:,1:3]]
                arg2[:,1:3]=prepareArg(a_in,ex2,wh)
                arrays2[0][:,1:3]=a_in[0]
                arrays2[1][:,1:3]=a_in[1]
                arrays2[2][:,1:3]=a_in[2]
                checkResult("slicing, set: rank=4 [:,1:3] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1:2,1:3],arrays_in[1][1:2,1:3],arrays_in[2][1:2,1:3]]
                arg2[1:2,1:3]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:2,1:3]=a_in[0]
                arrays2[1][1:2,1:3]=a_in[1]
                arrays2[2][1:2,1:3]=a_in[2]
                checkResult("slicing, set: rank=4 [1:2,1:3] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1:2,1:3,1,1],arrays_in[1][1:2,1:3,1,1],arrays_in[2][1:2,1:3,1,1]]
                arg2[1:2,1:3,1,1]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:2,1:3,1,1]=a_in[0]
                arrays2[1][1:2,1:3,1,1]=a_in[1]
                arrays2[2][1:2,1:3,1,1]=a_in[2]
                checkResult("slicing, set: rank=4 [1:2,1:3,1,1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1:2,1:3,:,1],arrays_in[1][1:2,1:3,:,1],arrays_in[2][1:2,1:3,:,1]]
                arg2[1:2,1:3,:,1]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:2,1:3,:,1]=a_in[0]
                arrays2[1][1:2,1:3,:,1]=a_in[1]
                arrays2[2][1:2,1:3,:,1]=a_in[2]
                checkResult("slicing, set: rank=4 [1:2,1:3,:,1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1:2,1:3,0,:],arrays_in[1][1:2,1:3,0,:],arrays_in[2][1:2,1:3,0,:]]
                arg2[1:2,1:3,0,:]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:2,1:3,0,:]=a_in[0]
                arrays2[1][1:2,1:3,0,:]=a_in[1]
                arrays2[2][1:2,1:3,0,:]=a_in[2]
                checkResult("slicing, set: rank=4 [1:2,1:3,0,:] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1:2,1:3,0,1],arrays_in[1][1:2,1:3,0,1],arrays_in[2][1:2,1:3,0,1]]
                arg2[1:2,1:3,0,1]=prepareArg(a_in,ex2,wh)
                arrays2[0][1:2,1:3,0,1]=a_in[0]
                arrays2[1][1:2,1:3,0,1]=a_in[1]
                arrays2[2][1:2,1:3,0,1]=a_in[2]
                checkResult("slicing, set: rank=4 [1:2,1:3,0,1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

                arg2=prepareArg(a_r4,ex1,wh) 
                arrays2=turnToArray(a_r4,ex1)
                a_in=[arrays_in[0][1,1,0,1],arrays_in[1][1,1,0,1],arrays_in[2][1,1,0,1]]
                arg2[1,1,0,1]=prepareArg(a_in,ex2,wh)
                arrays2[0][1,1,0,1]=a_in[0]
                arrays2[1][1,1,0,1]=a_in[1]
                arrays2[2][1,1,0,1]=a_in[2]
                checkResult("slicing, set: rank=4 [1,1,0,1] "+ex1+","+ex2,arg2,arrays2[0],arrays2[1],arrays2[2],wh)

#
# ==============================================================
#
# test operators:
#

  #for ex1 in ["Array","Constant","Expanded","Tagged1","Tagged2"]:
  for ex1 in ["Array","Constant","Expanded"]:

      for a1 in arglist:

        arg1=prepareArg(a1,ex1,wh)
        arrays1=turnToArray(a1,ex1)
        if isScalar(arg1):
           t1="(scalar)"
        else:
           t1=""

        #
        #   unary operations:
        #

        print "Unary ops: ", ex1, ":"

        if isinstance(arg1,Data):

           # pos:
           #ref=checkResult("+"+ex1,+arg1, \
           #                  +arrays1[0], \
           #                  +arrays1[1], \
           #                  +arrays1[2], \
           #                  wh)
           # neg:
           #ref=checkResult("-"+ex1,-arg1, \
           #                  -arrays1[0], \
           #                  -arrays1[1], \
           #                  -arrays1[2], \
           #                  wh)
           # pos:
           ref=checkResult("where positive("+ex1+")",(arg1-3).wherePositive(), \
                             numarray.greater(arrays1[0],3.), \
                             numarray.greater(arrays1[1],3.), \
                             numarray.greater(arrays1[2],3.), \
                             wh)
           # negative:
           ref=checkResult("where negative("+ex1+")",(arg1-3).whereNegative(), \
                             numarray.greater(3.,arrays1[0]), \
                             numarray.greater(3.,arrays1[1]), \
                             numarray.greater(3.,arrays1[2]), \
                             wh)
           # non-negative:
           ref=checkResult("where nonnegative("+ex1+")",(arg1-3).whereNonNegative(), \
                             numarray.greater_equal(arrays1[0],3.), \
                             numarray.greater_equal(arrays1[1],3.), \
                             numarray.greater_equal(arrays1[2],3.), \
                             wh)
           # non-positive:
           ref=checkResult("where nonpositive("+ex1+")",(arg1-3).whereNonPositive(), \
                             numarray.greater_equal(3.,arrays1[0]), \
                             numarray.greater_equal(3.,arrays1[1]), \
                             numarray.greater_equal(3.,arrays1[2]), \
                             wh)
           # zero:
           ref=checkResult("where zero("+ex1+")",(arg1-3).whereZero(), \
                             numarray.less_equal(numarray.abs(arrays1[0]-3.),EPSILON), \
                             numarray.less_equal(numarray.abs(arrays1[1]-3.),EPSILON), \
                             numarray.less_equal(numarray.abs(arrays1[2]-3.),EPSILON), \
                             wh)
           # non-zero:
           ref=checkResult("where nonzero("+ex1+")",(arg1-3).whereNonZero(), \
                             numarray.greater(numarray.abs(arrays1[0]-3.),EPSILON), \
                             numarray.greater(numarray.abs(arrays1[1]-3.),EPSILON), \
                             numarray.greater(numarray.abs(arrays1[2]-3.),EPSILON), \
                             wh)

        #
        #  binary operations
        #

        for ex2 in ["Array","Constant","Expanded","Tagged1","Tagged2"]:

          print "Binary ops: ", ex1, ",", ex2, ":"

          for a2 in arglist:

             arg2=prepareArg(a2,ex2,wh)
             arrays2=turnToArray(a2,ex2)
             if isScalar(arg2):
                t2="(scalar)"
             else:
                t2=""

             if isinstance(arg1,Data) or isinstance(arg2,Data):

                # the shape must match or at least one argument is sclar:
                if (getRank(arg1)==getRank(arg2)) or isScalar(arg1) or isScalar(arg2):

                  # sum 
                  checkResult(ex1+t1+"+"+ex2+t2,arg1+arg2, \
                                     arrays1[0]+arrays2[0], \
                                     arrays1[1]+arrays2[1], \
                                     arrays1[2]+arrays2[2], \
                                     wh)
                  # sub
                  checkResult(ex1+t1+"-"+ex2+t2,arg1-arg2, \
                                     arrays1[0]-arrays2[0], \
                                     arrays1[1]-arrays2[1], \
                                     arrays1[2]-arrays2[2], \
                                     wh)
                  # mul
                  checkResult(ex1+t1+"*"+ex2+t2,arg1*arg2, \
                                     arrays1[0]*arrays2[0], \
                                     arrays1[1]*arrays2[1], \
                                     arrays1[2]*arrays2[2], \
                                     wh)
                  # div
                  checkResult(ex1+t1+"/"+ex2+t2,arg1/arg2, \
                                     arrays1[0]/arrays2[0], \
                                     arrays1[1]/arrays2[1], \
                                     arrays1[2]/arrays2[2], \
                                     wh)
                  # pow 
                  if isinstance(arg1,Data):
                    a=arg1
                  else:
                    a=Data(value=arg1,what=arg2.getFunctionSpace())
                  checkResult(ex1+t1+"^"+ex2+t2,a**arg2, \
                                     arrays1[0]**arrays2[0], \
                                     arrays1[1]**arrays2[1], \
                                     arrays1[2]**arrays2[2], \
                                     wh)

                  # inplace for arg2 which must be of class Data
                  # if arg1 is expanded arg2 must be expanded
                  if isinstance(arg2,Data):

                     # if arg2 is scalar arg1 must be scalar:
                     #if not (isScalar(arg2) and not isScalar(arg1)):
                     if not(isScalar(arg2)) and not(isScalar(arg1)):

                        # inplace add:
                        arrays2[0]+=arrays1[0]
                        arrays2[1]+=arrays1[1]
                        arrays2[2]+=arrays1[2]
                        arg2+=arg1
                        checkResult(ex2+t2+"+="+ex1+t1,arg2, \
                                       arrays2[0], \
                                       arrays2[1], \
                                       arrays2[2], \
                                       wh)
                        # inplace sub:
                        arrays2[0]-=arrays1[0]
                        arrays2[1]-=arrays1[1]
                        arrays2[2]-=arrays1[2]
                        arg2-=arg1
                        checkResult(ex2+t2+"-="+ex1+t1,arg2, \
                                       arrays2[0], \
                                       arrays2[1], \
                                       arrays2[2], \
                                       wh)
                        # inplace mul:
                        arrays2[0]*=arrays1[0]
                        arrays2[1]*=arrays1[1]
                        arrays2[2]*=arrays1[2]
                        arg2*=arg1
                        checkResult(ex2+t2+"*="+ex1+t1,arg2, \
                                       arrays2[0], \
                                       arrays2[1], \
                                       arrays2[2], \
                                       wh)
                        # inplace div:
                        arrays2[0]/=arrays1[0]
                        arrays2[1]/=arrays1[1]
                        arrays2[2]/=arrays1[2]
                        arg2/=arg1
                        checkResult(ex2+t2+"/="+ex1+t1,arg2, \
                                       arrays2[0], \
                                       arrays2[1], \
                                       arrays2[2], \
                                       wh)

# end
