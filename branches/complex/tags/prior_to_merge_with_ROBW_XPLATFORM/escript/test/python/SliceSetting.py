"""

Test slice setting for Data objects.

Version $Id$

"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
import sys
import unittest
import os
                                                                                                                 
from esys.escript import *
from esys import bruce

import numarray

from numarray import array,Float64,ones,greater

Tag1=1
Tag2=2

tol=1.E-15

#  rank 1 test arguments
a_r1=[ [1,2,3], [-1,-2,-3], [100,200,300] ]
a_r1_in=[ [1./1,2,3], [-1./1,-1./2,-1./3], [1./100,1./200,1./300] ]

#  rank 4 test arguments
a_r4=[ \
[ [ [[ 1,2,3],[11,12,13]], [[21,22,23],[31,32,33]], [[41,42,43],[51,52,53]]  ], [ [[101,102,103],[111,112,113]], [[121,122,123],[131,132,133]], [[141,142,143],[151,152,153]]  ], [ [[201,202,203],[211,212,213]], [[221,222,223],[231,232,233]], [[241,242,243],[251,252,253]]  ] ], \
[ [ [[ -1,-2,-3],[-11,-12,-13]], [[-21,-22,-23],[-31,-32,-33]], [[-41,-42,-43],[-51,-52,-53]]  ], [ [[-101,-102,-103],[-111,-112,-113]], [[-121,-122,-123],[-131,-132,-133]], [[-141,-142,-143],[-151,-152,-153]]  ], [ [[-201,-202,-203],[-211,-212,-213]], [[-221,-222,-223],[-231,-232,-233]], [[-241,-242,-243],[-251,-252,-253]]  ] ], \
[ [[[ 11,12,13],[111,112,113]], [[121,122,123],[131,132,133]], [[141,142,143],[151,152,153]]  ], [ [[1101,1102,1103],[1111,1112,1113]], [[1121,1122,1123],[1131,1132,1133]], [[1141,1142,1143],[1151,1152,1153]]  ], [ [[1201,1202,1203],[1211,1212,1213]], [[1221,1222,1223],[1231,1232,1233]], [[1241,1242,1243],[1251,1252,1253]]  ] ] ]
a_r4_in=[ \
[ [ [[ 1./1,1./2,1./3],[1./11,1./12,1./13]], [[1./21,1./22,1./23],[1./31,1./32,1./33]], [[1./41,1./42,1./43],[1./51,1./52,1./53]]  ], [ [[1./101,1./102,1./103],[1./111,1./112,1./113]], [[1./121,1./122,1./123],[1./131,1./132,1./133]], [[1./141,1./142,1./143],[1./151,1./152,1./153]]  ], [ [[1./201,1./202,1./203],[1./211,1./212,1./213]], [[1./221,1./222,1./223],[1./231,1./232,1./233]], [[1./241,1./242,1./243],[1./251,1./252,1./253]]  ] ], \
[ [ [[ -1./1,-1./2,-1./3],[-1./11,-1./12,-1./13]], [[-1./21,-1./22,-1./23],[-1./31,-1./32,-1./33]], [[-1./41,-1./42,-1./43],[-1./51,-1./52,-1./53]]  ], [ [[-1./101,-1./102,-1./103],[-1./111,-1./112,-1./113]], [[-1./121,-1./122,-1./123],[-1./131,-1./132,-1./133]], [[-1./141,-1./142,-1./143],[1./-151,-1./152,-1./153]]  ], [ [[-1./201,-1./202,-1./203],[-1./211,-1./212,-1./213]], [[-1./221,-1./222,-1./223],[-1./231,-1./232,-1./233]], [[-1./241,-1./242,-1./243],[-1./251,-1./252,-1./253]]  ] ], \
[ [[[ 1./11,1./12,1./13],[1./111,1./112,1./113]], [[1./121,1./122,1./123],[1./131,1./132,1./133]], [[1./141,1./142,1./143],[1./151,1./152,1./153]]  ], [ [[1./1101,1./1102,1./1103],[1./1111,1./1112,1./1113]], [[1./1121,1./1122,1./1123],[1./1131,1./1132,1./1133]], [[1./1141,1./1142,1./1143],[1./1151,1./1152,1./1153]]  ], [ [[1./1201,1./1202,1./1203],[1./1211,1./1212,1./1213]], [[1./1221,1./1222,1./1223],[1./1231,1./1232,1./1233]], [[1./1241,1./1242,1./1243],[1./1251,1./1252,1./1253]]  ] ] ]

def prepareArg(val,ex,wh):
    if ex=="Array":
        out=val[0]
    else:
        if ex=="Expanded":
            exx=True
        else:
            exx=False
        out=Data(val[0],what=wh,expand=exx)
        if ex=="Tagged0":
            out.tag()
        elif ex=="Tagged1":
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
        raise SystemError,"@@ %s: error is too large"%(text)
        sys.exit(1)

#
# test slice setting:
#

msh=bruce.Rectangle(20,6)
for wh in [ContinuousFunction(msh),Function(msh)]:

  print wh

  for ex1 in ["Constant","Expanded","Tagged0","Tagged1","Tagged2"]:

    for ex2 in ["Array","Constant","Expanded","Tagged0","Tagged1","Tagged2"]:

          print "Slice setting: ", ex1, ",", ex2, ":"

          # rank 1

          arrays_in=turnToArray(a_r1_in,ex2)

          arg2=prepareArg(a_r1,ex1,wh) 
          arrays2=turnToArray(a_r1,ex1)
          a_in=[arrays_in[0][:], \
                arrays_in[1][:], \
                arrays_in[2][:]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[:]=expArg
          arrays2[0][:]=a_in[0]
          arrays2[1][:]=a_in[1]
          arrays2[2][:]=a_in[2]
          checkResult("slicing, set: rank=1 [:] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r1,ex1,wh) 
          arrays2=turnToArray(a_r1,ex1)
          a_in=[arrays_in[0][1], \
                arrays_in[1][1], \
                arrays_in[2][1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1]=expArg
          arrays2[0][1]=a_in[0]
          arrays2[1][1]=a_in[1]
          arrays2[2][1]=a_in[2]
          checkResult("slicing, set: rank=1 [1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r1,ex1,wh) 
          arrays2=turnToArray(a_r1,ex1)
          a_in=[arrays_in[0][1:3], \
                arrays_in[1][1:3], \
                arrays_in[2][1:3]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:3]=expArg
          arrays2[0][1:3]=a_in[0]
          arrays2[1][1:3]=a_in[1]
          arrays2[2][1:3]=a_in[2]
          checkResult("slicing, set: rank=1 [1:3] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r1,ex1,wh) 
          arrays2=turnToArray(a_r1,ex1)
          a_in=[arrays_in[0][:2], \
                arrays_in[1][:2], \
                arrays_in[2][:2]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[:2]=expArg
          arrays2[0][:2]=a_in[0]
          arrays2[1][:2]=a_in[1]
          arrays2[2][:2]=a_in[2]
          checkResult("slicing, set: rank=1 [:2] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r1,ex1,wh) 
          arrays2=turnToArray(a_r1,ex1)
          a_in=[arrays_in[0][2:], \
                arrays_in[1][2:], \
                arrays_in[2][2:]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[2:]=expArg
          arrays2[0][2:]=a_in[0]
          arrays2[1][2:]=a_in[1]
          arrays2[2][2:]=a_in[2]
          checkResult("slicing, set: rank=1 [2:] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          # rank 4

          arrays_in=turnToArray(a_r4_in,ex2)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][:], \
                arrays_in[1][:], \
                arrays_in[2][:]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[:]=expArg
          arrays2[0][:]=a_in[0]
          arrays2[1][:]=a_in[1]
          arrays2[2][:]=a_in[2]
          checkResult("slicing, set: rank=4 [:] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1], \
                arrays_in[1][1], \
                arrays_in[2][1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1]=expArg
          arrays2[0][1]=a_in[0]
          arrays2[1][1]=a_in[1]
          arrays2[2][1]=a_in[2]
          checkResult("slicing, set: rank=4 [1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1:3], \
                arrays_in[1][1:3], \
                arrays_in[2][1:3]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:3]=expArg
          arrays2[0][1:3]=a_in[0]
          arrays2[1][1:3]=a_in[1]
          arrays2[2][1:3]=a_in[2]
          checkResult("slicing, set: rank=4 [1:3] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][:,:], \
                arrays_in[1][:,:], \
                arrays_in[2][:,:]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[:,:]=expArg
          arrays2[0][:,:]=a_in[0]
          arrays2[1][:,:]=a_in[1]
          arrays2[2][:,:]=a_in[2]
          checkResult("slicing, set: rank=4 [:,:] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][:,1], \
                arrays_in[1][:,1], \
                arrays_in[2][:,1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[:,1]=expArg
          arrays2[0][:,1]=a_in[0]
          arrays2[1][:,1]=a_in[1]
          arrays2[2][:,1]=a_in[2]
          checkResult("slicing, set: rank=4 [:,1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][:,1:3], \
                arrays_in[1][:,1:3], \
                arrays_in[2][:,1:3]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[:,1:3]=expArg
          arrays2[0][:,1:3]=a_in[0]
          arrays2[1][:,1:3]=a_in[1]
          arrays2[2][:,1:3]=a_in[2]
          checkResult("slicing, set: rank=4 [:,1:3] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1:2,1:3], \
                arrays_in[1][1:2,1:3], \
                arrays_in[2][1:2,1:3]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:2,1:3]=expArg
          arrays2[0][1:2,1:3]=a_in[0]
          arrays2[1][1:2,1:3]=a_in[1]
          arrays2[2][1:2,1:3]=a_in[2]
          checkResult("slicing, set: rank=4 [1:2,1:3] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1:2,1:3,1,1], \
                arrays_in[1][1:2,1:3,1,1], \
                arrays_in[2][1:2,1:3,1,1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:2,1:3,1,1]=expArg
          arrays2[0][1:2,1:3,1,1]=a_in[0]
          arrays2[1][1:2,1:3,1,1]=a_in[1]
          arrays2[2][1:2,1:3,1,1]=a_in[2]
          checkResult("slicing, set: rank=4 [1:2,1:3,1,1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1:2,1:3,:,1], \
                arrays_in[1][1:2,1:3,:,1], \
                arrays_in[2][1:2,1:3,:,1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:2,1:3,:,1]=expArg
          arrays2[0][1:2,1:3,:,1]=a_in[0]
          arrays2[1][1:2,1:3,:,1]=a_in[1]
          arrays2[2][1:2,1:3,:,1]=a_in[2]
          checkResult("slicing, set: rank=4 [1:2,1:3,:,1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1:2,1:3,0,:], \
                arrays_in[1][1:2,1:3,0,:], \
                arrays_in[2][1:2,1:3,0,:]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:2,1:3,0,:]=expArg
          arrays2[0][1:2,1:3,0,:]=a_in[0]
          arrays2[1][1:2,1:3,0,:]=a_in[1]
          arrays2[2][1:2,1:3,0,:]=a_in[2]
          checkResult("slicing, set: rank=4 [1:2,1:3,0,:] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1:2,1:3,0,1], \
                arrays_in[1][1:2,1:3,0,1], \
                arrays_in[2][1:2,1:3,0,1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1:2,1:3,0,1]=expArg
          arrays2[0][1:2,1:3,0,1]=a_in[0]
          arrays2[1][1:2,1:3,0,1]=a_in[1]
          arrays2[2][1:2,1:3,0,1]=a_in[2]
          checkResult("slicing, set: rank=4 [1:2,1:3,0,1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

          arg2=prepareArg(a_r4,ex1,wh) 
          arrays2=turnToArray(a_r4,ex1)
          a_in=[arrays_in[0][1,1,0,1], \
                arrays_in[1][1,1,0,1], \
                arrays_in[2][1,1,0,1]]
          expArg=prepareArg(a_in,ex2,wh)
          arg2[1,1,0,1]=expArg
          arrays2[0][1,1,0,1]=a_in[0]
          arrays2[1][1,1,0,1]=a_in[1]
          arrays2[2][1,1,0,1]=a_in[2]
          checkResult("slicing, set: rank=4 [1,1,0,1] "+ex1+","+ex2, \
                      arg2, \
                      arrays2[0], \
                      arrays2[1], \
                      arrays2[2], \
                      wh)

sys.exit(0)
# end
