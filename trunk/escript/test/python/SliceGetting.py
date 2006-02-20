"""

Tests slice getting for Data objects.

Version $Id$

"""

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

# a rank 1 test argument
a_r1=[ [1,2,3], [-1,-2,-3], [100,200,300] ]

# a rank 4 test argument
a_r4=[ \
[ [ [[ 1,2,3],[11,12,13]], [[21,22,23],[31,32,33]], [[41,42,43],[51,52,53]]  ], [ [[101,102,103],[111,112,113]], [[121,122,123],[131,132,133]], [[141,142,143],[151,152,153]]  ], [ [[201,202,203],[211,212,213]], [[221,222,223],[231,232,233]], [[241,242,243],[251,252,253]]  ] ], \
[ [ [[ -1,-2,-3],[-11,-12,-13]], [[-21,-22,-23],[-31,-32,-33]], [[-41,-42,-43],[-51,-52,-53]]  ], [ [[-101,-102,-103],[-111,-112,-113]], [[-121,-122,-123],[-131,-132,-133]], [[-141,-142,-143],[-151,-152,-153]]  ], [ [[-201,-202,-203],[-211,-212,-213]], [[-221,-222,-223],[-231,-232,-233]], [[-241,-242,-243],[-251,-252,-253]]  ] ], \
[ [[[ 11,12,13],[111,112,113]], [[121,122,123],[131,132,133]], [[141,142,143],[151,152,153]]  ], [ [[1101,1102,1103],[1111,1112,1113]], [[1121,1122,1123],[1131,1132,1133]], [[1141,1142,1143],[1151,1152,1153]]  ], [ [[1201,1202,1203],[1211,1212,1213]], [[1221,1222,1223],[1231,1232,1233]], [[1241,1242,1243],[1251,1252,1253]]  ] ] ]

def prepareArg(val,ex,wh):
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

#
# ==============================================================
#
# test slice getting:
#

msh=bruce.Rectangle(20,6)

for wh in [ContinuousFunction(msh), Function(msh)]:

  print wh

  for ex1 in ["Constant", "Expanded", "Tagged1", "Tagged2"]:

    print "Slice getting: ", ex1

    # rank 1

    arg=prepareArg(a_r1,ex1,wh)
    arrays=turnToArray(a_r1,ex1)

    checkResult("slicing: rank=1 [:] "+ex1, \
                  arg[:], \
                  arrays[0][:], \
                  arrays[1][:], \
                  arrays[2][:], \
                  wh)

    checkResult("slicing: rank=1 [1] "+ex1, \
                  arg[1], \
                  arrays[0][1], \
                  arrays[1][1], \
                  arrays[2][1], \
                  wh)

    checkResult("slicing: rank=1 [1:3] "+ex1, \
                  arg[1:3], \
                  arrays[0][1:3], \
                  arrays[1][1:3], \
                  arrays[2][1:3], \
                  wh)

    checkResult("slicing: rank=1 [2:] "+ex1, \
                  arg[2:], \
                  arrays[0][2:], \
                  arrays[1][2:], \
                  arrays[2][2:], \
                  wh)

    checkResult("slicing: rank=1 [:1] "+ex1, \
                  arg[:1], \
                  arrays[0][:1], \
                  arrays[1][:1], \
                  arrays[2][:1], \
                  wh)

    # rank 4

    arg=prepareArg(a_r4,ex1,wh)
    arrays=turnToArray(a_r4,ex1)

    checkResult("slicing: rank=4 [:] "+ex1, \
                  arg[:], \
                  arrays[0][:], \
                  arrays[1][:], \
                  arrays[2][:], \
                  wh)

    checkResult("slicing: rank=4 [1] "+ex1, \
                  arg[1], \
                  arrays[0][1], \
                  arrays[1][1], \
                  arrays[2][1], \
                  wh)

    checkResult("slicing: rank=4 [1,1,1,1] "+ex1, \
                  arg[1,1,1,1], \
                  arrays[0][1,1,1,1], \
                  arrays[1][1,1,1,1], \
                  arrays[2][1,1,1,1], \
                  wh)

    checkResult("slicing: rank=4 [1:3] "+ex1, \
                  arg[1:3], \
                  arrays[0][1:3], \
                  arrays[1][1:3], \
                  arrays[2][1:3], \
                  wh)

    checkResult("slicing: rank=4 [:2] "+ex1, \
                  arg[:2], \
                  arrays[0][:2], \
                  arrays[1][:2], \
                  arrays[2][:2], \
                  wh)

    checkResult("slicing: rank=4 [1:] "+ex1, \
                  arg[1:], \
                  arrays[0][1:], \
                  arrays[1][1:], \
                  arrays[2][1:], \
                  wh)

    checkResult("slicing: rank=4 [:,:] "+ex1, \
                  arg[:,:], \
                  arrays[0][:,:], \
                  arrays[1][:,:], \
                  arrays[2][:,:], \
                  wh)

    checkResult("slicing: rank=4 [:,1] "+ex1, \
                  arg[:,1], \
                  arrays[0][:,1], \
                  arrays[1][:,1], \
                  arrays[2][:,1], \
                  wh)

    checkResult("slicing: rank=4 [:,1:3] "+ex1, \
                  arg[:,1:3], \
                  arrays[0][:,1:3], \
                  arrays[1][:,1:3], \
                  arrays[2][:,1:3], \
                  wh)

    checkResult("slicing: rank=4 [1:2,1:3] "+ex1, \
                  arg[1:2,1:3], \
                  arrays[0][1:2,1:3], \
                  arrays[1][1:2,1:3], \
                  arrays[2][1:2,1:3], \
                  wh)

    checkResult("slicing: rank=4 [:,1:2,1:2] "+ex1, \
                  arg[:,1:2,1:2], \
                  arrays[0][:,1:2,1:2], \
                  arrays[1][:,1:2,1:2], \
                  arrays[2][:,1:2,1:2], \
                  wh)

    checkResult("slicing: rank=4 [:,:,1:2,1:3] "+ex1, \
                  arg[:,:,1:2,1:3], \
                  arrays[0][:,:,1:2,1:3], \
                  arrays[1][:,:,1:2,1:3], \
                  arrays[2][:,:,1:2,1:3], \
                  wh)

    checkResult("slicing: rank=4 [1:2,1:3,1,1] "+ex1, \
                  arg[1:2,1:3,1,1], \
                  arrays[0][1:2,1:3,1,1], \
                  arrays[1][1:2,1:3,1,1], \
                  arrays[2][1:2,1:3,1,1], \
                  wh)

    checkResult("slicing: rank=4 [1:2,1:3,:,1] "+ex1, \
                  arg[1:2,1:3,:,1], \
                  arrays[0][1:2,1:3,:,1], \
                  arrays[1][1:2,1:3,:,1], \
                  arrays[2][1:2,1:3,:,1], \
                  wh)

    checkResult("slicing: rank=4 [1:2,1:3,0,:] "+ex1, \
                  arg[1:2,1:3,0,:], \
                  arrays[0][1:2,1:3,0,:], \
                  arrays[1][1:2,1:3,0,:], \
                  arrays[2][1:2,1:3,0,:], \
                  wh)

    checkResult("slicing: rank=4 [1:2,1:3,0,1] "+ex1, \
                  arg[1:2,1:3,0,1], \
                  arrays[0][1:2,1:3,0,1], \
                  arrays[1][1:2,1:3,0,1], \
                  arrays[2][1:2,1:3,0,1], \
                  wh)

    checkResult("slicing: rank=4 [1,1,0,1] "+ex1, \
                  arg[1,1,0,1], \
                  arrays[0][1,1,0,1], \
                  arrays[1][1,1,0,1], \
                  arrays[2][1,1,0,1], \
                  wh)

sys.exit(0)
# end
