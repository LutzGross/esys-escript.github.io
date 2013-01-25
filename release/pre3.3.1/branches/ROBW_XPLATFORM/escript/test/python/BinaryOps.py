"""

Test binary operations.

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

Tag1=10
Tag2=20

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
         print "**** %s : error is too large"%(text)
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
# test binary operators:
#

msh=bruce.Rectangle(20,6)
for wh in [ContinuousFunction(msh),Function(msh)]:

  print wh

  for ex1 in ["Array","Constant","Expanded","Tagged0","Tagged1","Tagged2"]:

    for ex2 in ["Array","Constant","Expanded","Tagged0","Tagged1","Tagged2"]:

      if ex1=="Array" and ex2=="Array":
        continue

      print "Binary ops: ", ex1, ",", ex2, ":"

      for a1 in arglist:

        arg1=prepareArg(a1,ex1,wh)
        arrays1=turnToArray(a1,ex1)
        if isScalar(arg1):
           t1="(scalar)"
        else:
           t1=""

        for a2 in arglist:

          arg2=prepareArg(a2,ex2,wh)
          arrays2=turnToArray(a2,ex2)
          if isScalar(arg2):
             t2="(scalar)"
          else:
             t2=""

          # the shape must match or at least one argument is scalar:
          if (getRank(arg1)==getRank(arg2)) or isScalar(arg1) or isScalar(arg2):

            # sum 
            checkResult(ex1+t1+"+"+ex2+t2, \
                               arg1+arg2, \
                               arrays1[0]+arrays2[0], \
                               arrays1[1]+arrays2[1], \
                               arrays1[2]+arrays2[2], \
                               wh)

            # sub
            checkResult(ex1+t1+"-"+ex2+t2, \
                               arg1-arg2, \
                               arrays1[0]-arrays2[0], \
                               arrays1[1]-arrays2[1], \
                               arrays1[2]-arrays2[2], \
                               wh)

            # mul
            checkResult(ex1+t1+"*"+ex2+t2, \
                               arg1*arg2, \
                               arrays1[0]*arrays2[0], \
                               arrays1[1]*arrays2[1], \
                               arrays1[2]*arrays2[2], \
                               wh)

            # div
            checkResult(ex1+t1+"/"+ex2+t2, \
                               arg1/arg2, \
                               arrays1[0]/arrays2[0], \
                               arrays1[1]/arrays2[1], \
                               arrays1[2]/arrays2[2], \
                               wh)

            # pow 
            if isinstance(arg1,Data):
              a=arg1
            else:
              a=Data(value=arg1,what=arg2.getFunctionSpace())
            checkResult(ex1+t1+"^"+ex2+t2, \
                               a**arg2, \
                               arrays1[0]**arrays2[0], \
                               arrays1[1]**arrays2[1], \
                               arrays1[2]**arrays2[2], \
                               wh)

            # test inplace operations

            # arg2 must be class Data
            if isinstance(arg2,Data):

              # if arg1 is expanded arg2 must be expanded
              if not (not ex2=="Expanded" and ex1=="Expanded") :

                # if arg2 is scalar arg1 must be scalar:
                if not (isScalar(arg2) and not isScalar(arg1)):

                  # inplace add:
                  arrays2[0]+=arrays1[0]
                  arrays2[1]+=arrays1[1]
                  arrays2[2]+=arrays1[2]
                  arg2+=arg1
                  checkResult(ex2+t2+"+="+ex1+t1, \
                                 arg2, \
                                 arrays2[0], \
                                 arrays2[1], \
                                 arrays2[2], \
                                 wh)

                  # inplace sub:
                  arrays2[0]-=arrays1[0]
                  arrays2[1]-=arrays1[1]
                  arrays2[2]-=arrays1[2]
                  arg2-=arg1
                  checkResult(ex2+t2+"-="+ex1+t1, \
                                 arg2, \
                                 arrays2[0], \
                                 arrays2[1], \
                                 arrays2[2], \
                                 wh)

                  # inplace mul:
                  arrays2[0]*=arrays1[0]
                  arrays2[1]*=arrays1[1]
                  arrays2[2]*=arrays1[2]
                  arg2*=arg1
                  checkResult(ex2+t2+"*="+ex1+t1, \
                                 arg2, \
                                 arrays2[0], \
                                 arrays2[1], \
                                 arrays2[2], \
                                 wh)

                  # inplace div:
                  arrays2[0]/=arrays1[0]
                  arrays2[1]/=arrays1[1]
                  arrays2[2]/=arrays1[2]
                  arg2/=arg1
                  checkResult(ex2+t2+"/="+ex1+t1, \
                                 arg2, \
                                 arrays2[0], \
                                 arrays2[1], \
                                 arrays2[2], \
                                 wh)

sys.exit(0)
# end
