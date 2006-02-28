"""

Load tests for escript/Data.

Version $Id$

"""

import sys
import unittest
import os

import time

from esys.escript import *
from esys import bruce

# ==============================================================

arglist = [ \
3.0, \
[3,4], \
[[1,2],[3,4]], \
[[[15,8],[12,8]],[[-9,9],[13,8]]], \
[[[[35,9],[82,1]],[[-3,5],[33,8]]],[[[95,8],[-18,2]],[[-2,7],[113,-8]]]] \
]

# ==============================================================

def prepareArg(val,ex,wh):
  if ex=="Expanded":
      exx=True
  else:
      exx=False
  out=Data(val,what=wh,expand=exx)
  if ex=="Tagged":
      out.tag()
  return out

def doTest(arg,ex,wh):

  arg._wherePositive()
  arg._whereZero()
  arg._trace()
  arg._log()
  arg._sin()
  arg._acosh()
  arg._Lsup()
  arg._maxval()
  arg._sign()

  arg+=arg
  arg-=arg
  arg*=arg
  arg/=arg

# ==============================================================

testNum = 0

totalTime = 0

for x0 in [10, 100, 1000]:
  for x1 in [10, 100, 1000]:

    print "#### x0:", x0, "#### x1:", x1, "####"

    msh=bruce.Rectangle(x0,x1)
    for wh in [ContinuousFunction(msh),Function(msh)]:

      print wh

      for ex in ["Constant","Tagged","Expanded"]:

        for a in arglist:

          arg=prepareArg(a,ex,wh)

          testNum+=1
          print "Test", testNum, ": ----------------------------------------------"
          print ex, "Rank", arg.getRank()

          starttime = time.clock()

          for i in range(1000):

            doTest(arg,ex,wh)

          stoptime = time.clock()
          testElapsed = stoptime - starttime
          totalTime += testElapsed

          print "Test elapsed time: ", testElapsed

print "Total elapsed time: ", totalTime

# end
