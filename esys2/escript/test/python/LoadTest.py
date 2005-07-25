"""

Load tests for escript/Data.

Version $Id$

"""

import sys
import unittest
import os

import time

from escript.escript import *
from finley import finley

# ==============================================================

arglist = [ \
[3,4], \
[[1,2],[3,4]], \
[[15,8],[12,8]], \
[[[15,8],[12,8]],[[-9,9],[13,8]]], \
3.0 \
]

# ==============================================================

def prepareArg(val,ex,wh):
  if ex=="Expanded":
      exx=True
  else:
      exx=False
  out=Data(val,what=wh,expand=exx)
  return out

def doTest(a,ex,wh):

  arg=prepareArg(a,ex,wh)

  arg.wherePositive()
  arg.whereZero()
  arg.trace()
  arg.log()
  arg.Linf()
  arg.maxval()
  arg.sign()

  arg+=arg
  arg-=arg
  arg*=arg
  arg/=arg

# ==============================================================

testNum = 0

totalTime = 0

for x0 in [1, 10, 100]:
  for x1 in [1, 10, 100]:

    print "#### x0:", x0, "#### x1:", x1, "####"
    msh=finley.Rectangle(x0,x1)

    for wh in [ContinuousFunction(msh),Function(msh)]:

      for ex in ["Constant","Expanded"]:

        for a in arglist:

          testNum+=1
          print testNum, ": ----------------------------------------------"

          testElapsed = 0

          for j in range(10):

            starttime = time.clock()

            for i in range(1000):

              doTest(a,ex,wh)

            stoptime = time.clock()
            elapsed = stoptime - starttime
            testElapsed += elapsed
            totalTime += elapsed

            print elapsed

          print "Test elapsed time: ", testElapsed

print "Total elapsed time: ", totalTime

# end
