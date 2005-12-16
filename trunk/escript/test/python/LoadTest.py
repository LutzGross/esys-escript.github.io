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

  arg._wherePositive()
  arg._whereZero()
  arg._trace()
  arg._log()
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

for x0 in [10, 100]:
  for x1 in [10, 100]:

    print "#### x0:", x0, "#### x1:", x1, "####"
    msh=bruce.Rectangle(x0,x1)

    for wh in [ContinuousFunction(msh),Function(msh)]:

      for ex in ["Constant","Expanded"]:

        for a in arglist:

          testNum+=1
          print testNum, ": ----------------------------------------------"
          print a, ex, wh

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
