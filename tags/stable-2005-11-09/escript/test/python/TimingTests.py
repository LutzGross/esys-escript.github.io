"""
Miscellaneous escript/Data timing tests.

Version $Id$
"""

import sys
import os
import time

from esys.escript import *
from esys import bruce

#
# ================== data values to test with =========================

arglist = [ \
3.0, \
[3.0,4.0], \
[[1.0,2.0],[3.0,4.0]], \
[[[15.0,8.0],[12.0,8.0]],[[-9.0,9.0],[13.0,8.0]]] \
]

testlist = [
"abs",
"maxval",
"minval",
"mindp",
"length",
"trace",
"sign",
"exp",
"sqrt",
"neg",
"pos",
"sin",
"cos",
"tan",
"log",
"ln",
"Lsup",
"Linf",
"sup",
"inf"
]

#
# ================== method definitions =========================

def prepareArg(val,ex,wh):
    if ex=="Expanded":
        exx=True
    else:
        exx=False
    out=Data(val,what=wh,expand=exx)
    return out

def getStartTime():
    return time.clock()

def calcElapsedTime(starttime):
    stoptime = time.clock()
    elapsed = stoptime - starttime
    print "\t\t", elapsed

def runTest(arg,test):
    print "\t\t", test,
    result = arg.__getattribute__(test)()
    del result

#
# ===================== main ==============================

msh=bruce.Rectangle(1000,1000)

for wh in [Function(msh),ContinuousFunction(msh)]:

  print "\n", wh, ":"

  for ex in ["Expanded"]:

    for a in arglist:

      print "\n\t", ex, a, "==>"
      print "\n\t\tFunction\tElapsed time"
      print "\t\t--------\t------------"

      arg=prepareArg(a,ex,wh)

      for test in testlist:

        starttime = getStartTime()

        runTest(arg,test)

        calcElapsedTime(starttime)

sys.exit(0)
# end