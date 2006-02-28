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
[[[15.0,8.0],[12.0,8.0]],[[-9.0,9.0],[13.0,8.0]]], \
[[[[14.0,7.0],[11.0,8.5]],[[-970,9.2],[18.0,8.0]]],[[[-4.4,7.0],[93.0,8.0]],[[-1.0,9.4],[12.0,9.0]]]] \
]

testlist = [
"_trace           ",
"_maxval          ",
"_minval          ",
"_wherePositive   ",
"_whereNegative   ",
"_whereNonNegative",
"_whereNonPositive",
"_whereZero       ",
"_whereNonZero    ",
"_sin             ",
"_cos             ",
"_tan             ",
"_asin            ",
"_acos            ",
"_atan            ",
"_sinh            ",
"_cosh            ",
"_tanh            ",
"_asinh           ",
"_acosh           ",
"_atanh           ",
"_exp             ",
"_sqrt            ",
"_log10           ",
"_log             ",
"_sign            ",
"_Lsup            ",
"_sup             ",
"_inf             "
]

#
# ================== method definitions =========================

def prepareArg(val,ex,wh):
    if ex=="Expanded":
        exx=True
    else:
        exx=False
    out=Data(val,what=wh,expand=exx)
    if ex=="Tagged":
        out.tag()
    return out

def getStartTime():
    return time.clock()

def calcElapsedTime(starttime):
    stoptime = time.clock()
    elapsed = stoptime - starttime
    print "\t\t", elapsed

def runTest(arg,test):
    print "\t\t", test,
    test_name = test.rstrip()
    result = arg.__getattribute__(test_name)()
    del result

#
# ===================== main ==============================

msh=bruce.Rectangle(1000,1000)

for wh in [Function(msh),ContinuousFunction(msh)]:

  print "\n", wh, ":"

  for ex in ["Constant", "Tagged", "Expanded"]:

    for a in arglist:

      arg=prepareArg(a,ex,wh)

      print "\n\t", ex, "Rank", arg.getRank(), "==>"
      print "\n\t\tFunction\t\tElapsed time"
      print "\t\t--------\t\t------------"

      for test in testlist:

        starttime = getStartTime()

        runTest(arg,test)

        calcElapsedTime(starttime)

sys.exit(0)
# end
