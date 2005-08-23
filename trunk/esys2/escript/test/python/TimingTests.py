import sys
import unittest
import os

import time

from escript.escript import *
from finley import finley

import numarray
from numarray import array,Float64,ones,greater

"""
Miscellaneous escript/Data timing tests.

Version $Id$
"""

arglist = [ \
3.0, \
[3,4], \
[[1,2],[3,4]], \
[[15,8],[12,8]], \
[[[15,8],[12,8]],[[-9,9],[13,8]]] \
]

def turnToArray(val):
     out=array(val,Float64)
     return out

def prepareArg(val,ex,wh):
     if ex=="Expanded":
         exx=True
     else:
         exx=False
     out=Data(val,what=wh,expand=exx)
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
       print "**** %s : error is too large"%(text)

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

print "\n\n"

msh=finley.Rectangle(1000,1000,1)

for wh in [Function(msh)]:

  print wh

  for ex in ["Expanded"]:

    for a in arglist:

      print "\n", ex, a, "==>"

      arg=prepareArg(a,ex,wh)

      starttime = time.clock()

      print "\nabs:",
      arg_abs = arg.abs()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nmaxval:",
      arg_maxval = arg.maxval()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nminval:",
      arg_minval = arg.minval()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nmindp:",
      arg_mindp = arg.mindp()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nlength",
      arg_length = arg.length()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\ntrace",
      arg_trace = arg.trace()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nsign",
      arg_sin = arg.sign()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nexp",
      arg_exp =  arg.exp()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nsqrt",
      arg_sqrt = arg.sqrt()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nneg",
      arg_neg = arg.neg()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\npos",
      arg_pos = arg.pos()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nsin",
      arg_sin = arg.sin()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\ncos",
      arg_cos = arg.cos()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\ntan",
      arg_tab = arg.tan()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nlog",
      arg_log = arg.log()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nln",
      arg_ln = arg.ln()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nLsup",
      arg_Lsup = arg.Lsup()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nLinf",
      arg_Linf = arg.Linf()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\nsup",
      arg_sup = arg.sup()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

      print "\ninf",
      arg_inf = arg.inf()

      stoptime = time.clock()
      elapsed = stoptime - starttime
      starttime = time.clock()
      print elapsed

sys.exit(0)
# end
