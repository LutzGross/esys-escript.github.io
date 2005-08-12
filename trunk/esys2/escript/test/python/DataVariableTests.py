import sys
import unittest
import os

from escript.escript import *
from finley import finley

import numarray

"""

Tests for DataVariable

Version $Id$

"""

from numarray import array,Float64,ones,greater


#  list of arguments: a list item has the form [a0,a1,a2]
#  what a0 is the default value and a1 is used for tag Tag1
#  and a2 for tag2. a0,a1,a2 are converted into numarrays.

arglist = [ \
[3,4], \
[[1,2],[3,4]], \
[[15,8],[12,8]], \
[[[15,8],[12,8]],[[-9,9],[13,8]]], \
3.0 \
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

msh=finley.Rectangle(2,3,1)

print "\n\nTests of DataVariable"

d1 = Data(1.0,Function(msh),True)
d2 = Data(2.0,Function(msh),True)
d3 = Data(3.0,Function(msh),True)
d4 = Data(4.0,Function(msh),True)

dv0 = DataVariable()
dv1 = DataVariable(d1)
dv2 = DataVariable(d2)
dv3 = DataVariable(d3)
dv4 = DataVariable(d4)

dv3.sum(dv4)
dv2.diff(dv3)
dv1.sum(dv2)

d5 = dv1.evaluate()

print d5.toString()

sys.exit(0)
# end
