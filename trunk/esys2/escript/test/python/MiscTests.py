import sys
import unittest
import os

from esys.escript import *
from esys import finley

import numarray

"""

Miscellaneous escript/Data tests.

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

msh=finley.Rectangle(1,1,1)

for wh in [ContinuousFunction(msh),Function(msh)]:

  print wh.toString()

  for ex in ["Constant","Expanded"]:

    for a in arglist:

      print "\n", ex, a, "==>"

      arg=prepareArg(a,ex,wh)

      narry1 = arg.convertToNumArray()
      narry2 = arg.convertToNumArrayFromSampleNo(0)
      narry3 = arg.convertToNumArrayFromDPNo(0,0)

      print "\n\nTests of conversion to numarray:"

      print "arg.convertToNumArray()"
      print narry1

      print "arg.convertToNumArrayFromSampleNo(0)"
      print narry2

      print "arg.convertToNumArrayFromDPNo(0,0)"
      print narry3

      print "\n\nTests of misc python functions:"

      print "\nabs:"
      print arg.abs()

      print "\nmaxval:"
      print arg.maxval()

      print "\nminval:"
      print arg.minval()

      print "\nmindp:"
      print arg.mindp()

      print "\nlength"
      print arg.length()

      print "\ntrace"
      print arg.trace()

      print "\nsign"
      print arg.sign()

      print "\nexp"
      print arg.exp()

      print "\nsqrt"
      print arg.sqrt()

      print "\nneg"
      print arg.neg()

      print "\npos"
      print arg.pos()

      print "\nsin"
      print arg.sin()

      print "\ncos"
      print arg.cos()

      print "\ntan"
      print arg.tan()

      print "\nlog"
      print arg.log()

      print "\nln"
      print arg.ln()

      print "\nLsup"
      print arg.Lsup()

      print "\nLinf"
      print arg.Linf()

      print "\nsup"
      print arg.sup()

      print "\ninf"
      print arg.inf()

print "\n\nTests of archiveData and extractData:"

print "\nDataExpanded:"
archDataE=Data([[1.00001],[2.00001]],Function(msh),True)
archDataE.archiveData("data-archiveE")
exDataE=Data()
exDataE.extractData("data-archiveE",Function(msh))
exDataE.archiveData("data-archive2E");

print "\nDataTagged:"
archDataT=Data([[1.00001],[2.00001]],Function(msh))
archDataT.tag()
archDataT.archiveData("data-archiveT")
exDataT=Data()
exDataT.extractData("data-archiveT",Function(msh))
exDataT.archiveData("data-archive2T");

print "\nDataConstant:"
archDataC=Data([1.00001], Function(msh))
archDataC.archiveData("data-archiveC")
exDataC=Data()
exDataC.extractData("data-archiveC",Function(msh))
exDataC.archiveData("data-archive2C");

print "\nDataEmpty:"
archDataM=Data()
archDataM.archiveData("data-archiveE")
exDataM=Data()
exDataM.extractData("data-archiveE",FunctionSpace())
exDataM.archiveData("data-archive2E");

print "\n\nDo some Tagged data tests:"

print "\nCreate a Tagged data:"
tagData=Data([ [1.0,1.1],[2.0,2.1] ],Function(msh))
tagData.setTaggedValue(1,[[3.0,3.1],[4.0,4.1]])
tagData.setTaggedValue(2,[[5.0,5.1],[6.0,6.1]])
print tagData
print "\nSlice it [0:1,:]"
print tagData[0:1]
print "\nSlice it [0,0:1]"
print tagData[0]
print "\nSlice it [1:2]"
print tagData[1:2]
print "\nSlice it [1,1]"
print tagData[1,1]

print "\nSlice to it [1,1]"
tagData[0,0] = tagData[1,1]
print tagData

print "\nSlice to it [1:2,0:1]"
tagData[1:2,0:1] = tagData[0,1]
print tagData

print "\nSlice to it [0,1]"
tagData[0,1] = 9.0
print tagData

# end
