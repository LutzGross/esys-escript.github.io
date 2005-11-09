import sys
import unittest
import os

from esys.escript import *
from esys import bruce

import numarray
from numarray import array,Float64,ones,greater

"""

Miscellaneous escript/Data tests.

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

msh=bruce.Rectangle()

for wh in [ContinuousFunction(msh),Function(msh)]:

  print wh

  for ex in ["Constant","Expanded"]:

    for a in arglist:

      print "\n", ex, a, "==>"

      arg=prepareArg(a,ex,wh)

      print "\n\nTests of copy method:"

      arg_copy = Data()
      arg_copy.copy(arg)
      print arg
      print arg_copy

      print "\n\nTests of conversion to numarray:"

      narray1 = arg.convertToNumArray()
      narray2 = arg.convertToNumArrayFromSampleNo(0)
      narray3 = arg.convertToNumArrayFromDPNo(0,0)

      print "arg.convertToNumArray()"
      print narray1

      print "arg.convertToNumArrayFromSampleNo(0)"
      print narray2

      print "arg.convertToNumArrayFromDPNo(0,0)"
      print narray3

      if (ex == "Constant"):
        print "\n\nTest of getTagNumber:"
        arg_copy = Data()
        arg_copy.copy(arg)
        arg_copy.tag()
        for dpno in range(narray1.shape[0]):
          print arg_copy.getTagNumber(dpno)
          print wh.getTagFromDataPointNo(dpno)

      if (ex == "Expanded"):
        print "\n\nTests of get/setRefValue functions:"
        result_array = numarray.array(a)
        print result_array
        result_array+=1
        arg.setRefValue(0,result_array)
        arg.getRefValue(0,result_array)
        print result_array

      print "\n\nTests of misc python functions:"

      print "\nmindp:"
      print arg.mindp()

      print "\nabs:"
      print arg.abs()

      print "\nmaxval:"
      print arg.maxval()

      print "\nminval:"
      print arg.minval()

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

      print "\nasin"
      print arg.asin()

      print "\nacos"
      print arg.acos()

      print "\natan"
      print arg.atan()

      print "\nsinh"
      print arg.sinh()

      print "\ncosh"
      print arg.cosh()

      print "\ntanh"
      print arg.tanh()

      print "\nasinh"
      print arg.asinh()

      print "\nacosh"
      print arg.acosh()

      print "\natanh"
      print arg.atanh()

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

sys.exit(0)
# end
