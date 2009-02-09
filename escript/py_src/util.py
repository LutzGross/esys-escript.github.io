
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
Utility functions for escript

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
@var EPSILON: smallest positive value with 1.<1.+EPSILON
@var DBLE_MAX: largest positive float
"""

__author__="Lutz Gross, l.gross@uq.edu.au"


import math
import numarray
import escript
import os
from esys.escript import C_GeneralTensorProduct
from esys.escript import getVersion
from esys.escript import printParallelThreadCounts
from esys.escript import listEscriptParams

#=========================================================
#   some helpers:
#=========================================================
def getEpsilon():
     return escript.getMachinePrecision()
EPSILON=getEpsilon()

def getMaxFloat():
     return escript.getMaxFloat()
DBLE_MAX=getMaxFloat()


def getTagNames(domain):
    """
    Returns a list of tag names used by the domain.

    @param domain: a domain object
    @type domain: L{escript.Domain}
    @return: a list of tag names used by the domain
    @rtype: C{list} of C{str}
    """
    return [n.strip() for n in domain.showTagNames().split(",") ]

def insertTagNames(domain,**kwargs):
    """
    Inserts tag names into the domain.

    @param domain: a domain object
    @type domain: C{escript.Domain}
    @keyword <tag_name>: tag key assigned to <tag_name>
    @type <tag_name>: C{int}
    """
    for  k in kwargs:
         domain.setTagMap(k,kwargs[k])

def insertTaggedValues(target,**kwargs):
    """
    Inserts tagged values into the target using tag names.

    @param target: data to be filled by tagged values
    @type target: L{escript.Data}
    @keyword <tag_name>: value to be used for <tag_name>
    @type <tag_name>: C{float} or C{numarray.NumArray}
    @return: C{target}
    @rtype: L{escript.Data}
    """
    for k in kwargs:
        target.setTaggedValue(k,kwargs[k])
    return target

def saveVTK(filename,domain=None,**data):
    """
    Writes L{Data} objects and their mesh into a file using the VTK XML file
    format.

    Example::

        tmp=Scalar(..)
        v=Vector(..)
        saveVTK("solution.vtu", temperature=tmp, velocity=v)

    C{tmp} and C{v} are written into "solution.vtu" where C{tmp} is named
    "temperature" and C{v} is named "velocity".

    @param filename: file name of the output file
    @type filename: C{str}
    @param domain: domain of the L{Data} objects. If not specified, the domain
                   of the given L{Data} objects is used.
    @type domain: L{escript.Domain}
    @keyword <name>: writes the assigned value to the VTK file using <name> as
                     identifier
    @type <name>: L{Data} object
    @note: The data objects have to be defined on the same domain. They may not
           be in the same L{FunctionSpace} but one cannot expect that all
           L{FunctionSpace}s can be mixed. Typically, data on the boundary and
           data on the interior cannot be mixed.
    """
    new_data={}
    for n,d in data.items():
          if not d.isEmpty():
            fs=d.getFunctionSpace()
            domain2=fs.getDomain()
            if fs == escript.Solution(domain2):
               new_data[n]=interpolate(d,escript.ContinuousFunction(domain2))
            elif fs == escript.ReducedSolution(domain2):
               new_data[n]=interpolate(d,escript.ReducedContinuousFunction(domain2))
            else:
               new_data[n]=d
            if domain==None: domain=domain2
    if domain==None:
        raise ValueError,"saveVTK: no domain detected."
    domain.saveVTK(filename,new_data)

def saveDX(filename,domain=None,**data):
    """
    Writes L{Data} objects into a file using the OpenDX file format.

    Example::

        tmp=Scalar(..)
        v=Vector(..)
        saveDX("solution.dx", temperature=tmp, velocity=v)

    C{tmp} and C{v} are written into "solution.dx" where C{tmp} is named
    "temperature" and C{v} is named "velocity".

    @param filename: file name of the output file
    @type filename: C{str}
    @param domain: domain of the L{Data} objects. If not specified, the domain
                   of the given L{Data} objects is used.
    @type domain: L{escript.Domain}
    @keyword <name>: writes the assigned value to the DX file using <name> as
                     identifier. The identifier can be used to select the data
                     set when data are imported into DX.
    @type <name>: L{Data} object
    @note: The data objects have to be defined on the same domain. They may not
           be in the same L{FunctionSpace} but one cannot expect that all
           L{FunctionSpace}s can be mixed. Typically, data on the boundary and
           data on the interior cannot be mixed.
    """
    new_data={}
    for n,d in data.items():
          if not d.isEmpty():
            fs=d.getFunctionSpace()
            domain2=fs.getDomain()
            if fs == escript.Solution(domain2):
               new_data[n]=interpolate(d,escript.ReducedContinuousFunction(domain2))
            elif fs == escript.ReducedSolution(domain2):
               new_data[n]=interpolate(d,escript.ReducedContinuousFunction(domain2))
            elif fs == escript.ContinuousFunction(domain2):
               new_data[n]=interpolate(d,escript.ReducedContinuousFunction(domain2))
            else:
               new_data[n]=d
            if domain==None: domain=domain2
    if domain==None:
        raise ValueError,"saveDX: no domain detected."
    domain.saveDX(filename,new_data)

def saveESD(datasetName, dataDir=".", domain=None, **data):
    """
    Saves L{Data} objects to files and creates an I{escript dataset} (ESD) file
    for convenient processing/visualisation.

    Example::

        tmp = Scalar(..)
        v = Vector(..)
        saveESD("solution", "data", temperature=tmp, velocity=v)

    tmp, v and the domain are saved in native format in the "data"
    directory and the file "solution.esd" is created that refers to tmp by
    the name "temperature" and to v by the name "velocity".

    @param datasetName: name of the dataset, used to name the ESD file
    @type datasetName: C{str}
    @param dataDir: optional directory where the data files should be saved
    @type dataDir: C{str}
    @param domain: domain of the L{Data} object(s). If not specified, the
                   domain of the given L{Data} objects is used.
    @type domain: L{escript.Domain}
    @keyword <name>: writes the assigned value to the file using <name> as
                     identifier
    @type <name>: L{Data} object.
    @note: The data objects have to be defined on the same domain. They may not
           be in the same L{FunctionSpace} but one cannot expect that all
           L{FunctionSpace}s can be mixed. Typically, data on the boundary and
           data on the interior cannot be mixed.
    """
    new_data = {}
    for n,d in data.items():
          if not d.isEmpty(): 
            fs = d.getFunctionSpace() 
            domain2 = fs.getDomain()
            if fs == escript.Solution(domain2):
               new_data[n]=interpolate(d,escript.ContinuousFunction(domain2))
            elif fs == escript.ReducedSolution(domain2):
               new_data[n]=interpolate(d,escript.ReducedContinuousFunction(domain2))
            else:
               new_data[n]=d
            if domain==None: domain=domain2
    if domain==None:
        raise ValueError, "saveESD: no domain detected."

    if domain.onMasterProcessor() and not os.path.isdir(dataDir):
        os.mkdir(dataDir)

    meshFile = os.path.join(dataDir, datasetName+"_mesh")
    domain.dump(meshFile + ".nc")
    outputString = ""
    if domain.onMasterProcessor():
        outputString += "#escript datafile V1.0\n"
        # number of timesteps (currently only 1 is supported)
        outputString += "T=1\n"
        # name of the mesh file
        outputString += "M=%s\n" % meshFile
        # number of blocks (MPI size)
        outputString += "N=%d\n" % domain.getMPISize()

    # now add the variables
    for varName, d in new_data.items():
        varFile = os.path.join(dataDir, datasetName+"_"+varName)
        d.dump(varFile + ".nc")
        if domain.onMasterProcessor():
            outputString += "V=%s:%s\n" % (varFile, varName)

    if domain.onMasterProcessor():
        esdfile = open(datasetName+".esd", "w")
        esdfile.write(outputString)
        esdfile.close()

def kronecker(d=3):
   """
   Returns the kronecker S{delta}-symbol.

   @param d: dimension or an object that has the C{getDim} method defining the
             dimension
   @type d: C{int}, L{escript.Domain} or L{escript.FunctionSpace}
   @return: the object u of rank 2 with M{u[i,j]=1} for M{i=j} and M{u[i,j]=0}
            otherwise
   @rtype: C{numarray.NumArray} or L{escript.Data} of rank 2
   """
   return identityTensor(d)

def identity(shape=()):
   """
   Returns the C{shape} x C{shape} identity tensor.

   @param shape: input shape for the identity tensor
   @type shape: C{tuple} of C{int}
   @return: array whose shape is shape x shape where M{u[i,k]=1} for M{i=k} and
            M{u[i,k]=0} otherwise for len(shape)=1. If len(shape)=2:
            M{u[i,j,k,l]=1} for M{i=k and j=l} and M{u[i,j,k,l]=0} otherwise.
   @rtype: C{numarray.NumArray} of rank 1, rank 2 or rank 4
   @raise ValueError: if len(shape)>2
   """
   if len(shape)>0:
      out=numarray.zeros(shape+shape,numarray.Float64)
      if len(shape)==1:
          for i0 in range(shape[0]):
             out[i0,i0]=1.
      elif len(shape)==2:
          for i0 in range(shape[0]):
             for i1 in range(shape[1]):
                out[i0,i1,i0,i1]=1.
      else:
          raise ValueError,"identity: length of shape is restricted to 2."
   else:
      out=1.
   return out

def identityTensor(d=3):
   """
   Returns the C{d} x C{d} identity matrix.

   @param d: dimension or an object that has the C{getDim} method defining the
             dimension
   @type d: C{int}, L{escript.Domain} or L{escript.FunctionSpace}
   @return: the object u of rank 2 with M{u[i,j]=1} for M{i=j} and M{u[i,j]=0}
            otherwise
   @rtype: C{numarray.NumArray} or L{escript.Data} of rank 2
   """
   if isinstance(d,escript.FunctionSpace):
       return escript.Data(identity((d.getDim(),)),d)
   elif isinstance(d,escript.Domain):
       return identity((d.getDim(),))
   else:
       return identity((d,))

def identityTensor4(d=3):
   """
   Returns the C{d} x C{d} x C{d} x C{d} identity tensor.

   @param d: dimension or an object that has the C{getDim} method defining the
             dimension
   @type d: C{int} or any object with a C{getDim} method
   @return: the object u of rank 4 with M{u[i,j,k,l]=1} for M{i=k and j=l} and
            M{u[i,j,k,l]=0} otherwise
   @rtype: C{numarray.NumArray} or L{escript.Data} of rank 4
   """
   if isinstance(d,escript.FunctionSpace):
       return escript.Data(identity((d.getDim(),d.getDim())),d)
   elif isinstance(d,escript.Domain):
       return identity((d.getDim(),d.getDim()))
   else:
       return identity((d,d))

def unitVector(i=0,d=3):
   """
   Returns a unit vector u of dimension d whose non-zero element is at index i.

   @param i: index for non-zero element
   @type i: C{int}
   @param d: dimension or an object that has the C{getDim} method defining the
             dimension
   @type d: C{int}, L{escript.Domain} or L{escript.FunctionSpace}
   @return: the object u of rank 1 with M{u[j]=1} for M{j=index} and M{u[j]=0}
            otherwise
   @rtype: C{numarray.NumArray} or L{escript.Data} of rank 1
   """
   return kronecker(d)[i]

#=========================================================================
#   global reduction operations (these functions have no symbolic version)
#=========================================================================
def Lsup(arg):
    """
    Returns the Lsup-norm of argument C{arg}. This is the maximum absolute value
    over all data points. This function is equivalent to C{sup(abs(arg))}.

    @param arg: argument
    @type arg: C{float}, C{int}, L{escript.Data}, C{numarray.NumArray}
    @return: maximum value of the absolute value of C{arg}over all components
             and all data points
    @rtype: C{float}
    @raise TypeError: if type of C{arg}cannot be processed
    """
    if isinstance(arg,numarray.NumArray):
        return sup(abs(arg))
    elif isinstance(arg,escript.Data):
        return arg._Lsup()
    elif isinstance(arg,float):
        return abs(arg)
    elif isinstance(arg,int):
        return abs(float(arg))
    else:
      raise TypeError,"Lsup: Unknown argument type."

def sup(arg):
    """
    Returns the maximum value over all data points.

    @param arg: argument
    @type arg: C{float}, C{int}, L{escript.Data}, C{numarray.NumArray}
    @return: maximum value of C{arg}over all components and all data points
    @rtype: C{float}
    @raise TypeError: if type of C{arg}cannot be processed
    """
    if isinstance(arg,numarray.NumArray):
        return arg.max()
    elif isinstance(arg,escript.Data):
        return arg._sup()
    elif isinstance(arg,float):
        return arg
    elif isinstance(arg,int):
        return float(arg)
    else:
      raise TypeError,"sup: Unknown argument type."

def inf(arg):
    """
    Returns the minimum value over all data points.

    @param arg: argument
    @type arg: C{float}, C{int}, L{escript.Data}, C{numarray.NumArray}
    @return: minimum value of C{arg}over all components and all data points
    @rtype: C{float}
    @raise TypeError: if type of C{arg}cannot be processed
    """
    if isinstance(arg,numarray.NumArray):
        return arg.min()
    elif isinstance(arg,escript.Data):
        return arg._inf()
    elif isinstance(arg,float):
        return arg
    elif isinstance(arg,int):
        return float(arg)
    else:
      raise TypeError,"inf: Unknown argument type."


#=========================================================================
#   some little helpers
#=========================================================================
def getRank(arg):
    """
    Identifies the rank of the argument.

    @param arg: an object whose rank is to be returned
    @type arg: C{numarray.NumArray}, L{escript.Data}, C{float}, C{int},
               C{Symbol}
    @return: the rank of the argument
    @rtype: C{int}
    @raise TypeError: if type of C{arg}cannot be processed
    """

    if isinstance(arg,numarray.NumArray):
        return arg.rank
    elif isinstance(arg,escript.Data):
        return arg.getRank()
    elif isinstance(arg,float):
        return 0
    elif isinstance(arg,int):
        return 0
    elif isinstance(arg,Symbol):
        return arg.getRank()
    else:
      raise TypeError,"getRank: Unknown argument type."

def getShape(arg):
    """
    Identifies the shape of the argument.

    @param arg: an object whose shape is to be returned
    @type arg: C{numarray.NumArray}, L{escript.Data}, C{float}, C{int},
               C{Symbol}
    @return: the shape of the argument
    @rtype: C{tuple} of C{int}
    @raise TypeError: if type of C{arg}cannot be processed
    """

    if isinstance(arg,numarray.NumArray):
        return arg.shape
    elif isinstance(arg,escript.Data):
        return arg.getShape()
    elif isinstance(arg,float):
        return ()
    elif isinstance(arg,int):
        return ()
    elif isinstance(arg,Symbol):
        return arg.getShape()
    else:
      raise TypeError,"getShape: Cannot identify shape"

def pokeDim(arg):
    """
    Identifies the spatial dimension of the argument.

    @param arg: an object whose spatial dimension is to be returned
    @type arg: any
    @return: the spatial dimension of the argument, if available, or C{None}
    @rtype: C{int} or C{None}
    """

    if isinstance(arg,escript.Data):
        return arg.getFunctionSpace().getDim()
    elif isinstance(arg,Symbol):
        return arg.getDim()
    else:
        return None

def commonShape(arg0, arg1):
    """
    Returns a shape to which C{arg0} can be extended from the right and C{arg1}
    can be extended from the left.

    @param arg0: an object with a shape (see L{getShape})
    @param arg1: an object with a shape (see L{getShape})
    @return: the shape of C{arg0} or C{arg1} such that the left part equals the
             shape of C{arg0} and the right end equals the shape of C{arg1}
    @rtype: C{tuple} of C{int}
    @raise ValueError: if no shape can be found
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)<len(sh1):
       if not sh0==sh1[:len(sh0)]:
             raise ValueError,"argument 0 cannot be extended to the shape of argument 1"
       return sh1
    elif len(sh0)>len(sh1):
       if not sh1==sh0[:len(sh1)]:
             raise ValueError,"argument 1 cannot be extended to the shape of argument 0"
       return sh0
    else:
       if not sh0==sh1:
             raise ValueError,"argument 1 and argument 0 have not the same shape."
       return sh0

def commonDim(*args):
    """
    Identifies, if possible, the spatial dimension across a set of objects
    which may or may not have a spatial dimension.

    @param args: given objects
    @return: the spatial dimension of the objects with identifiable dimension
             (see L{pokeDim}). If none of the objects has a spatial dimension
             C{None} is returned.
    @rtype: C{int} or C{None}
    @raise ValueError: if the objects with identifiable dimension don't have
                       the same spatial dimension.
    """
    out=None
    for a in args:
       d=pokeDim(a)
       if not out==None:
          if not (d==None or out==d):
             raise ValueError,"dimension of arguments don't match"
       else:
          out=d
    return out

def testForZero(arg):
    """
    Tests if the argument is identical to zero.

    @param arg: the object to test for zero
    @type arg: typically C{numarray.NumArray}, L{escript.Data}, C{float}, C{int}
    @return: True if the argument is identical to zero, False otherwise
    @rtype: C{bool}
    """
    if isinstance(arg,numarray.NumArray):
       return not Lsup(arg)>0.
    elif isinstance(arg,escript.Data):
       return False
    elif isinstance(arg,float):
       return not Lsup(arg)>0.
    elif isinstance(arg,int):
       return not Lsup(arg)>0.
    elif isinstance(arg,Symbol):
       return False
    else:
       return False

def matchType(arg0=0.,arg1=0.):
    """
    Converts C{arg0} and C{arg1} both to the same type C{numarray.NumArray} or
    L{escript.Data} or, if one of C{arg0} or C{arg1} is of type L{Symbol}, the
    other one to be of type C{numarray.NumArray} or L{escript.Data}.

    @param arg0: first argument
    @type arg0: C{numarray.NumArray},L{escript.Data},C{float}, C{int}, C{Symbol}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray},L{escript.Data},C{float}, C{int}, C{Symbol}
    @return: a tuple representing C{arg0} and C{arg1} with the same type or
             with one of them being a L{Symbol}
    @rtype: C{tuple} of two C{numarray.NumArray}, two L{escript.Data},
            a C{Symbol} and one of the types C{numarray.NumArray} or
            L{escript.Data}
    @raise TypeError: if type of C{arg0} or C{arg1} cannot be processed
    """
    if isinstance(arg0,numarray.NumArray):
       if isinstance(arg1,numarray.NumArray):
          pass
       elif isinstance(arg1,escript.Data):
          arg0=escript.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg1=numarray.array(arg1,type=numarray.Float64)
       elif isinstance(arg1,int):
          arg1=numarray.array(float(arg1),type=numarray.Float64)
       elif isinstance(arg1,Symbol):
          pass
       else:
          raise TypeError,"function: Unknown type of second argument."
    elif isinstance(arg0,escript.Data):
       if isinstance(arg1,numarray.NumArray):
          arg1=escript.Data(arg1,arg0.getFunctionSpace())
       elif isinstance(arg1,escript.Data):
          pass
       elif isinstance(arg1,float):
          arg1=escript.Data(arg1,(),arg0.getFunctionSpace())
       elif isinstance(arg1,int):
          arg1=escript.Data(float(arg1),(),arg0.getFunctionSpace())
       elif isinstance(arg1,Symbol):
          pass
       else:
          raise TypeError,"function: Unknown type of second argument."
    elif isinstance(arg0,Symbol):
       if isinstance(arg1,numarray.NumArray):
          pass
       elif isinstance(arg1,escript.Data):
          pass
       elif isinstance(arg1,float):
          arg1=numarray.array(arg1,type=numarray.Float64)
       elif isinstance(arg1,int):
          arg1=numarray.array(float(arg1),type=numarray.Float64)
       elif isinstance(arg1,Symbol):
          pass
       else:
          raise TypeError,"function: Unknown type of second argument."
    elif isinstance(arg0,float):
       if isinstance(arg1,numarray.NumArray):
          arg0=numarray.array(arg0,type=numarray.Float64)
       elif isinstance(arg1,escript.Data):
          arg0=escript.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numarray.array(arg0,type=numarray.Float64)
          arg1=numarray.array(arg1,type=numarray.Float64)
       elif isinstance(arg1,int):
          arg0=numarray.array(arg0,type=numarray.Float64)
          arg1=numarray.array(float(arg1),type=numarray.Float64)
       elif isinstance(arg1,Symbol):
          arg0=numarray.array(arg0,type=numarray.Float64)
       else:
          raise TypeError,"function: Unknown type of second argument."
    elif isinstance(arg0,int):
       if isinstance(arg1,numarray.NumArray):
          arg0=numarray.array(float(arg0),type=numarray.Float64)
       elif isinstance(arg1,escript.Data):
          arg0=escript.Data(float(arg0),arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numarray.array(float(arg0),type=numarray.Float64)
          arg1=numarray.array(arg1,type=numarray.Float64)
       elif isinstance(arg1,int):
          arg0=numarray.array(float(arg0),type=numarray.Float64)
          arg1=numarray.array(float(arg1),type=numarray.Float64)
       elif isinstance(arg1,Symbol):
          arg0=numarray.array(float(arg0),type=numarray.Float64)
       else:
          raise TypeError,"function: Unknown type of second argument."
    else:
      raise TypeError,"function: Unknown type of first argument."

    return arg0,arg1

def matchShape(arg0,arg1):
    """
    Returns a representation of C{arg0} and C{arg1} which have the same shape.

    @param arg0: first argument
    @type arg0: C{numarray.NumArray},L{escript.Data},C{float}, C{int}, L{Symbol}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray},L{escript.Data},C{float}, C{int}, L{Symbol}
    @return: C{arg0} and C{arg1} where copies are returned when the shape has
             to be changed
    @rtype: C{tuple}
    """
    sh=commonShape(arg0,arg1)
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)<len(sh):
       return outer(arg0,numarray.ones(sh[len(sh0):],numarray.Float64)),arg1
    elif len(sh1)<len(sh):
       return arg0,outer(arg1,numarray.ones(sh[len(sh1):],numarray.Float64))
    else:
       return arg0,arg1

#=========================================================
#   symbolic tool box starts here:
#=========================================================
class Symbol(object):
   """
   Symbol class objects provide the same functionality as C{numarray.NumArray}
   and L{escript.Data} objects but they do not have a value and therefore
   cannot be plotted or visualized. The main purpose is the possibility to
   calculate derivatives with respect to other Symbols used to define a Symbol.

   """
   def __init__(self,shape=(),args=[],dim=None):
       """
       Creates an instance of a symbol of a given shape. The symbol may depend
       on a list of arguments C{args} which may be symbols or any other object.

       @param args: the arguments of the symbol
       @type args: C{list}
       @param shape: the shape of the symbol
       @type shape: C{tuple} of C{int}
       @param dim: spatial dimension of the symbol. If dim=C{None} the spatial
                   dimension is undefined.
       @type dim: C{None} or C{int}

       """
       if len(shape)>4:
           raise ValueError,"Symbol only supports tensors up to order 4"
       self.__args=args
       self.__shape=shape
       self.__dim=dim

   def getArgument(self,i=None):
       """
       Returns the i-th argument of the symbol.

       @param i: index of the argument requested
       @type i: C{int} or C{None}
       @raise IndexError: if the requested index does not exist
       @return: the value of the i-th argument or if i is not specified the
                list of all arguments
       @rtype: a single object or a list of objects
       """
       if i==None:
          return self.__args
       else:
          if i<0 or i>=len(self.__args):
             raise IndexError,"there are only %s arguments"%len(self.__args)
          return self.__args[i]

   def getRank(self):
       """
       Returns the rank of the symbol.

       @return: the rank of the symbol. This is length of the shape.
       @rtype: C{int}
       """
       return len(self.getShape())

   def getShape(self):
       """
       Returns the shape of the symbol.

       @return: the shape of the symbol
       @rtype: C{tuple} of C{int}
       """
       return self.__shape

   def getDim(self):
       """
       Returns the spatial dimension.

       @return: the symbol's spatial dimension
       @rtype: C{int} if the dimension is defined, C{None} otherwise
       """
       return self.__dim

   def __str__(self):
       """
       Returns a string representation of the symbol.

       @return: a string representation of the object
       @rtype: C{str}
       """
       args=[]
       for arg in self.getArgument():
          args.append(str(arg))
       try:
           out=self.getMyCode(args,format="str")
       except NotImplementedError:
           out="<Symbol %s>"%id(self)
       return out

   def getSubstitutedArguments(self,argvals):
       """
       Substitutes symbols in the arguments of this object and returns the
       result as a list.

       @param argvals: L{Symbol} and their substitutes. The L{Symbol} u in the
                       expression defining this object is replaced by
                       argvals[u].
       @type argvals: C{dict} with keywords of type L{Symbol}
       @rtype: C{list} of objects
       @return: list of the object assigned to the arguments through
                substitution or for the arguments which are not L{Symbol}s the
                value assigned to the argument at instantiation.
       """
       out=[]
       for a in self.getArgument():
          if isinstance(a,Symbol):
             out.append(a.substitute(argvals))
          else:
             out.append(a)
       return out

   def getDifferentiatedArguments(self,arg):
       """
       Applies differentials to the arguments of this object and returns the
       result as a list.

       @param arg: the derivative is calculated with respect to C{arg}
       @type arg: typically L{escript.Symbol} but can also be C{float},
                  L{escript.Data}, C{numarray.NumArray} depending on the
                  involved functions and data
       @rtype: C{list} of objects
       @return: list of object obtained by calculating the derivatives of the
                arguments with respect to C{arg}
       """
       out=[]
       for a in self.getArgument():
          if isinstance(a,Symbol):
             out.append(a.substitute(argvals))
          else:
              s=getShape(s)+arg.getShape()
              if len(s)>0:
                 out.append(numarray.zeros(s),numarray.Float64)
              else:
                 out.append(a)
       return out

   def isAppropriateValue(self,arg):
      """
      Checks if the given argument C{arg} can be used as a substitution for
      this object. The method checks the shape of C{arg} and, if the spatial
      dimension is defined, the spatial dimension of C{arg}.

      @param arg: object to be checked
      @type arg: C{numarray.NumArray}, L{escript.Data}, C{float}, C{int},
                 C{Symbol}
      @return: True if C{arg} is a suitable object to be used for substitution,
               False otherwise
      @rtype: C{bool}
      """
      if isinstance(arg,numarray.NumArray):
          return arg.shape==self.getShape()
      elif isinstance(arg,escript.Data):
          if self.getDim()==None:
              return arg.getShape()==self.getShape()
          elif self.getDim()==arg.getFunctionSpace().getDim():
              return arg.getShape()==self.getShape()
          else:
              return False
      elif isinstance(arg,Symbol):
          if self.getDim()==None:
              return arg.getShape()==self.getShape()
          elif self.getDim()==arg.getDim():
              return arg.getShape()==self.getShape()
          else:
              return False
      elif isinstance(arg,float):
          return ()==self.getShape()
      elif isinstance(arg,int):
          return ()==self.getShape()
      else:
         return False

   def getMyCode(self,argstrs,format="escript"):
       """
       Returns program code that can be used to evaluate the symbol.

       @param argstrs: a string for each argument representing the argument
                       for the evaluation
       @type argstrs: C{list} of C{str}
       @param format: specifies the format to be used. At the moment only
                      "escript", "str" and "text" are supported.
       @type format: C{str}
       @return: a piece of program code which can be used to evaluate the
                expression assuming the values for the arguments are available
       @rtype: C{str}
       @raise NotImplementedError: if no implementation for the given format
                                   is available
       @note: This method has to be overwritten by subclasses.
       """
       raise NotImplementedError,"no code for %s representation available"%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.

      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @note: this method has to be overwritten by a particular L{Symbol}
      @raise NotImplementedError: if no implementation for the given format is
                                  available
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"Symbol: new value is not appropriate."
      else:
         raise NotImplementedError,"no substitution in %s avialable"%str(self)

   def diff(self,arg):
       """
       Returns the derivative of the symbol with respect to L{Symbol} C{arg}.

       @param arg: the derivative is calculated with respect to C{arg}
       @type arg: typically L{escript.Symbol} but can also be C{float},
                  L{escript.Data}, C{numarray.NumArray} depending on the
                  involved functions and data
       @return: derivative with respect to C{arg}
       @rtype: typically L{escript.Symbol} but other types such as C{float},
               L{escript.Data}, C{numarray.NumArray} are possible
       @note: this method is overwritten by a particular L{Symbol}.
       """
       if arg==self:
          return identity(self.getShape())
       else:
          s=self.getShape()+arg.getShape()
          if len(s)>0:
             return numarray.zeros(s,numarray.Float64)
          else:
             return 0.

   def __neg__(self):
       """
       Returns -self.

       @return: a L{Symbol} representing the negative of the object
       @rtype: L{DependendSymbol}
       """
       return self*(-1.)

   def __pos__(self):
       """
       Returns +self.

       @return: a L{Symbol} representing the positive of the object
       @rtype: L{DependendSymbol}
       """
       return self*(1.)

   def __abs__(self):
       """
       Returns a L{Symbol} representing the absolute value of the object.
       """
       return Abs_Symbol(self)

   def __add__(self,other):
       """
       Adds another object to this object.

       @param other: object to be added to this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}.
       @return: a L{Symbol} representing the sum of this object and C{other}
       @rtype: L{DependendSymbol}
       """
       return add(self,other)

   def __radd__(self,other):
       """
       Adds this object to another object.

       @param other: object to add this object to
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the sum of C{other} and this object
       @rtype: L{DependendSymbol}
       """
       return add(other,self)

   def __sub__(self,other):
       """
       Subtracts another object from this object.

       @param other: object to be subtracted from this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the difference of C{other} and this
                object
       @rtype: L{DependendSymbol}
       """
       return add(self,-other)

   def __rsub__(self,other):
       """
       Subtracts this object from another object.

       @param other: object this object is to be subtracted from
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the difference of this object and
                C{other}.
       @rtype: L{DependendSymbol}
       """
       return add(-self,other)

   def __mul__(self,other):
       """
       Multiplies this object with another object.

       @param other: object to be mutiplied by this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the product of the object and C{other}
       @rtype: L{DependendSymbol} or 0 if other is identical to zero.
       """
       return mult(self,other)

   def __rmul__(self,other):
       """
       Multiplies another object by this object.

       @param other: object this object is multiplied with
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the product of C{other} and the object
       @rtype: L{DependendSymbol} or 0 if other is identical to zero
       """
       return mult(other,self)

   def __div__(self,other):
       """
       Divides this object by another object.

       @param other: object dividing this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the quotient of this object and
                C{other}
       @rtype: L{DependendSymbol}
       """
       return quotient(self,other)

   def __rdiv__(self,other):
       """
       Divides another object by this object.

       @param other: object to be divided by this object
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the quotient of C{other} and this
                object
       @rtype: L{DependendSymbol} or 0 if C{other} is identical to zero
       """
       return quotient(other,self)

   def __pow__(self,other):
       """
       Raises this object to the power of C{other}.

       @param other: exponent
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the power of this object to C{other}
       @rtype: L{DependendSymbol} or 1 if C{other} is identical to zero
       """
       return power(self,other)

   def __rpow__(self,other):
       """
       Raises an object to the power of this object.

       @param other: basis
       @type other: L{escript.Symbol}, C{float}, L{escript.Data},
                    C{numarray.NumArray}
       @return: a L{Symbol} representing the power of C{other} to this object
       @rtype: L{DependendSymbol} or 0 if C{other} is identical to zero
       """
       return power(other,self)

   def __getitem__(self,index):
       """
       Returns the slice defined by C{index}.

       @param index: the slice index
       @type index: C{slice} or C{int} or a C{tuple} of them
       @return: a L{Symbol} representing the slice defined by index
       @rtype: L{DependendSymbol}
       """
       return GetSlice_Symbol(self,index)

class DependendSymbol(Symbol):
   """
   DependendSymbol extents L{Symbol} by modifying the == operator to allow two
   instances to be equal. Two C{DependendSymbol}s are equal if they have the
   same shape, the same arguments and one of them has an unspecified spatial
   dimension or the spatial dimension is identical.

   Example::

     u1=Symbol(shape=(3,4),dim=2,args=[4.])
     u2=Symbol(shape=(3,4),dim=2,args=[4.])
     print u1==u2
     False

   but::

     u1=DependendSymbol(shape=(3,4),dim=2,args=[4.])
     u2=DependendSymbol(shape=(3,4),dim=2,args=[4.])
     u3=DependendSymbol(shape=(2,),dim=2,args=[4.])
     print u1==u2, u1==u3
     True False

   @note: DependendSymbol should be used as return value of functions with
          L{Symbol} arguments. This will allow the optimizer to remove
          redundant function calls.
   """
   def __eq__(self,other):
      """
      Checks if C{other} equals self.

      @param other: any object
      @return: True if other has the same class as self and the shape, the
               spatial dimension and the arguments are equal, False otherwise
      @rtype: C{bool}
      """
      if isinstance(other,DependendSymbol):
         if self.__class__==other.__class__:
            if self.getShape()==other.getShape():
               if self.getArgument()==other.getArgument():
                  if self.getDim()==None or other.getDim()==None or self.getDim()==other.getDim():
                     return True
      return False

   def __ne__(self,other):
      """
      Checks if C{other} is not equal to self.

      @param other: any object
      @return: False if other has the same class as self and the shape, the
               spatial dimension and the arguments are equal, True otherwise
      @rtype: C{bool}
      """
      return not self==other
#=========================================================
#  Unary operations preserving the shape
#========================================================
class GetSlice_Symbol(DependendSymbol):
   """
   L{Symbol} representing getting a slice for a L{Symbol}.
   """
   def __init__(self,arg,index):
      """
      Initialization of the L{Symbol} with argument C{arg}.

      @param arg: argument
      @type arg: L{Symbol}
      @param index: defines index
      @type index: C{slice} or C{int} or a C{tuple} of them
      @raises IndexError: if length of index is larger than rank of arg or
                          index start or stop is out of range
      @raises ValueError: if a step is given
      """
      if not isinstance(index,tuple): index=(index,)
      if len(index)>arg.getRank():
           raise IndexError,"GetSlice_Symbol: index out of range."
      sh=()
      index2=()
      for i in range(len(index)):
         ix=index[i]
         if isinstance(ix,int):
            if ix<0 or ix>=arg.getShape()[i]:
               raise IndexError,"GetSlice_Symbol: index out of range."
            index2=index2+(ix,)
         else:
           if not ix.step==None:
             raise ValueError,"GetSlice_Symbol: stepping is not supported."
           if ix.start==None:
              s=0
           else:
              s=ix.start
           if ix.stop==None:
              e=arg.getShape()[i]
           else:
              e=ix.stop
              if e>arg.getShape()[i]:
                 raise IndexError,"GetSlice_Symbol: index out of range."
           index2=index2+(slice(s,e),)
           if e>s:
               sh=sh+(e-s,)
           elif s>e:
               raise IndexError,"GetSlice_Symbol: slice start must be less or equal slice end"
      for i in range(len(index),arg.getRank()):
          index2=index2+(slice(0,arg.getShape()[i]),)
          sh=sh+(arg.getShape()[i],)
      super(GetSlice_Symbol, self).__init__(args=[arg,index2],shape=sh,dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "%s.__getitem__(%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"GetItem_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         arg=args[0]
         index=args[1]
         return arg.__getitem__(index)

def log10(arg):
   """
   Returns base-10 logarithm of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.log10(arg)
   elif isinstance(arg,escript.Data):
      return arg._log10()
   elif isinstance(arg,float):
      return math.log10(arg)
   elif isinstance(arg,int):
      return math.log10(float(arg))
   elif isinstance(arg,Symbol):
      return log(arg)/log(10.)
   else:
      raise TypeError,"log10: Unknown argument type."

def wherePositive(arg):
   """
   Returns mask of positive values of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}.
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.greater(arg,numarray.zeros(arg.shape,numarray.Float64))*1.
      if isinstance(out,float): out=numarray.array(out,type=numarray.Float64)
      return out
   elif isinstance(arg,escript.Data):
      return arg._wherePositive()
   elif isinstance(arg,float):
      if arg>0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if arg>0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return WherePositive_Symbol(arg)
   else:
      raise TypeError,"wherePositive: Unknown argument type."

class WherePositive_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the mask of positive values function.
   """
   def __init__(self,arg):
      """
      Initialization of wherePositive L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "wherePositive(%s)"%argstrs
      else:
         raise NotImplementedError,"WherePositive_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return wherePositive(arg)

def whereNegative(arg):
   """
   Returns mask of negative values of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.less(arg,numarray.zeros(arg.shape,numarray.Float64))*1.
      if isinstance(out,float): out=numarray.array(out,type=numarray.Float64)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNegative()
   elif isinstance(arg,float):
      if arg<0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if arg<0:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return WhereNegative_Symbol(arg)
   else:
      raise TypeError,"whereNegative: Unknown argument type."

class WhereNegative_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the mask of negative values function.
   """
   def __init__(self,arg):
      """
      Initialization of whereNegative L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "whereNegative(%s)"%argstrs
      else:
         raise NotImplementedError,"WhereNegative_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return whereNegative(arg)

def whereNonNegative(arg):
   """
   Returns mask of non-negative values of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.greater_equal(arg,numarray.zeros(arg.shape,numarray.Float64))*1.
      if isinstance(out,float): out=numarray.array(out,type=numarray.Float64)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNonNegative()
   elif isinstance(arg,float):
      if arg<0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,int):
      if arg<0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,Symbol):
      return 1.-whereNegative(arg)
   else:
      raise TypeError,"whereNonNegative: Unknown argument type."

def whereNonPositive(arg):
   """
   Returns mask of non-positive values of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.less_equal(arg,numarray.zeros(arg.shape,numarray.Float64))*1.
      if isinstance(out,float): out=numarray.array(out,type=numarray.Float64)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNonPositive()
   elif isinstance(arg,float):
      if arg>0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,int):
      if arg>0:
        return 0.
      else:
        return 1.
   elif isinstance(arg,Symbol):
      return 1.-wherePositive(arg)
   else:
      raise TypeError,"whereNonPositive: Unknown argument type."

def whereZero(arg,tol=0.):
   """
   Returns mask of zero entries of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @param tol: tolerance. Values with absolute value less than tol are accepted
               as zero.
   @type tol: C{float}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.less_equal(abs(arg)-tol,numarray.zeros(arg.shape,numarray.Float64))*1.
      if isinstance(out,float): out=numarray.array(out,type=numarray.Float64)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereZero(tol)
   elif isinstance(arg,float):
      if abs(arg)<=tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if abs(float(arg))<=tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return WhereZero_Symbol(arg,tol)
   else:
      raise TypeError,"whereZero: Unknown argument type."

class WhereZero_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the mask of zero entries function.
   """
   def __init__(self,arg,tol=0.):
      """
      Initialization of whereZero L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg,tol],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "whereZero(%s,tol=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"WhereZero_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return whereZero(arg[0],arg[1])

def whereNonZero(arg,tol=0.):
   """
   Returns mask of values different from zero of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}.
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      out=numarray.greater(abs(arg)-tol,numarray.zeros(arg.shape,numarray.Float64))*1.
      if isinstance(out,float): out=numarray.array(out,type=numarray.Float64)
      return out
   elif isinstance(arg,escript.Data):
      return arg._whereNonZero(tol)
   elif isinstance(arg,float):
      if abs(arg)>tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,int):
      if abs(float(arg))>tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return 1.-whereZero(arg,tol)
   else:
      raise TypeError,"whereNonZero: Unknown argument type."

def erf(arg):
   """
   Returns the error function M{erf} of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}.
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,escript.Data):
      return arg._erf()
   else:
      raise TypeError,"erf: Unknown argument type."

def sin(arg):
   """
   Returns sine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}.
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.sin(arg)
   elif isinstance(arg,escript.Data):
      return arg._sin()
   elif isinstance(arg,float):
      return math.sin(arg)
   elif isinstance(arg,int):
      return math.sin(arg)
   elif isinstance(arg,Symbol):
      return Sin_Symbol(arg)
   else:
      raise TypeError,"sin: Unknown argument type."

class Sin_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the sine function.
   """
   def __init__(self,arg):
      """
      Initialization of sin L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "sin(%s)"%argstrs
      else:
         raise NotImplementedError,"Sin_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return sin(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(cos(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def cos(arg):
   """
   Returns cosine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.cos(arg)
   elif isinstance(arg,escript.Data):
      return arg._cos()
   elif isinstance(arg,float):
      return math.cos(arg)
   elif isinstance(arg,int):
      return math.cos(arg)
   elif isinstance(arg,Symbol):
      return Cos_Symbol(arg)
   else:
      raise TypeError,"cos: Unknown argument type."

class Cos_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the cosine function.
   """
   def __init__(self,arg):
      """
      Initialization of cos L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "cos(%s)"%argstrs
      else:
         raise NotImplementedError,"Cos_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return cos(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(-sin(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def tan(arg):
   """
   Returns tangent of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.tan(arg)
   elif isinstance(arg,escript.Data):
      return arg._tan()
   elif isinstance(arg,float):
      return math.tan(arg)
   elif isinstance(arg,int):
      return math.tan(arg)
   elif isinstance(arg,Symbol):
      return Tan_Symbol(arg)
   else:
      raise TypeError,"tan: Unknown argument type."

class Tan_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the tangent function.
   """
   def __init__(self,arg):
      """
      Initialization of tan L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "tan(%s)"%argstrs
      else:
         raise NotImplementedError,"Tan_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return tan(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./cos(myarg)**2,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def asin(arg):
   """
   Returns the inverse sine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arcsin(arg)
   elif isinstance(arg,escript.Data):
      return arg._asin()
   elif isinstance(arg,float):
      return math.asin(arg)
   elif isinstance(arg,int):
      return math.asin(arg)
   elif isinstance(arg,Symbol):
      return Asin_Symbol(arg)
   else:
      raise TypeError,"asin: Unknown argument type."

class Asin_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse sine function.
   """
   def __init__(self,arg):
      """
      Initialization of asin L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "asin(%s)"%argstrs
      else:
         raise NotImplementedError,"Asin_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return asin(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./sqrt(1.-myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def acos(arg):
   """
   Returns the inverse cosine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arccos(arg)
   elif isinstance(arg,escript.Data):
      return arg._acos()
   elif isinstance(arg,float):
      return math.acos(arg)
   elif isinstance(arg,int):
      return math.acos(arg)
   elif isinstance(arg,Symbol):
      return Acos_Symbol(arg)
   else:
      raise TypeError,"acos: Unknown argument type."

class Acos_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse cosine function.
   """
   def __init__(self,arg):
      """
      Initialization of acos L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "acos(%s)"%argstrs
      else:
         raise NotImplementedError,"Acos_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return acos(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(-1./sqrt(1.-myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def atan(arg):
   """
   Returns inverse tangent of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arctan(arg)
   elif isinstance(arg,escript.Data):
      return arg._atan()
   elif isinstance(arg,float):
      return math.atan(arg)
   elif isinstance(arg,int):
      return math.atan(arg)
   elif isinstance(arg,Symbol):
      return Atan_Symbol(arg)
   else:
      raise TypeError,"atan: Unknown argument type."

class Atan_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse tangent function.
   """
   def __init__(self,arg):
      """
      Initialization of atan L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "atan(%s)"%argstrs
      else:
         raise NotImplementedError,"Atan_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return atan(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./(1+myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def sinh(arg):
   """
   Returns the hyperbolic sine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.sinh(arg)
   elif isinstance(arg,escript.Data):
      return arg._sinh()
   elif isinstance(arg,float):
      return math.sinh(arg)
   elif isinstance(arg,int):
      return math.sinh(arg)
   elif isinstance(arg,Symbol):
      return Sinh_Symbol(arg)
   else:
      raise TypeError,"sinh: Unknown argument type."

class Sinh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the hyperbolic sine function.
   """
   def __init__(self,arg):
      """
      Initialization of sinh L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "sinh(%s)"%argstrs
      else:
         raise NotImplementedError,"Sinh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return sinh(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(cosh(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def cosh(arg):
   """
   Returns the hyperbolic cosine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.cosh(arg)
   elif isinstance(arg,escript.Data):
      return arg._cosh()
   elif isinstance(arg,float):
      return math.cosh(arg)
   elif isinstance(arg,int):
      return math.cosh(arg)
   elif isinstance(arg,Symbol):
      return Cosh_Symbol(arg)
   else:
      raise TypeError,"cosh: Unknown argument type."

class Cosh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the hyperbolic cosine function.
   """
   def __init__(self,arg):
      """
      Initialization of cosh L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "cosh(%s)"%argstrs
      else:
         raise NotImplementedError,"Cosh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return cosh(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(sinh(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def tanh(arg):
   """
   Returns the hyperbolic tangent of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.tanh(arg)
   elif isinstance(arg,escript.Data):
      return arg._tanh()
   elif isinstance(arg,float):
      return math.tanh(arg)
   elif isinstance(arg,int):
      return math.tanh(arg)
   elif isinstance(arg,Symbol):
      return Tanh_Symbol(arg)
   else:
      raise TypeError,"tanh: Unknown argument type."

class Tanh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the hyperbolic tangent function.
   """
   def __init__(self,arg):
      """
      Initialization of tanh L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "tanh(%s)"%argstrs
      else:
         raise NotImplementedError,"Tanh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return tanh(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./cosh(myarg)**2,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def asinh(arg):
   """
   Returns the inverse hyperbolic sine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arcsinh(arg)
   elif isinstance(arg,escript.Data):
      return arg._asinh()
   elif isinstance(arg,float):
      return numarray.arcsinh(arg)
   elif isinstance(arg,int):
      return numarray.arcsinh(float(arg))
   elif isinstance(arg,Symbol):
      return Asinh_Symbol(arg)
   else:
      raise TypeError,"asinh: Unknown argument type."

class Asinh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse hyperbolic sine function.
   """
   def __init__(self,arg):
      """
      Initialization of asinh L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "asinh(%s)"%argstrs
      else:
         raise NotImplementedError,"Asinh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return asinh(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./sqrt(myarg**2+1),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def acosh(arg):
   """
   Returns the inverse hyperbolic cosine of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arccosh(arg)
   elif isinstance(arg,escript.Data):
      return arg._acosh()
   elif isinstance(arg,float):
      return numarray.arccosh(arg)
   elif isinstance(arg,int):
      return numarray.arccosh(float(arg))
   elif isinstance(arg,Symbol):
      return Acosh_Symbol(arg)
   else:
      raise TypeError,"acosh: Unknown argument type."

class Acosh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse hyperbolic cosine function.
   """
   def __init__(self,arg):
      """
      Initialization of acosh L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "acosh(%s)"%argstrs
      else:
         raise NotImplementedError,"Acosh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return acosh(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./sqrt(myarg**2-1),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def atanh(arg):
   """
   Returns the inverse hyperbolic tangent of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.arctanh(arg)
   elif isinstance(arg,escript.Data):
      return arg._atanh()
   elif isinstance(arg,float):
      return numarray.arctanh(arg)
   elif isinstance(arg,int):
      return numarray.arctanh(float(arg))
   elif isinstance(arg,Symbol):
      return Atanh_Symbol(arg)
   else:
      raise TypeError,"atanh: Unknown argument type."

class Atanh_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse hyperbolic tangent function.
   """
   def __init__(self,arg):
      """
      Initialization of atanh L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "atanh(%s)"%argstrs
      else:
         raise NotImplementedError,"Atanh_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return atanh(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./(1.-myarg**2),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def exp(arg):
   """
   Returns M{e} to the power of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}.
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of arg
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.exp(arg)
   elif isinstance(arg,escript.Data):
      return arg._exp()
   elif isinstance(arg,float):
      return math.exp(arg)
   elif isinstance(arg,int):
      return math.exp(arg)
   elif isinstance(arg,Symbol):
      return Exp_Symbol(arg)
   else:
      raise TypeError,"exp: Unknown argument type."

class Exp_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the exponential function.
   """
   def __init__(self,arg):
      """
      Initialization of exp L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "exp(%s)"%argstrs
      else:
         raise NotImplementedError,"Exp_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return exp(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(self,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def sqrt(arg):
   """
   Returns the square root of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
           depending on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.sqrt(arg)
   elif isinstance(arg,escript.Data):
      return arg._sqrt()
   elif isinstance(arg,float):
      return math.sqrt(arg)
   elif isinstance(arg,int):
      return math.sqrt(arg)
   elif isinstance(arg,Symbol):
      return Sqrt_Symbol(arg)
   else:
      raise TypeError,"sqrt: Unknown argument type."

class Sqrt_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the square root function.
   """
   def __init__(self,arg):
      """
      Initialization of sqrt L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "sqrt(%s)"%argstrs
      else:
         raise NotImplementedError,"Sqrt_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return sqrt(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(0.5/self,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def log(arg):
   """
   Returns the natural logarithm of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}.
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return numarray.log(arg)
   elif isinstance(arg,escript.Data):
      return arg._log()
   elif isinstance(arg,float):
      return math.log(arg)
   elif isinstance(arg,int):
      return math.log(arg)
   elif isinstance(arg,Symbol):
      return Log_Symbol(arg)
   else:
      raise TypeError,"log: Unknown argument type."

class Log_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the natural logarithm function.
   """
   def __init__(self,arg):
      """
      Initialization of log L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "log(%s)"%argstrs
      else:
         raise NotImplementedError,"Log_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return log(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(1./arg,self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def sign(arg):
   """
   Returns the sign of argument C{arg}.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending
           on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      return wherePositive(arg)-whereNegative(arg)
   elif isinstance(arg,escript.Data):
      return arg._sign()
   elif isinstance(arg,float):
      if arg>0:
        return 1.
      elif arg<0:
        return -1.
      else:
        return 0.
   elif isinstance(arg,int):
      if float(arg)>0:
        return 1.
      elif float(arg)<0:
        return -1.
      else:
        return 0.
   elif isinstance(arg,Symbol):
      return wherePositive(arg)-whereNegative(arg)
   else:
      raise TypeError,"sign: Unknown argument type."

class Abs_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the absolute value function.
   """
   def __init__(self,arg):
      """
      Initialization of abs L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "abs(%s)"%argstrs
      else:
         raise NotImplementedError,"Abs_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return abs(arg)

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myarg=self.getArgument()[0]
         val=matchShape(sign(myarg),self.getDifferentiatedArguments(arg)[0])
         return val[0]*val[1]

def minval(arg):
   """
   Returns the minimum value over all components of C{arg} at each data point.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol} depending on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      if arg.rank==0:
         return float(arg)
      else:
         return arg.min()
   elif isinstance(arg,escript.Data):
      return arg._minval()
   elif isinstance(arg,float):
      return arg
   elif isinstance(arg,int):
      return float(arg)
   elif isinstance(arg,Symbol):
      return Minval_Symbol(arg)
   else:
      raise TypeError,"minval: Unknown argument type."

class Minval_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the minimum value function.
   """
   def __init__(self,arg):
      """
      Initialization of minimum value L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "minval(%s)"%argstrs
      else:
         raise NotImplementedError,"Minval_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return minval(arg)

def maxval(arg):
   """
   Returns the maximum value over all components of C{arg} at each data point.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol} depending on the type of C{arg}
   @raises TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numarray.NumArray):
      if arg.rank==0:
         return float(arg)
      else:
         return arg.max()
   elif isinstance(arg,escript.Data):
      return arg._maxval()
   elif isinstance(arg,float):
      return arg
   elif isinstance(arg,int):
      return float(arg)
   elif isinstance(arg,Symbol):
      return Maxval_Symbol(arg)
   else:
      raise TypeError,"maxval: Unknown argument type."

class Maxval_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the maximum value function.
   """
   def __init__(self,arg):
      """
      Initialization of maximum value L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: typically L{Symbol}
      """
      DependendSymbol.__init__(self,args=[arg],shape=(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if isinstance(argstrs,list):
          argstrs=argstrs[0]
      if format=="escript" or format=="str"  or format=="text":
         return "maxval(%s)"%argstrs
      else:
         raise NotImplementedError,"Maxval_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)[0]
         return maxval(arg)

def length(arg):
   """
   Returns the length (Euclidean norm) of argument C{arg} at each data point.

   @param arg: argument
   @type arg: C{float}, L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @rtype: C{float}, L{escript.Data}, L{Symbol} depending on the type of C{arg}
   """
   return sqrt(inner(arg,arg))

def trace(arg,axis_offset=0):
   """
   Returns the trace of C{arg} which is the sum of C{arg[k,k]} over k.

   @param arg: argument
   @type arg: L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @param axis_offset: C{axis_offset} to components to sum over. C{axis_offset}
                       must be non-negative and less than the rank of C{arg} +1.
                       The dimensions of component C{axis_offset} and
                       axis_offset+1 must be equal.
   @type axis_offset: C{int}
   @return: trace of arg. The rank of the returned object is rank of C{arg}
            minus 2.
   @rtype: L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending on the
           type of C{arg}
   """
   if isinstance(arg,numarray.NumArray):
      sh=arg.shape
      if len(sh)<2:
        raise ValueError,"rank of argument must be greater than 1"
      if axis_offset<0 or axis_offset>len(sh)-2:
        raise ValueError,"axis_offset must be between 0 and %d"%(len(sh)-2)
      s1=1
      for i in range(axis_offset): s1*=sh[i]
      s2=1
      for i in range(axis_offset+2,len(sh)): s2*=sh[i]
      if not sh[axis_offset] == sh[axis_offset+1]:
        raise ValueError,"dimensions of component %d and %d must match."%(axis_offset,axis_offset+1)
      arg_reshaped=numarray.reshape(arg,(s1,sh[axis_offset],sh[axis_offset],s2))
      out=numarray.zeros([s1,s2],numarray.Float64)
      for i1 in range(s1):
        for i2 in range(s2):
            for j in range(sh[axis_offset]): out[i1,i2]+=arg_reshaped[i1,j,j,i2]
      out.resize(sh[:axis_offset]+sh[axis_offset+2:])
      return out
   elif isinstance(arg,escript.Data):
      if arg.getRank()<2:
        raise ValueError,"rank of argument must be greater than 1"
      if axis_offset<0 or axis_offset>arg.getRank()-2:
        raise ValueError,"axis_offset must be between 0 and %d"%(arg.getRank()-2)
      s=list(arg.getShape())
      if not s[axis_offset] == s[axis_offset+1]:
        raise ValueError,"dimensions of component %d and %d must match."%(axis_offset,axis_offset+1)
      return arg._trace(axis_offset)
   elif isinstance(arg,float):
      raise TypeError,"illegal argument type float."
   elif isinstance(arg,int):
      raise TypeError,"illegal argument type int."
   elif isinstance(arg,Symbol):
      return Trace_Symbol(arg,axis_offset)
   else:
      raise TypeError,"Unknown argument type."

class Trace_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the trace function.
   """
   def __init__(self,arg,axis_offset=0):
      """
      Initialization of trace L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: L{Symbol}
      @param axis_offset: C{axis_offset} to components to sum over.
                          C{axis_offset} must be non-negative and less than the
                          rank of C{arg} +1. The dimensions on component
                          C{axis_offset} and axis_offset+1 must be equal.
      @type axis_offset: C{int}
      """
      if arg.getRank()<2:
        raise ValueError,"rank of argument must be greater than 1"
      if axis_offset<0 or axis_offset>arg.getRank()-2:
        raise ValueError,"axis_offset must be between 0 and %d"%(arg.getRank()-2)
      s=list(arg.getShape())
      if not s[axis_offset] == s[axis_offset+1]:
        raise ValueError,"dimensions of component %d and %d must match."%(axis_offset,axis_offset+1)
      super(Trace_Symbol,self).__init__(args=[arg,axis_offset],shape=tuple(s[0:axis_offset]+s[axis_offset+2:]),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "trace(%s,axis_offset=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"Trace_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return trace(arg[0],axis_offset=arg[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return trace(self.getDifferentiatedArguments(arg)[0],axis_offset=self.getArgument()[1])

def transpose(arg,axis_offset=None):
   """
   Returns the transpose of C{arg} by swapping the first C{axis_offset} and the
   last C{rank-axis_offset} components.

   @param arg: argument
   @type arg: L{escript.Data}, L{Symbol}, C{numarray.NumArray}, C{float}, C{int}
   @param axis_offset: the first C{axis_offset} components are swapped with the
                       rest. C{axis_offset} must be non-negative and less or
                       equal to the rank of C{arg}. If C{axis_offset} is not
                       present C{int(r/2)} where r is the rank of C{arg} is
                       used.
   @type axis_offset: C{int}
   @return: transpose of C{arg}
   @rtype: L{escript.Data}, L{Symbol}, C{numarray.NumArray}, C{float}, C{int}
           depending on the type of C{arg}
   """
   if isinstance(arg,numarray.NumArray):
      if axis_offset==None: axis_offset=int(arg.rank/2)
      return numarray.transpose(arg,axes=range(axis_offset,arg.rank)+range(0,axis_offset))
   elif isinstance(arg,escript.Data):
      r=arg.getRank()
      if axis_offset==None: axis_offset=int(r/2)
      if axis_offset<0 or axis_offset>r:
        raise ValueError,"axis_offset must be between 0 and %s"%r
      return arg._transpose(axis_offset)
   elif isinstance(arg,float):
      if not ( axis_offset==0 or axis_offset==None):
        raise ValueError,"axis_offset must be 0 for float argument"
      return arg
   elif isinstance(arg,int):
      if not ( axis_offset==0 or axis_offset==None):
        raise ValueError,"axis_offset must be 0 for int argument"
      return float(arg)
   elif isinstance(arg,Symbol):
      if axis_offset==None: axis_offset=int(arg.getRank()/2)
      return Transpose_Symbol(arg,axis_offset)
   else:
      raise TypeError,"Unknown argument type."

class Transpose_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the transpose function.
   """
   def __init__(self,arg,axis_offset=None):
      """
      Initialization of transpose L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: L{Symbol}
      @param axis_offset: the first C{axis_offset} components are swapped with
                          the rest. C{axis_offset} must be non-negative and
                          less than or equal to the rank of C{arg}. If
                          C{axis_offset} is not present C{int(r/2)} where r is
                          the rank of C{arg} is used.
      @type axis_offset: C{int}
      """
      if axis_offset==None: axis_offset=int(arg.getRank()/2)
      if axis_offset<0 or axis_offset>arg.getRank():
        raise ValueError,"axis_offset must be between 0 and %s"%arg.getRank()
      s=arg.getShape()
      super(Transpose_Symbol,self).__init__(args=[arg,axis_offset],shape=s[axis_offset:]+s[:axis_offset],dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "transpose(%s,axis_offset=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"Transpose_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return transpose(arg[0],axis_offset=arg[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return transpose(self.getDifferentiatedArguments(arg)[0],axis_offset=self.getArgument()[1])

def swap_axes(arg,axis0=0,axis1=1):
   """
   Returns the swap of C{arg} by swapping the components C{axis0} and C{axis1}.

   @param arg: argument
   @type arg: L{escript.Data}, L{Symbol}, C{numarray.NumArray}
   @param axis0: first axis. C{axis0} must be non-negative and less than the
                 rank of C{arg}.
   @type axis0: C{int}
   @param axis1: second axis. C{axis1} must be non-negative and less than the
                 rank of C{arg}.
   @type axis1: C{int}
   @return: C{arg} with swapped components
   @rtype: L{escript.Data}, L{Symbol}, C{numarray.NumArray} depending on the
           type of C{arg}
   """
   if axis0 > axis1:
      axis0,axis1=axis1,axis0
   if isinstance(arg,numarray.NumArray):
      return numarray.swapaxes(arg,axis0,axis1)
   elif isinstance(arg,escript.Data):
      return arg._swap_axes(axis0,axis1)
   elif isinstance(arg,float):
      raise TypeError,"float argument is not supported."
   elif isinstance(arg,int):
      raise TypeError,"int argument is not supported."
   elif isinstance(arg,Symbol):
      return SwapAxes_Symbol(arg,axis0,axis1)
   else:
      raise TypeError,"Unknown argument type."

class SwapAxes_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the swap function.
   """
   def __init__(self,arg,axis0=0,axis1=1):
      """
      Initialization of swap L{Symbol} with argument C{arg}.

      @param arg: argument
      @type arg: L{Symbol}
      @param axis0: first axis. C{axis0} must be non-negative and less than the
                    rank of C{arg}.
      @type axis0: C{int}
      @param axis1: second axis. C{axis1} must be non-negative and less than
                    the rank of C{arg}.
      @type axis1: C{int}
      """
      if arg.getRank()<2:
         raise ValueError,"argument must have at least rank 2."
      if axis0<0 or axis0>arg.getRank()-1:
         raise ValueError,"axis0 must be between 0 and %s"%arg.getRank()-1
      if axis1<0 or axis1>arg.getRank()-1:
         raise ValueError,"axis1 must be between 0 and %s"%arg.getRank()-1
      if axis0 == axis1:
         raise ValueError,"axis indices must be different."
      if axis0 > axis1:
         axis0,axis1=axis1,axis0
      s=arg.getShape()
      s_out=[]
      for i in range(len(s)):
         if i == axis0:
            s_out.append(s[axis1])
         elif i == axis1:
            s_out.append(s[axis0])
         else:
            s_out.append(s[i])
      super(SwapAxes_Symbol,self).__init__(args=[arg,axis0,axis1],shape=tuple(s_out),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "swap(%s,axis_offset=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"SwapAxes_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return swap_axes(arg[0],axis0=arg[1],axis1=arg[2])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return swap_axes(self.getDifferentiatedArguments(arg)[0],axis0=self.getArgument()[1],axis1=self.getArgument()[2])

def symmetric(arg):
    """
    Returns the symmetric part of the square matrix C{arg}. That is,
    M{(arg+transpose(arg))/2}.

    @param arg: input matrix. Must have rank 2 or 4 and be square.
    @type arg: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: symmetric part of C{arg}
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    if isinstance(arg,numarray.NumArray):
      if arg.rank==2:
        if not (arg.shape[0]==arg.shape[1]):
           raise ValueError,"argument must be square."
      elif arg.rank==4:
        if not (arg.shape[0]==arg.shape[2] and arg.shape[1]==arg.shape[3]):
           raise ValueError,"argument must be square."
      else:
        raise ValueError,"rank 2 or 4 is required."
      return (arg+transpose(arg))/2
    elif isinstance(arg,escript.Data):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError,"argument must be square."
        return arg._symmetric()
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError,"argument must be square."
        return arg._symmetric()
      else:
        raise ValueError,"rank 2 or 4 is required."
    elif isinstance(arg,float):
      return arg
    elif isinstance(arg,int):
      return float(arg)
    elif isinstance(arg,Symbol):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError,"argument must be square."
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError,"argument must be square."
      else:
        raise ValueError,"rank 2 or 4 is required."
      return (arg+transpose(arg))/2
    else:
      raise TypeError,"symmetric: Unknown argument type."

def nonsymmetric(arg):
    """
    Returns the non-symmetric part of the square matrix C{arg}. That is,
    M{(arg-transpose(arg))/2}.

    @param arg: input matrix. Must have rank 2 or 4 and be square.
    @type arg: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: non-symmetric part of C{arg}
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    if isinstance(arg,numarray.NumArray):
      if arg.rank==2:
        if not (arg.shape[0]==arg.shape[1]):
           raise ValueError,"nonsymmetric: argument must be square."
      elif arg.rank==4:
        if not (arg.shape[0]==arg.shape[2] and arg.shape[1]==arg.shape[3]):
           raise ValueError,"nonsymmetric: argument must be square."
      else:
        raise ValueError,"nonsymmetric: rank 2 or 4 is required."
      return (arg-transpose(arg))/2
    elif isinstance(arg,escript.Data):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError,"argument must be square."
        return arg._nonsymmetric()
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError,"argument must be square."
        return arg._nonsymmetric()
      else:
        raise ValueError,"rank 2 or 4 is required."
    elif isinstance(arg,float):
      return arg
    elif isinstance(arg,int):
      return float(arg)
    elif isinstance(arg,Symbol):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError,"nonsymmetric: argument must be square."
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError,"nonsymmetric: argument must be square."
      else:
        raise ValueError,"nonsymmetric: rank 2 or 4 is required."
      return (arg-transpose(arg))/2
    else:
      raise TypeError,"nonsymmetric: Unknown argument type."

def inverse(arg):
    """
    Returns the inverse of the square matrix C{arg}.

    @param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal.
    @type arg: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: inverse of the argument. C{matrix_mult(inverse(arg),arg)} will be
             almost equal to C{kronecker(arg.getShape()[0])}
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    @note: for L{escript.Data} objects the dimension is restricted to 3.
    """
    import numarray.linear_algebra # This statement should be after the next statement but then somehow numarray is gone.
    if isinstance(arg,numarray.NumArray):
      return numarray.linear_algebra.inverse(arg)
    elif isinstance(arg,escript.Data):
      return escript_inverse(arg)
    elif isinstance(arg,float):
      return 1./arg
    elif isinstance(arg,int):
      return 1./float(arg)
    elif isinstance(arg,Symbol):
      return Inverse_Symbol(arg)
    else:
      raise TypeError,"inverse: Unknown argument type."

def escript_inverse(arg): # this should be escript._inverse and use LAPACK
      "arg is a Data object!!!"
      if not arg.getRank()==2:
        raise ValueError,"escript_inverse: argument must have rank 2"
      s=arg.getShape()
      if not s[0] == s[1]:
        raise ValueError,"escript_inverse: argument must be a square matrix."
      out=escript.Data(0.,s,arg.getFunctionSpace())
      if s[0]==1:
          if inf(abs(arg[0,0]))==0: # in c this should be done point wise as abs(arg[0,0](i))<=0.
              raise ZeroDivisionError,"escript_inverse: argument not invertible"
          out[0,0]=1./arg[0,0]
      elif s[0]==2:
          A11=arg[0,0]
          A12=arg[0,1]
          A21=arg[1,0]
          A22=arg[1,1]
          D = A11*A22-A12*A21
          if inf(abs(D))==0: # in c this should be done point wise as abs(D(i))<=0.
              raise ZeroDivisionError,"escript_inverse: argument not invertible"
          D=1./D
          out[0,0]= A22*D
          out[1,0]=-A21*D
          out[0,1]=-A12*D
          out[1,1]= A11*D
      elif s[0]==3:
          A11=arg[0,0]
          A21=arg[1,0]
          A31=arg[2,0]
          A12=arg[0,1]
          A22=arg[1,1]
          A32=arg[2,1]
          A13=arg[0,2]
          A23=arg[1,2]
          A33=arg[2,2]
          D = A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22)
          if inf(abs(D))==0: # in c this should be done point wise as abs(D(i))<=0.
              raise ZeroDivisionError,"escript_inverse: argument not invertible"
          D=1./D
          out[0,0]=(A22*A33-A23*A32)*D
          out[1,0]=(A31*A23-A21*A33)*D
          out[2,0]=(A21*A32-A31*A22)*D
          out[0,1]=(A13*A32-A12*A33)*D
          out[1,1]=(A11*A33-A31*A13)*D
          out[2,1]=(A12*A31-A11*A32)*D
          out[0,2]=(A12*A23-A13*A22)*D
          out[1,2]=(A13*A21-A11*A23)*D
          out[2,2]=(A11*A22-A12*A21)*D
      else:
         raise TypeError,"escript_inverse: only matrix dimensions 1,2,3 are supported right now."
      return out

class Inverse_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the inverse function.
   """
   def __init__(self,arg):
      """
      Initialization of inverse L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: L{Symbol}
      """
      if not arg.getRank()==2:
        raise ValueError,"Inverse_Symbol:: argument must have rank 2"
      s=arg.getShape()
      if not s[0] == s[1]:
        raise ValueError,"Inverse_Symbol:: argument must be a square matrix."
      super(Inverse_Symbol,self).__init__(args=[arg],shape=s,dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 1 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "inverse(%s)"%argstrs[0]
      else:
         raise NotImplementedError,"Inverse_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return inverse(arg[0])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return -matrix_mult(matrix_mult(self,self.getDifferentiatedArguments(arg)[0]),self)

def eigenvalues(arg):
    """
    Returns the eigenvalues of the square matrix C{arg}.

    @param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal. It must also be symmetric, ie.
                C{transpose(arg)==arg} (this is not checked).
    @type arg: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the eigenvalues in increasing order
    @rtype: C{numarray.NumArray},L{escript.Data}, L{Symbol} depending on the
            input
    @note: for L{escript.Data} and L{Symbol} objects the dimension is
           restricted to 3.
    """
    if isinstance(arg,numarray.NumArray):
      out=numarray.linear_algebra.eigenvalues((arg+numarray.transpose(arg))/2.)
      out.sort()
      return out
    elif isinstance(arg,escript.Data):
      return arg._eigenvalues()
    elif isinstance(arg,Symbol):
      if not arg.getRank()==2:
        raise ValueError,"eigenvalues: argument must have rank 2"
      s=arg.getShape()
      if not s[0] == s[1]:
        raise ValueError,"eigenvalues: argument must be a square matrix."
      if s[0]==1:
          return arg[0]
      elif s[0]==2:
          arg1=symmetric(arg)
          A11=arg1[0,0]
          A12=arg1[0,1]
          A22=arg1[1,1]
          trA=(A11+A22)/2.
          A11-=trA
          A22-=trA
          s=sqrt(A12**2-A11*A22)
          return trA+s*numarray.array([-1.,1.],type=numarray.Float64)
      elif s[0]==3:
          arg1=symmetric(arg)
          A11=arg1[0,0]
          A12=arg1[0,1]
          A22=arg1[1,1]
          A13=arg1[0,2]
          A23=arg1[1,2]
          A33=arg1[2,2]
          trA=(A11+A22+A33)/3.
          A11-=trA
          A22-=trA
          A33-=trA
          A13_2=A13**2
          A23_2=A23**2
          A12_2=A12**2
          p=A13_2+A23_2+A12_2+(A11**2+A22**2+A33**2)/2.
          q=A13_2*A22+A23_2*A11+A12_2*A33-A11*A22*A33-2*A12*A23*A13
          sq_p=sqrt(p/3.)
          alpha_3=acos(clip(-q*(sq_p+whereZero(p,0.)*1.e-15)**(-3.)/2.,-1.,1.))/3.  # whereZero is protection against divison by zero
          sq_p*=2.
          f=cos(alpha_3)               *numarray.array([0.,0.,1.],type=numarray.Float64) \
           -cos(alpha_3+numarray.pi/3.)*numarray.array([0.,1.,0.],type=numarray.Float64) \
           -cos(alpha_3-numarray.pi/3.)*numarray.array([1.,0.,0.],type=numarray.Float64)
          return trA+sq_p*f
      else:
         raise TypeError,"eigenvalues: only matrix dimensions 1,2,3 are supported right now."
    elif isinstance(arg,float):
      return arg
    elif isinstance(arg,int):
      return float(arg)
    else:
      raise TypeError,"eigenvalues: Unknown argument type."

def eigenvalues_and_eigenvectors(arg):
    """
    Returns the eigenvalues and eigenvectors of the square matrix C{arg}.

    @param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal. It must also be symmetric, ie.
                C{transpose(arg)==arg} (this is not checked).
    @type arg: L{escript.Data}
    @return: the eigenvalues and eigenvectors. The eigenvalues are ordered by
             increasing value. The eigenvectors are orthogonal and normalized.
             If V are the eigenvectors then V[:,i] is the eigenvector
             corresponding to the i-th eigenvalue.
    @rtype: L{tuple} of L{escript.Data}
    @note: The dimension is restricted to 3.
    """
    if isinstance(arg,numarray.NumArray):
      raise TypeError,"eigenvalues_and_eigenvectors does not support numarray arguments"
    elif isinstance(arg,escript.Data):
      return arg._eigenvalues_and_eigenvectors()
    elif isinstance(arg,Symbol):
      raise TypeError,"eigenvalues_and_eigenvectors does not support Symbol arguments"
    elif isinstance(arg,float):
      return (numarray.array([[arg]],numarray.Float),numarray.ones((1,1),numarray.Float))
    elif isinstance(arg,int):
      return (numarray.array([[arg]],numarray.Float),numarray.ones((1,1),numarray.Float))
    else:
      raise TypeError,"eigenvalues: Unknown argument type."

#=======================================================
#  Binary operations:
#=======================================================
def add(arg0,arg1):
       """
       Adds C{arg0} and C{arg1} together.

       @param arg0: first term
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second term
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @return: the sum of C{arg0} and C{arg1}
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
               C{numarray.NumArray}
       @note: The shape of both arguments is matched according to the rules
              used in L{matchShape}.
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]):
          return args[1]
       elif testForZero(args[1]):
          return args[0]
       else:
          if isinstance(args[0],Symbol) or isinstance(args[1],Symbol) :
              return Add_Symbol(args[0],args[1])
          elif isinstance(args[0],numarray.NumArray):
              return args[1]+args[0]
          else:
              return args[0]+args[1]

class Add_Symbol(DependendSymbol):
   """
   Symbol representing the sum of two arguments.
   """
   def __init__(self,arg0,arg1):
       """
       Initialization of the L{Symbol} representing the sum of two arguments.

       @param arg0: first term in the sum
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second term in the sum
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: if both arguments do not have the same shape
       @note: if both arguments have a spatial dimension, they must be equal.
       """
       sh0=getShape(arg0)
       sh1=getShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Add_Symbol: shape of arguments must match"
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="str" or format=="text":
         return "(%s)+(%s)"%(argstrs[0],argstrs[1])
      elif format=="escript":
         return "add(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         out=add(args[0],args[1])
         return out

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         dargs=self.getDifferentiatedArguments(arg)
         return add(dargs[0],dargs[1])

def mult(arg0,arg1):
       """
       Product of C{arg0} and C{arg1}.

       @param arg0: first term
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second term
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @return: the product of C{arg0} and C{arg1}
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
               C{numarray.NumArray}
       @note: The shape of both arguments is matched according to the rules
              used in L{matchShape}.
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]) or testForZero(args[1]):
          return numarray.zeros(getShape(args[0]),numarray.Float64)
       else:
          if isinstance(args[0],Symbol) or isinstance(args[1],Symbol) :
              return Mult_Symbol(args[0],args[1])
          elif isinstance(args[0],numarray.NumArray):
              return args[1]*args[0]
          else:
              return args[0]*args[1]

class Mult_Symbol(DependendSymbol):
   """
   Symbol representing the product of two arguments.
   """
   def __init__(self,arg0,arg1):
       """
       Initialization of the L{Symbol} representing the product of two
       arguments.

       @param arg0: first factor
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second factor
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: if both arguments do not have the same shape
       @note: if both arguments have a spatial dimension, they must be equal.
       """
       sh0=getShape(arg0)
       sh1=getShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Mult_Symbol: shape of arguments must match"
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="str" or format=="text":
         return "(%s)*(%s)"%(argstrs[0],argstrs[1])
      elif format=="escript":
         return "mult(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}.
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return mult(args[0],args[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myargs=self.getArgument()
         dargs=self.getDifferentiatedArguments(arg)
         return add(mult(myargs[0],dargs[1]),mult(myargs[1],dargs[0]))

def quotient(arg0,arg1):
       """
       Quotient of C{arg0} and C{arg1}.

       @param arg0: numerator
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: denominator
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @return: the quotient of C{arg0} and C{arg1}
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
               C{numarray.NumArray}
       @note: The shape of both arguments is matched according to the rules
              used in L{matchShape}.
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]):
          return numarray.zeros(getShape(args[0]),numarray.Float64)
       elif isinstance(args[0],Symbol):
          if isinstance(args[1],Symbol):
             return Quotient_Symbol(args[0],args[1])
          else:
             return mult(args[0],1./args[1])
       else:
          if isinstance(args[1],Symbol):
             return Quotient_Symbol(args[0],args[1])
          elif isinstance(args[0],numarray.NumArray) and not isinstance(args[1],numarray.NumArray):
             return 1./args[1]*args[0]
          else:
             return args[0]/args[1]

class Quotient_Symbol(DependendSymbol):
   """
   Symbol representing the quotient of two arguments.
   """
   def __init__(self,arg0,arg1):
       """
       Initialization of L{Symbol} representing the quotient of two arguments.

       @param arg0: numerator
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: denominator
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: if both arguments do not have the same shape
       @note: if both arguments have a spatial dimension, they must be equal.
       """
       sh0=getShape(arg0)
       sh1=getShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Quotient_Symbol: shape of arguments must match"
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="str" or format=="text":
         return "(%s)/(%s)"%(argstrs[0],argstrs[1])
      if format=="escript":
         return "quotient(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return quotient(args[0],args[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myargs=self.getArgument()
         dargs=self.getDifferentiatedArguments(arg)
         return quotient(add(mult(myargs[1],dargs[0]),mult(-myargs[0],dargs[1])),myargs[1]*myargs[1])


def power(arg0,arg1):
       """
       Raises C{arg0} to the power of C{arg1}.

       @param arg0: basis
       @type arg0: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: exponent
       @type arg1: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
                   C{numarray.NumArray}
       @return: C{arg0} to the power of C{arg1}
       @rtype: L{escript.Symbol}, C{float}, C{int}, L{escript.Data},
               C{numarray.NumArray}
       @note: The shape of both arguments is matched according to the rules
              used in L{matchShape}
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]):
          return numarray.zeros(getShape(args[0]),numarray.Float64)
       elif testForZero(args[1]):
          return numarray.ones(getShape(args[1]),numarray.Float64)
       elif isinstance(args[0],Symbol) or isinstance(args[1],Symbol):
          return Power_Symbol(args[0],args[1])
       elif isinstance(args[0],numarray.NumArray) and not isinstance(args[1],numarray.NumArray):
          return exp(args[1]*log(args[0]))
       else:
           return args[0]**args[1]

class Power_Symbol(DependendSymbol):
   """
   Symbol representing the first argument to the power of the second argument.
   """
   def __init__(self,arg0,arg1):
       """
       Initialization of the L{Symbol} representing raising the first argument
       to the power of the second.

       @param arg0: basis
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: exponent
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: if both arguments do not have the same shape
       @note: if both arguments have a spatial dimension, they must be equal
       """
       sh0=getShape(arg0)
       sh1=getShape(arg1)
       if not sh0==sh1:
          raise ValueError,"Power_Symbol: shape of arguments must match"
       d0=pokeDim(arg0)
       d1=pokeDim(arg1)
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0,args=[arg0,arg1])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str" or format=="text":
         return "(%s)**(%s)"%(argstrs[0],argstrs[1])
      elif format=="escript":
         return "power(%s,%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return power(args[0],args[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         myargs=self.getArgument()
         dargs=self.getDifferentiatedArguments(arg)
         return mult(self,add(mult(log(myargs[0]),dargs[1]),mult(quotient(myargs[1],myargs[0]),dargs[0])))

def maximum(*args):
    """
    The maximum over arguments C{args}.

    @param args: arguments
    @type args: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or
                C{float}
    @return: an object which in each entry gives the maximum of the
             corresponding values in C{args}
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or
            C{float} depending on the input
    """
    out=None
    for a in args:
       if out==None:
          out=a*1.
       else:
          if isinstance(out,escript.Data) and isinstance(a,escript.Data):
             out.copyWithMask(a,wherePositive(a-out))
          else:
             diff=add(a,-out)
             out=add(out,mult(wherePositive(diff),diff))
    return out

def minimum(*args):
    """
    The minimum over arguments C{args}.

    @param args: arguments
    @type args: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or
                C{float}
    @return: an object which gives in each entry the minimum of the
             corresponding values in C{args}
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or
            C{float} depending on the input
    """
    out=None
    for a in args:
       if out==None:
          out=a*1.
       else:
          if isinstance(out,escript.Data) and isinstance(a,escript.Data):
             out.copyWithMask(a,whereNegative(a-out))
          else:
             diff=add(a,-out)
             out=add(out,mult(whereNegative(diff),diff))
    return out

def clip(arg,minval=None,maxval=None):
    """
    Cuts the values of C{arg} between C{minval} and C{maxval}.

    @param arg: argument
    @type arg: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or
               C{float}
    @param minval: lower range. If None no lower range is applied
    @type minval: C{float} or C{None}
    @param maxval: upper range. If None no upper range is applied
    @type maxval: C{float} or C{None}
    @return: an object that contains all values from C{arg} between C{minval}
             and C{maxval}
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{int} or
            C{float} depending on the input
    @raise ValueError: if C{minval>maxval}
    """
    if not minval==None and not maxval==None:
       if minval>maxval:
          raise ValueError,"minval = %s must be less than maxval %s"%(minval,maxval)
    if minval == None:
        tmp=arg
    else:
        tmp=maximum(minval,arg)
    if maxval == None:
        return tmp
    else:
        return minimum(tmp,maxval)


def inner(arg0,arg1):
    """
    Inner product of the two arguments. The inner product is defined as:

    C{out=S{Sigma}_s arg0[s]*arg1[s]}

    where s runs through C{arg0.Shape}.

    C{arg0} and C{arg1} must have the same shape.

    @param arg0: first argument
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @return: the inner product of C{arg0} and C{arg1} at each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float}
            depending on the input
    @raise ValueError: if the shapes of the arguments are not identical
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not sh0==sh1:
        raise ValueError,"inner: shape of arguments does not match"
    return generalTensorProduct(arg0,arg1,axis_offset=len(sh0))

def outer(arg0,arg1):
    """
    The outer product of the two arguments. The outer product is defined as:

    C{out[t,s]=arg0[t]*arg1[s]}

    where
        - s runs through C{arg0.Shape}
        - t runs through C{arg1.Shape}

    @param arg0: first argument
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @return: the outer product of C{arg0} and C{arg1} at each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    return generalTensorProduct(arg0,arg1,axis_offset=0)

def matrixmult(arg0,arg1):
    """
    See L{matrix_mult}.
    """
    return matrix_mult(arg0,arg1)

def matrix_mult(arg0,arg1):
    """
    matrix-matrix or matrix-vector product of the two arguments.

    C{out[s0]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0]}

    or

    C{out[s0,s1]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0,s1]}

    The second dimension of C{arg0} and the first dimension of C{arg1} must
    match.

    @param arg0: first argument of rank 2
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of at least rank 1
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the matrix-matrix or matrix-vector product of C{arg0} and C{arg1}
             at each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    @raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not len(sh0)==2 :
        raise ValueError,"first argument must have rank 2"
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError,"second argument must have rank 1 or 2"
    return generalTensorProduct(arg0,arg1,axis_offset=1)

def tensormult(arg0,arg1):
    """
    See L{tensor_mult}.
    """
    return tensor_mult(arg0,arg1)

def tensor_mult(arg0,arg1):
    """
    The tensor product of the two arguments.

    For C{arg0} of rank 2 this is

    C{out[s0]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0]}

    or

    C{out[s0,s1]=S{Sigma}_{r0} arg0[s0,r0]*arg1[r0,s1]}

    and for C{arg0} of rank 4 this is

    C{out[s0,s1,s2,s3]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2,s3]}

    or

    C{out[s0,s1,s2]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2]}

    or

    C{out[s0,s1]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1]}

    In the first case the second dimension of C{arg0} and the last dimension of
    C{arg1} must match and in the second case the two last dimensions of C{arg0}
    must match the two first dimensions of C{arg1}.

    @param arg0: first argument of rank 2 or 4
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of shape greater than 1 or 2 depending on the
                 rank of C{arg0}
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the tensor product of C{arg0} and C{arg1} at each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
       return generalTensorProduct(arg0,arg1,axis_offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
       return generalTensorProduct(arg0,arg1,axis_offset=2)
    else:
        raise ValueError,"tensor_mult: first argument must have rank 2 or 4"

def generalTensorProduct(arg0,arg1,axis_offset=0):
    """
    Generalized tensor product.

    C{out[s,t]=S{Sigma}_r arg0[s,r]*arg1[r,t]}

    where
        - s runs through C{arg0.Shape[:arg0.Rank-axis_offset]}
        - r runs through C{arg0.Shape[:axis_offset]}
        - t runs through C{arg1.Shape[axis_offset:]}

    @param arg0: first argument
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @return: the general tensor product of C{arg0} and C{arg1} at each data
             point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    if isinstance(arg0,float) and isinstance(arg1,float): return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg1 are both numarray.NumArray or escript.Data or
    # Symbols
    if isinstance(arg0,numarray.NumArray):
       if isinstance(arg1,Symbol):
           return GeneralTensorProduct_Symbol(arg0,arg1,axis_offset)
       else:
           if not arg0.shape[arg0.rank-axis_offset:]==arg1.shape[:axis_offset]:
               raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset)
           arg0_c=arg0.copy()
           arg1_c=arg1.copy()
           sh0,sh1=arg0.shape,arg1.shape
           d0,d1,d01=1,1,1
           for i in sh0[:arg0.rank-axis_offset]: d0*=i
           for i in sh1[axis_offset:]: d1*=i
           for i in sh1[:axis_offset]: d01*=i
           arg0_c.resize((d0,d01))
           arg1_c.resize((d01,d1))
           out=numarray.zeros((d0,d1),numarray.Float64)
           for i0 in range(d0):
                    for i1 in range(d1):
                         out[i0,i1]=numarray.sum(arg0_c[i0,:]*arg1_c[:,i1])
           out.resize(sh0[:arg0.rank-axis_offset]+sh1[axis_offset:])
           return out
    elif isinstance(arg0,escript.Data):
       if isinstance(arg1,Symbol):
           return GeneralTensorProduct_Symbol(arg0,arg1,axis_offset)
       else:
           return escript_generalTensorProduct(arg0,arg1,axis_offset) # this calls has to be replaced by escript._generalTensorProduct(arg0,arg1,axis_offset)
    else:
       return GeneralTensorProduct_Symbol(arg0,arg1,axis_offset)

class GeneralTensorProduct_Symbol(DependendSymbol):
   """
   Symbol representing the general tensor product of two arguments.
   """
   def __init__(self,arg0,arg1,axis_offset=0):
       """
       Initialization of L{Symbol} representing the general tensor product of
       two arguments.

       @param arg0: first argument
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second argument
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: illegal dimension
       @note: if both arguments have a spatial dimension, they must be equal.
       """
       sh_arg0=getShape(arg0)
       sh_arg1=getShape(arg1)
       sh0=sh_arg0[:len(sh_arg0)-axis_offset]
       sh01=sh_arg0[len(sh_arg0)-axis_offset:]
       sh10=sh_arg1[:axis_offset]
       sh1=sh_arg1[axis_offset:]
       if not sh01==sh10:
           raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset)
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0+sh1,args=[arg0,arg1,axis_offset])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str" or format=="text":
         return "generalTensorProduct(%s,%s,axis_offset=%s)"%(argstrs[0],argstrs[1],argstrs[2])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return generalTensorProduct(args[0],args[1],args[2])

def escript_generalTensorProduct(arg0,arg1,axis_offset,transpose=0):
    "arg0 and arg1 are both Data objects but not necessarily on the same function space. They could be identical!!!"
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, transpose)

def transposed_matrix_mult(arg0,arg1):
    """
    transposed(matrix)-matrix or transposed(matrix)-vector product of the two
    arguments.

    C{out[s0]=S{Sigma}_{r0} arg0[r0,s0]*arg1[r0]}

    or

    C{out[s0,s1]=S{Sigma}_{r0} arg0[r0,s0]*arg1[r0,s1]}

    The function call C{transposed_matrix_mult(arg0,arg1)} is equivalent to
    C{matrix_mult(transpose(arg0),arg1)}.

    The first dimension of C{arg0} and C{arg1} must match.

    @param arg0: first argument of rank 2
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of at least rank 1
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the product of the transpose of C{arg0} and C{arg1} at each data
             point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    @raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not len(sh0)==2 :
        raise ValueError,"first argument must have rank 2"
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError,"second argument must have rank 1 or 2"
    return generalTransposedTensorProduct(arg0,arg1,axis_offset=1)

def transposed_tensor_mult(arg0,arg1):
    """
    The tensor product of the transpose of the first and the second argument.

    For C{arg0} of rank 2 this is

    C{out[s0]=S{Sigma}_{r0} arg0[r0,s0]*arg1[r0]}

    or

    C{out[s0,s1]=S{Sigma}_{r0} arg0[r0,s0]*arg1[r0,s1]}

    and for C{arg0} of rank 4 this is

    C{out[s0,s1,s2,s3]=S{Sigma}_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1,s2,s3]}

    or

    C{out[s0,s1,s2]=S{Sigma}_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1,s2]}

    or

    C{out[s0,s1]=S{Sigma}_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1]}

    In the first case the first dimension of C{arg0} and the first dimension of
    C{arg1} must match and in the second case the two first dimensions of
    C{arg0} must match the two first dimensions of C{arg1}.

    The function call C{transposed_tensor_mult(arg0,arg1)} is equivalent to
    C{tensor_mult(transpose(arg0),arg1)}.

    @param arg0: first argument of rank 2 or 4
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of shape greater of 1 or 2 depending on the
                 rank of C{arg0}
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the tensor product of transpose of arg0 and arg1 at each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
       return generalTransposedTensorProduct(arg0,arg1,axis_offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
       return generalTransposedTensorProduct(arg0,arg1,axis_offset=2)
    else:
        raise ValueError,"first argument must have rank 2 or 4"

def generalTransposedTensorProduct(arg0,arg1,axis_offset=0):
    """
    Generalized tensor product of transposed of C{arg0} and C{arg1}.

    C{out[s,t]=S{Sigma}_r arg0[r,s]*arg1[r,t]}

    where
        - s runs through C{arg0.Shape[axis_offset:]}
        - r runs through C{arg0.Shape[:axis_offset]}
        - t runs through C{arg1.Shape[axis_offset:]}

    The function call C{generalTransposedTensorProduct(arg0,arg1,axis_offset)}
    is equivalent to
    C{generalTensorProduct(transpose(arg0,arg0.rank-axis_offset),arg1,axis_offset)}.

    @param arg0: first argument
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @return: the general tensor product of C{transpose(arg0)} and C{arg1} at
             each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    if isinstance(arg0,float) and isinstance(arg1,float): return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg1 are both numarray.NumArray or escript.Data or
    # Symbols
    if isinstance(arg0,numarray.NumArray):
       if isinstance(arg1,Symbol):
           return GeneralTransposedTensorProduct_Symbol(arg0,arg1,axis_offset)
       else:
           if not arg0.shape[:axis_offset]==arg1.shape[:axis_offset]:
               raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset)
           arg0_c=arg0.copy()
           arg1_c=arg1.copy()
           sh0,sh1=arg0.shape,arg1.shape
           d0,d1,d01=1,1,1
           for i in sh0[axis_offset:]: d0*=i
           for i in sh1[axis_offset:]: d1*=i
           for i in sh0[:axis_offset]: d01*=i
           arg0_c.resize((d01,d0))
           arg1_c.resize((d01,d1))
           out=numarray.zeros((d0,d1),numarray.Float64)
           for i0 in range(d0):
                    for i1 in range(d1):
                         out[i0,i1]=numarray.sum(arg0_c[:,i0]*arg1_c[:,i1])
           out.resize(sh0[axis_offset:]+sh1[axis_offset:])
           return out
    elif isinstance(arg0,escript.Data):
       if isinstance(arg1,Symbol):
           return GeneralTransposedTensorProduct_Symbol(arg0,arg1,axis_offset)
       else:
           # this call has to be replaced by escript._generalTensorProduct(arg0,arg1,axis_offset)
           return escript_generalTransposedTensorProduct(arg0,arg1,axis_offset)
    else:
       return GeneralTransposedTensorProduct_Symbol(arg0,arg1,axis_offset)

class GeneralTransposedTensorProduct_Symbol(DependendSymbol):
   """
   Symbol representing the general tensor product of the transpose of C{arg0}
   and C{arg1}
   """
   def __init__(self,arg0,arg1,axis_offset=0):
       """
       Initialization of L{Symbol} representing tensor product of the
       transpose of C{arg0} and C{arg1}.

       @param arg0: first argument
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second argument
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: inconsistent dimensions of arguments
       @note: if both arguments have a spatial dimension, they must be equal.
       """
       sh_arg0=getShape(arg0)
       sh_arg1=getShape(arg1)
       sh01=sh_arg0[:axis_offset]
       sh10=sh_arg1[:axis_offset]
       sh0=sh_arg0[axis_offset:]
       sh1=sh_arg1[axis_offset:]
       if not sh01==sh10:
           raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset)
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0+sh1,args=[arg0,arg1,axis_offset])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str" or format=="text":
         return "generalTransposedTensorProduct(%s,%s,axis_offset=%s)"%(argstrs[0],argstrs[1],argstrs[2])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return generalTransposedTensorProduct(args[0],args[1],args[2])

# this should be escript._generalTransposedTensorProduct
def escript_generalTransposedTensorProduct(arg0,arg1,axis_offset):
    "arg0 and arg1 are both Data objects but not necessarily on the same function space. They could be identical!!!"
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, 1)

def matrix_transposed_mult(arg0,arg1):
    """
    matrix-transposed(matrix) product of the two arguments.

    C{out[s0,s1]=S{Sigma}_{r0} arg0[s0,r0]*arg1[s1,r0]}

    The function call C{matrix_transposed_mult(arg0,arg1)} is equivalent to
    C{matrix_mult(arg0,transpose(arg1))}.

    The last dimensions of C{arg0} and C{arg1} must match.

    @param arg0: first argument of rank 2
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of rank 1 or 2
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the product of C{arg0} and the transposed of C{arg1} at each data
             point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    @raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not len(sh0)==2 :
        raise ValueError,"first argument must have rank 2"
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError,"second argument must have rank 1 or 2"
    return generalTensorTransposedProduct(arg0,arg1,axis_offset=1)

def tensor_transposed_mult(arg0,arg1):
    """
    The tensor product of the first and the transpose of the second argument.

    For C{arg0} of rank 2 this is

    C{out[s0,s1]=S{Sigma}_{r0} arg0[s0,r0]*arg1[s1,r0]}

    and for C{arg0} of rank 4 this is

    C{out[s0,s1,s2,s3]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[s2,s3,r0,r1]}

    or

    C{out[s0,s1,s2]=S{Sigma}_{r0,r1} arg0[s0,s1,r0,r1]*arg1[s2,r0,r1]}

    In the first case the the second dimension of C{arg0} and C{arg1} must
    match and in the second case the two last dimensions of C{arg0} must match
    the two last dimensions of C{arg1}.

    The function call C{tensor_transpose_mult(arg0,arg1)} is equivalent to
    C{tensor_mult(arg0,transpose(arg1))}.

    @param arg0: first argument of rank 2 or 4
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @param arg1: second argument of shape greater of 1 or 2 depending on rank
                 of C{arg0}
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}
    @return: the tensor product of the transposed of C{arg0} and C{arg1} at
             each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
       return generalTensorTransposedProduct(arg0,arg1,axis_offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
       return generalTensorTransposedProduct(arg0,arg1,axis_offset=2)
    else:
        raise ValueError,"first argument must have rank 2 or 4"

def generalTensorTransposedProduct(arg0,arg1,axis_offset=0):
    """
    Generalized tensor product of C{arg0} and transpose of C{arg1}.

    C{out[s,t]=S{Sigma}_r arg0[s,r]*arg1[t,r]}

    where
        - s runs through C{arg0.Shape[:arg0.Rank-axis_offset]}
        - r runs through C{arg0.Shape[arg1.Rank-axis_offset:]}
        - t runs through C{arg1.Shape[arg1.Rank-axis_offset:]}

    The function call C{generalTensorTransposedProduct(arg0,arg1,axis_offset)}
    is equivalent to
    C{generalTensorProduct(arg0,transpose(arg1,arg1.Rank-axis_offset),axis_offset)}.

    @param arg0: first argument
    @type arg0: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @param arg1: second argument
    @type arg1: C{numarray.NumArray}, L{escript.Data}, L{Symbol}, C{float},
                C{int}
    @return: the general tensor product of C{arg0} and C{transpose(arg1)} at
             each data point
    @rtype: C{numarray.NumArray}, L{escript.Data}, L{Symbol} depending on the
            input
    """
    if isinstance(arg0,float) and isinstance(arg1,float): return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg1 are both numarray.NumArray or escript.Data or
    # Symbols
    if isinstance(arg0,numarray.NumArray):
       if isinstance(arg1,Symbol):
           return GeneralTensorTransposedProduct_Symbol(arg0,arg1,axis_offset)
       else:
           if not arg0.shape[arg0.rank-axis_offset:]==arg1.shape[arg1.rank-axis_offset:]:
               raise ValueError,"dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset)
           arg0_c=arg0.copy()
           arg1_c=arg1.copy()
           sh0,sh1=arg0.shape,arg1.shape
           d0,d1,d01=1,1,1
           for i in sh0[:arg0.rank-axis_offset]: d0*=i
           for i in sh1[:arg1.rank-axis_offset]: d1*=i
           for i in sh1[arg1.rank-axis_offset:]: d01*=i
           arg0_c.resize((d0,d01))
           arg1_c.resize((d1,d01))
           out=numarray.zeros((d0,d1),numarray.Float64)
           for i0 in range(d0):
                    for i1 in range(d1):
                         out[i0,i1]=numarray.sum(arg0_c[i0,:]*arg1_c[i1,:])
           out.resize(sh0[:arg0.rank-axis_offset]+sh1[:arg1.rank-axis_offset])
           return out
    elif isinstance(arg0,escript.Data):
       if isinstance(arg1,Symbol):
           return GeneralTensorTransposedProduct_Symbol(arg0,arg1,axis_offset)
       else:
           # this call has to be replaced by escript._generalTensorProduct(arg0,arg1,axis_offset)
           return escript_generalTensorTransposedProduct(arg0,arg1,axis_offset)
    else:
       return GeneralTensorTransposedProduct_Symbol(arg0,arg1,axis_offset)

class GeneralTensorTransposedProduct_Symbol(DependendSymbol):
   """
   Symbol representing the general tensor product of C{arg0} and the transpose
   of C{arg1}.
   """
   def __init__(self,arg0,arg1,axis_offset=0):
       """
       Initialization of L{Symbol} representing the general tensor product of
       C{arg0} and the transpose of C{arg1}.

       @param arg0: first argument
       @type arg0: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @param arg1: second argument
       @type arg1: L{escript.Symbol}, C{float}, L{escript.Data},
                   C{numarray.NumArray}
       @raise ValueError: inconsistent dimensions of arguments
       @note: if both arguments have a spatial dimension, they must be equal.
       """
       sh_arg0=getShape(arg0)
       sh_arg1=getShape(arg1)
       sh0=sh_arg0[:len(sh_arg0)-axis_offset]
       sh01=sh_arg0[len(sh_arg0)-axis_offset:]
       sh10=sh_arg1[len(sh_arg1)-axis_offset:]
       sh1=sh_arg1[:len(sh_arg1)-axis_offset]
       if not sh01==sh10:
           raise ValueError,"dimensions of last %s components in left argument don't match the last %s components in the right argument."%(axis_offset,axis_offset)
       DependendSymbol.__init__(self,dim=commonDim(arg0,arg1),shape=sh0+sh1,args=[arg0,arg1,axis_offset])

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str" or format=="text":
         return "generalTensorTransposedProduct(%s,%s,axis_offset=%s)"%(argstrs[0],argstrs[1],argstrs[2])
      else:
         raise NotImplementedError,"%s does not provide program code for format %s."%(str(self),format)

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         args=self.getSubstitutedArguments(argvals)
         return generalTensorTransposedProduct(args[0],args[1],args[2])

# this should be escript._generalTensorTransposedProduct
def escript_generalTensorTransposedProduct(arg0,arg1,axis_offset):
    "arg0 and arg1 are both Data objects but not necessarily on the same function space. They could be identical!!!"
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, 2)

#=========================================================
#  functions dealing with spatial dependency
#=========================================================
def grad(arg,where=None):
    """
    Returns the spatial gradient of C{arg} at C{where}.

    If C{g} is the returned object, then

      - if C{arg} is rank 0 C{g[s]} is the derivative of C{arg} with respect to
        the C{s}-th spatial dimension
      - if C{arg} is rank 1 C{g[i,s]} is the derivative of C{arg[i]} with
        respect to the C{s}-th spatial dimension
      - if C{arg} is rank 2 C{g[i,j,s]} is the derivative of C{arg[i,j]} with
        respect to the C{s}-th spatial dimension
      - if C{arg} is rank 3 C{g[i,j,k,s]} is the derivative of C{arg[i,j,k]}
        with respect to the C{s}-th spatial dimension.

    @param arg: function of which the gradient is to be calculated. Its rank
                has to be less than 3.
    @type arg: L{escript.Data} or L{Symbol}
    @param where: FunctionSpace in which the gradient is calculated.
                  If not present or C{None} an appropriate default is used.
    @type where: C{None} or L{escript.FunctionSpace}
    @return: gradient of C{arg}
    @rtype: L{escript.Data} or L{Symbol}
    """
    if isinstance(arg,Symbol):
       return Grad_Symbol(arg,where)
    elif isinstance(arg,escript.Data):
       if where==None:
          return arg._grad()
       else:
          return arg._grad(where)
    else:
       raise TypeError,"grad: Unknown argument type."

class Grad_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the gradient operator.
   """
   def __init__(self,arg,where=None):
      """
      Initialization of gradient L{Symbol} with argument C{arg}.

      @param arg: argument of function
      @type arg: L{Symbol}
      @param where: FunctionSpace in which the gradient will be calculated.
                    If not present or C{None} an appropriate default is used.
      @type where: C{None} or L{escript.FunctionSpace}
      """
      d=arg.getDim()
      if d==None:
         raise ValueError,"argument must have a spatial dimension"
      super(Grad_Symbol,self).__init__(args=[arg,where],shape=arg.getShape()+(d,),dim=d)

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "grad(%s,where=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"Grad_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return grad(arg[0],where=arg[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to arg
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return grad(self.getDifferentiatedArguments(arg)[0],where=self.getArgument()[1])

def integrate(arg,where=None):
    """
    Returns the integral of the function C{arg} over its domain. If C{where} is
    present C{arg} is interpolated to C{where} before integration.

    @param arg: the function which is integrated
    @type arg: L{escript.Data} or L{Symbol}
    @param where: FunctionSpace in which the integral is calculated.
                  If not present or C{None} an appropriate default is used.
    @type where: C{None} or L{escript.FunctionSpace}
    @return: integral of C{arg}
    @rtype: C{float}, C{numarray.NumArray} or L{Symbol}
    """
    if isinstance(arg,Symbol):
       return Integrate_Symbol(arg,where)
    elif isinstance(arg,escript.Data):
       if not where==None: arg=escript.Data(arg,where)
       if arg.getRank()==0:
          return arg._integrate()[0]
       else:
          return arg._integrate()
    else:
       arg2=escript.Data(arg,where)
       if arg2.getRank()==0:
          return arg2._integrate()[0]
       else:
          return arg2._integrate()

class Integrate_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the spatial integration operator.
   """
   def __init__(self,arg,where=None):
      """
      Initialization of integration L{Symbol} with argument C{arg}.

      @param arg: argument of the integration
      @type arg: L{Symbol}
      @param where: FunctionSpace in which the integration will be performed.
                    If not present or C{None} an appropriate default is used.
      @type where: C{None} or L{escript.FunctionSpace}
      """
      super(Integrate_Symbol,self).__init__(args=[arg,where],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "integrate(%s,where=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"Integrate_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return integrate(arg[0],where=arg[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: typically L{Symbol} but other types such as C{float},
              L{escript.Data}, C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return integrate(self.getDifferentiatedArguments(arg)[0],where=self.getArgument()[1])


def interpolate(arg,where):
    """
    Interpolates the function into the L{FunctionSpace} C{where}. If the
    argument C{arg} has the requested function space C{where} no interpolation
    is performed and C{arg} is returned.

    @param arg: interpolant
    @type arg: L{escript.Data} or L{Symbol}
    @param where: L{FunctionSpace} to be interpolated to
    @type where: L{escript.FunctionSpace}
    @return: interpolated argument
    @rtype: C{escript.Data} or L{Symbol}
    """
    if isinstance(arg,Symbol):
       return Interpolate_Symbol(arg,where)
    elif isinstance(arg,escript.Data):
       if arg.isEmpty():
          return arg
       elif where == arg.getFunctionSpace():
          return arg
       else:
          return escript.Data(arg,where)
    else:
       return escript.Data(arg,where)

class Interpolate_Symbol(DependendSymbol):
   """
   L{Symbol} representing the result of the interpolation operator.
   """
   def __init__(self,arg,where):
      """
      Initialization of interpolation L{Symbol} with argument C{arg}.

      @param arg: argument of the interpolation
      @type arg: L{Symbol}
      @param where: L{FunctionSpace} into which the argument is interpolated
      @type where: L{escript.FunctionSpace}
      """
      super(Interpolate_Symbol,self).__init__(args=[arg,where],shape=arg.getShape(),dim=arg.getDim())

   def getMyCode(self,argstrs,format="escript"):
      """
      Returns program code that can be used to evaluate the symbol.

      @param argstrs: a string for each argument representing the argument
                      for the evaluation
      @type argstrs: C{str} or a C{list} of length 2 of C{str}
      @param format: specifies the format to be used. At the moment only
                     "escript", "str" and "text" are supported.
      @type format: C{str}
      @return: a piece of program code which can be used to evaluate the
               expression assuming the values for the arguments are available
      @rtype: C{str}
      @raise NotImplementedError: if no implementation for the given format
                                  is available
      """
      if format=="escript" or format=="str"  or format=="text":
         return "interpolate(%s,where=%s)"%(argstrs[0],argstrs[1])
      else:
         raise NotImplementedError,"Interpolate_Symbol does not provide program code for format %s."%format

   def substitute(self,argvals):
      """
      Assigns new values to symbols in the definition of the symbol.
      The method replaces the L{Symbol} u by argvals[u] in the expression
      defining this object.

      @param argvals: new values assigned to symbols
      @type argvals: C{dict} with keywords of type L{Symbol}
      @return: result of the substitution process. Operations are executed as
               much as possible.
      @rtype: L{escript.Symbol}, C{float}, L{escript.Data}, C{numarray.NumArray}
              depending on the degree of substitution
      @raise TypeError: if a value for a L{Symbol} cannot be substituted
      """
      if argvals.has_key(self):
         arg=argvals[self]
         if self.isAppropriateValue(arg):
            return arg
         else:
            raise TypeError,"%s: new value is not appropriate."%str(self)
      else:
         arg=self.getSubstitutedArguments(argvals)
         return interpolate(arg[0],where=arg[1])

   def diff(self,arg):
      """
      Differential of this object.

      @param arg: the derivative is calculated with respect to C{arg}
      @type arg: L{escript.Symbol}
      @return: derivative with respect to C{arg}
      @rtype: L{Symbol} but other types such as L{escript.Data},
              C{numarray.NumArray} are possible
      """
      if arg==self:
         return identity(self.getShape())
      else:
         return interpolate(self.getDifferentiatedArguments(arg)[0],where=self.getArgument()[1])


def div(arg,where=None):
    """
    Returns the divergence of C{arg} at C{where}.

    @param arg: function of which the divergence is to be calculated. Its
                shape has to be (d,) where d is the spatial dimension.
    @type arg: L{escript.Data} or L{Symbol}
    @param where: L{FunctionSpace} in which the divergence will be calculated.
                  If not present or C{None} an appropriate default is used.
    @type where: C{None} or L{escript.FunctionSpace}
    @return: divergence of C{arg}
    @rtype: L{escript.Data} or L{Symbol}
    """
    if isinstance(arg,Symbol):
        dim=arg.getDim()
    elif isinstance(arg,escript.Data):
        dim=arg.getDomain().getDim()
    else:
        raise TypeError,"div: argument type not supported"
    if not arg.getShape()==(dim,):
      raise ValueError,"div: expected shape is (%s,)"%dim
    return trace(grad(arg,where))

def jump(arg,domain=None):
    """
    Returns the jump of C{arg} across the continuity of the domain.

    @param arg: argument
    @type arg: L{escript.Data} or L{Symbol}
    @param domain: the domain where the discontinuity is located. If domain is
                   not present or equal to C{None} the domain of C{arg} is used.
                   If C{arg} is a L{Symbol} the domain must be present.
    @type domain: C{None} or L{escript.Domain}
    @return: jump of C{arg}
    @rtype: L{escript.Data} or L{Symbol}
    """
    if domain==None: domain=arg.getDomain()
    return interpolate(arg,escript.FunctionOnContactOne(domain))-interpolate(arg,escript.FunctionOnContactZero(domain))

def L2(arg):
    """
    Returns the L2 norm of C{arg} at C{where}.

    @param arg: function of which the L2 norm is to be calculated
    @type arg: L{escript.Data} or L{Symbol}
    @return: L2 norm of C{arg}
    @rtype: L{float} or L{Symbol}
    @note: L2(arg) is equivalent to C{sqrt(integrate(inner(arg,arg)))}
    """
    return sqrt(integrate(inner(arg,arg)))

def getClosestValue(arg,origin=0):
    """
    Returns the value in C{arg} which is closest to origin.

    @param arg: function
    @type arg: L{escript.Data} or L{Symbol}
    @param origin: reference value
    @type origin: C{float} or L{escript.Data}
    @return: value in C{arg} closest to origin
    @rtype: C{numarray.NumArray} or L{Symbol}
    """
    return arg.getValueOfGlobalDataPoint(*(length(arg-origin).minGlobalDataPoint()))

def normalize(arg,zerolength=0):
    """
    Returns the normalized version of C{arg} (=C{arg/length(arg)}).

    @param arg: function
    @type arg: L{escript.Data} or L{Symbol}
    @param zerolength: relative tolerance for arg == 0
    @type zerolength: C{float}
    @return: normalized C{arg} where C{arg} is non-zero, and zero elsewhere
    @rtype: L{escript.Data} or L{Symbol}
    """
    l=length(arg)
    m=whereZero(l,zerolength*Lsup(l))
    mm=1-m
    return arg*(mm/(l*mm+m))

def deviatoric(arg):
    """
    Returns the deviatoric version of C{arg}.
    """
    return arg-(trace(arg)/trace(kronecker(arg.getDomain())))*kronecker(arg.getDomain())

def diameter(domain):
    """
    Returns the diameter of a domain.

    @param domain: a domain
    @type domain: L{escript.Domain}
    @rtype: C{float}
    """
    x=domain.getX()
    out=0.
    for i in xrange(domain.getDim()):
       x_i=x[i]
       out=max(out,sup(x_i)-inf(x_i))
    return out
#=============================
#

def reorderComponents(arg,index):
    """
    Resorts the components of C{arg} according to index.

    """
    raise NotImplementedError

def showEscriptParams():
    """
    Displays the parameters escript recognises with an explanation.
    """
    p=listEscriptParams()
    for name,desc in p:
	print name+':\t'+desc

