# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file in the repository root for contributors and development history
# https://github.com/LutzGross/esys-escript.github.io/blob/master/CREDITS
#
##############################################################################

"""
Utility functions for escript.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
:var EPSILON: smallest positive value with 1.<1.+EPSILON
:var DBLE_MAX: largest positive float
"""


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"
__author__="Lutz Gross, l.gross@uq.edu.au"


import math
import cmath
import os
import warnings
import numpy
import numbers
warnings.simplefilter('default', category=DeprecationWarning)
# suppress the following which comes from sympy with python 3.5
warnings.filterwarnings('ignore', category=DeprecationWarning, message='inspect.getargspec.*')

from . import escriptcpp as escore
from .escriptcpp import C_GeneralTensorProduct, Data
from .escriptcpp import getVersion, getMPIRankWorld, getMPIWorldMax, hasFeature
#from .escriptcpp import printParallelThreadCounts
from .escriptcpp import listEscriptParams
from .start import HAVE_SYMBOLS

if HAVE_SYMBOLS:
    from . import symboliccore as sym
else:
    class sym:
        Symbol=type(memoryview(b'1')) #This type should never be passed to util.py




#=========================================================
#   some helpers:
#=========================================================
def getEpsilon():
     """
     Returns the machine epsilon (smallest positive value where 1.0 + epsilon != 1.0).

     :return: machine precision epsilon
     :rtype: ``float``
     """
     return escore.getMachinePrecision()
EPSILON=getEpsilon()

def getMaxFloat():
     """
     Returns the largest representable positive floating point number.

     :return: maximum float value
     :rtype: ``float``
     """
     return escore.getMaxFloat()
DBLE_MAX=getMaxFloat()

def getTagNames(domain):
    """
    Returns a list of tag names used by the domain.

    :param domain: a domain object
    :type domain: `escript.Domain`
    :return: a list of tag names used by the domain
    :rtype: ``list`` of ``str``
    """
    return [n.strip() for n in domain.showTagNames().split(",") ]

def insertTagNames(domain,**kwargs):
    """
    Inserts tag names into the domain.

    :param domain: a domain object
    :type domain: ``escript.Domain``
    :keyword <tag_name>: tag key assigned to <tag_name>
    :type <tag_name>: ``int``
    """
    for  k in kwargs:
         domain.setTagMap(k,kwargs[k])

def insertTaggedValues(target,**kwargs):
    """
    Inserts tagged values into the target using tag names.

    :param target: data to be filled by tagged values
    :type target: `escript.Data`
    :keyword <tag_name>: value to be used for <tag_name>
    :type <tag_name>: ``float`` or ``numpy.ndarray``
    :return: ``target``
    :rtype: `escript.Data`
    """
    for k in kwargs:
        target.setTaggedValue(k,kwargs[k])
    return target


def interpolateTable(tab, dat, start, step, undef=1.e50, check_boundaries=False):
    """
    Interpolates values from a lookup table.

    :param tab: the lookup table array
    :type tab: ``numpy.ndarray``
    :param dat: input data to interpolate
    :type dat: `Data` or ``numpy.ndarray``
    :param start: starting coordinate(s) for the table
    :type start: ``float`` or ``tuple`` of ``float``
    :param step: step size(s) for the table
    :type step: ``float`` or ``tuple`` of ``float``
    :param undef: value to use for undefined/out-of-bounds results
    :type undef: ``float``
    :param check_boundaries: if ``True``, check that input is within table bounds
    :type check_boundaries: ``bool``
    :return: interpolated values
    :rtype: `Data` or ``numpy.ndarray``

    .. deprecated::
       This function is deprecated and is known to contain bugs.
    """
    print("WARNING: This function is deprecated and is known to contain bugs.")
    try:
        dim=len(start)
    except TypeError:
        start=(start,)
        dim=1
    try:
        slen=len(step)
    except TypeError:
        step=(step,)
        slen=1
    if dim<1 or dim>3:
        raise ValueError("Length of start list must be between 1 and 3.")
    if dim!=slen:
        raise ValueError("Length of start and step must be the same.")
    dshape=dat.getShape()
    if len(dshape)==0:
        datdim=0
        firstdim=dat
    else:
        datdim=dshape[0]
        firstdim=dat[0]
    #So now we know firstdim is a scalar
    if (dim==1 and datdim>1) or (dim>1 and datdim!=dim):
        print((dim, datdim))
        raise ValueError("The dimension of dat must be equal to the length of start.")
    if dim==3:
        d1=dat[1]
        d2=dat[2]
        return firstdim._interpolateTable3d(tab, start[0], step[0], d1, start[1], step[1], d2, start[2], step[2], undef, check_boundaries)
    if dim==2:
        d1=dat[1]
        return firstdim.interpolateTable(tab, start[0], step[0], d1, start[1], step[1], undef, check_boundaries)
#   return d1.interpolateTable(tab, start[1], step[1], firstdim, start[0], step[0], undef, check_boundaries)
    else:
        return firstdim.interpolateTable(tab, start[0], step[0], undef, check_boundaries)


def saveDataCSV(filename, append=False, refid=False, sep=", ", csep="_", **data):
    """
    Writes `Data` objects to a CSV file.
    These objects must have compatible FunctionSpaces, i.e. it must be possible
    to interpolate all data to one `FunctionSpace`. Note, that with more than
    one MPI rank this function  will fail for some function spaces on some
    domains.

    :param filename: file to save data to.
    :type filename: ``string``
    :param append: If ``True``, then open file at end rather than beginning
    :type append: ``bool``
    :param refid: If ``True``, then a list of reference ids will be printed in the first column
    :type refid: ``bool``
    :param sep: separator between fields
    :type sep: ``string``
    :param csep: separator for components of rank 2 and above (e.g. '_' -> c0_1)

    The keyword args are Data objects to save.
    If a scalar `Data` object is passed with the name ``mask``, then only
    samples which correspond to positive values in ``mask`` will be output.
    Example::

        s=Scalar(..)
        v=Vector(..)
        t=Tensor(..)
        f=float()
        saveDataCSV("f.csv", a=s, b=v, c=t, d=f)
    
    Will result in a file
    
    a, b0, b1, c0_0, c0_1, .., c1_1, d
    1.0, 1.5, 2.7, 3.1, 3.4, .., 0.89,  0.0
    0.9, 8.7, 1.9, 3.4, 7.8, .., 1.21,  0.0
    
    The first line is a header, the remaining lines give the values.

    """
    # find a function space:
    fs = None
    for n,d in sorted(data.items(), key=lambda x: x[0]):
        if isinstance(d, Data): fs=d.getFunctionSpace()
    if fs is None:
        raise ValueError("saveDataCSV: there must be at least one Data object in the argument list.")
    
    new_data={}
    for n,d in sorted(data.items(), key=lambda x: x[0]):
        if isinstance(d, Data):
            new_data[n]=d
        else:
            try:
               new_data[n]=Data(d,fs)
            except:
               raise ValueError("saveDataCSV: unknown non-data argument type for %s"%(str(n)))
    escore._saveDataCSV(filename, new_data, sep, csep, refid, append)


def getNumpy(**data):
    """
    Converts `Data` objects to numpy arrays.

    The keyword args are Data objects to convert.
    If a scalar `Data` object is passed with the name ``mask``, then only
    samples which correspond to positive values in ``mask`` will be output.

    :keyword <name>: Data object or scalar to convert
    :type <name>: `Data` or ``float``
    :return: dictionary mapping names to numpy record arrays
    :rtype: ``dict`` of ``numpy.ndarray``
    :raise ValueError: if no Data object is provided

    Example usage::

        s=Scalar(..)
        v=Vector(..)
        t=Tensor(..)
        f=float()
        array = getNumpy(a=s, b=v, c=t, d=f)
    """
    # find a function space:
    fs = None
    for n,d in sorted(data.items(), key=lambda x: x[0]):
        if isinstance(d, Data): fs=d.getFunctionSpace()
    if fs is None:
        raise ValueError("getNumpy: there must be at least one Data object in the argument list.")
    
    new_data={}
    for n,d in sorted(data.items(), key=lambda x: x[0]):
        if isinstance(d, Data):
            new_data[n]=d
        else:
            try:
                new_data[n]=Data(d,fs)
            except:
                raise ValueError("getNumpy: unknown non-data argument type for %s"%(str(n)))

    numpy_arrays = escore._getNumpy(new_data)

    # Get some information about the shape of the data
    counter = -1
    shape=[]
    names=[]
    index=[]
    for i in range(0, len(numpy_arrays)):
        tmp = numpy_arrays[i]
        # New variable?
        if isinstance(tmp, str) is True:
            name = tmp
            counter+=1
            shape.append(0)
            names.append(name)
            index.append(i+1)
        else:
            shape[counter]+=1

    # Dictionary to store the information in
    answer = {}
    for i in range(0, len(names)):
        if(shape[i] == 1):
            answer.update({names[i] : numpy.rec.fromarrays(numpy_arrays[index[i]],names=names[i])})
        else:
            tmp=numpy.asarray(numpy.rec.fromarrays(numpy_arrays[index[i]],names=names[i]))
            for j in range(1, shape[i]):
                tmp=numpy.hstack([tmp,numpy.rec.fromarrays(numpy_arrays[index[i]]+j,names=names[i])])                
            answer.update({names[i] : tmp})

    return answer

def convertToNumpy(data):
    """
    Converts a `Data` object to a numpy array.

    :param data: the Data object to convert
    :type data: `Data`
    :return: numpy array containing the data values
    :rtype: ``numpy.ndarray``
    """
    if hasFeature("boostnumpy"):
        return escore._convertToNumpy(Data(data,data.getFunctionSpace()))
    else:
        N=data.getNumberOfDataPoints()
        if data.getShape() == ():
            S = ( 1 , N )
        else:
            S = data.getShape() + (N,)
        if data.isComplex():
            out = numpy.empty(S, dtype = numpy.complex128)
        else:
            out = numpy.empty( S, dtype = numpy.float64)

        for k in range(data.getNumberOfDataPoints()):
            out[:,k] = numpy.array( data.getTupleForDataPoint(k), dtype = out.dtype, subok = False, ndmin = 1 )
        return out


def NumpyToData(array, isComplex, functionspace):
    """
    Creates a `Data` object from a numpy ndarray.

    :param array: the numpy array to convert
    :type array: ``numpy.ndarray``
    :param isComplex: whether the data contains complex values
    :type isComplex: ``bool``
    :param functionspace: the function space for the new Data object
    :type functionspace: `FunctionSpace`
    :return: Data object created from the numpy array
    :rtype: `Data`
    :raise ValueError: if arguments are invalid or boost numpy is not available

    Example usage::

        NewDataObject = NumpyToData(ndarray, isComplex, FunctionSpace)
    """
    if hasFeature("boostnumpy"):
      if(not isinstance(array, (numpy.ndarray, numpy.generic))):
        raise ValueError("NumpyToData: Invalid argument for array.")
      if(not isinstance(functionspace, escore.FunctionSpace)):
        raise ValueError("NumpyToData: Invalid argument for functionspace.")
      # escore.internal_initBoostNumpy();
      return escore._numpyToData(array, isComplex, functionspace)
    else:
      raise ValueError("NumpyToData: Please recompile escript with boost numpy.")
      

def saveESD(datasetName, dataDir=".", domain=None, timeStep=0, deltaT=1, dynamicMesh=0, timeStepFormat="%04d", **data):
    """
    Saves `Data` objects to files and creates an `escript dataset` (ESD) file
    for convenient processing/visualisation.

    Single timestep example::

        tmp = Scalar(..)
        v = Vector(..)
        saveESD("solution", "data", temperature=tmp, velocity=v)

    Time series example::

        while t < t_end:
            tmp = Scalar(..)
            v = Vector(..)
            # save every 10 timesteps
            if t % 10 == 0:
                saveESD("solution", "data", timeStep=t, deltaT=10, temperature=tmp, velocity=v)
            t = t + 1

    tmp, v and the domain are saved in native format in the "data"
    directory and the file "solution.esd" is created that refers to tmp by
    the name "temperature" and to v by the name "velocity".

    :param datasetName: name of the dataset, used to name the ESD file
    :type datasetName: ``str``
    :param dataDir: optional directory where the data files should be saved
    :type dataDir: ``str``
    :param domain: domain of the `Data` object(s). If not specified, the
                   domain of the given `Data` objects is used.
    :type domain: `escript.Domain`
    :param timeStep: current timestep or sequence number - first one must be 0
    :type timeStep: `int`
    :param deltaT: timestep or sequence increment, see example above
    :type deltaT: `int`
    :param dynamicMesh: by default the mesh is assumed to be static and thus
                        only saved once at timestep 0 to save disk space.
                        Setting this to 1 changes the behaviour and the mesh
                        is saved at each timestep.
    :type dynamicMesh: `int`
    :param timeStepFormat: timestep format string (defaults to "%04d")
    :type timeStepFormat: ``str``
    :keyword <name>: writes the assigned value to the file using <name> as
                     identifier
    :type <name>: `Data` object
    :note: The ESD concept is experimental and the file format likely to
           change so use this function with caution.
    :note: The data objects have to be defined on the same domain (but not
           necessarily on the same `FunctionSpace`).
    :note: When saving a time series the first timestep must be 0 and it is
           assumed that data from all timesteps share the domain. The dataset
           file is updated in each iteration.
    """
    new_data = {}
    for n,d in sorted(data.items(), key=lambda x: x[0]):
          if not d.isEmpty(): 
            fs = d.getFunctionSpace() 
            domain2 = fs.getDomain()
            if fs == escore.Solution(domain2):
               new_data[n]=interpolate(d,escore.ContinuousFunction(domain2))
            elif fs == escore.ReducedSolution(domain2):
               new_data[n]=interpolate(d,escore.ReducedContinuousFunction(domain2))
            else:
               new_data[n]=d
            if domain is None: domain=domain2
    if domain is None:
        raise ValueError("saveESD: no domain detected.")

    if domain.isRootRank() and not os.path.isdir(dataDir):
        os.mkdir(dataDir)

    meshFile = datasetName+"_mesh"
    fileNumber = timeStep / deltaT

    if dynamicMesh == 0:
        # later timesteps reuse mesh from t=0
        if timeStep == 0:
            domain.dump(os.path.join(dataDir, meshFile))
    else:
        meshFile += ".%s"%timeStepFormat
        domain.dump(os.path.join(dataDir, (meshFile) % fileNumber))

    outputString = ""

    if domain.isRootRank():
        outputString += "#escript datafile V1.0\n"
        # number of timesteps (currently only 1 is supported)
        outputString += "T=%d\n" % (fileNumber+1)
        # timestep increment
        outputString += "DT=1\n"
        # name of the mesh file
        outputString += "M=%s\n" % meshFile
        # number of blocks (MPI size)
        outputString += "N=%d\n" % domain.getMPISize()

    # now add the variables
    for varName,d in sorted(new_data.items(), key=lambda x: x[0]):
        varFile = datasetName+"_"+varName+".%s"%timeStepFormat
        d.dump(os.path.join(dataDir, (varFile ) % fileNumber))
        if domain.isRootRank():
            outputString += "V=%s:%s\n" % (varFile, varName)

    if domain.isRootRank():
        esdfile = open(os.path.join(dataDir, datasetName+".esd"), "w")
        esdfile.write(outputString)
        esdfile.close()

def kronecker(d=3):
   """
   Returns the kronecker delta-symbol.

   :param d: dimension or an object that has the ``getDim`` method defining the
             dimension
   :type d: ``int``, `escript.Domain` or `escript.FunctionSpace`
   :return: the object u of rank 2 with *u[i,j]=1* for *i=j* and *u[i,j]=0*
            otherwise
   :rtype: ``numpy.ndarray`` or `escript.Data` of rank 2
   """
   return identityTensor(d)

def identity(shape=()):
   """
   Returns the ``shape`` x ``shape`` identity tensor.

   :param shape: input shape for the identity tensor
   :type shape: ``tuple`` of ``int``
   :return: array whose shape is shape x shape where *u[i,k]=1* for *i=k* and
            *u[i,k]=0* otherwise for len(shape)=1. If len(shape)=2:
            *u[i,j,k,l]=1* for *i=k and j=l* and *u[i,j,k,l]=0* otherwise.
   :rtype: ``numpy.ndarray`` of rank 1, rank 2 or rank 4
   :raise ValueError: if len(shape)>2
   """
   if len(shape)>0:
      out=numpy.zeros(shape+shape,numpy.float64)
      if len(shape)==1:
          for i0 in range(shape[0]):
             out[i0,i0]=1.
      elif len(shape)==2:
          for i0 in range(shape[0]):
             for i1 in range(shape[1]):
                out[i0,i1,i0,i1]=1.
      else:
          raise ValueError("identity: length of shape is restricted to 2.")
   else:
      out=1.
   return out

def zeros(shape=()):
   """
   Returns the ``shape`` zero tensor.

   :param shape: input shape for the identity tensor
   :type shape: ``tuple`` of ``int``
   :return: array of shape filled with zeros
   :rtype: ``numpy.ndarray``
   """
   if len(shape)>0:
      out=numpy.zeros(shape,numpy.float64)
   else:
      out=0.
   return out

def identityTensor(d=3):
   """
   Returns the ``d`` x ``d`` identity matrix.

   :param d: dimension or an object that has the ``getDim`` method defining the
             dimension
   :type d: ``int``, `escript.Domain` or `escript.FunctionSpace`
   :return: the object u of rank 2 with *u[i,j]=1* for *i=j* and *u[i,j]=0*
            otherwise
   :rtype: ``numpy.ndarray`` or `escript.Data` of rank 2
   """
   if isinstance(d,escore.FunctionSpace):
       return escore.Data(identity((d.getDim(),)),d)
   elif isinstance(d,escore.Domain):
       return identity((d.getDim(),))
   else:
       return identity((d,))

def identityTensor4(d=3):
   """
   Returns the ``d`` x ``d`` x ``d`` x ``d`` identity tensor.

   :param d: dimension or an object that has the ``getDim`` method defining the
             dimension
   :type d: ``int`` or any object with a ``getDim`` method
   :return: the object u of rank 4 with *u[i,j,k,l]=1* for *i=k and j=l* and
            *u[i,j,k,l]=0* otherwise
   :rtype: ``numpy.ndarray`` or `escript.Data` of rank 4
   """
   if isinstance(d,escore.FunctionSpace):
       return escore.Data(identity((d.getDim(),d.getDim())),d)
   elif isinstance(d,escore.Domain):
       return identity((d.getDim(),d.getDim()))
   else:
       return identity((d,d))

def unitVector(i=0,d=3):
   """
   Returns a unit vector u of dimension d whose non-zero element is at index i.

   :param i: index for non-zero element
   :type i: ``int``
   :param d: dimension or an object that has the ``getDim`` method defining the
             dimension
   :type d: ``int``, `escript.Domain` or `escript.FunctionSpace`
   :return: the object u of rank 1 with *u[j]=1* for *j=index* and *u[j]=0*
            otherwise
   :rtype: ``numpy.ndarray`` or `escript.Data` of rank 1
   """
   return kronecker(d)[i]

#=========================================================================
#   global reduction operations
#=========================================================================
def Lsup(arg):
    """
    Returns the Lsup-norm of argument ``arg``. This is the maximum absolute value
    over all data points. This function is equivalent to ``sup(abs(arg))``.

    :param arg: argument
    :type arg: ``float``, ``int``, `escript.Data`, ``numpy.ndarray``
    :return: maximum value of the absolute value of ``arg`` over all components
             and all data points
    :rtype: ``float``
    :raise TypeError: if type of ``arg`` cannot be processed
    """
    if isinstance(arg,numpy.ndarray):
        return sup(abs(arg))
    elif isinstance(arg,escore.Data):
        return arg._Lsup()
    elif isinstance(arg,float) or isinstance(arg, complex):
        return abs(arg)
    elif isinstance(arg,int):
        return abs(float(arg))
    else:
        raise TypeError("Lsup: Unknown argument type ("+str(type(arg))+").")

def sup(arg):
    """
    Returns the maximum value over all data points.

    :param arg: argument
    :type arg: ``float``, ``int``, `escript.Data`, ``numpy.ndarray``
    :return: maximum value of ``arg`` over all components and all data points
    :rtype: ``float``
    :raise TypeError: if type of ``arg`` cannot be processed
    """
    if isinstance(arg,numpy.ndarray):
        if arg.dtype.kind=='c':
            raise TypeError("sup: operation not supported for complex");         
        return arg.max()
    elif isinstance(arg,escore.Data):
        return arg._sup()
    elif isinstance(arg,float):
        return arg
    elif isinstance(arg,int):
        return float(arg)
    elif isinstance(arg,complex):
        raise TypeError("sup:  Operation not supported for complex.")
    else:
        raise TypeError("sup: Unknown argument type.")

def inf(arg):
    """
    Returns the minimum value over all data points.

    :param arg: argument
    :type arg: ``float``, ``int``, `escript.Data`, ``numpy.ndarray``
    :return: minimum value of ``arg`` over all components and all data points
    :rtype: ``float``
    :raise TypeError: if type of ``arg`` cannot be processed
    """
    if isinstance(arg,numpy.ndarray):
        if arg.dtype.kind=='c':
            raise TypeError("inf: operation not supported for complex");         
        return arg.min()
    elif isinstance(arg,escore.Data):
        return arg._inf()
    elif isinstance(arg,float):
        return arg
    elif isinstance(arg,int):
        return float(arg)
    elif isinstance(arg,complex):
        raise TypeError("inf:  Operation not supported for complex.")
    else:
      raise TypeError("inf: Unknown argument type.")

#=========================================================================
#   some little helpers
#=========================================================================
def getRank(arg):
    """
    Identifies the rank of the argument.

    :param arg: an object whose rank is to be returned
    :type arg: ``numpy.ndarray``, `escript.Data`, ``float``, ``int``,
               ``Symbol``
    :return: the rank of the argument
    :rtype: ``int``
    :raise TypeError: if type of ``arg`` cannot be processed
    """
    if isinstance(arg,list) or isinstance(arg,tuple):
        return numpy.array(arg).ndim
    elif isinstance(arg,numpy.ndarray):
        return arg.ndim
    elif isinstance(arg,escore.Data):
        return arg.getRank()
    elif isinstance(arg,float) or isinstance(arg,int) or isinstance(arg,complex):
        return 0
    elif isinstance(arg,sym.Symbol):
        return arg.getRank()
    else:
      raise TypeError("getRank: Unknown argument type.")

def getShape(arg):
    """
    Identifies the shape of the argument.

    :param arg: an object whose shape is to be returned
    :type arg: ``numpy.ndarray``, `escript.Data`, ``float``, ``int``,
               ``Symbol``
    :return: the shape of the argument
    :rtype: ``tuple`` of ``int``
    :raise TypeError: if type of ``arg`` cannot be processed
    """

    if isinstance(arg,numpy.ndarray):
        return arg.shape
    elif isinstance(arg,list) or isinstance(arg,tuple):
        return numpy.array(arg).shape
    elif isinstance(arg,escore.Data):
        return arg.getShape()
    elif isinstance(arg,float) or isinstance(arg,int) or isinstance(arg,complex):
        return ()
    elif isinstance(arg,sym.Symbol):
        return arg.getShape()
    else:
      raise TypeError("getShape: Cannot identify shape")

def pokeDim(arg):
    """
    Identifies the spatial dimension of the argument.

    :param arg: an object whose spatial dimension is to be returned
    :type arg: any
    :return: the spatial dimension of the argument, if available, or ``None``
    :rtype: ``int`` or ``None``
    """

    if isinstance(arg,escore.Data):
        return arg.getFunctionSpace().getDim()
    else:
        return None

def commonShape(arg0, arg1):
    """
    Returns a shape to which ``arg0`` can be extended from the right and ``arg1``
    can be extended from the left.

    :param arg0: an object with a shape (see `getShape`)
    :param arg1: an object with a shape (see `getShape`)
    :return: the shape of ``arg0`` or ``arg1`` such that the left part equals the
             shape of ``arg0`` and the right end equals the shape of ``arg1``
    :rtype: ``tuple`` of ``int``
    :raise ValueError: if no shape can be found
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)<len(sh1):
       if not sh0==sh1[:len(sh0)]:
             raise ValueError("argument 0 cannot be extended to the shape of argument 1")
       return sh1
    elif len(sh0)>len(sh1):
       if not sh1==sh0[:len(sh1)]:
             raise ValueError("argument 1 cannot be extended to the shape of argument 0")
       return sh0
    else:
       if not sh0==sh1:
             raise ValueError("argument 1 and argument 0 have not the same shape.")
       return sh0

def commonDim(*args):
    """
    Identifies, if possible, the spatial dimension across a set of objects
    which may or may not have a spatial dimension.

    :param args: given objects
    :return: the spatial dimension of the objects with identifiable dimension
             (see `pokeDim`). If none of the objects has a spatial dimension
             ``None`` is returned.
    :rtype: ``int`` or ``None``
    :raise ValueError: if the objects with identifiable dimension don't have
                       the same spatial dimension.
    """
    out=None
    for a in args:
       d=pokeDim(a)
       if not out is None:
          if not (d is None or out==d):
             raise ValueError("dimension of arguments don't match")
       else:
          out=d
    return out

def testForZero(arg):
    """
    Tests if the argument is identical to zero.

    :param arg: the object to test for zero
    :type arg: typically ``numpy.ndarray``, `escript.Data`, ``float``, ``int``
    :return: True if the argument is identical to zero, False otherwise
    :rtype: ``bool``
    """
    if isinstance(arg,numpy.ndarray):
       return not Lsup(arg)>0.
    elif isinstance(arg,escore.Data):
       return False
    elif isinstance(arg,float) or isinstance(arg,int) or isinstance(arg,complex):
       return not Lsup(arg)>0.
    else:
       return False

def matchType(arg0=0.,arg1=0.):
    """
    Converts ``arg0`` and ``arg1`` both to the same type ``numpy.ndarray`` or
    `escript.Data`

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``,`escript.Data`,``float``, ``int``, ``Symbol``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``,`escript.Data`,``float``, ``int``, ``Symbol``
    :return: a tuple representing ``arg0`` and ``arg1`` with the same type or
             with at least one of them being a `Symbol`
    :rtype: ``tuple`` of two ``numpy.ndarray`` or two `escript.Data`
    :raise TypeError: if type of ``arg0`` or ``arg1`` cannot be processed
    """
    if isinstance(arg0,numpy.ndarray):
       if isinstance(arg1,numpy.ndarray):
          pass
       elif isinstance(arg1,escore.Data):
          arg0=escore.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg1=numpy.array(arg1,dtype=numpy.double)
       elif isinstance(arg1,int):
          arg1=numpy.array(float(arg1),dtype=numpy.double)
       elif isinstance(arg1,complex):
          arg1=numpy.array(arg1, dtype=numpy.cdouble)
       elif isinstance(arg1,sym.Symbol):
          pass
       else:
          raise TypeError("function: Unknown type of second argument.")
    elif isinstance(arg0,escore.Data):
       if isinstance(arg1,numpy.ndarray):
          arg1=escore.Data(arg1,arg0.getFunctionSpace())
       elif isinstance(arg1,escore.Data):
          pass
       elif isinstance(arg1,float) or isinstance(arg1,complex):
          arg1=escore.Data(arg1,(),arg0.getFunctionSpace())
       elif isinstance(arg1,int):
          arg1=escore.Data(float(arg1),(),arg0.getFunctionSpace())
       elif isinstance(arg1,sym.Symbol):
          pass
       else:
          raise TypeError("function: Unknown type of second argument.")
    elif isinstance(arg0,sym.Symbol):
       if isinstance(arg1,numpy.ndarray):
          pass
       elif isinstance(arg1,escore.Data):
          pass
       elif isinstance(arg1,complex):  
          pass
       elif isinstance(arg1,float):
          pass
       elif isinstance(arg1,int):
          pass
       elif isinstance(arg1,sym.Symbol):
          pass
       else:
          raise TypeError("function: Unknown type of second argument.")
    elif isinstance(arg0,complex):
       if isinstance(arg1,numpy.ndarray):
          arg0=numpy.array(arg0,dtype=numpy.cdouble)
       elif isinstance(arg1,escore.Data):
          arg0=escore.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numpy.array(arg0,dtype=numpy.cdouble)
          arg1=numpy.array(arg1,dtype=numpy.cdouble)
       elif isinstance(arg1,int):
          arg0=numpy.array(arg0,dtype=numpy.cdouble)
          arg1=numpy.array(float(arg1),dtype=numpy.cdouble)
       elif isinstance(arg1,sym.Symbol):
          pass
       elif isinstance(arg1,complex):
          pass
       else:
          raise TypeError("function: Unknown type of second argument.") 
    elif isinstance(arg0,float):
       if isinstance(arg1,numpy.ndarray):
          arg0=numpy.array(arg0,dtype=numpy.double)
       elif isinstance(arg1,escore.Data):
          arg0=escore.Data(arg0,arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numpy.array(arg0,dtype=numpy.double)
          arg1=numpy.array(arg1,dtype=numpy.double)
       elif isinstance(arg1,int):
          arg0=numpy.array(arg0,dtype=numpy.double)
          arg1=numpy.array(float(arg1),dtype=numpy.double)
       elif isinstance(arg1,complex):
          arg0=numpy.array(complex(arg0),dtype=numpy.cdouble)
          arg1=numpy.array(complex(arg1),dtype=numpy.cdouble)    
       elif isinstance(arg1,sym.Symbol):
          pass
       else:
          raise TypeError("function: Unknown type of second argument.")
    elif isinstance(arg0,int):
       if isinstance(arg1,numpy.ndarray):
          arg0=numpy.array(float(arg0),dtype=numpy.double)
       elif isinstance(arg1,escore.Data):
          arg0=escore.Data(float(arg0),arg1.getFunctionSpace())
       elif isinstance(arg1,float):
          arg0=numpy.array(float(arg0),dtype=numpy.double)
          arg1=numpy.array(arg1,dtype=numpy.double)
       elif isinstance(arg1,complex):
          arg0=numpy.array(complex(arg0),dtype=numpy.cdouble)
          arg1=numpy.array(complex(arg1),dtype=numpy.cdouble)          
       elif isinstance(arg1,int):
          arg0=numpy.array(float(arg0),dtype=numpy.double)
          arg1=numpy.array(float(arg1),dtype=numpy.double)
       elif isinstance(arg1,sym.Symbol):
          pass
       else:
          raise TypeError("function: Unknown type of second argument.")
    else:
      raise TypeError("function: Unknown type of first argument.")

    return arg0,arg1

def matchShape(arg0,arg1):
    """
    Returns a representation of ``arg0`` and ``arg1`` which have the same shape.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``,`escript.Data`,``float``, ``int``, `Symbol`
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``,`escript.Data`,``float``, ``int``, `Symbol`
    :return: ``arg0`` and ``arg1`` where copies are returned when the shape has
             to be changed
    :rtype: ``tuple``
    """
    sh=commonShape(arg0,arg1)
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)<len(sh):
       return outer(arg0,numpy.ones(sh[len(sh0):],numpy.double)),arg1
    elif len(sh1)<len(sh):
       return arg0,outer(arg1,numpy.ones(sh[len(sh1):],numpy.double))
    else:
       return arg0,arg1

def log10(arg):
   """
   Returns base-10 logarithm of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.log10(arg)
   elif isinstance(arg,escore.Data):
      return arg._log10()
   elif isinstance(arg,complex):
      return cmath.log10(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.log10(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.log10)
   else:
      raise TypeError("log10: Unknown argument type.")

def wherePositive(arg):
   """
   Returns mask of positive values of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("wherePositive: operation not supported for complex");
      out=numpy.greater(arg,numpy.zeros(arg.shape,numpy.double))*1.
      return out
   elif isinstance(arg,escore.Data):
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
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.wherePositive)
   elif isinstance(arg,complex):
      raise TypeError("wherePositive: operation not supported for complex");
   else:
      raise TypeError("wherePositive: Unknown argument type.")

def whereNegative(arg):
   """
   Returns mask of negative values of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("whereNegative: operation not supported for complex");
      out=numpy.less(arg,numpy.zeros(arg.shape,numpy.double))*1.
      return out
   elif isinstance(arg,escore.Data):
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
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.whereNegative)
   elif isinstance(arg,complex):
      raise TypeError("whereNegative: operation not supported for complex");
   else:
      raise TypeError("whereNegative: Unknown argument type.")

def whereNonNegative(arg):
   """
   Returns mask of non-negative values of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("whereNonNegative: operation not supported for complex");       
      out=numpy.greater_equal(arg,numpy.zeros(arg.shape,numpy.double))*1.
      return out
   elif isinstance(arg,escore.Data):
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
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.whereNonNegative)
   elif isinstance(arg,complex):
      raise TypeError("whereNonNegative: operation not supported for complex");
   else:
      raise TypeError("whereNonNegative: Unknown argument type.")

def whereNonPositive(arg):
   """
   Returns mask of non-positive values of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("whereNonPositive: operation not supported for complex");       
      out=numpy.less_equal(arg,numpy.zeros(arg.shape,numpy.double))*1.
      return out
   elif isinstance(arg,escore.Data):
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
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.whereNonPositive)
   elif isinstance(arg,complex):
      raise TypeError("whereNonPositive: operation not supported for complex");
   else:
      raise TypeError("whereNonPositive: Unknown argument type.")

def whereZero(arg,tol=None,rtol=math.sqrt(EPSILON)):
   """
   Returns mask of zero entries of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :param tol: absolute tolerance. Values with absolute value less than tol are accepted
               as zero. If ``tol`` is not present ``rtol``*```Lsup` (arg)`` is used. 
   :type tol: ``float``
   :param rtol: relative tolerance used to define the absolute tolerance if ``tol`` is not present.
   :type rtol: non-negative ``float``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise ValueError: if ``rtol`` is non-negative.
   :raise TypeError: if the type of the argument is not expected
   """
   if tol is None and not isinstance(arg, sym.Symbol):
      if rtol<0: raise ValueError("rtol must be non-negative.")
      tol = Lsup(arg)*rtol
   if isinstance(arg,numpy.ndarray):
      out=numpy.less_equal(abs(arg)-tol,numpy.zeros(arg.shape,numpy.double))*1.
      if isinstance(out,float): out=numpy.array(out,dtype=numpy.double)
      return out
   elif isinstance(arg,escore.Data):
      return arg._whereZero(tol)
   elif isinstance(arg,float) or isinstance(arg,complex) or isinstance(arg, int):
      if abs(arg)<=tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.whereZero)
   else:
      raise TypeError("whereZero: Unknown argument type.")

def whereNonZero(arg,tol=0.):
   """
   Returns mask of values different from zero of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :param tol: absolute tolerance. Values with absolute value less than tol are accepted
               as zero. If ``tol`` is not present ``rtol``*```Lsup` (arg)`` is used. 
   :type tol: ``float``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise ValueError: if ``rtol`` is non-negative.
   :raise TypeError: if the type of the argument is not expected
   """
   if tol is None:
      if rtol<=0: raise ValueError("rtol must be non-negative.")
      tol = Lsup(arg)*rtol
   if isinstance(arg,numpy.ndarray):
      out=numpy.greater(abs(arg)-tol,numpy.zeros(arg.shape,numpy.double))*1.
      if isinstance(out,float): out=numpy.array(out,dtype=numpy.double)
      return out
   elif isinstance(arg,escore.Data):
      return arg._whereNonZero(tol)
   elif isinstance(arg,float) or isinstance(arg,complex) or isinstance(arg, int):
      if abs(arg)>tol:
        return 1.
      else:
        return 0.
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.whereNonZero)
   else:
      raise TypeError("whereNonZero: Unknown argument type.")

def Abs(arg):
   """
   Returns the absolute value of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.sym.symfn.abs)
   else:
      return abs(arg)

def erf(arg):
   """
   Returns the error function *erf* of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,escore.Data):
      return arg._erf()
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.erf)
   else:
      raise TypeError("erf: Unknown argument type.")

def sin(arg):
   """
   Returns sine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.sin(arg)
   elif isinstance(arg,escore.Data):
      return arg._sin()
   elif isinstance(arg,complex):
      return cmath.sin(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.sin(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.sin)
   else:
      raise TypeError("sin: Unknown argument type.")

def cos(arg):
   """
   Returns cosine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.cos(arg)
   elif isinstance(arg,escore.Data):
      return arg._cos()
   elif isinstance(arg,complex):
      return cmath.cos(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.cos(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.cos)
   else:
      raise TypeError("cos: Unknown argument type.")

def tan(arg):
   """
   Returns tangent of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.tan(arg)
   elif isinstance(arg,escore.Data):
      return arg._tan()
   elif isinstance(arg,complex):
      return cmath.tan(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.tan(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.tan)
   else:
      raise TypeError("tan: Unknown argument type.")

def asin(arg):
   """
   Returns the inverse sine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.arcsin(arg)
   elif isinstance(arg,escore.Data):
      return arg._asin()
   elif isinstance(arg,complex):
      return cmath.asin(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.asin(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.asin)
   else:
      raise TypeError("asin: Unknown argument type.")

def acos(arg):
   """
   Returns the inverse cosine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.arccos(arg)
   elif isinstance(arg,escore.Data):
      return arg._acos()
   elif isinstance(arg,complex):
      return cmath.acos(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.acos(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.acos)
   else:
      raise TypeError("acos: Unknown argument type.")

def atan(arg):
   """
   Returns inverse tangent of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.arctan(arg)
   elif isinstance(arg,escore.Data):
      return arg._atan()
   elif isinstance(arg,complex):
      return cmath.atan(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.atan(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.atan)
   else:
      raise TypeError("atan: Unknown argument type.")

def atan2(arg0, arg1):
   """
   Returns inverse tangent of argument ``arg0`` over ``arg1`` 
   """
   m=whereZero(arg1, rtol=EPSILON)
   m2=whereNegative(arg1*arg0)
   # return atan(arg0/(arg1+m))*(1-m)+(numpy.pi/2)*(1-2*m2)*m
   s=atan(arg0/(arg1+m))(1-m)+(np.pi/2)(1-2*m2)*m
   s+=(wherePositive(arg0)*whereNegative(arg1)-whereNegative(arg0)*whereNegative(arg1))*numpy.pi
   return s

      
def sinh(arg):
   """
   Returns the hyperbolic sine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.sinh(arg)
   elif isinstance(arg,escore.Data):
      return arg._sinh()
   elif isinstance(arg,complex):
      return cmath.sinh(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.sinh(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.sinh)
   else:
      raise TypeError("sinh: Unknown argument type.")

def cosh(arg):
   """
   Returns the hyperbolic cosine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.cosh(arg)
   elif isinstance(arg,escore.Data):
      return arg._cosh()
   elif isinstance(arg,complex):
      return cmath.cosh(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.cosh(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.cosh)
   else:
      raise TypeError("cosh: Unknown argument type.")

def tanh(arg):
   """
   Returns the hyperbolic tangent of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.tanh(arg)
   elif isinstance(arg,escore.Data):
      return arg._tanh()
   elif isinstance(arg,complex):
      return cmath.tanh(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.tanh(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.tanh)
   else:
      raise TypeError("tanh: Unknown argument type.")

def asinh(arg):
   """
   Returns the inverse hyperbolic sine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.arcsinh(arg)
   elif isinstance(arg,escore.Data):
      return arg._asinh()
   elif isinstance(arg,complex):
      return numpy.arcsinh(complex(arg))
   elif isinstance(arg,float) or isinstance(arg,int):
      return numpy.arcsinh(float(arg))
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.asinh)
   else:
      raise TypeError("asinh: Unknown argument type.")

def acosh(arg):
   """
   Returns the inverse hyperbolic cosine of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.arccosh(arg)
   elif isinstance(arg,escore.Data):
      return arg._acosh()
   elif isinstance(arg,complex):
      return numpy.arccosh(complex(arg))
   elif isinstance(arg,float) or isinstance(arg,int):
      return numpy.arccosh(float(arg))
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.acosh)
   else:
      raise TypeError("acosh: Unknown argument type.")

def atanh(arg):
   """
   Returns the inverse hyperbolic tangent of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.arctanh(arg)
   elif isinstance(arg,escore.Data):
      return arg._atanh()
   elif isinstance(arg,complex):
      return numpy.arctanh(complex(arg))
   elif isinstance(arg,float) or isinstance(arg,int):
      return numpy.arctanh(float(arg))
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.atanh)
   else:
      raise TypeError("atanh: Unknown argument type.")

def exp(arg):
   """
   Returns *e* to the power of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of arg
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.exp(arg)
   elif isinstance(arg,escore.Data):
      return arg._exp()
   elif isinstance(arg,complex):
      return cmath.exp(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.exp(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.exp)
   else:
      raise TypeError("exp: Unknown argument type.")

def sqrt(arg):
   """
   Returns the square root of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
           depending on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.sqrt(arg)
   elif isinstance(arg,escore.Data):
      return arg._sqrt()
   elif isinstance(arg,complex):
      return cmath.sqrt(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.sqrt(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.sqrt)
   else:
      raise TypeError("sqrt: Unknown argument type.")

def log(arg):
   """
   Returns the natural logarithm of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``.
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      return numpy.log(arg)
   elif isinstance(arg,escore.Data):
      return arg._log()
   elif isinstance(arg,complex):
      return cmath.log(arg)
   elif isinstance(arg,float) or isinstance(arg,int):
      return math.log(arg)
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.log)
   else:
      raise TypeError("log: Unknown argument type.")

def sign(arg):
   """
   Returns the sign of argument ``arg``.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray`` depending
           on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("sign: operation not supported for complex")        
      return wherePositive(arg)-whereNegative(arg)
   elif isinstance(arg,escore.Data):
      return arg._sign()
   elif isinstance(arg,complex):
       raise TypeError("sign: operation not supported for complex")
   elif isinstance(arg,float) or isinstance(arg,int):
      if arg>0:
        return 1.
      elif arg<0:
        return -1.
      else:
        return 0.
   elif isinstance(arg,sym.Symbol):
      return arg.applyfunc(sym.symfn.sign)
   else:
      raise TypeError("sign: Unknown argument type.")

def minval(arg):
   """
   Returns the minimum value over all components of ``arg`` at each data point.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol` depending on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("minval: operation not supported for complex");       
      if arg.ndim==0:
         return float(arg)
      else:
         return arg.min()
   elif isinstance(arg,escore.Data):
      return arg._minval()
   elif isinstance(arg,float):
      return arg
   elif isinstance(arg,int):
      return float(arg)
   elif isinstance(arg,sym.Symbol):
      return sym.symfn.minval(arg)
   else:
      raise TypeError("minval: Unknown argument type.")

def maxval(arg):
   """
   Returns the maximum value over all components of ``arg`` at each data point.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol` depending on the type of ``arg``
   :raise TypeError: if the type of the argument is not expected
   """
   if isinstance(arg,numpy.ndarray):
      if arg.dtype.kind=='c':
          raise TypeError("maxval: operation not supported for complex");
      if arg.ndim==0:
         return float(arg)
      else:
         return arg.max()
   elif isinstance(arg,escore.Data):
      return arg._maxval()
   elif isinstance(arg,float):
      return arg
   elif isinstance(arg,int):
      return float(arg)
   elif isinstance(arg,sym.Symbol):
      return sym.symfn.maxval(arg)
   else:
      raise TypeError("maxval: Unknown argument type.")

def length(arg):
   """
   Returns the length (Euclidean norm) of argument ``arg`` at each data point.

   :param arg: argument
   :type arg: ``float``, `escript.Data`, `Symbol`, ``numpy.ndarray``
   :rtype: ``float``, `escript.Data`, `Symbol` depending on the type of ``arg``
   """
   a=abs(arg)
   return sqrt(inner(a,a))

def trace(arg,axis_offset=0):
   """
   Returns the trace of ``arg`` which is the sum of ``arg[k,k]`` over k.

   :param arg: argument
   :type arg: `escript.Data`, `Symbol`, ``numpy.ndarray``
   :param axis_offset: ``axis_offset`` to components to sum over. ``axis_offset``
                       must be non-negative and less than the rank of ``arg`` +1.
                       The dimensions of component ``axis_offset`` and
                       axis_offset+1 must be equal.
   :type axis_offset: ``int``
   :return: trace of arg. The rank of the returned object is rank of ``arg``
            minus 2.
   :rtype: `escript.Data`, `Symbol` or ``numpy.ndarray`` depending on the type
           of ``arg``
   """
   if isinstance(arg,numpy.ndarray):
      sh=arg.shape
      if len(sh)<2:
        raise ValueError("rank of argument must be greater than 1")
      if axis_offset<0 or axis_offset>len(sh)-2:
        raise ValueError("axis_offset must be between 0 and %d"%(len(sh)-2))
      s1=1
      for i in range(axis_offset): s1*=sh[i]
      s2=1
      for i in range(axis_offset+2,len(sh)): s2*=sh[i]
      if not sh[axis_offset] == sh[axis_offset+1]:
        raise ValueError("dimensions of component %d and %d must match."%(axis_offset,axis_offset+1))
      arg_reshaped=numpy.reshape(arg,(s1,sh[axis_offset],sh[axis_offset],s2))
      out=numpy.zeros([s1,s2],numpy.double)
      for i1 in range(s1):
        for i2 in range(s2):
            for j in range(sh[axis_offset]): out[i1,i2]+=arg_reshaped[i1,j,j,i2]
      out.resize(sh[:axis_offset]+sh[axis_offset+2:])
      return out
   elif isinstance(arg,escore.Data):
      if arg.getRank()<2:
        raise ValueError("rank of argument must be greater than 1")
      if axis_offset<0 or axis_offset>arg.getRank()-2:
        raise ValueError("axis_offset must be between 0 and %d"%(arg.getRank()-2))
      s=list(arg.getShape())
      if not s[axis_offset] == s[axis_offset+1]:
        raise ValueError("dimensions of component %d and %d must match."%(axis_offset,axis_offset+1))
      return arg._trace(axis_offset)
   elif isinstance(arg,sym.Symbol):
      if arg.getRank()<2:
        raise ValueError("rank of argument must be greater than 1")
      if axis_offset<0 or axis_offset>arg.getRank()-2:
        raise ValueError("axis_offset must be between 0 and %d"%(arg.getRank()-2))
      s=list(arg.getShape())
      if not s[axis_offset] == s[axis_offset+1]:
        raise ValueError("dimensions of component %d and %d must match."%(axis_offset,axis_offset+1))
      return arg.trace(axis_offset)
   elif isinstance(arg,complex):
      raise TypeError("illegal argument type complex.")
   elif isinstance(arg,float):
      raise TypeError("illegal argument type float.")
   elif isinstance(arg,int):
      raise TypeError("illegal argument type int.")
   else:
      raise TypeError("Unknown argument type.")

def transpose(arg,axis_offset=None):
   """
   Returns the transpose of ``arg`` by swapping the first ``axis_offset`` and
   the last ``rank-axis_offset`` components.

   :param arg: argument
   :type arg: `escript.Data`, `Symbol`, ``numpy.ndarray``, ``float``, ``int``
   :param axis_offset: the first ``axis_offset`` components are swapped with the
                       rest. ``axis_offset`` must be non-negative and less or
                       equal to the rank of ``arg``. If ``axis_offset`` is not
                       present ``int(r/2)`` where r is the rank of ``arg`` is
                       used.
   :type axis_offset: ``int``
   :return: transpose of ``arg``
   :rtype: `escript.Data`, `Symbol`, ``numpy.ndarray``, ``float``, ``int``
           depending on the type of ``arg``
   """
   if isinstance(arg,numpy.ndarray):
      if axis_offset is None: axis_offset=int(arg.ndim/2)
      return numpy.transpose(arg,axes=list(range(axis_offset,arg.ndim))+list(range(0,axis_offset)))
   elif isinstance(arg,escore.Data):
      r=arg.getRank()
      if axis_offset is None: axis_offset=int(r/2)
      if axis_offset<0 or axis_offset>r:
        raise ValueError("axis_offset must be between 0 and %s"%r)
      return arg._transpose(axis_offset)
   elif isinstance(arg,complex):
      if not ( axis_offset==0 or axis_offset is None):
        raise ValueError("axis_offset must be 0 for complex argument")
      return arg
   elif isinstance(arg,float):
      if not ( axis_offset==0 or axis_offset is None):
        raise ValueError("axis_offset must be 0 for float argument")
      return arg
   elif isinstance(arg,int):
      if not ( axis_offset==0 or axis_offset is None):
        raise ValueError("axis_offset must be 0 for int argument")
      return float(arg)
   elif isinstance(arg,sym.Symbol):
      r=arg.getRank()
      if axis_offset is None: axis_offset=int(r/2)
      if axis_offset<0 or axis_offset>r:
        raise ValueError("axis_offset must be between 0 and %s"%r)
      return arg.transpose(axis_offset)
   else:
      raise TypeError("Unknown argument type.")

def swap_axes(arg,axis0=0,axis1=1):
   """
   Returns the swap of ``arg`` by swapping the components ``axis0`` and ``axis1``.

   :param arg: argument
   :type arg: `escript.Data`, `Symbol`, ``numpy.ndarray``
   :param axis0: first axis. ``axis0`` must be non-negative and less than the
                 rank of ``arg``.
   :type axis0: ``int``
   :param axis1: second axis. ``axis1`` must be non-negative and less than the
                 rank of ``arg``.
   :type axis1: ``int``
   :return: ``arg`` with swapped components
   :rtype: `escript.Data`, `Symbol` or ``numpy.ndarray`` depending on the type
           of ``arg``
   """
   if axis0 > axis1:
      axis0,axis1=axis1,axis0
   if isinstance(arg,numpy.ndarray):
      return numpy.swapaxes(arg,axis0,axis1)
   elif isinstance(arg,escore.Data):
      return arg._swap_axes(axis0,axis1)
   elif isinstance(arg,sym.Symbol):
      return arg.swap_axes(axis0,axis1)
   elif isinstance(arg,complex):
      raise TypeError("complex argument is not supported.")
   elif isinstance(arg,float):
      raise TypeError("float argument is not supported.")
   elif isinstance(arg,int):
      raise TypeError("int argument is not supported.")
   else:
      raise TypeError("Unknown argument type.")

def symmetric(arg):
    """
    Returns the symmetric part of the square matrix ``arg``. That is,
    *(arg+transpose(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: symmetric part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg,numpy.ndarray):
      if arg.ndim==2:
        if not (arg.shape[0]==arg.shape[1]):
           raise ValueError("argument must be square.")
      elif arg.ndim==4:
        if not (arg.shape[0]==arg.shape[2] and arg.shape[1]==arg.shape[3]):
           raise ValueError("argument must be square.")
      else:
        raise ValueError("rank 2 or 4 is required.")
      return (arg+transpose(arg))/2
    elif isinstance(arg,escore.Data):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError("argument must be square.")
        return arg._symmetric()
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError("argument must be square.")
        return arg._symmetric()
      else:
        raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank()==2:
            if arg.getShape()[0]!=arg.getShape()[1]:
                raise ValueError("symmetric: argument must be square.")
        elif arg.getRank()==4:
            if arg.getShape()[0]!=arg.getShape()[2] or arg.getShape()[1]!=arg.getShape()[3]:
                raise ValueError("symmetric: argument must be square.")
        else:
            raise ValueError("symmetric: rank 2 or 4 is required.")
        return (arg+transpose(arg))/2
    elif isinstance(arg,complex):
      return arg
    elif isinstance(arg,float):
      return arg
    elif isinstance(arg,int):
      return float(arg)
    else:
      raise TypeError("symmetric: Unknown argument type.")

def nonsymmetric(arg):
    """
    Deprecated alias for antisymmetric
    """
    return antisymmetric(arg)

def antisymmetric(arg):
    """
    Returns the anti-symmetric part of the square matrix ``arg``. That is,
    *(arg-transpose(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: anti-symmetric part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg,numpy.ndarray):
      if arg.ndim==2:
        if not (arg.shape[0]==arg.shape[1]):
           raise ValueError("antisymmetric: argument must be square.")
      elif arg.ndim==4:
        if not (arg.shape[0]==arg.shape[2] and arg.shape[1]==arg.shape[3]):
           raise ValueError("antisymmetric: argument must be square.")
      else:
        raise ValueError("antisymmetric: rank 2 or 4 is required.")
      return (arg-transpose(arg))/2
    elif isinstance(arg,escore.Data):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError("argument must be square.")
        return arg._antisymmetric()
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError("argument must be square.")
        return arg._antisymmetric()
      else:
        raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank()==2:
            if arg.getShape()[0]!=arg.getShape()[1]:
                raise ValueError("antisymmetric: argument must be square.")
        elif arg.getRank()==4:
            if arg.getShape()[0]!=arg.getShape()[2] or arg.getShape()[1]!=arg.getShape()[3]:
                raise ValueError("antisymmetric: argument must be square.")
        else:
            raise ValueError("antisymmetric: rank 2 or 4 is required.")
        return (arg-transpose(arg))/2
    elif isinstance(arg,complex):
        return complex(0)
    elif isinstance(arg,float):
        return float(0)
    elif isinstance(arg,int):
        return float(0)
    else:
        raise TypeError("antisymmetric: Unknown argument type.")

def hermitian(arg):
    """
    Returns the hermitian part of the square matrix ``arg``. That is,
    *(arg+adjoint(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: hermitian part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg,numpy.ndarray):
      if arg.ndim==2:
        if not (arg.shape[0]==arg.shape[1]):
           raise ValueError("argument must be square.")
      elif arg.ndim==4:
        if not (arg.shape[0]==arg.shape[2] and arg.shape[1]==arg.shape[3]):
           raise ValueError("argument must be square.")
      else:
        raise ValueError("rank 2 or 4 is required.")
      return (arg+transpose(arg).conj())/2
    elif isinstance(arg,escore.Data):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError("argument must be square.")
        return arg._hermitian()
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError("argument must be square.")
        return arg._hermitian()
      else:
        raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank()==2:
            if arg.getShape()[0]!=arg.getShape()[1]:
                raise ValueError("hermitian: argument must be square.")
        elif arg.getRank()==4:
            if arg.getShape()[0]!=arg.getShape()[2] or arg.getShape()[1]!=arg.getShape()[3]:
                raise ValueError("hermitian: argument must be square.")
        else:
            raise ValueError("hermitian: rank 2 or 4 is required.")
        return (arg+adjoint(arg))/2
    elif isinstance(arg,complex):
      return complex(arg.real)
    elif isinstance(arg,float):
      return arg
    elif isinstance(arg,int):
      return float(arg)
    else:
      raise TypeError("hermitian: Unknown argument type.")

def antihermitian(arg):
    """
    Returns the anti-hermitian part of the square matrix ``arg``. That is,
    *(arg-adjoint(arg))/2*.

    :param arg: input matrix. Must have rank 2 or 4 and be square.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: anti-hermitian part of ``arg``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if isinstance(arg,numpy.ndarray):
      if arg.ndim==2:
        if not (arg.shape[0]==arg.shape[1]):
           raise ValueError("antihermitian: argument must be square.")
      elif arg.ndim==4:
        if not (arg.shape[0]==arg.shape[2] and arg.shape[1]==arg.shape[3]):
           raise ValueError("antihermitian: argument must be square.")
      else:
        raise ValueError("antihermitian: rank 2 or 4 is required.")
      return (arg-transpose(arg).conj())/2
    elif isinstance(arg,escore.Data):
      if arg.getRank()==2:
        if not (arg.getShape()[0]==arg.getShape()[1]):
           raise ValueError("argument must be square.")
        return arg._antihermitian()
      elif arg.getRank()==4:
        if not (arg.getShape()[0]==arg.getShape()[2] and arg.getShape()[1]==arg.getShape()[3]):
           raise ValueError("argument must be square.")
        return arg._antihermitian()
      else:
        raise ValueError("rank 2 or 4 is required.")
    elif isinstance(arg, sym.Symbol):
        if arg.getRank()==2:
            if arg.getShape()[0]!=arg.getShape()[1]:
                raise ValueError("antihermitian: argument must be square.")
        elif arg.getRank()==4:
            if arg.getShape()[0]!=arg.getShape()[2] or arg.getShape()[1]!=arg.getShape()[3]:
                raise ValueError("antihermitian: argument must be square.")
        else:
            raise ValueError("antihermitian: rank 2 or 4 is required.")
        return (arg-hermitian(arg))/2
    elif isinstance(arg,complex):
        return complex(arg.imag*1j)
    elif isinstance(arg,float):
        return float(0)
    elif isinstance(arg,int):
        return float(0)
    else:
        raise TypeError("antihermitian: Unknown argument type.")        
        
        
def inverse(arg):
    """
    Returns the inverse of the square matrix ``arg``.

    :param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal.
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: inverse of the argument. ``matrix_mult(inverse(arg),arg)`` will be
             almost equal to ``kronecker(arg.getShape()[0])``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :note: for `escript.Data` objects the dimension is restricted to 3.
    """
    import numpy.linalg
    if isinstance(arg,numpy.ndarray):
      return numpy.linalg.tensorinv(arg,ind=1)
    elif isinstance(arg,escore.Data):
      return escript_inverse(arg)
    elif isinstance(arg,complex):
      return 1./arg
    elif isinstance(arg,float):
      return 1./arg
    elif isinstance(arg,int):
      return 1./float(arg)
    elif isinstance(arg,sym.Symbol):
      return arg.inverse()
    else:
      raise TypeError("inverse: Unknown argument type.")

def escript_inverse(arg): # this should be escript._inverse and use LAPACK
      "arg is a Data object!"
      return arg._inverse()

      #if not arg.getRank()==2:
        #raise ValueError,"escript_inverse: argument must have rank 2"
      #s=arg.getShape()
      #if not s[0] == s[1]:
        #raise ValueError,"escript_inverse: argument must be a square matrix."
      #out=escore.Data(0.,s,arg.getFunctionSpace())
      #if s[0]==1:
          #if inf(abs(arg[0,0]))==0: # in c this should be done point wise as abs(arg[0,0](i))<=0.
              #raise ZeroDivisionError,"escript_inverse: argument not invertible"
          #out[0,0]=1./arg[0,0]
      #elif s[0]==2:
          #A11=arg[0,0]
          #A12=arg[0,1]
          #A21=arg[1,0]
          #A22=arg[1,1]
          #D = A11*A22-A12*A21
          #if inf(abs(D))==0: # in c this should be done point wise as abs(D(i))<=0.
              #raise ZeroDivisionError,"escript_inverse: argument not invertible"
          #D=1./D
          #out[0,0]= A22*D
          #out[1,0]=-A21*D
          #out[0,1]=-A12*D
          #out[1,1]= A11*D
      #elif s[0]==3:
          #A11=arg[0,0]
          #A21=arg[1,0]
          #A31=arg[2,0]
          #A12=arg[0,1]
          #A22=arg[1,1]
          #A32=arg[2,1]
          #A13=arg[0,2]
          #A23=arg[1,2]
          #A33=arg[2,2]
          #D = A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22)
          #if inf(abs(D))==0: # in c this should be done point wise as abs(D(i))<=0.
              #raise ZeroDivisionError,"escript_inverse: argument not invertible"
          #D=1./D
          #out[0,0]=(A22*A33-A23*A32)*D
          #out[1,0]=(A31*A23-A21*A33)*D
          #out[2,0]=(A21*A32-A31*A22)*D
          #out[0,1]=(A13*A32-A12*A33)*D
          #out[1,1]=(A11*A33-A31*A13)*D
          #out[2,1]=(A12*A31-A11*A32)*D
          #out[0,2]=(A12*A23-A13*A22)*D
          #out[1,2]=(A13*A21-A11*A23)*D
          #out[2,2]=(A11*A22-A12*A21)*D
      #else:
         #raise TypeError,"escript_inverse: only matrix dimensions 1,2,3 are supported right now."
      #return out

def eigenvalues(arg):
    """
    Returns the eigenvalues of the square matrix ``arg``.

    :param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal. It must also be symmetric, ie.
                ``transpose(arg)==arg`` (this is not checked).
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the eigenvalues in increasing order
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :note: for `escript.Data` and `Symbol` objects the dimension is
           restricted to 3.
    """
    if isinstance(arg,numpy.ndarray):
      out=numpy.linalg.eigvals((arg+numpy.transpose(arg))/2.)
      out.sort()
      return out
    elif isinstance(arg,escore.Data):
      return arg._eigenvalues()
    elif isinstance(arg,complex):
      return arg
    elif isinstance(arg,float):
      return arg
    elif isinstance(arg,int):
      return float(arg)
    elif isinstance(arg,sym.Symbol):
      return sym.symfn.eigenvalues(arg)
    else:
      raise TypeError("eigenvalues: Unknown argument type.")

def eigenvalues_and_eigenvectors(arg):
    """
    Returns the eigenvalues and eigenvectors of the square matrix ``arg``.

    :param arg: square matrix. Must have rank 2 and the first and second
                dimension must be equal. It must also be symmetric, ie.
                ``transpose(arg)==arg`` (this is not checked).
    :type arg: `escript.Data`
    :return: the eigenvalues and eigenvectors. The eigenvalues are ordered by
             increasing value. The eigenvectors are orthogonal and normalized.
             If V are the eigenvectors then V[:,i] is the eigenvector
             corresponding to the i-th eigenvalue.
    :rtype: `tuple` of `escript.Data`
    :note: The dimension is restricted to 3.
    """
    if isinstance(arg,numpy.ndarray):
      raise TypeError("eigenvalues_and_eigenvectors does not support numpy.ndarray arguments")
    elif isinstance(arg,escore.Data):
      return arg._eigenvalues_and_eigenvectors()
    elif isinstance(arg,complex):
      return (numpy.array([[arg]],numpy.cdouble_),numpy.ones((1,1),numpy.cdouble_))
    elif isinstance(arg,float):
      return (numpy.array([[arg]],numpy.float_),numpy.ones((1,1),numpy.float_))
    elif isinstance(arg,int):
      return (numpy.array([[arg]],numpy.float_),numpy.ones((1,1),numpy.float_))
    elif isinstance(arg,sym.Symbol):
      return sym.symfn.eigenvalues_and_eigenvectors(arg)
    else:
      raise TypeError("eigenvalues: Unknown argument type.")

def mult(arg0,arg1):
       """
       Product of ``arg0`` and ``arg1``.

       :param arg0: first term
       :type arg0: `Symbol`, ``float``, ``int``, `escript.Data` or
                   ``numpy.ndarray``
       :param arg1: second term
       :type arg1: `Symbol`, ``float``, ``int``, `escript.Data` or
                   ``numpy.ndarray``
       :return: the product of ``arg0`` and ``arg1``
       :rtype: `Symbol`, ``float``, ``int``, `escript.Data` or
               ``numpy.ndarray``
       :note: The shape of both arguments is matched according to the rules
              used in `matchShape`.
       """
       args=matchShape(arg0,arg1)
       if testForZero(args[0]) or testForZero(args[1]):
          return numpy.zeros(getShape(args[0]),numpy.double)
       else:
          if isinstance(args[0],numpy.ndarray):
              return args[1]*args[0]
          else:
              return args[0]*args[1]

def maximum(*args):
    """
    The maximum over arguments ``args``.

    :param args: arguments
    :type args: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or
                ``float``
    :return: an object which in each entry gives the maximum of the
             corresponding values in ``args``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or
            ``float`` depending on the input
    """
    if max([isinstance(v,sym.Symbol) for v in args]):
        return sym.symfn.maximum(*args)
    out=None
    for a in args:
       if out is None:
          out=a*1.
       else:
          if isinstance(out,escore.Data) and isinstance(a,escore.Data):
             if out.getRank()==0 and a.getRank()>0:
                # We need to consider the case where we have scalars and
                # higher ranked objects mixed. If the scalar was first it will
                # get picked as the initial out and we have a problem, so we
                # swap the objects
                res=a.copy() #Deep copy of a
                res.copyWithMask(out,wherePositive(out-a))
                out=res
             else:
                out.copyWithMask(a,wherePositive(a-out))
          else:
             if isinstance(a, numpy.ndarray): 
                diff=-out+a
             else:
                diff=a-out
             temp=mult(whereNonPositive(diff),out)+mult(wherePositive(diff),a)
             if isinstance(out,numpy.ndarray) and isinstance(a,numpy.ndarray):
                # we need to convert the result to an array 
                temp=numpy.array(temp)             
             out=temp
    return out

def minimum(*args):
    """
    The minimum over arguments ``args``.

    :param args: arguments
    :type args: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or
                ``float``
    :return: an object which gives in each entry the minimum of the
             corresponding values in ``args``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or
            ``float`` depending on the input
    """
    if max([isinstance(v,sym.Symbol) for v in args]):
        return sym.symfn.minimum(*args)
    out=None
    for a in args:
       if out is None:
          if isinstance(a, numpy.ndarray):
              out=a.copy()
          else:
              out=a*1.
       else:
          if isinstance(out,escore.Data) and isinstance(a,escore.Data):
             if out.getRank()==0 and a.getRank()>0:
                # We need to consider the case where we have scalars and
                # higher ranked objects mixed. If the scalar was first it will
                # get picked as the initial out and we have a problem, so we
                # swap the objects
                res=a.copy() # Deep copy of a
                res.copyWithMask(out,whereNegative(out-a))
                out=res
             else:
                out.copyWithMask(a,whereNegative(a-out))
          else:
             if isinstance(a, numpy.ndarray): 
                diff=-out+a
             else:
                diff=a-out
             #out=add(out,mult(whereNegative(diff),diff))
             temp=mult(whereNonNegative(diff),out)+mult(whereNegative(diff),a)
             if isinstance(out,numpy.ndarray) and isinstance(a,numpy.ndarray):
                # we need to convert the result to an array 
                temp=numpy.array(temp)
             out=temp
    return out

def clip(arg,minval=None,maxval=None):
    """
    Cuts the values of ``arg`` between ``minval`` and ``maxval``.

    :param arg: argument
    :type arg: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or
               ``float``
    :param minval: lower range. If None no lower range is applied
    :type minval: ``float`` or ``None``
    :param maxval: upper range. If None no upper range is applied
    :type maxval: ``float`` or ``None``
    :return: an object that contains all values from ``arg`` between ``minval``
             and ``maxval``
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``int`` or
            ``float`` depending on the input
    :raise ValueError: if ``minval>maxval``
    """
    if isinstance(arg, sym.Symbol):
        clip_item=lambda item: sym.symfn.clip(item, minval, maxval)
        return arg.applyfunc(clip_item)
    if not minval is None and not maxval is None:
       if minval>maxval:
          raise ValueError("minval = %s must be less than maxval %s"%(minval,maxval))
    if minval is None:
        tmp=arg
    else:
        tmp=maximum(minval,arg)
    if maxval is None:
        return tmp
    else:
        return minimum(tmp,maxval)

def cross(arg0,arg1):
    """
    Cross product of the two arguments ``arg0`` and ``arg1`` which need to be shape (3,).

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol` 
    :return: the cross product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`
            depending on ``arg0``
    :raise ValueError: if the shapes of the arguments are not (3,)
    """
    if isinstance(arg0, numpy.ndarray) and isinstance(arg1, numpy.ndarray):
        out=numpy.cross(arg0, arg1)
    else:
        sh0=getShape(arg0)
        sh1=getShape(arg1)
        if not sh0 == (3,):
            raise ValueError("cross: arg0 needs to be of shape (3,)")
        if not sh1 == (3,):
            raise ValueError("cross: arg1 needs to be of shape (3,)")
        
        if isinstance(arg0, sym.Symbol):
            if isinstance(arg1, sym.Symbol):
                out=Symbol(arg0.name+"x"+arg1.name, (3,))
            else: 
                out=Symbol(arg0.name+"x"+str(type(arg1)), (3,))
        elif isinstance(arg0, escore.Data):
            if isinstance(arg1, sym.Symbol):
                out=Symbol(str(type(arg0))+"x"+arg1.name, (3,))
            else:
                out=escore.Data(0.,(3,), arg0.getFunctionSpace())
        elif isinstance(arg1, escore.Data):
            if isinstance(arg0, sym.Symbol):
                out=Symbol(str(type(arg0))+"x"+arg1.name, (3,))
            else:
                out=escore.Data(0.,(3,), arg1.getFunctionSpace())
        else:
            raise TypeError("cross: argument type not supported")
        
        out[0]=arg0[1]*arg1[2]-arg0[2]*arg1[1]
        out[1]=arg0[2]*arg1[0]-arg0[0]*arg1[2]
        out[2]=arg0[0]*arg1[1]-arg0[1]*arg1[0]
    return out

def inner(arg0,arg1):
    """
    Inner product of the two arguments. The inner product is defined as:

    `out=Sigma_s arg0[s]*arg1[s]`

    where s runs through ``arg0.Shape``.

    ``arg0`` and ``arg1`` must have the same shape.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the inner product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``
            depending on the input
    :raise ValueError: if the shapes of the arguments are not identical
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not sh0==sh1:
        raise ValueError("inner: shape of arguments does not match")
    return generalTensorProduct(arg0,arg1,axis_offset=len(sh0))

def outer(arg0,arg1):
    """
    The outer product of the two arguments. The outer product is defined as:

    ``out[t,s]=arg0[t]*arg1[s]``

    where
        - s runs through ``arg0.Shape``
        - t runs through ``arg1.Shape``

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the outer product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    return generalTensorProduct(arg0,arg1,axis_offset=0)

def matrixmult(arg0,arg1):
    """
    See `matrix_mult`.
    """
    return matrix_mult(arg0,arg1)

def matrix_mult(arg0,arg1):
    """
    matrix-matrix or matrix-vector product of the two arguments.

    `out[s0]=Sigma_{r0} arg0[s0,r0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[r0,s1]`

    The second dimension of ``arg0`` and the first dimension of ``arg1`` must
    match.

    :param arg0: first argument of rank 2
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of at least rank 1
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the matrix-matrix or matrix-vector product of ``arg0`` and ``arg1``
             at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not len(sh0)==2 :
        raise ValueError("first argument must have rank 2")
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError("second argument must have rank 1 or 2")
    return generalTensorProduct(arg0,arg1,axis_offset=1)

def tensormult(arg0,arg1):
    """
    See `tensor_mult`.
    """
    return tensor_mult(arg0,arg1)

def tensor_mult(arg0,arg1):
    """
    The tensor product of the two arguments.

    For ``arg0`` of rank 2 this is

    `out[s0]=Sigma_{r0} arg0[s0,r0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[r0,s1]`

    and for ``arg0`` of rank 4 this is

    `out[s0,s1,s2,s3]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2,s3]`

    or

    `out[s0,s1,s2]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1,s2]`

    or

    `out[s0,s1]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[r0,r1]`

    In the first case the second dimension of ``arg0`` and the last dimension of
    ``arg1`` must match and in the second case the two last dimensions of ``arg0``
    must match the two first dimensions of ``arg1``.

    :param arg0: first argument of rank 2 or 4
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of shape greater than 1 or 2 depending on the
                 rank of ``arg0``
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the tensor product of ``arg0`` and ``arg1`` at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
        return generalTensorProduct(arg0,arg1,axis_offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
        return generalTensorProduct(arg0,arg1,axis_offset=2)
    else:
        raise ValueError("tensor_mult: first argument must have rank 2 or 4 and second rank must be in (1,2) or (2,3,4) respectively.")

def generalTensorProduct(arg0,arg1,axis_offset=0):
    """
    Generalized tensor product.

    `out[s,t]=Sigma_r arg0[s,r]*arg1[r,t]`

    where
        - s runs through ``arg0.Shape[:arg0.ndim-axis_offset]``
        - r runs through ``arg1.Shape[:axis_offset]``
        - t runs through ``arg1.Shape[axis_offset:]``

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the general tensor product of ``arg0`` and ``arg1`` at each data
             point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if (isinstance(arg0,float) or isinstance(arg0,complex)) and (isinstance(arg1,float) or isinstance(arg1,complex)):
         return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg1 are both numpy.ndarray or escript.Data,
    # or one is a Symbol and the other either of the allowed types
    if isinstance(arg0,sym.Symbol):
       sh0=arg0.getShape()
       sh1=getShape(arg1)
       if not sh0[arg0.getRank()-axis_offset:]==sh1[:axis_offset]:
          raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset))
       if isinstance(arg1,float):
          return arg0*arg1
       elif isinstance(arg1,numpy.ndarray) or isinstance(arg1, sym.Symbol):
          return arg0.tensorProduct(arg1, axis_offset)
       elif isinstance(arg1, escore.Data):
          raise TypeError("tensor product of Symbol and Data not supported yet")
    elif isinstance(arg0,numpy.ndarray):
       if not arg0.shape[arg0.ndim-axis_offset:]==arg1.shape[:axis_offset]:
          raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset))
       arg0_c=arg0.copy()
       arg1_c=arg1.copy()
       sh0,sh1=arg0.shape,arg1.shape
       d0,d1,d01=1,1,1
       for i in sh0[:arg0.ndim-axis_offset]: d0*=i
       for i in sh1[axis_offset:]: d1*=i
       for i in sh1[:axis_offset]: d01*=i
       arg0_c.resize((d0,d01))
       arg1_c.resize((d01,d1))
       if arg0_c.dtype.kind=='c':
           restype=arg0_c.dtype
       else:
           restype=arg1_c.dtype
       out=numpy.zeros((d0,d1),restype)
       for i0 in range(d0):
          for i1 in range(d1):
             out[i0,i1]=numpy.sum(arg0_c[i0,:]*arg1_c[:,i1])
       out.resize(sh0[:arg0.ndim-axis_offset]+sh1[axis_offset:])
       return out
    elif isinstance(arg0,escore.Data):
       if isinstance(arg1, sym.Symbol):
          raise TypeError("tensor product of Data and Symbol not supported yet")
       return escript_generalTensorProduct(arg0,arg1,axis_offset) # this call has to be replaced by escript._generalTensorProduct(arg0,arg1,axis_offset)
    raise TypeError("generalTensorProduct: Unsupported argument type")

def escript_generalTensorProduct(arg0,arg1,axis_offset,transpose=0):
    "arg0 and arg1 are both Data objects but not necessarily on the same function space. They could be identical!!!"
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, transpose)

def transposed_matrix_mult(arg0,arg1):
    """
    transposed(matrix)-matrix or transposed(matrix)-vector product of the two
    arguments.

    `out[s0]=Sigma_{r0} arg0[r0,s0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[r0,s0]*arg1[r0,s1]`

    The function call ``transposed_matrix_mult(arg0,arg1)`` is equivalent to
    ``matrix_mult(transpose(arg0),arg1)``.

    The first dimension of ``arg0`` and ``arg1`` must match.

    :param arg0: first argument of rank 2
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of at least rank 1
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the product of the transpose of ``arg0`` and ``arg1`` at each data
             point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not len(sh0)==2 :
        raise ValueError("first argument must have rank 2")
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError("second argument must have rank 1 or 2")
    return generalTransposedTensorProduct(arg0,arg1,axis_offset=1)

def transposed_tensor_mult(arg0,arg1):
    """
    The tensor product of the transpose of the first and the second argument.

    For ``arg0`` of rank 2 this is

    `out[s0]=Sigma_{r0} arg0[r0,s0]*arg1[r0]`

    or

    `out[s0,s1]=Sigma_{r0} arg0[r0,s0]*arg1[r0,s1]`

    and for ``arg0`` of rank 4 this is

    `out[s0,s1,s2,s3]=Sigma_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1,s2,s3]`

    or

    `out[s0,s1,s2]=Sigma_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1,s2]`

    or

    `out[s0,s1]=Sigma_{r0,r1} arg0[r0,r1,s0,s1]*arg1[r0,r1]`

    In the first case the first dimension of ``arg0`` and the first dimension of
    ``arg1`` must match and in the second case the two first dimensions of
    ``arg0`` must match the two first dimensions of ``arg1``.

    The function call ``transposed_tensor_mult(arg0,arg1)`` is equivalent to
    ``tensor_mult(transpose(arg0),arg1)``.

    :param arg0: first argument of rank 2 or 4
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of shape greater of 1 or 2 depending on the
                 rank of ``arg0``
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the tensor product of transpose of arg0 and arg1 at each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
       return generalTransposedTensorProduct(arg0,arg1,axis_offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
       return generalTransposedTensorProduct(arg0,arg1,axis_offset=2)
    else:
        raise ValueError("first argument must have rank 2 or 4")

def generalTransposedTensorProduct(arg0,arg1,axis_offset=0):
    """
    Generalized tensor product of transposed of ``arg0`` and ``arg1``.

    `out[s,t]=Sigma_r arg0[r,s]*arg1[r,t]`

    where
        - s runs through ``arg0.Shape[axis_offset:]``
        - r runs through ``arg0.Shape[:axis_offset]``
        - t runs through ``arg1.Shape[axis_offset:]``

    The function call ``generalTransposedTensorProduct(arg0,arg1,axis_offset)``
    is equivalent to
    ``generalTensorProduct(transpose(arg0,arg0.ndim-axis_offset),arg1,axis_offset)``.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the general tensor product of ``transpose(arg0)`` and ``arg1`` at
             each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if (isinstance(arg0,float) and isinstance(arg1,float)) or (isinstance(arg0,complex) and isinstance(arg1,complex)): return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg1 are both numpy.ndarray or escript.Data,
    # or one is a Symbol and the other either of the allowed types
    if isinstance(arg0,sym.Symbol):
       sh0=arg0.getShape()
       sh1=getShape(arg1)
       if not sh0[:axis_offset]==sh1[:axis_offset]:
          raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset))
       if isinstance(arg1,float):
          return arg0*arg1
       elif isinstance(arg1,numpy.ndarray) or isinstance(arg1, sym.Symbol):
          return arg0.transposedTensorProduct(arg1, axis_offset)
       elif isinstance(arg1, escore.Data):
          raise TypeError("tensor product of Symbol and Data not supported yet")
    elif isinstance(arg0,numpy.ndarray):
       if not arg0.shape[:axis_offset]==arg1.shape[:axis_offset]:
           raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset))
       arg0_c=arg0.copy()
       arg1_c=arg1.copy()
       sh0,sh1=arg0.shape,arg1.shape
       d0,d1,d01=1,1,1
       for i in sh0[axis_offset:]: d0*=i
       for i in sh1[axis_offset:]: d1*=i
       for i in sh0[:axis_offset]: d01*=i
       arg0_c.resize((d01,d0))
       arg1_c.resize((d01,d1))
       target_type=arg0.dtype if arg0.dtype.kind=='c' else arg1.dtype
       out=numpy.zeros((d0,d1), target_type)
       for i0 in range(d0):
                for i1 in range(d1):
                     out[i0,i1]=numpy.sum(arg0_c[:,i0]*arg1_c[:,i1])
       out.resize(sh0[axis_offset:]+sh1[axis_offset:])
       return out
    elif isinstance(arg0,escore.Data):
       if isinstance(arg1, sym.Symbol):
          raise TypeError("tensor product of Data and Symbol not supported yet")
       # this call has to be replaced by escript._generalTensorProduct(arg0,arg1,axis_offset)
       return escript_generalTransposedTensorProduct(arg0,arg1,axis_offset)

# this should be escript._generalTransposedTensorProduct
def escript_generalTransposedTensorProduct(arg0,arg1,axis_offset):
    "arg0 and arg1 are both Data objects but not necessarily on the same function space. They could be identical!!!"
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, 1)

def matrix_transposed_mult(arg0,arg1):
    """
    matrix-transposed(matrix) product of the two arguments.

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[s1,r0]`

    The function call ``matrix_transposed_mult(arg0,arg1)`` is equivalent to
    ``matrix_mult(arg0,transpose(arg1))``.

    The last dimensions of ``arg0`` and ``arg1`` must match.

    :param arg0: first argument of rank 2
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of rank 1 or 2
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the product of ``arg0`` and the transposed of ``arg1`` at each data
             point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    :raise ValueError: if the shapes of the arguments are not appropriate
    """
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if not len(sh0)==2 :
        raise ValueError("first argument must have rank 2")
    if not len(sh1)==2 and not len(sh1)==1:
        raise ValueError("second argument must have rank 1 or 2")
    return generalTensorTransposedProduct(arg0,arg1,axis_offset=1)

def tensor_transposed_mult(arg0,arg1):
    """
    The tensor product of the first and the transpose of the second argument.

    For ``arg0`` of rank 2 this is

    `out[s0,s1]=Sigma_{r0} arg0[s0,r0]*arg1[s1,r0]`

    and for ``arg0`` of rank 4 this is

    `out[s0,s1,s2,s3]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[s2,s3,r0,r1]`

    or

    `out[s0,s1,s2]=Sigma_{r0,r1} arg0[s0,s1,r0,r1]*arg1[s2,r0,r1]`

    In the first case the second dimension of ``arg0`` and ``arg1`` must
    match and in the second case the two last dimensions of ``arg0`` must match
    the two last dimensions of ``arg1``.

    The function call ``tensor_transpose_mult(arg0,arg1)`` is equivalent to
    ``tensor_mult(arg0,transpose(arg1))``.

    :param arg0: first argument of rank 2 or 4
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :param arg1: second argument of shape greater of 1 or 2 depending on rank
                 of ``arg0``
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`
    :return: the tensor product of the transposed of ``arg0`` and ``arg1`` at
             each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    #return tensor_mult(arg0, transpose(arg1))
    # The code below has a bug, for now we use the less efficient call above
    sh0=getShape(arg0)
    sh1=getShape(arg1)
    if len(sh0)==2 and ( len(sh1)==2 or len(sh1)==1 ):
       return generalTensorTransposedProduct(arg0,arg1,axis_offset=1)
    elif len(sh0)==4 and (len(sh1)==2 or len(sh1)==3 or len(sh1)==4):
        if len(sh1)==2:
            return generalTensorTransposedProduct(arg0,transpose(arg1),axis_offset=2)
        else:
            return generalTensorTransposedProduct(arg0,arg1,axis_offset=2)
    else:
        raise ValueError("first argument must have rank 2 or 4")

def generalTensorTransposedProduct(arg0,arg1,axis_offset=0):
    """
    Generalized tensor product of ``arg0`` and transpose of ``arg1``.

    `out[s,t]=Sigma_r arg0[s,r]*arg1[t,r]`

    where
        - s runs through ``arg0.Shape[:arg0.ndim-axis_offset]``
        - r runs through ``arg0.Shape[arg1.ndim-axis_offset:]``
        - t runs through ``arg1.Shape[arg1.ndim-axis_offset:]``

    The function call ``generalTensorTransposedProduct(arg0,arg1,axis_offset)``
    is equivalent to
    ``generalTensorProduct(arg0,transpose(arg1,arg1.ndim-axis_offset),axis_offset)``.

    :param arg0: first argument
    :type arg0: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :param arg1: second argument
    :type arg1: ``numpy.ndarray``, `escript.Data`, `Symbol`, ``float``, ``int``
    :return: the general tensor product of ``arg0`` and ``transpose(arg1)`` at
             each data point
    :rtype: ``numpy.ndarray``, `escript.Data`, `Symbol` depending on the input
    """
    if ((isinstance(arg0,float) or isinstance(arg0,complex)) and 
        (isinstance(arg1,float) or isinstance(arg1,complex))):
            return arg1*arg0
    arg0,arg1=matchType(arg0,arg1)
    # at this stage arg0 and arg1 are both numpy.ndarray or escript.Data,
    # or one is a Symbol and the other either of the allowed types
    if isinstance(arg0,sym.Symbol):
       sh0=arg0.getShape()
       sh1=getShape(arg1)
       r1=getRank(arg1)
       if not sh0[arg0.getRank()-axis_offset:]==sh1[r1-axis_offset:]:
          raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset))
       if isinstance(arg1,float) or isinstance(arg1,complex):
          return arg0*arg1
       elif isinstance(arg1,numpy.ndarray) or isinstance(arg1, sym.Symbol):
          return arg0.tensorTransposedProduct(arg1, axis_offset)
       elif isinstance(arg1, escore.Data):
          raise TypeError("tensor product of Symbol and Data not supported yet")
    elif isinstance(arg0,numpy.ndarray):
       if not (arg0.shape[arg0.ndim-axis_offset:]==arg1.shape[arg1.ndim-axis_offset:] or 
            arg0.shape[arg0.ndim-axis_offset:],tuple(reversed(arg1.shape[arg1.ndim-axis_offset:]))):
          raise ValueError("dimensions of last %s components in left argument don't match the first %s components in the right argument."%(axis_offset,axis_offset))
       arg0_c=arg0.copy()
       arg1_c=arg1.copy()
       sh0,sh1=arg0.shape,arg1.shape
       d0,d1,d01=1,1,1
       for i in sh0[:arg0.ndim-axis_offset]: d0*=i
       for i in sh1[:arg1.ndim-axis_offset]: d1*=i
       for i in sh1[arg1.ndim-axis_offset:]: d01*=i
       arg0_c.resize((d0,d01))
       arg1_c.resize((d1,d01))
       if arg0_c.dtype!=numpy.double:
           out=numpy.zeros((d0,d1),arg0_c.dtype)
       else:
           out=numpy.zeros((d0,d1),numpy.double)       
       for i0 in range(d0):
          for i1 in range(d1):
             out[i0,i1]=numpy.sum(arg0_c[i0,:]*arg1_c[i1,:])
       out.resize(sh0[:arg0.ndim-axis_offset]+sh1[:arg1.ndim-axis_offset])
       return out
    elif isinstance(arg0,escore.Data):
       if isinstance(arg1, sym.Symbol):
          raise TypeError("tensor product of Data and Symbol not supported yet")
       # this call has to be replaced by escript._generalTensorProduct(arg0,arg1,axis_offset)
       return escript_generalTensorTransposedProduct(arg0,arg1,axis_offset)

# this should be escript._generalTensorTransposedProduct
def escript_generalTensorTransposedProduct(arg0,arg1,axis_offset):
    "arg0 and arg1 are both Data objects but not necessarily on the same function space. They could be identical!!!"
    return C_GeneralTensorProduct(arg0, arg1, axis_offset, 2)

#=========================================================
#  functions dealing with spatial dependency
#=========================================================
def grad(arg,where=None):
    """
    Returns the spatial gradient of ``arg`` at ``where``.

    If ``g`` is the returned object, then

      - if ``arg`` is rank 0 ``g[s]`` is the derivative of ``arg`` with respect to
        the ``s``-th spatial dimension
      - if ``arg`` is rank 1 ``g[i,s]`` is the derivative of ``arg[i]`` with
        respect to the ``s``-th spatial dimension
      - if ``arg`` is rank 2 ``g[i,j,s]`` is the derivative of ``arg[i,j]`` with
        respect to the ``s``-th spatial dimension
      - if ``arg`` is rank 3 ``g[i,j,k,s]`` is the derivative of ``arg[i,j,k]``
        with respect to the ``s``-th spatial dimension.

    :param arg: function of which the gradient is to be calculated. Its rank
                has to be less than 3.
    :type arg: `escript.Data` or `Symbol`
    :param where: FunctionSpace in which the gradient is calculated.
                  If not present or ``None`` an appropriate default is used.
    :type where: ``None`` or `escript.FunctionSpace`
    :return: gradient of ``arg``
    :rtype: `escript.Data` or `Symbol`
    """
    if isinstance(arg,sym.Symbol):
       if where is None:
           return arg.grad()
       else:
           return arg.grad(where)
    elif isinstance(arg,escore.Data):
       if where is None:
          return arg._grad()
       else:
          return arg._grad(where)
    else:
       raise TypeError("grad: Unknown argument type.")

def grad_n(arg, n, where=None):
    return grad(arg, where)[n]

def curl(arg,where=None):
    """
    Returns the (spatial) curl of ``arg`` at ``where``..

    :param arg: function of which the curl is to be calculated. Its shape
                has to be (), (2,) or (3,)
    :type arg: `escript.Data` or `Symbol`
    :param where: FunctionSpace in which the gradient is calculated.
                  If not present or ``None`` an appropriate default is used.
    :type where: ``None`` or `escript.FunctionSpace`
    :return: curl of ``arg``
    :rtype: `escript.Data` or `Symbol`
    """
    if not getShape(arg) in [ (), (2,) , (3,) ]:
        raise TypeError("curl: illegal argument rank.")
    g=grad(arg, where)
    if isinstance(arg,sym.Symbol):
        
        if arg.getShape() == ():
            out=sym.Symbol("curl("+arg.name+")", (2,), g.getFunctionSpace())
            out[0]=g[1]
            out[1]=-g[0]
        elif arg.getShape() == (2,):
            out=g[1,0]-g[0,1]
        else:
            out=sym.Symbol("curl("+arg.name+")", (3,), g.getFunctionSpace())
            out[0]=g[2,1]-g[1,2]
            out[1]=g[0,2]-g[2,0]
            out[2]=g[1,0]-g[0,1]
            
    elif isinstance(arg,escore.Data):
        if arg.getShape() == ():
            if g.isComplex():
                out=escore.Data(0j, (2,), g.getFunctionSpace())
            else:
                out=escore.Data(0, (2,), g.getFunctionSpace())
            out[0]=g[1]
            out[1]=-g[0]
        elif arg.getShape() == (2,):
            out=g[1,0]-g[0,1]
        else:
            if g.isComplex():
                out=escore.Data(0j, (3,), g.getFunctionSpace())
            else:
                out=escore.Data(0, (3,), g.getFunctionSpace())
            out[0]=g[2,1]-g[1,2]
            out[1]=g[0,2]-g[2,0]
            out[2]=g[1,0]-g[0,1]
    else:
       raise TypeError("curl: Unknown argument type.")
    return out

def integrate(arg,where=None):
    """
    Returns the integral of the function ``arg`` over its domain. If ``where`` is
    present ``arg`` is interpolated to ``where`` before integration.

    :param arg: the function which is integrated
    :type arg: `escript.Data` or `Symbol`
    :param where: FunctionSpace in which the integral is calculated.
                  If not present or ``None`` an appropriate default is used.
    :type where: ``None`` or `escript.FunctionSpace`
    :return: integral of ``arg``
    :rtype: ``float``, ``numpy.ndarray`` or `Symbol`
    """
    if isinstance(arg,escore.Data):
       if not where is None: arg=escore.Data(arg,where)
       if arg.getRank()==0:
          return arg._integrateToTuple()[0]
       else:
          return numpy.array(arg._integrateToTuple())
    elif isinstance(arg,sym.Symbol):
       return sym.symfn.integrate(arg, where)
    else:
       arg2=escore.Data(arg,where)
       if arg2.getRank()==0:
          return arg2._integrateToTuple()[0]
       else:
          return numpy.array(arg2._integrateToTuple())

def interpolate(arg,where):
    """
    Interpolates the function into the `FunctionSpace` ``where``. If the
    argument ``arg`` has the requested function space ``where`` no interpolation
    is performed and ``arg`` is returned.

    :param arg: interpolant
    :type arg: `escript.Data` or `Symbol`
    :param where: `FunctionSpace` to be interpolated to
    :type where: `escript.FunctionSpace`
    :return: interpolated argument
    :rtype: ``escript.Data`` or `Symbol`
    """
    if isinstance(arg,escore.Data):
       if arg.isEmpty():
          return arg
       elif where == arg.getFunctionSpace():
          return arg
       else:
          # work around for Bug #390
          if arg.isComplex():
              return interpolate(arg.real(), where)+1j*interpolate(arg.imag(), where)
          else:
              if arg.isConstant():
                  arg.expand()
              return escore.Data(arg,where)
    elif isinstance(arg,sym.Symbol):
       return sym.symfn.interpolate(arg, where)
    else:
       return escore.Data(arg, where)

def div(arg,where=None):
    """
    Returns the divergence of ``arg`` at ``where``.

    :param arg: function of which the divergence is to be calculated. Its
                shape has to be (d,) where d is the spatial dimension.
    :type arg: `escript.Data` or `Symbol`
    :param where: `FunctionSpace` in which the divergence will be calculated.
                  If not present or ``None`` an appropriate default is used.
    :type where: ``None`` or `escript.FunctionSpace`
    :return: divergence of ``arg``
    :rtype: `escript.Data` or `Symbol`
    """
    if isinstance(arg,escore.Data):
        dim=arg.getDomain().getDim()
        if not arg.getShape()==(dim,):
            raise ValueError("div: expected shape is (%s,)"%dim)
    elif not isinstance(arg, sym.Symbol):
        raise TypeError("div: argument type not supported")
    return trace(grad(arg,where))

def jump(arg,domain=None):
    """
    Returns the jump of ``arg`` across the continuity of the domain.

    :param arg: argument
    :type arg: `escript.Data` or `Symbol`
    :param domain: the domain where the discontinuity is located. If domain is
                   not present or equal to ``None`` the domain of ``arg`` is used.
    :type domain: ``None`` or `escript.Domain`
    :return: jump of ``arg``
    :rtype: `escript.Data` or `Symbol`
    """
    if domain is None: domain=arg.getDomain()
    return interpolate(arg,escore.FunctionOnContactOne(domain))-interpolate(arg,escore.FunctionOnContactZero(domain))

def L2(arg):
    """
    Returns the L2 norm of ``arg`` at ``where``.

    :param arg: function of which the L2 norm is to be calculated
    :type arg: `escript.Data` or `Symbol`
    :return: L2 norm of ``arg``
    :rtype: `float` or `Symbol`
    :note: L2(arg) is equivalent to ``sqrt(integrate(inner(arg,arg)))``
    """
    if isinstance(arg,sym.Symbol):
        return sym.symfn.L2(arg)
    return sqrt(integrate(inner(arg,arg)))

def getClosestValue(arg,origin=0):
    """
    Returns the value in ``arg`` which is closest to origin.

    :param arg: function
    :type arg: `escript.Data`
    :param origin: reference value
    :type origin: ``float`` or `escript.Data`
    :return: value in ``arg`` closest to origin
    :rtype: ``numpy.ndarray``
    """
    return arg.getValueOfGlobalDataPoint(*(length(arg-origin).minGlobalDataPoint()))

def normalize(arg,zerolength=0):
    """
    Returns the normalized version of ``arg`` (=``arg/length(arg)``).

    :param arg: function
    :type arg: `escript.Data` or `Symbol`
    :param zerolength: relative tolerance for arg == 0
    :type zerolength: ``float``
    :return: normalized ``arg`` where ``arg`` is non-zero, and zero elsewhere
    :rtype: `escript.Data` or `Symbol`
    """
    l=length(arg)
    m=whereZero(l,zerolength*Lsup(l))
    mm=1-m
    return arg*(mm/(l*mm+m))

def deviatoric(arg):
    """
    Returns the deviatoric version of ``arg``.
    """
    return arg-(trace(arg)/trace(kronecker(arg.getDomain())))*kronecker(arg.getDomain())

def vol(arg):
    """
    Returns the volume or area of the oject ``arg``

    :param arg: a geometrical object
    :type arg: `escript.FunctionSpace` or `escript.Domain`
    :rtype: ``float``
    """
    if isinstance(arg,escore.Domain): arg=escore.Function(arg)
    return integrate(escore.Scalar(1.,arg))

def meanValue(arg):
    """
    return the mean value of the argument over its domain

    :param arg: function
    :type arg: `escript.Data`
    :return: mean value
    :rtype: ``float`` or ``numpy.ndarray``
    """
    fs=arg.getFunctionSpace()
    d=fs.getDomain()
    if fs == escore.Solution(d) or fs == escore.ContinuousFunction(d):
       fs=escore.Function(d)
    if fs == escore.ReducedSolution(d) or fs == escore.ReducedContinuousFunction(d):
       fs=escore.ReducedFunction(d)
    a=vol(fs)
    if a == 0:
        raise ValueError("FunctionSpace %s with zero volume."%str(fs))
    return integrate(arg,fs)/a
 
def diameter(domain):
    """
    Returns the diameter of a domain.

    :param domain: a domain
    :type domain: `escript.Domain`
    :rtype: ``float``
    """
    return sqrt(sum( [ v**2 for v in boundingBoxEdgeLengths(domain) ] ))

def boundingBoxEdgeLengths(domain):
    """
    Returns the edge lengths of the bounding box of a domain

    :param domain: a domain
    :type domain: `escript.Domain`
    :rtype: ``list`` of ``float``
    """
    return  [ v[1]-v[0] for v in boundingBox(domain) ] 

    
def boundingBox(domain):
    """
    Returns the bounding box of a domain

    :param domain: a domain
    :type domain: `escript.Domain`
    :return: bounding box of the domain
    :rtype: ``list`` of pairs of ``float``
    """
    x=domain.getX()
    out=[]
    for i in range(domain.getDim()):
       x_i=x[i]
       out.append((inf(x_i),sup(x_i)))
    return out

class BBxInterval(object):
    def __init__(self, minval=0, maxval=0):
        self.min=minval
        self.max=maxval
    def __str__(self):
        return f"[{self.min}, {self.max}]"
    def __repr__(self):
        return f"[{self.min}, {self.max}]"
def getBoundingBox(domain):
    """
    Returns the bounding box of a domain as a list of intervals

    :param domain: a domain
    :type domain: `escript.Domain`
    :return: bounding box of the domain
    :rtype: ``list`` of ``BBxInterval``
    """
    return [ BBxInterval(minval, maxval) for minval, maxval in boundingBox(domain)]

def longestEdge(domain):
    """
    Returns the length of the longest edge of the domain

    :param domain: a domain
    :type domain: `escript.Domain`
    :return: longest edge of the domain parallel to the Cartesian axis 
    :rtype: ``float``
    """
    return max(boundingBoxEdgeLengths(domain))

def mkDir(*pathname):
    """
    creates a directory of name ``pathname`` if the directory does not exist.

    :param pathname: valid path name
    :type pathname: ``str`` or ``sequence of strings``
    :note: The method is MPI safe.
    """
    errno = 0
    p_fail = None
    ex = None
    if getMPIRankWorld() == 0:
      for p in pathname:
       if os.path.exists(p):
          if not os.path.isdir(p):
                errno = 2
                p_fail = p
       else:
          try:
              os.makedirs(p)
          except Exception as e:
              ex = e
              errno = 1
              p_fail = p
    
    errno = getMPIWorldMax(errno)
    if errno > 0:
         if errno==2:
            if p_fail is None:
               raise IOError("Unable to create directory.")
            else:
               raise IOError("Unable to create directory %s. It already exists and is not a directory."%p_fail)
         elif ex is None:
            if p_fail is None:
               raise IOError("Unable to create directory.")
            else:
               raise IOError("Unable to create directory %s."%p_fail)
         else:
            if len(str(ex)) > 0:
               raise IOError(str(ex))
            else:
               if p_fail is None:
                  raise IOError("Unable to create directory.")
               else:
                  raise IOError("Unable to create directory %s."%p_fail)

def getRegionTags(function_space):
    """
    Returns the tag distribution of a function_space as a Data object.

    :param function_space: function_space to be used
    :type function_space: `escript.FunctionSpace`
    :return: a data object which contains the tag distribution
    :rtype: `escript.Data` of rank 0 with ReducedFunction attributes.
    """

    out = escore.Scalar(0., function_space)
    for tag in function_space.getListOfTags():
        out.setTaggedValue(tag,float(tag))
    out.expand()
    return out

class FileWriter(object):
    """
    Interface to write data to a file. In essence this class wraps the standard ``file`` object to write data that are global in MPI
    to a file. In fact, data are written on the processor with MPI rank 0 only. It is recommended to use ``FileWriter`` rather than ``open`` in order to write
    code that is running with as well as without MPI. It is safe to use ``open`` under MPI to read data which are global under MPI.

    :var name: name of file
    :var mode: access mode (='w' or ='a')
    :var closed: True to indicate closed file
    :var newlines: line separator
    """
    def __init__(self,fn,append=False,createLocalFiles=False):
         """
         Opens a file of name ``fn`` for writing. If running under MPI only the first processor with rank==0
         will open the file and write to it. If ``createLocalFiles`` each individual processor will create a file
         where for any processor with rank>0 the file name is extended by its rank. This option is normally only used for 
         debug purposes.

         :param fn: filename.
         :type fn: ``str``
         :param append: if True, open file for appending rather than overwriting
         :type append: ``bool``
         :param createLocalFiles: switches on the creation of local files.
         :type createLocalFiles: ``bool``
         """
         error=None
         errno=0
         if len(fn)==0:
             errno=1
             error="No filename provided"
         else:
             self.name=fn
             if append:
                 self.mode='a'
             else:
                 self.mode='w'
             self.__file=None
             self.closed=False
             self.newlines=os.linesep
             # if not the master:
             if getMPIRankWorld()>0:
                  if createLocalFiles:
                      fn2=fn+".%s"%getMPIRankWorld()
                      try:
                         self.__file=open(fn2,self.mode)
                      except Exception as e:
                         errno=1
                         error=e
             else:
                  try:
                      self.__file=open(fn,self.mode)
                  except Exception as e:
                      errno=1
                      error=e
         self.__handelerror(errno, error, "opening")

    def __handelerror(self,errno,e,operation):
         errno=getMPIWorldMax(errno)
         if errno>0:
            if e is None:
               raise IOError("Unable to access file %s in mode %s for %s."%(self.name,self.mode,operation))
            else:
               raise IOError(str(e))
         
    def close(self):
        """
        Closes the file
        """
        errno=0
        e=None
        try:
           if not self.__file is None:
               self.__file.close()
        except Exception as e:
           errno=1
        self.__handelerror(errno,e,"closing")
        self.closed=True

    def flush(self):
        """
        Flush the internal I/O buffer.
        """
        errno=0
        e=None
        try:
           if not self.__file is None:
               self.__file.flush()
        except Exception as e:
           errno=1
        self.__handelerror(errno,e,"flushing")

    def write(self,txt):
        """
        Write string ``txt`` to file.

        :param txt: string ``txt`` to be written to file
        :type txt: ``str``
        """
        errno=0
        e=None
        try:
           if not self.__file is None:
               self.__file.write(txt)
        except Exception as e:
           errno=1
        self.__handelerror(errno,e,"writing")

    def writelines(self, txts):
        """
        Write the list ``txts`` of strings to the file.

        :param txts: sequence of strings to be written to file
        :type txts: any iterable object producing strings
        :note: Note that newlines are not added. This method is equivalent to calling write() for each string.
        """
        errno=0
        e=None
        try:
           if not self.__file is None:
               self.__file.writelines(txts)
        except Exception as e:
           errno=1
        self.__handelerror(errno,e,"writing strings")

def reorderComponents(arg,index):
    """
    Resorts the components of ``arg`` according to index.

    """
    raise NotImplementedError

def showEscriptParams():
    """
    Displays the parameters escript recognises with an explanation and their
    current value.
    """
    p=listEscriptParams()
    for name,value,desc in p:
       print('%s (=%s): %s'%(name, value, desc))

#Lazy related things
#These are just wrappers
def resolve(arg):
   """
   Returns the value of arg resolved.
   """
   if not isinstance(arg,Data):
        raise TypeError("Can only resolve Data.")
   if arg.isLazy():
        arg.resolve()
   return arg
   
def delay(arg):
   """
   Returns a lazy version of arg
   """
   if not isinstance(arg,Data):
         raise TypeError("Can only delay Data.")
   return arg.delay()

def positive(arg):
   """
   returns the positive part of arg
   """
   return (abs(arg)+arg)/2.

def negative(arg):
   """
   returns the negative part of arg
   """
   return (arg-abs(arg))/2.

def safeDiv(arg0, arg1, rtol=None):
    """
    returns arg0/arg1 but return 0 where arg1 is (almost) zero
    """
    if rtol is None:
      m1=whereZero(arg1,tol=0)
    else:
      m1=whereZero(arg1,tol=None, rtol=rtol)
    return arg0/(arg1+m1)*whereNonPositive(m1)

def condEval(f, tval, fval):
    """
    Wrapper to allow non-data objects to be used.
    """
    if not isinstance(tval,Data) and not isinstance(fval,Data):
        raise TypeError("At least one of the alternatives must be a Data object.")
    if isinstance(tval,Data) and isinstance(fval, Data):
        return escore._condEval(f,tval,fval)
    if not isinstance(fval, Data):
        return escore._condEval(f, tval, Data(fval, tval.getShape(), tval.getFunctionSpace()))
    return escore._condEval(f, Data(fval, fval.getShape(), fval.getFunctionSpace()), fval )

def polarToCart(r, phase):
    """
    conversion from cartesian to polar coordinates
    
    :param r: length
    :type r: any float type object
    :param phase: the phase angle in rad
    :type phase: any float type object
    :return: cartesian representation as complex number
    :rtype: appropriate complex
    """
    return r*exp(1*j*phase)

def phase(arg):
    """
    return the "phase"/"arg"/"angle" of a number
    """
    if isinstance(arg, numbers.Number):
        return cmath.phase(arg)
    if isinstance(arg, Data):
        return arg.phase()
    if isinstance(arg, numpy.ndarray):
        return numpy.phase(arg)
    return arg.phase()

def makeTagMap(fs):
    """
    Produce an expanded Data over the function space where
    the value is the tag associated with the sample
    """
    out=escore.Scalar(0, fs)    # The default tag is zero anyway
    for t in fs.getListOfTags():
        out.setTaggedValue(t,t)
    out.expand()
    return out

def real(arg):
    """
    returns the real part of arg
    """
    return arg.real()

def imag(arg):
    """
    returns the imaginary part of arg
    """
    return arg.imag()

def conjugate(arg):
    """
    returns the complex conjugate of arg
    """
    return arg.conjugate()


#=========================================================
#   Re-exports from submodules for backward compatibility
#=========================================================
# These imports make functions from the submodules available through util
# This enables gradual migration without breaking existing code
from .util_math import *
from .util_tensor import *
from .util_io import *
