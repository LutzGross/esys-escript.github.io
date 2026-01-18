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
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""
I/O utilities and domain helper functions for escript.

This module provides:
- Constants: EPSILON, DBLE_MAX, getEpsilon, getMaxFloat
- Data I/O: interpolateTable, saveDataCSV, getNumpy, convertToNumpy, NumpyToData, saveESD
- Domain utilities: getTagNames, insertTagNames, insertTaggedValues, diameter,
  boundingBoxEdgeLengths, boundingBox, getBoundingBox, longestEdge, getRegionTags, makeTagMap
- Filesystem: mkDir, FileWriter class, BBxInterval class

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var EPSILON: smallest positive value with 1.<1.+EPSILON
:var DBLE_MAX: largest positive float
"""


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"
__author__="Lutz Gross, l.gross@uq.edu.au"

import math
import os
import numpy

from . import escriptcpp as escore
from .escriptcpp import Data, getMPIRankWorld, getMPIWorldMax, hasFeature


#=========================================================
#   Constants and helpers
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


#=========================================================
#   Tag utilities
#=========================================================
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

def makeTagMap(fs):
    """
    Produces an expanded Data over the function space where
    the value is the tag associated with the sample.

    :param fs: the function space
    :type fs: `escript.FunctionSpace`
    :return: Data object with tag values at each sample point
    :rtype: `escript.Data`
    """
    out=escore.Scalar(0, fs)    # The default tag is zero anyway
    for t in fs.getListOfTags():
        out.setTaggedValue(t,t)
    out.expand()
    return out


#=========================================================
#   Data I/O functions
#=========================================================
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
    # Local import to avoid circular dependency
    from . import util

    new_data = {}
    for n,d in sorted(data.items(), key=lambda x: x[0]):
          if not d.isEmpty():
            fs = d.getFunctionSpace()
            domain2 = fs.getDomain()
            if fs == escore.Solution(domain2):
               new_data[n]=util.interpolate(d,escore.ContinuousFunction(domain2))
            elif fs == escore.ReducedSolution(domain2):
               new_data[n]=util.interpolate(d,escore.ReducedContinuousFunction(domain2))
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


#=========================================================
#   Domain utilities
#=========================================================
def diameter(domain):
    """
    Returns the diameter of a domain.

    :param domain: a domain
    :type domain: `escript.Domain`
    :rtype: ``float``
    """
    return math.sqrt(sum( [ v**2 for v in boundingBoxEdgeLengths(domain) ] ))

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
    # Local import to avoid circular dependency
    from . import util

    x=domain.getX()
    out=[]
    for i in range(domain.getDim()):
       x_i=x[i]
       out.append((util.inf(x_i), util.sup(x_i)))
    return out

class BBxInterval(object):
    """
    Represents an interval for bounding box coordinates.
    """
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


#=========================================================
#   Filesystem utilities
#=========================================================
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
