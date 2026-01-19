
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


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

from esys.escript import convertToNumpy, hasFeature
from .weipacpp import visitInitialize, visitPublishData

__nodocorecursion=['weipacpp']

import numpy as np
import warnings

try:
    # to suppress some warnings when importing VTK
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        from tvtk.api import tvtk
    HAVE_TVTK = True
except:
    HAVE_TVTK = False


def interpolateEscriptData(domain, data):
    """
    esys.weipa does not support the function spaces Solution and
    ReducedSolution. This function interpolates Data defined on those function
    spaces to compatible alternatives.
    """
    from esys.escript import Solution, ReducedSolution
    from esys.escript import ContinuousFunction, ReducedContinuousFunction
    from esys.escript.util import interpolate
    
    new_data={}
    for n,d in sorted(list(data.items()), key=lambda x: x[0]):
        if not d.isEmpty():
            fs=d.getFunctionSpace()
            if domain is None:
                domain=fs.getDomain()
            elif domain != fs.getDomain():
                raise ValueError("weipa: All Data must be on the same domain!")
            new_data[n]=d
            try:
                if fs == Solution(domain):
                    new_data[n]=interpolate(d, ContinuousFunction(domain))
                elif domain.getDescription().startswith("speckley"):
                    new_data[n]=interpolate(d, ContinuousFunction(domain))
                elif fs == ReducedSolution(domain):
                    new_data[n]=interpolate(d, ReducedContinuousFunction(domain))
            except RuntimeError as e:
                if str(e).startswith("FunctionSpaceException"):
                    pass
                else:
                    raise e

    return domain,new_data

def createDataset(domain=None, **data):
    """
    Creates and returns an esys.weipa dataset consisting of a Domain and Data
    objects. The returned object provides methods to access and export data.
    """
    from .weipacpp import EscriptDataset
    dataset=EscriptDataset()
    domain,new_data=interpolateEscriptData(domain, data)
    dataset.setDomain(domain)
    for n,d in sorted(new_data.items()):
        #TODO: data units are not supported here yet
        dataset.addData(d, n, "")
    return dataset

def saveSilo(filename, domain=None, write_meshdata=False, time=0., cycle=0,
        **data):
    """
    Writes `Data` objects and their mesh to a file using the SILO file format.

    Example::

        temp=Scalar(..)
        v=Vector(..)
        saveSilo("solution.silo", temperature=temp, velocity=v)

    ``temp`` and ``v`` are written to "solution.silo" where ``temp`` is named
    "temperature" and ``v`` is named "velocity".

    :param filename: name of the output file ('.silo' is added if required)
    :type filename: ``str``
    :param domain: domain of the `Data` objects. If not specified, the domain
                   of the given `Data` objects is used.
    :type domain: `escript.Domain`
    :param write_meshdata: whether to save mesh-related data such as element
                           identifiers, ownership etc. This is mainly useful
                           for debugging.
    :type write_meshdata: ``bool``
    :param time: the timestamp to save within the file
    :type time: ``float``
    :param cycle: the cycle (or timestep) of the data
    :type cycle: ``int``
    :keyword <name>: writes the assigned value to the Silo file using <name> as
                     identifier
    :note: All data objects have to be defined on the same domain but they may
           be defined on separate `FunctionSpace` s.
    """

    dataset = createDataset(domain, **data)
    dataset.setCycleAndTime(cycle, time)
    dataset.setSaveMeshData(write_meshdata)
    return dataset.saveSilo(filename)

def saveVTK(filename, domain=None, metadata='', metadata_schema=None,
        write_meshdata=False, time=0., cycle=0, **data):
    """
    Writes `Data` objects and their mesh to a file using the VTK XML file
    format.

    Example::

        temp=Scalar(..)
        v=Vector(..)
        saveVTK("solution.vtu", temperature=temp, velocity=v)

    ``temp`` and ``v`` are written to "solution.vtu" where ``temp`` is named
    "temperature" and ``v`` is named "velocity".

    Meta tags, e.g. a timeStamp, can be added to the file, for instance::

        tmp=Scalar(..)
        v=Vector(..)
        saveVTK("solution.vtu", temperature=tmp, velocity=v,
                metadata="<timeStamp>1.234</timeStamp>",
                metadata_schema={"gml":"http://www.opengis.net/gml"})

    The argument ``metadata_schema`` allows the definition of name spaces with
    a schema used in the definition of meta tags.

    :param filename: name of the output file ('.vtu' is added if required)
    :type filename: ``str``
    :param domain: domain of the `Data` objects. If not specified, the domain
                   of the given `Data` objects is used.
    :type domain: `escript.Domain`
    :keyword <name>: writes the assigned value to the VTK file using <name> as
                     identifier
    :param metadata: additional XML meta data which are inserted into the VTK
                     file. The meta data are marked by the tag ``<MetaData>``.
    :type metadata: ``str``
    :param metadata_schema: assigns schemas to namespaces which have been used
                            to define meta data.
    :type metadata_schema: ``dict`` with ``metadata_schema[<namespace>]=<URI>``
                           to assign the scheme ``<URI>`` to the name space
                           ``<namespace>``.
    :param write_meshdata: whether to save mesh-related data such as element
                           identifiers, ownership etc. This is mainly useful
                           for debugging.
    :type write_meshdata: ``bool``
    :param time: the timestamp to save within the file, seperate to metadata
    :type time: ``float``
    :param cycle: the cycle (or timestep) of the data
    :type cycle: ``int``
    :note: All data objects have to be defined on the same domain. They may not
           be in the same `FunctionSpace` but not all combinations of
           `FunctionSpace` s can be written to a single VTK file.
           Typically, data on the boundary and on the interior cannot be mixed.
    """

    dataset = createDataset(domain, **data)
    dataset.setCycleAndTime(cycle, time)
    ss=''
    ms=''
    if not metadata is None:
        ms=metadata
    if not metadata_schema is None:
        if hasattr(metadata_schema, 'items'):
            for i,p in sorted(metadata_schema.items(), key=lambda x: x[0]):
                ss="%s xmlns:%s=\"%s\""%(ss, i, p)
        else:
            ss=metadata_schema
    dataset.setMetadataSchemaString(ss.strip(), ms.strip())
    dataset.setSaveMeshData(write_meshdata)
    return dataset.saveVTK(filename)

def saveVoxet(filename, **data):
    """
    Writes `Data` objects to a file using the GOCAD Voxet file format as
    separate properties on the same grid.
    At the moment only Data on a `ripley` domain can be saved in this format.
    Note that this function will produce one header file (ending in .vo) and
    a separate property file for each `Data` object.

    :param filename: name of the output file ('.vo' is added if required)
    :type filename: ``str``
    :note: All data objects have to be defined on the same ripley domain and
           either defined on reduced Function or on a `FunctionSpace` that
           allows interpolation to reduced Function.
    """

    from esys.escript import ReducedFunction
    from esys.escript.util import interpolate
    from esys.ripley.ripleycpp import DATATYPE_FLOAT32, BYTEORDER_BIG_ENDIAN

    new_data={}
    domain=None
    for n,d in sorted(data.items(), key=lambda x: x[0]):
        if d.isEmpty():
            continue
        fs=d.getFunctionSpace()
        if domain is None:
            domain=fs.getDomain()
        elif domain != fs.getDomain():
            raise ValueError("saveVoxet: All Data must be on the same domain!")

        try:
            nd=interpolate(d, ReducedFunction(domain))
        except:
            raise ValueError("saveVoxet: Unable to interpolate all Data to reduced Function!")
        new_data[n]=nd

    if filename[-3:]=='.vo':
        fileprefix=filename[:-3]+"_"
    else:
        fileprefix=filename+"_"
        filename=filename+'.vo'

    origin, spacing, NE = domain.getGridParameters()

    # Voxet "origin" actually refers to centre of first cell so shift:
    origin = tuple([ origin[i] + spacing[i]/2. for i in range(len(origin)) ])

    # flip vertical origin
    origin=origin[:-1]+(-origin[-1],)
    axis_max=NE[:-1]+(-NE[-1],)
    midpoint=tuple([n/2 for n in NE])

    if domain.getDim() == 2:
        origin=origin+(0.,)
        spacing=spacing+(1.,)
        NE=NE+(1,)
        midpoint=midpoint+(0,)
        axis_max=axis_max+(0,)

    mainvar=list(new_data.keys())[0]
    f=open(filename,'w')
    f.write("GOCAD Voxet 1\nHEADER {\nname: escriptdata\n")
    f.write("sections: 3 1 1 %d 2 1 %d 3 1 %d\n"%midpoint)
    f.write("painted: on\nascii: off\n*painted*variable: %s\n}"%mainvar)
    f.write("""
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME "gocad Local"
AXIS_NAME X Y Z
AXIS_UNIT m m m
ZPOSITIVE Depth
END_ORIGINAL_COORDINATE_SYSTEM\n""")

    f.write("AXIS_O %0.2f %0.2f %0.2f\n"%origin)
    f.write("AXIS_U %0.2f 0 0\n"%spacing[0])
    f.write("AXIS_V 0 %0.2f 0\n"%spacing[1])
    f.write("AXIS_W 0 0 %0.2f\n"%spacing[2])
    f.write("AXIS_MIN 0 0 0\n")
    f.write("AXIS_MAX %d %d %d\n"%axis_max)
    f.write("AXIS_N %d %d %d\n"%NE)
    f.write("\n")

    num=0
    for n,d in sorted(new_data.items(), key=lambda x: x[0]):
        num=num+1
        propfile=fileprefix+n
        domain.writeBinaryGrid(d, propfile, BYTEORDER_BIG_ENDIAN, DATATYPE_FLOAT32)
        f.write("\nPROPERTY %d %s\n"%(num, n))
        f.write("PROPERTY_SUBCLASS %d QUANTITY Float\n"%num)
        f.write("PROP_ESIZE %d 4\n"%num)
        f.write("PROP_ETYPE %d IEEE\n"%num)
        f.write("PROP_FORMAT %d RAW\n"%num)
        f.write("PROP_OFFSET %d 0\n"%num)
        f.write("PROP_FILE %d %s\n"%(num,propfile))
    f.close()

class EscriptToTVTK(object):
    """
    a simple interface from escript to TVTK for rendering with mayavi.mlab
    """
    def __init__(self, domain=None):
        """
        sets up driver for translating Data objects in domain to TVTK object.
        """
        if HAVE_TVTK == False:
            raise Exception("Failed to load the TVTK module")
        if hasFeature("boostnumpy") == False:
            raise Exception("This feature requires boost version 1.64 or higher")
        if domain is None:
            self.domain=None
        else:
            self.setDomain(domain)
    def setDomain(self, domain):
        """
        resets the domain
        """
        x=domain.getNumpyX()
        cells=domain.getConnectivityInfo()
        cell_type=domain.getVTKElementType()
        points=np.zeros(shape=(x.shape[1],3),dtype=np.float32)
        for i in range(0, x.shape[1]):
            if domain.getDim() == 2:
                points[i]=[x[0][i],x[1][i],0.0]
            else:
                points[i]=[x[0][i],x[1][i],x[2][i]]
        self.tvtk = tvtk.UnstructuredGrid(points=points)
        self.tvtk.set_cells(cell_type, cells)
        self.domain = domain
        return self
    
    def getDomain(self):
        return self.domain()
    
    def setData(self, **kwargs):
        """
        set the scalar data set:
        """
        assert self.domain is not None, "no domain set."
        for n in kwargs:
            d=kwargs[n]

            # Check that the domains match
            assert kwargs[n].getDomain() == self.domain, "domain of argument %s does not match."%n

            # Check data rank
            assert d.getRank() < 2, "data %s is of rank greater than 2"%n

            # Convert to a numpy array. This requires boost v. 1.64 or higher
            data=convertToNumpy(d)

            # Work out if we have point-centered or cell-centered data
            data_type=''
            if self.domain.isCellOriented(d.getFunctionSpace().getTypeCode()):
                data_type = 'cell'
            else:
                data_type = 'point'

            # Add a third component of zeros if the array is two dimensional
            if d.getShape() == (): 
                data=np.stack((data[0]),axis=-1)
            else:
                if data.shape[0] == 2:
                    padding=np.zeros((1,data.shape[1]))
                    data=np.append(data,padding,axis=0)
                data=np.stack((data[0],data[1],data[2]),axis=-1)


            if d.getShape() == (): # scalar data
                if data_type == 'point':
                    self.tvtk.point_data.scalars = data
                    self.tvtk.point_data.scalars.name = n
                else:
                    self.tvtk.cell_data.scalars = data
                    self.tvtk.cell_data.scalars.name = n
            elif d.getShape() == (self.domain.getDim(),): # vector data
                if data_type == 'point':
                    self.tvtk.point_data.vectors = data
                    self.tvtk.point_data.vectors.name = n
                else:
                    self.tvtk.cell_data.vectors = data
                    self.tvtk.cell_data.vectors.name = n

        return self

    def getTVTK(self):
        return self.tvtk
    
