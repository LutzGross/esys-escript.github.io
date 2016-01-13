
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from .weipacpp import visitInitialize, visitPublishData

__nodocorecursion=['weipacpp']

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
    for n,d in list(data.items()):
        if not d.isEmpty():
            fs=d.getFunctionSpace()
            if domain is None:
                domain=fs.getDomain()
            elif domain != fs.getDomain():
                raise ValueError("weipa: All Data must be on the same domain!")
            if fs == Solution(domain):
                new_data[n]=interpolate(d, ContinuousFunction(domain))
            elif fs == ReducedSolution(domain):
                new_data[n]=interpolate(d, ReducedContinuousFunction(domain))
            else:
                new_data[n]=d
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

def saveSilo(filename, domain=None, write_meshdata=False, **data):
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
    :keyword <name>: writes the assigned value to the Silo file using <name> as
                     identifier
    :note: All data objects have to be defined on the same domain but they may
           be defined on separate `FunctionSpace` s.
    """

    dataset = createDataset(domain, **data)
    dataset.setSaveMeshData(write_meshdata)
    return dataset.saveSilo(filename)

def saveVTK(filename, domain=None, metadata='', metadata_schema=None, write_meshdata=False, **data):
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
    :note: All data objects have to be defined on the same domain. They may not
           be in the same `FunctionSpace` but not all combinations of
           `FunctionSpace` s can be written to a single VTK file.
           Typically, data on the boundary and on the interior cannot be mixed.
    """

    dataset = createDataset(domain, **data)
    ss=''
    ms=''
    if not metadata is None:
        ms=metadata
    if not metadata_schema is None:
        if hasattr(metadata_schema, 'items'):
            for i,p in list(metadata_schema.items()):
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
    for n,d in list(data.items()):
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
    for n,d in list(new_data.items()):
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

