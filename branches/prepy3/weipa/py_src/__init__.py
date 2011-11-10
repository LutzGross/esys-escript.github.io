
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from .weipacpp import visitInitialize, visitPublishData

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
    for n,d in data.items():
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
    for n,d in sorted(new_data.iteritems()):
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
            for i,p in metadata_schema.items():
                ss="%s xmlns:%s=\"%s\""%(ss, i, p)
        else:
            ss=metadata_schema
    dataset.setMetadataSchemaString(ss.strip(), ms.strip())
    dataset.setSaveMeshData(write_meshdata)
    return dataset.saveVTK(filename)

def _saveVTK(filename, domain=None, metadata='', metadata_schema=None, data={}):
    """
    This is only here to support the deprecated domain C++ member saveVTK().
    """
    return saveVTK(filename, domain, metadata, metadata_schema, **data)

