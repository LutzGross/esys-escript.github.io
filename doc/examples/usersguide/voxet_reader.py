##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

import os
from esys.downunder import CartesianReferenceSystem
from esys.escript import ReducedFunction, getMPIRankWorld

def readVoxet(domain, filename, voproperty=1, origin=None, fillValue=0.,
              referenceSystem=CartesianReferenceSystem()):
    """
    Reads a single property from a GOCAD Voxet file and returns a data
    object on the given domain with the property data.
    Restrictions:
    - Voxet origin in UVW space (i.e. AXIS_MIN) needs to be [0,0,0]
    - samples size must be 4 (float32) or 8 (float64)
    - data type must be IEEE
    - format must be RAW
    - domain resolution must be (approximately) a multiple of voxet resolution

    :param domain: the domain to use for data (must be a ripley domain)
    :type domain: `Domain`
    :param filename: Voxet header filename (usually ends in .vo)
    :type filename: ``string``
    :param voproperty: identifier of the property to read. Either the numeric
                       property ID, the property name, or the filename of the
                       property data.
    :type voproperty: ``int`` or ``string``
    :param origin: if supplied will override the Voxet origin as read from the
                   file.
    :type origin: ``list`` or ``tuple`` or ``None``
    :param fillValue: value to use for cells that are not covered by property
                      data (if applicable)
    :type fillValue: ``float``
    :param referenceSystem: coordinate system of domain. Used to scale vertical
                            axis accordingly
    :type referenceSystem: `ReferenceSystem`
    """
    from esys.ripley import readBinaryGrid, BYTEORDER_BIG_ENDIAN, DATATYPE_FLOAT32, DATATYPE_FLOAT64
    header=open(filename).readlines()
    if not header[0].startswith('GOCAD Voxet'):
        raise ValueError("Voxet header not found. Invalid Voxet file?!")
    NE=None
    axis_uvw=[None,None,None]
    axis_min=[0.,0.,0.]
    axis_max=[1.,1.,1.]
    # props[id]=[name,file,datatype]
    props={}
    for line in header:
        if line.startswith('AXIS_O '):
            if origin is None:
                origin=[float(i) for i in line.split()[1:4]]
        elif line.startswith('AXIS_U '):
            u=[float(i) for i in line.split()[1:4]]
            if (u[1] != 0) or (u[2] != 0):
                raise ValueError('This coordinate system is not supported')
            axis_uvw[0]=u[0]
        elif line.startswith('AXIS_V '):
            v=[float(i) for i in line.split()[1:4]]
            if (v[0] != 0) or (v[2] != 0):
                raise ValueError('This coordinate system is not supported')
            axis_uvw[1]=v[1]
        elif line.startswith('AXIS_W '):
            w=[float(i) for i in line.split()[1:4]]
            if (w[0] != 0) or (w[1] != 0):
                raise ValueError('This coordinate system is not supported')
            axis_uvw[2]=w[2]
        elif line.startswith('AXIS_MIN '):
            axis_min=[float(i) for i in line.split()[1:4]]
            if axis_min != [0,0,0]:
                raise ValueError('AXIS_MIN != [0,0,0] is not supported')
        elif line.startswith('AXIS_MAX '):
            axis_max=[float(i) for i in line.split()[1:4]]
        elif line.startswith('AXIS_N '):
            NE=[int(i) for i in line.split()[1:4]]
        elif line.startswith('PROPERTY '):
            propid=int(line.split()[1])
            if not propid in props:
                props[propid]=[None,None,None]
            props[propid][0]=line.split()[2].strip()
        elif line.startswith('PROP_ESIZE '):
            propid=int(line.split()[1])
            t=int(line.split()[2])
            if t==4:
                props[propid][2]=DATATYPE_FLOAT32
            elif t==8:
                props[propid][2]=DATATYPE_FLOAT64
            else:
                raise ValueError('Unsupported data size '+t)
        elif line.startswith('PROP_ETYPE '):
            t=line.split()[2].strip()
            if t != 'IEEE':
                raise ValueError('Unsupported data type '+t)
        elif line.startswith('PROP_FORMAT '):
            t=line.split()[2].strip()
            if t != 'RAW':
                raise ValueError('Unsupported data format '+t)
        elif line.startswith('PROP_OFFSET '):
            dataoffset=int(line.split()[2])
            if dataoffset != 0:
                raise ValueError('data offset != 0 not supported')
        elif line.startswith('PROP_FILE '):
            propid=int(line.split()[1])
            props[propid][1]=line.split()[2].strip()

    if (axis_uvw[0] is None) or (axis_uvw[1] is None) or (axis_uvw[2] is None)\
            or (NE is None) or (origin is None):
        raise ValueError('Could not determine data configuration. Invalid file?!')
    if len(props)==0:
        raise ValueError('No properties found.')

    # voxets have these conventions:
    # AXIS_N = number of samples (=cells!) in each dimension
    # AXIS_UVW * AXIS_MAX = voxet length in each dimension
    # AXIS_O = origin of voxet (cell centres!)
    # see also http://paulbourke.net/dataformats/gocad/gocad.pdf

    length = [axis_uvw[i]*axis_max[i] for i in range(3)]

    # modify length and origin to account for the fact that Voxet cells are
    # centred at the data points, i.e.:
    # BEFORE:                   AFTER:
    #
    #       O----length---->|        O------length------>|
    #      ___________________        ___________________
    #     | * | * | * | * | * |      | * | * | * | * | * |
    #      -------------------        -------------------

    for i in range(3):
        dz = length[i] / (NE[i]-1)
        origin[i] -= dz/2.
        length[i] += dz

    if referenceSystem.isCartesian():
        v_scale=1.
    else:
        v_scale=1./referenceSystem.getHeightUnit()

    origin[-1] = origin[-1]*v_scale

    # retrieve domain configuration so we know where to place the voxet data
    gridorigin, gridspacing, gridNE = domain.getGridParameters()

    # determine base location of this dataset within the domain
    first=[int((origin[i]-gridorigin[i])/gridspacing[i]) for i in range(domain.getDim())]

    # determine the resolution difference between domain and data.
    # If domain has twice the resolution we can double up the data etc.
    multiplier=[int(round((abs(length[i])/NE[i])/gridspacing[i])) for i in range(domain.getDim())]

    # NOTE: Depending on your data you might have to multiply your vertical
    # multiplier by 1000. to convert km in meters.
    #multiplier[-1] = int(multiplier[-1] * v_scale * 1000.)
    multiplier[-1] = int(multiplier[-1] * v_scale)

    datatype=None
    propfile=None
    for pid in props.keys():
        p=props[pid]
        if (isinstance(voproperty, int) and pid == voproperty) or \
           (isinstance(voproperty, str) and (p[0]==voproperty or p[1]==voproperty)):
            datatype=p[2]
            name=p[1]
            #remove quotes which GoCAD introduces for filenames with spaces
            if name.startswith('"') and name.endswith('"'):
                name=name[1:-1]
            propfile=os.path.join(os.path.dirname(filename), name)
            print("Voxet property file: %s"%propfile)
            break

    if propfile is None or datatype is None:
        raise ValueError("Invalid property "+str(voproperty))

    reverse=[0]*domain.getDim()
    if axis_uvw[-1] < 0:
        reverse[-1]=1

    print("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s, reverse=%s"%(str(first),str(NE),str(multiplier),str(reverse)))
    data=readBinaryGrid(propfile, ReducedFunction(domain), shape=(),
            fill=fillValue, byteOrder=BYTEORDER_BIG_ENDIAN,
            dataType=p[2], first=first, numValues=NE, multiplier=multiplier,
            reverse=reverse)

    return data


if __name__ == "__main__":
    try:
        from esys.ripley import Brick
        HAVE_RIPLEY = True
    except ImportError:
        HAVE_RIPLEY = False
        print("Ripley module not available")

    if HAVE_RIPLEY:
        from esys.escript import *
        from esys.escript.linearPDEs import Poisson
        from esys.weipa import saveSilo, saveVoxet
        
        import tempfile

        dom = Brick(l0=1.,l1=1.,n0=9, n1=9, n2=9)
        x = dom.getX()
        gammaD = whereZero(x[0])+whereZero(x[1])
        pde = Poisson(dom)
        q = gammaD
        pde.setValue(f=1, q=q)
        u = pde.getSolution()
        u=interpolate(u+dom.getX()[2], ReducedFunction(dom))
        print(u)
        if os.name == "nt":
            filename = os.environ["TEMP"]+os.path.sep+"temp.vo"
        else:
            filename = "/tmp/temp.vo"
        saveVoxet(filename, u=u)

        print("-------")
        dom = Brick(l0=1.,l1=1.,l2=4.,n0=18, n1=18, n2=36)
        v=readVoxet(dom, filename, 'u', fillValue=0.5)
        print(v)
        if getMPIRankWorld() == 0:
            os.remove(filename)
            os.remove(filename[:-3] + '_u')
        #saveSilo('/tmp/poisson', v=v)

