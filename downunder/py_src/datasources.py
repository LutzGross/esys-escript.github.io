
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['DataSource','UBCDataSource','ERSDataSource','SyntheticDataSource','SmoothAnomaly']

import logging
import numpy as np
import struct
from esys.escript import *
from esys.escript.linearPDEs import *
import esys.escript.unitsSI as U
try:
    from scipy.io.netcdf import netcdf_file
    __all__ += ['NetCDFDataSource']
except:
    pass

def LatLonToUTM(lon, lat, wkt_string=None):
    """
    Converts one or more longitude,latitude pairs to the corresponding x,y
    coordinates in the Universal Transverse Mercator projection.
    If wkt_string is not given or invalid or the gdal module is not available
    to convert the string, then the input values are assumed to be given in the
    Clarke 1866 projection.
    """

    # not really optimal: if pyproj is not installed we return the input
    # values without modification.
    try:
        import pyproj
    except:
        print("Warning, pyproj not available. Domain extents will be wrong")
        return lon,lat

    # determine UTM zone from the input data
    zone=int(np.median((np.floor((np.array(lon) + 180)/6) + 1) % 60))
    try:
        import osgeo.osr
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(wkt_string)
        p_src = pyproj.Proj(srs.ExportToProj4())
    except:
        p_src = pyproj.Proj('+proj=longlat +ellps=clrk66 +no_defs')
    p_dest = pyproj.Proj(proj='utm', zone=zone) # ellps?
    x,y=pyproj.transform(p_src, p_dest, lon, lat)
    return x,y

class DataSource(object):
    """
    A class that provides survey data for the inversion process.
    """
    # this is currently specific to gravity inversion and should be generalised

    def __init__(self):
        """
        """
        self._constrainBottom=False
        self._domain=None
        self._pad_l=0.1
        self._pad_h=0.1
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)

    def _addPadding(self, pad_l, pad_h, NE, l, origin):
        """
        Helper method that computes new number of elements, length and origin
        after adding padding to the input values.
        """
        DIM=len(NE)
        frac=[0.]*DIM
        # padding is applied to each side so multiply by 2 to get the total
        # amount of padding per dimension
        if pad_l>0 and pad_l<1:
            for i in xrange(DIM-1):
                frac[i]=2*pad_l
        elif pad_l>=1:
            for i in xrange(DIM-1):
                frac[i]=2*pad_l/float(NE[i])
        if pad_h>0 and pad_h<1:
            frac[DIM-1]=2*pad_h
        elif pad_h>=1:
            frac[DIM-1]=2*pad_h/(float(NE[DIM-1]))
        # calculate new number of elements
        NE_new=[int(NE[i]*(1+frac[i])) for i in xrange(DIM)]
        NEdiff=[NE_new[i]-NE[i] for i in xrange(DIM)]
        spacing=[l[i]/NE[i] for i in xrange(DIM)]
        l_new=[NE_new[i]*spacing[i] for i in xrange(DIM)]
        origin_new=[origin[i]-NEdiff[i]/2.*spacing[i] for i in xrange(DIM)]
        return NE_new, l_new, origin_new

    def _interpolateOnDomain(self, data, shape, origin, spacing, length):
        """
        Helper method that interpolates data arrays onto the domain.
        """
        dim=len(shape)
        arrays=np.zeros(((len(data[0])-dim),)+tuple(shape))
        for entry in data:
            index=()
            for i in xrange(dim):
                index=((entry[i]-origin[i])/spacing[i],)+index
            for i in xrange(arrays.shape[0]):
                arrays[i][index]=entry[dim+i]
        dom=self.getDomain()
        x=dom.getX()
        delta=[length[i]/(shape[dim-i-1]-1) for i in xrange(dim)]
        realorigin=[inf(x[i]) for i in xrange(dim)]
        res=[]
        for i in xrange(arrays.shape[0]):
            res.append(interpolateTable(arrays[i], x[:dim], realorigin, delta, 1e9))
        return res

    def setPadding(self, pad_l=0.1, pad_h=0.1):
        """
        Sets the amount of padding around the dataset. If pad_l/pad_h is >=1
        they are treated as number of elements to be added to the domain.
        If 0 < pad_l;pad_h < 1, the padding amount is relative.
        """
        self._pad_l=pad_l
        self._pad_h=pad_h

    def setConstrainBottom(self, constrain):
        """
        If `constrain` is True, then the density mask will be set to 1 in the
        padding area at the bottom of the domain. By default this area is
        unconstrained.
        """
        self._constrainBottom=constrain

    def getDomain(self):
        """
        Returns a domain that spans the data area plus padding.
        """
        if self._domain is None:
            self._domain=self._createDomain(self._pad_l, self._pad_h)
        return self._domain

    def getDensityMask(self):
        """
        Returns the density mask data object, where mask has value 1 on the
        padding area, 0 elsewhere.
        """
        raise NotImplementedError

    def getGravityAndStdDev(self):
        """
        Returns the gravity anomaly and standard deviation data objects as a
        tuple.
        """
        raise NotImplementedError

    def _createDomain(self, padding_l, padding_h):
        """
        creates and returns an escript domain that spans the entire area of
        available data plus a buffer zone.
        """
        raise NotImplementedError


##############################################################################
class UBCDataSource(DataSource):
    def __init__(self, domainclass, meshfile, gravfile, topofile=None):
        super(UBCDataSource,self).__init__()
        self._meshfile=meshfile
        self._gravfile=gravfile
        self._topofile=topofile
        self._domainclass=domainclass
        self._readMesh()

    def getDensityMask(self):
        #topodata=self._readTopography()
        #shape=[self.NE[1]+1, self.NE[0]+1]
        #mask=self._interpolateOnDomain(topodata, shape, self._origin, self._spacing, self._meshlen)
        #mask=wherePositive(self.getDomain().getX()[2]-mask[0])
        return self._mask

    def getGravityAndStdDev(self):
        gravlist=self._readGravity() # x,y,z,g,s
        shape=[self.NE[2]+1, self.NE[1]+1, self.NE[0]+1]
        g_and_sigma=self._interpolateOnDomain(gravlist, shape, self._origin, self._spacing, self._meshlen)
        return g_and_sigma[0]*[0,0,1], g_and_sigma[1]

    def _readMesh(self):
        meshdata=open(self._meshfile).readlines()
        numDataPoints=meshdata[0].split()
        origin=meshdata[1].split()
        self._numDataPoints=[int(x) for x in numDataPoints]
        self._origin=[float(x) for x in origin]
        self._spacing=[float(X.split('*')[1]) for X in meshdata[2:]]
        # vertical data is upside down
        self._origin[2]-=(self._numDataPoints[2]-1)*self._spacing[2]

    def _createDomain(self, padding_l, padding_h):
        NE=[self._numDataPoints[i]-1 for i in xrange(3)]
        l=[NE[i]*self._spacing[i] for i in xrange(3)]
        NE_new, l_new, origin_new = self._addPadding(padding_l, padding_h, \
                NE, l, self._origin)

        self._meshlen=l_new
        self.NE=NE_new
        self._origin=origin_new
        lo=[(self._origin[i], self._origin[i]+l_new[i]) for i in xrange(3)]
        NEdiff=[NE_new[i]-NE[i] for i in xrange(3)]
        try:
            dom=self._domainclass(*self.NE, l0=lo[0], l1=lo[1], l2=lo[2])
            x=dom.getX()-[self._origin[i]+NEdiff[i]/2.*self._spacing[i] for i in xrange(3)]
            mask=wherePositive(dom.getX()[2])

        except TypeError:
            dom=self._domainclass(*self.NE, l0=l_new[0], l1=l_new[1], l2=l_new[2])
            x=dom.getX()-[NEdiff[i]/2.*self._spacing[i] for i in xrange(3)]
            mask=wherePositive(x[2]+self._origin[2])

        if self._constrainBottom:
            M=3 # constrain bottom
        else:
            M=2 # do not constrain bottom
        for i in xrange(M):
            mask=mask + whereNegative(x[i]) + \
                    wherePositive(x[i]-l_new[i]+NEdiff[i]*self._spacing[i])
        self._mask=wherePositive(mask)

        self.logger.debug("Domain size: %d x %d x %d elements"%(self.NE[0],self.NE[1],self.NE[2]))
        self.logger.debug("     length: %g x %g x %g"%(l_new[0],l_new[1],l_new[2]))
        self.logger.debug("     origin: %g x %g x %g"%(origin_new[0],origin_new[1],origin_new[2]))

        return dom

    def _readTopography(self):
        f=open(self._topofile)
        n=int(f.readline())
        topodata=np.zeros((n,3))
        for i in xrange(n):
            x=f.readline().split()
            x[0]=float(x[0])
            x[1]=float(x[1])
            x[2]=float(x[2])
            topodata[i]=x
        f.close()
        return topodata

    def _readGravity(self):
        f=open(self._gravfile)
        n=int(f.readline())
        gravdata=np.zeros((n,5))
        for i in xrange(n):
            x=f.readline().split()
            x[0]=float(x[0]) # x
            x[1]=float(x[1]) # y
            x[2]=float(x[2]) # z
            x[3]=float(x[3]) # gravity anomaly in mGal
            x[4]=float(x[4]) # stddev
            # convert gravity anomaly units to m/s^2 and rescale error
            x[3]*=-1e-5
            x[4]*=1e-5
            gravdata[i]=x
        f.close()
        return gravdata

##############################################################################
class NetCDFDataSource(DataSource):
    def __init__(self, domainclass, gravfile, topofile=None, vertical_extents=(-40000,10000,26), alt_of_data=0.):
        """
        vertical_extents - (alt_min, alt_max, num_elements)
        alt_of_data - altitude of measurements
        """
        super(NetCDFDataSource,self).__init__()
        self._topofile=topofile
        self._gravfile=gravfile
        self._domainclass=domainclass
        self._determineExtents(vertical_extents)
        self._altOfData=alt_of_data

    def getDensityMask(self):
        #topodata=self._readTopography()
        #shape=self._numDataPoints[1::-1]
        #mask=self._interpolateOnDomain(topodata, shape, self._origin, self._spacing, self._meshlen)
        #mask=wherePositive(self.getDomain().getX()[2]-mask[0])
        #rho=fill*(1.-mask) + RHO_AIR*mask
        return self._mask

    def getGravityAndStdDev(self):
        gravlist=self._readGravity() # x,y,z,g,s
        shape=[self.NE[2]+1, self.NE[1]+1, self.NE[0]+1]
        g_and_sigma=self._interpolateOnDomain(gravlist, shape, self._origin, self._spacing, self._meshlen)
        return g_and_sigma[0]*[0,0,1], g_and_sigma[1]

    def _determineExtents(self, ve):
        f=netcdf_file(self._gravfile, 'r')
        NX=0
        for n in ['lon','longitude','x']:
            if n in f.dimensions:
                NX=f.dimensions[n]
                break
        if NX==0:
            raise RuntimeError("Could not determine extents of data")
        NY=0
        for n in ['lat','latitude','y']:
            if n in f.dimensions:
                NY=f.dimensions[n]
                break
        if NY==0:
            raise RuntimeError("Could not determine extents of data")

        # find longitude and latitude variables
        lon_name=None
        for n in ['lon','longitude']:
            if n in f.variables:
                lon_name=n
                longitude=f.variables.pop(n)
                break
        if lon_name is None:
            raise RuntimeError("Could not determine longitude variable")
        lat_name=None
        for n in ['lat','latitude']:
            if n in f.variables:
                lat_name=n
                latitude=f.variables.pop(n)
                break
        if lat_name is None:
            raise RuntimeError("Could not determine latitude variable")

        # try to figure out gravity variable name
        grav_name=None
        if len(f.variables)==1:
            grav_name=f.variables.keys()[0]
        else:
            for n in f.variables.keys():
                dims=f.variables[n].dimensions
                if (lat_name in dims) and (lon_name in dims):
                    grav_name=n
                    break
        if grav_name is None:
            raise RuntimeError("Could not determine gravity variable")

        # see if there is a wkt string to convert coordinates
        try:
            wkt_string=f.variables[grav_name].esri_pe_string
        except:
            wkt_string=None

        # we don't trust actual_range & geospatial_lon_min/max since subset
        # data does not seem to have these fields updated it seems.
        # Getting min/max from the arrays is obviously not very efficient...
        #lon_range=longitude.actual_range
        #lat_range=latitude.actual_range
        #lon_range=[f.geospatial_lon_min,f.geospatial_lon_max]
        #lat_range=[f.geospatial_lat_min,f.geospatial_lat_max]
        lon_range=longitude.data.min(),longitude.data.max()
        lat_range=latitude.data.min(),latitude.data.max()
        lon_range,lat_range=LatLonToUTM(lon_range, lat_range, wkt_string)
        origin=[lon_range[0],lat_range[0],ve[0]]
        lengths=[lon_range[1]-lon_range[0], lat_range[1]-lat_range[0],ve[1]-ve[0]]

        f.close()

        self._numDataPoints=[NX, NY, ve[2]]
        self._origin=origin
        # we are rounding to avoid interpolateOnDomain issues
        self._spacing=[np.round(lengths[i]/(self._numDataPoints[i]-1)) for i in xrange(3)]
        self._meshlen=[(self._numDataPoints[i]-1)*self._spacing[i] for i in xrange(3)]
        self._wkt_string=wkt_string
        self._lon=lon_name
        self._lat=lat_name
        self._grv=grav_name

    def _createDomain(self, padding_l, padding_h):
        NE=[self._numDataPoints[i]-1 for i in xrange(3)]
        l=self._meshlen
        NE_new, l_new, origin_new = self._addPadding(padding_l, padding_h, \
                NE, l, self._origin)

        self._meshlen=l_new
        self.NE=NE_new
        self._dataorigin=self._origin
        self._origin=origin_new
        lo=[(self._origin[i], self._origin[i]+l_new[i]) for i in xrange(3)]
        NEdiff=[NE_new[i]-NE[i] for i in xrange(3)]
        try:
            dom=self._domainclass(*self.NE, l0=lo[0], l1=lo[1], l2=lo[2])
            # ripley may adjust NE and length
            self._meshlen=[sup(dom.getX()[i])-inf(dom.getX()[i]) for i in xrange(3)]
            self.NE=[self._meshlen[i]/self._spacing[i] for i in xrange(3)]
            x=dom.getX()-[self._origin[i]+NEdiff[i]/2.*self._spacing[i] for i in xrange(3)]
            mask=wherePositive(dom.getX()[2])

        except TypeError:
            dom=self._domainclass(*self.NE, l0=l_new[0], l1=l_new[1], l2=l_new[2])
            x=dom.getX()-[NEdiff[i]/2.*self._spacing[i] for i in xrange(3)]
            mask=wherePositive(x[2]+self._origin[2])

        if self._constrainBottom:
            M=3 # constrain bottom
        else:
            M=2 # do not constrain bottom
        for i in xrange(M):
            mask=mask + whereNegative(x[i]) + \
                    wherePositive(x[i]-l_new[i]+NEdiff[i]*self._spacing[i])
        self._mask=wherePositive(mask)

        self.logger.debug("Domain size: %d x %d x %d elements"%(self.NE[0],self.NE[1],self.NE[2]))
        self.logger.debug("     length: %s x %s x %s"%(self._meshlen[0],self._meshlen[1],self._meshlen[2]))
        self.logger.debug("     origin: %s x %s x %s"%(self._origin[0],self._origin[1],self._origin[2]))

        return dom

    def _readTopography(self):
        f=netcdf_file(self._topofile, 'r')
        lon=None
        for n in ['lon','longitude']:
            if n in f.variables:
                lon=f.variables[n][:]
                break
        if lon is None:
            raise RuntimeError("Could not determine longitude variable")
        lat=None
        for n in ['lat','latitude']:
            if n in f.variables:
                lat=f.variables[n][:]
                break
        if lat is None:
            raise RuntimeError("Could not determine latitude variable")
        alt=None
        for n in ['altitude','alt']:
            if n in f.variables:
                alt=f.variables[n][:]
                break
        if alt is None:
            raise RuntimeError("Could not determine altitude variable")

        topodata=np.column_stack((lon,lat,alt))
        f.close()
        return topodata

    def _readGravity(self):
        f=netcdf_file(self._gravfile, 'r')
        #lon=f.variables[self._lon][:]
        #lat=f.variables[self._lat][:]
        NE=[self._numDataPoints[i]-1 for i in xrange(2)]
        lon=np.linspace(self._dataorigin[0], self._dataorigin[0]+NE[0]*self._spacing[0], NE[0]+1)
        lat=np.linspace(self._dataorigin[1], self._dataorigin[1]+NE[1]*self._spacing[1], NE[1]+1)
        lon,lat=np.meshgrid(lon,lat)
        grav=f.variables[self._grv][:]
        f.close()
        lon=lon.flatten()
        lat=lat.flatten()
        grav=grav.flatten()
        alt=self._altOfData*np.ones(grav.shape)
        # error value is an assumption
        try:
            missing=f.variables[self._grv].missing_value
            err=np.where(grav==missing, 0.0, 20.0)
        except:
            err=20.0*np.ones(lon.shape)
        # convert units
        err=1e-6*err
        grav=1e-6*grav
        gravdata=np.column_stack((lon,lat,alt,grav,err))
        return gravdata

##############################################################################
class ERSDataSource(DataSource):
    """
    Data Source for ER Mapper raster data.
    Note that this class only accepts a very specific type of ER Mapper data
    input and will raise an exception if other data is found.
    """
    def __init__(self, domainclass, headerfile, datafile=None, vertical_extents=(-40000,10000,26), alt_of_data=0.):
        """
        headerfile - usually ends in .ers
        datafile - usually has the same name as the headerfile without '.ers'
        """
        super(ERSDataSource,self).__init__()
        self._headerfile=headerfile
        if datafile is None:
            self._datafile=headerfile[:-4]
        else:
            self._datafile=datafile
        self._domainclass=domainclass
        self._readHeader(vertical_extents)
        self._altOfData=alt_of_data

    def _readHeader(self, ve):
        metadata=open(self._headerfile, 'r').readlines()
        # parse metadata
        start=-1
        for i in xrange(len(metadata)):
            if metadata[i].strip() == 'DatasetHeader Begin':
                start=i+1
        if start==-1:
            raise RuntimeError('Invalid ERS file (DatasetHeader not found)')

        md_dict={}
        section=[]
        for i in xrange(start, len(metadata)):
            line=metadata[i]
            if line[-6:].strip() == 'Begin':
                section.append(line[:-6].strip())
            elif line[-4:].strip() == 'End':
                if len(section)>0:
                    section.pop()
            else:
                vals=line.split('=')
                if len(vals)==2:
                    key = vals[0].strip()
                    value = vals[1].strip()
                    fullkey='.'.join(section+[key])
                    md_dict[fullkey]=value

        try:
            if md_dict['RasterInfo.CellType'] != 'IEEE4ByteReal':
                raise RuntimeError('Unsupported data type '+md_dict['RasterInfo.CellType'])
        except KeyError:
            self.logger.warn("Cell type not specified. Assuming IEEE4ByteReal.")

        try:
            NX = int(md_dict['RasterInfo.NrOfCellsPerLine'])
            NY = int(md_dict['RasterInfo.NrOfLines'])
        except:
            raise RuntimeError("Could not determine extents of data")

        try:
            maskval = float(md_dict['RasterInfo.NullCellValue'])
        except:
            raise RuntimeError("Could not determine NullCellValue")

        try:
            spacingX = float(md_dict['RasterInfo.CellInfo.Xdimension'])
            spacingY = float(md_dict['RasterInfo.CellInfo.Ydimension'])
        except:
            raise RuntimeError("Could not determine cell dimensions")

        try:
            originX = float(md_dict['RasterInfo.RegistrationCoord.Eastings'])
            originY = float(md_dict['RasterInfo.RegistrationCoord.Northings'])
        except:
            self.logger.warn("Could not determine coordinate origin")
            originX,originY = 0.0, 0.0

        # test data sets have origin in top-left corner so y runs top-down
        self._dataorigin=[originX, originY]
        originY-=(NY-1)*spacingY
        self._maskval=maskval
        self._spacing=[spacingX, spacingY, (ve[1]-ve[0])/(ve[2]-1)]
        self._numDataPoints = [NX, NY, ve[2]]
        self._origin = [originX, originY, ve[0]]
        self._meshlen=[self._numDataPoints[i]*self._spacing[i] for i in xrange(3)]

    def getDensityMask(self):
        return self._mask

    def getGravityAndStdDev(self):
        gravlist=self._readGravity() # x,y,z,g,s
        shape=[self.NE[2]+1, self.NE[1]+1, self.NE[0]+1]
        g_and_sigma=self._interpolateOnDomain(gravlist, shape, self._origin, self._spacing, self._meshlen)
        return g_and_sigma[0]*[0,0,1], g_and_sigma[1]

    def _createDomain(self, padding_l, padding_h):
        NE=[self._numDataPoints[i]-1 for i in xrange(3)]
        l=[NE[i]*self._spacing[i] for i in xrange(3)]
        NE_new, l_new, origin_new = self._addPadding(padding_l, padding_h, \
                NE, l, self._origin)

        self._meshlen=l_new
        self.NE=NE_new
        self._origin=origin_new
        lo=[(self._origin[i], self._origin[i]+l_new[i]) for i in xrange(3)]
        NEdiff=[NE_new[i]-NE[i] for i in xrange(3)]
        try:
            dom=self._domainclass(*self.NE, l0=lo[0], l1=lo[1], l2=lo[2])
            x=dom.getX()-[self._origin[i]+NEdiff[i]/2.*self._spacing[i] for i in xrange(3)]
            mask=wherePositive(dom.getX()[2])

        except TypeError:
            dom=self._domainclass(*self.NE, l0=l_new[0], l1=l_new[1], l2=l_new[2])
            x=dom.getX()-[NEdiff[i]/2.*self._spacing[i] for i in xrange(3)]
            mask=wherePositive(x[2]+self._origin[2])

        if self._constrainBottom:
            M=3 # constrain bottom
        else:
            M=2 # do not constrain bottom
        for i in xrange(M):
            mask=mask + whereNegative(x[i]) + \
                    wherePositive(x[i]-l_new[i]+NEdiff[i]*self._spacing[i])
        self._mask=wherePositive(mask)

        self.logger.debug("Domain size: %d x %d x %d elements"%(self.NE[0],self.NE[1],self.NE[2]))
        self.logger.debug("     length: %g x %g x %g"%(l_new[0],l_new[1],l_new[2]))
        self.logger.debug("     origin: %g x %g x %g"%(origin_new[0],origin_new[1],origin_new[2]))

        return dom

    def _readGravity(self):
        f=open(self._datafile, 'r')
        n=self._numDataPoints[0]*self._numDataPoints[1]
        grav=[]
        err=[]
        for i in xrange(n):
            v=struct.unpack('f',f.read(4))[0]
            grav.append(v)
            if abs((self._maskval-v)/v) < 1e-6:
                err.append(0.)
            else:
                err.append(1.)
        f.close()

        NE=[self._numDataPoints[i] for i in xrange(2)]
        x=self._dataorigin[0]+np.arange(start=0, stop=NE[0]*self._spacing[0], step=self._spacing[0])
        # test data sets have origin in top-left corner so y runs top-down
        y=self._dataorigin[1]-np.arange(start=0, stop=NE[1]*self._spacing[1], step=self._spacing[1])
        x,y=np.meshgrid(x,y)
        x,y=x.flatten(),y.flatten()
        alt=self._altOfData*np.ones((n,))

        # convert units
        err=2e-5*np.array(err)
        grav=1e-5*np.array(grav)
        gravdata=np.column_stack((x,y,alt,grav,err))
        return gravdata

##############################################################################
class SourceFeature(object):
    """
    A feature adds a density distribution to (parts of) a domain of a synthetic
    data source, for example a layer of a specific rock type or a simulated
    ore body.
    """
    def getDensity(self):
        """
        Returns the density for the area covered by mask. It can be constant
        or a data object with spatial dependency.
        """
        raise NotImplementedError
    def getMask(self, x):
        """
        Returns the mask of the area of interest for this feature. That is,
        mask is non-zero where the density returned by getDensity() should be
        applied, zero elsewhere.
        """
        raise NotImplementedError

class SmoothAnomaly(SourceFeature):
    def __init__(self, lx, ly, lz, x, y, depth, rho_inner, rho_outer):
        self.x=x
        self.y=y
        self.lx=lx
        self.ly=ly
        self.lz=lz
        self.depth=depth
        self.rho_inner=rho_inner
        self.rho_outer=rho_outer
        self.rho=None
    
    def getDensity(self):
        return self.rho

    def getMask(self, x):
        DIM=x.getDomain().getDim()
        m=whereNonNegative(x[DIM-1]-(sup(x[DIM-1])-self.depth-self.lz/2)) * whereNonPositive(x[DIM-1]-(sup(x[DIM-1])-self.depth+self.lz/2)) \
            *whereNonNegative(x[0]-(self.x-self.lx/2)) * whereNonPositive(x[0]-(self.x+self.lx/2)) 
        if DIM>2:
            m*=whereNonNegative(x[1]-(self.y-self.ly/2)) * whereNonPositive(x[1]-(self.y+self.ly/2))
        if self.rho is None:
            alpha=-log(abs(self.rho_outer/self.rho_inner))*4
            rho=exp(-alpha*((x[0]-self.x)/self.lx)**2)
            rho=rho*exp(-alpha*((x[DIM-1]-(sup(x[DIM-1])-self.depth))/self.lz)**2)
            self.rho=maximum(abs(self.rho_outer), abs(self.rho_inner*rho))
            if self.rho_inner<0: self.rho=-self.rho
        return m	 

##############################################################################
class SyntheticDataSource(DataSource):
    def __init__(self, DIM, NE, l, h, features):
        super(SyntheticDataSource,self).__init__()
        self._features = features
        self.DIM=DIM
        self.NE=NE
        self.l=l
        self.h=h

    def _createDomain(self, padding_l, padding_h):
        NE_H=self.NE
        NE_L=int((self.l/self.h)*NE_H+0.5)
        l=[self.l]*(self.DIM-1)+[self.h]
        NE=[NE_L]*(self.DIM-1)+[NE_H]
        origin=[0.]*self.DIM
        NE_new, l_new, origin_new = self._addPadding(padding_l, padding_h, \
                NE, l, origin)

        self.NE=NE_new
        self.l=l_new[0]
        self.h=l_new[self.DIM-1]

        if self.DIM==2:
            from esys.finley import Rectangle
            dom = Rectangle(n0=NE_new[0], n1=NE_new[1], l0=l_new[0], l1=l_new[1])
            self._x = dom.getX() + origin_new
            self.logger.debug("Domain size: %d x %d elements"%(NE_new[0], NE_new[1]))
            self.logger.debug("     length: %g x %g"%(l_new[0],l_new[1]))
            self.logger.debug("     origin: %g x %g"%(origin_new[0],origin_new[1]))
        else:   
            from esys.finley import Brick
            dom = Brick(n0=NE_new[0], n1=NE_new[1], n2=NE_new[2], l0=l_new[0], l1=l_new[1], l2=l_new[2])
            self._x = dom.getX() + origin_new
            self.logger.debug("Domain size: %d x %d x %d elements"%(self.NE[0],self.NE[1],self.NE[2]))
            self.logger.debug("     length: %g x %g x %g"%(l_new[0],l_new[1],l_new[2]))
            self.logger.debug("     origin: %g x %g x %g"%(origin_new[0],origin_new[1],origin_new[2]))

        dz=l_new[self.DIM-1]/NE_new[self.DIM-1]
        self._g_mask=wherePositive(dom.getX()[0]-origin_new[0]) \
                * whereNegative(dom.getX()[0]-(l_new[0]-origin_new[0])) \
                * whereNonNegative(dom.getX()[self.DIM-1]-(l_new[self.DIM-1]+origin_new[self.DIM-1])) \
                * whereNonPositive(dom.getX()[self.DIM-1]-(l_new[self.DIM-1]+(origin_new[self.DIM-1]+dz)))
        self._mask=whereNegative(self._x[self.DIM-1]) + \
                wherePositive(self._x[self.DIM-1]-l[self.DIM-1])
        for i in xrange(self.DIM-1):
            self._mask=self._mask + whereNegative(self._x[i]) + \
                    wherePositive(self._x[i]-l[i])
        self._mask=wherePositive(self._mask)

        rho_ref=0.
        for f in self._features:
            m=f.getMask(self._x)
            rho_ref = rho_ref * (1-m) + f.getDensity() * m
        self._rho=rho_ref

        return dom

    def getDensityMask(self):
        return self._mask

    def getReferenceDensity(self):
        return self._rho

    def getGravityAndStdDev(self):
        pde=LinearSinglePDE(self.getDomain())
        G=6.6742e-11*U.m**3/(U.kg*U.sec**2)
        m_psi_ref=0.
        for i in range(self.DIM):
            m_psi_ref=m_psi_ref + whereZero(self._x[i]-inf(self._x[i])) \
                    + whereZero(self._x[i]-sup(self._x[i]))

        pde.setValue(A=kronecker(self.getDomain()), Y=-4*np.pi*G*self._rho, q=m_psi_ref)
        pde.setSymmetryOn()
        psi_ref=pde.getSolution()
        del pde
        g=-grad(psi_ref)
        sigma=self._g_mask
        return g,sigma

