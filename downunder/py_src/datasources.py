
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

"""Data readers/providers for inversions"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['DataSource', 'ErMapperData', 'NumpyData', \
        'SyntheticDataBase', 'SyntheticFeatureData', 'SyntheticData',
        'SmoothAnomaly', 'SeismicSource']

import logging
import numpy as np
import tempfile
import esys.escript as es
from esys.escript import ReducedFunction, FunctionOnBoundary, Scalar
from esys.escript import unitsSI as U
from esys.escript.linearPDEs import LinearSinglePDE
import esys.escript.util as esu
from .coordinates import ReferenceSystem,  CartesianReferenceSystem

HAS_RIPLEY = True
try:
    from esys.ripley import *
except ImportError as e:
    HAS_RIPLEY = False

try:
    from scipy.io.netcdf import netcdf_file
    __all__ += ['NetCdfData']
    HAVE_SCIPY_NETCDF = True
except:
    HAVE_SCIPY_NETCDF = False
    pass

try:
    import pyproj
    HAVE_PYPROJ=True
except:
    HAVE_PYPROJ=False

try:
    import warnings
    with warnings.catch_warnings(record=True):
       warnings.simplefilter("ignore", category=DeprecationWarning)
       import osgeo.osr
    HAVE_GDAL=True
except ImportError:
    HAVE_GDAL=False


def getUTMZone(lon, lat, wkt_string=None):
    """
    """

    logger = logging.getLogger('inv.datasources.getUTMZone')
    zone = 0

    nplon=np.array(lon)
    nplat=np.array(lat)
    if np.abs(nplon).max()>360.0 or np.abs(nplat).max()>180.0:
        if HAVE_GDAL and (wkt_string is not None):
            srs = osgeo.osr.SpatialReference()
            result=srs.ImportFromWkt(wkt_string)
            if result==0:
                zone = srs.GetUTMZone()
    else:
        # determine UTM zone from the input data
        zone = int(np.median((np.floor((nplon + 180)/6) + 1) % 60))

    logger.debug("Determined UTM zone %d."%zone)
    return zone

def LatLonToUTM(lon, lat, wkt_string=None):
    """
    Converts one or more longitude,latitude pairs to the corresponding x,y
    coordinates in the Universal Transverse Mercator projection.

    :note: The ``pyproj`` module is required unless the input coordinates are
           determined to be already projected. If it is not found an exception
           is raised.
    :note: If `wkt_string` is not given or invalid or the ``gdal`` module is
           not available to convert the string, then the input values are
           assumed to be using the Clarke 1866 ellipsoid.

    :param lon: longitude value(s)
    :type lon: ``float``, ``list``, ``tuple``, or ``numpy.array``
    :param lat: latitude value(s)
    :type lat: ``float``, ``list``, ``tuple``, or ``numpy.array``
    :param wkt_string: Well-known text (WKT) string describing the coordinate
                       system used. The ``gdal`` module is used to convert
                       the string to the corresponding Proj4 string.
    :type wkt_string: ``str``
    :rtype: ``tuple``
    """

    logger = logging.getLogger('inv.datasources.LatLonToUTM')

    nplon=np.array(lon)
    nplat=np.array(lat)
    zone = getUTMZone(nplon, nplat, wkt_string)

    if np.abs(nplon).max()>360.0 or np.abs(nplat).max()>180.0:
        logger.debug('Coordinates appear to be projected. Passing through.')
        return lon,lat,zone

    logger.debug('Need to project coordinates.')
    if not HAVE_PYPROJ:
        logger.error("Cannot import pyproj! Exiting.")
        raise ImportError("In order to perform coordinate transformations on "
           "the data you are using the 'pyproj' Python module is required but "
           "was not found. Please install the module and try again.")

    p_src=None
    if HAVE_GDAL and (wkt_string is not None) and len(wkt_string)>0:
        srs = osgeo.osr.SpatialReference()
        result=srs.ImportFromWkt(wkt_string)
        try:
            p_src = pyproj.Proj(srs.ExportToProj4())
        except RuntimeError as e:
            logger.warning('pyproj returned exception: %s [wkt=%s]'%(e,wkt_string))

    if p_src is None:
        if HAVE_GDAL:
            reason="no wkt string provided."
        else:
            reason="the gdal python module not available."
        logger.warning("Assuming lon/lat coordinates on Clarke 1866 ellipsoid since "+reason)
        p_src = pyproj.Proj('+proj=longlat +ellps=clrk66 +no_defs')

    # check for hemisphere
    if np.median(nplat) < 0.:
        south='+south '
    else:
        south=''
    p_dest = pyproj.Proj('+proj=utm +zone=%d %s+units=m +ellps=WGS84'%(zone,south))
    x,y=pyproj.transform(p_src, p_dest, lon, lat)
    return x,y,zone

class DataSource(object):
    """
    A class that provides survey data for the inversion process.
    This is an abstract base class that implements common functionality.
    Methods to be overwritten by subclasses are marked as such.
    This class assumes 2D data which is mapped to a slice of a 3D domain.
    For other setups override the methods as required.
    """

    GRAVITY, MAGNETIC, ACOUSTIC, MT = list(range(4))

    def __init__(self, reference_system=None, tags=[]):
        """
        Constructor. Sets some defaults and initializes logger.
        
        :param tags: a list of tags associated with the data set.
        :type tags: ``list`` of almost any type (typically `str`) 
        :param reference_system: the reference coordinate system
        :type reference_system: ``None`` or `ReferenceSystem`
        """
        if not isinstance(tags ,list):
            raise ValueError("tags argument must be a list.")
        self.__tags=tags
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)
        self.__subsampling_factor=1
        if not reference_system:
             self.__reference_system = CartesianReferenceSystem()
        else:
             self.__reference_system = reference_system

        if self.__reference_system.isCartesian():
            self.__v_scale=1.
        else:
            self.__v_scale=1./self.getReferenceSystem().getHeightUnit()

    def getTags(self):
        """
        returns the list of tags
        
        :rtype: ``list``
        """
        return self.__tags

    def hasTag(self, tag):
        """
        returns true if the data set has tag ``tag``
        
        :rtype: ``bool``
        """
        return tag in self.__tags
                    
             
    def getReferenceSystem(self):
        """
        returns the reference coordinate system
        
        :rtype: `ReferenceSystem`
        """
        return self.__reference_system

    def getHeightScale(self):
        """
        returns the height scale factor to convert from meters to the
        appropriate units of the reference system used.

        :rtype: ``float``
        """
        return self.__v_scale

    def getDataExtents(self):
        """
        returns a tuple of tuples ``( (x0, y0), (nx, ny), (dx, dy) )``, where

        - ``x0``, ``y0`` = coordinates of data origin
        - ``nx``, ``ny`` = number of data points in x and y
        - ``dx``, ``dy`` = spacing of data points in x and y

        This method must be implemented in subclasses.
        """
        raise NotImplementedError

    def getDataType(self):
        """
        Returns the type of survey data managed by this source.
        Subclasses must return `GRAVITY` or `MAGNETIC` or `ACOUSTIC` as appropriate.
        """
        raise NotImplementedError

    def getSurveyData(self, domain, origin, NE, spacing):
        """
        This method is called by the `DomainBuilder` to retrieve the survey
        data as `Data` objects on the given domain.

        Subclasses should return one or more `Data` objects with survey data
        interpolated on the given `escript` domain. The exact return type
        depends on the type of data.

        :param domain: the escript domain to use
        :type domain: `esys.escript.Domain`
        :param origin: the origin coordinates of the domain
        :type origin: ``tuple`` or ``list``
        :param NE: the number of domain elements in each dimension
        :type NE: ``tuple`` or ``list``
        :param spacing: the cell sizes (node spacing) in the domain
        :type spacing: ``tuple`` or ``list``
        """
        raise NotImplementedError

    def getUtmZone(self):
        """
        All data source coordinates are converted to UTM (Universal Transverse
        Mercator) in order to have useful domain extents. Subclasses should
        implement this method and return the UTM zone number of the projected
        coordinates.
        """
        raise NotImplementedError

    def setSubsamplingFactor(self, f):
        """
        Sets the data subsampling factor (default=1).

        The factor is applied in all dimensions. For example a 2D dataset
        with 300 x 150 data points will be reduced to 150 x 75 when a
        subsampling factor of 2 is used.
        This becomes important when adding data of varying resolution to
        a `DomainBuilder`.
        """
        self.__subsampling_factor=f

    def getSubsamplingFactor(self):
        """
        Returns the subsampling factor that was set via `setSubsamplingFactor`
        (see there).
        """
        return self.__subsampling_factor


##############################################################################
class ErMapperData(DataSource):
    """
    Data Source for ER Mapper raster data.
    Note that this class only accepts a very specific type of ER Mapper data
    input and will raise an exception if other data is found.
    """
    def __init__(self, data_type, headerfile, datafile=None, altitude=0.,
                 error=None, scale_factor=None, null_value=None,
                 reference_system=None):
        """
        :param data_type: type of data, must be `GRAVITY` or `MAGNETIC`
        :type data_type: ``int``
        :param headerfile: ER Mapper header file (usually ends in .ers)
        :type headerfile: ``str``
        :param datafile: ER Mapper binary data file name. If not supplied the
                         name of the header file without '.ers' is assumed
        :type datafile: ``str``
        :param altitude: altitude of measurements above ground in meters
        :type altitude: ``float``
        :param error: constant value to use for the data uncertainties.
                      If a value is supplied, it is scaled by the same factor
                      as the measurements. If not provided the error is
                      assumed to be 2 units for all measurements (i.e. 0.2
                      mGal and 2 nT for gravity and magnetic, respectively)
        :type error: ``float``
        :param scale_factor: the measurements and error values are scaled by
                             this factor. By default, gravity data is assumed
                             to be given in 1e-6 m/s^2 (0.1 mGal), while
                             magnetic data is assumed to be in 1e-9 T (1 nT).
        :type scale_factor: ``float``
        :param null_value: value that is used in the file to mark undefined
                           areas. This information is usually included in the
                           file.
        :type null_value: ``float``
        :param reference_system: reference coordinate system to be used.
                                 For a Cartesian reference (default) the
                                 appropriate UTM transformation is applied.
        :type reference_system: `ReferenceSystem`
        :note: consistence in the reference coordinate system and the reference
               coordinate system used in the data source is not checked.
        """
        super(ErMapperData, self).__init__(reference_system, [ headerfile ] )
        self.__headerfile=headerfile
        if datafile is None:
            self.__datafile=headerfile[:-4]
        else:
            self.__datafile=datafile
        self.__altitude=altitude
        self.__data_type=data_type
        self.__utm_zone = None
        self.__scale_factor = scale_factor
        self.__null_value = null_value
        self.__error_value = error
        self.__readHeader()

    def __readHeader(self):
        self.logger.debug("Checking Data Source: %s (header: %s)"%(self.__datafile, self.__headerfile))
        metadata=open(self.__headerfile, 'r').readlines()
        start=-1
        for i in range(len(metadata)):
            if metadata[i].strip() == 'DatasetHeader Begin':
                start=i+1
        if start==-1:
            raise RuntimeError('Invalid ER Mapper header file ("DatasetHeader" not found)')

        # parse header file filling dictionary of found values
        md_dict={}
        section=[]
        for i in range(start, len(metadata)):
            line=metadata[i].strip()
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

        # check that the data format/type is supported
        try:
            if md_dict['ByteOrder'] != 'LSBFirst':
                raise RuntimeError('Unsupported byte order '+md_dict['ByteOrder'])
        except KeyError:
            self.logger.warning("Byte order not specified. Assuming LSB first.")

        try:
            if md_dict['DataType'] != 'Raster':
                raise RuntimeError('Unsupported data type '+md_dict['DataType'])
        except KeyError:
            self.logger.warning("Data type not specified. Assuming raster data.")

        try:
            if md_dict['RasterInfo.CellType'] == 'IEEE4ByteReal':
                self.__celltype = DATATYPE_FLOAT32
            elif md_dict['RasterInfo.CellType'] == 'IEEE8ByteReal':
                self.__celltype = DATATYPE_FLOAT64
            elif md_dict['RasterInfo.CellType'] == 'Signed32BitInteger':
                self.__celltype = DATATYPE_INT32
            else:
                raise RuntimeError('Unsupported data type '+md_dict['RasterInfo.CellType'])
        except KeyError:
            self.logger.warning("Cell type not specified. Assuming IEEE4ByteReal.")
            self.__celltype = DATATYPE_FLOAT32

        try:
            fileOffset = int(md_dict['HeaderOffset'])
        except:
            fileOffset = 0
        if fileOffset > 0:
            raise RuntimeError("ER Mapper data with header offset >0 not supported.")

        # now extract required information
        try:
            NX = int(md_dict['RasterInfo.NrOfCellsPerLine'])
            NY = int(md_dict['RasterInfo.NrOfLines'])
        except:
            raise RuntimeError("Could not determine extents of data")

        ### mask/null value
        # note that NaN is always filtered out in ripley
        if self.__null_value is None:
            try:
                self.__null_value = float(md_dict['RasterInfo.NullCellValue'])
            except:
                self.logger.debug("Could not determine null value, using default.")
                self.__null_value = 99999
        elif not isinstance(self.__null_value,float) and not isinstance(self.__null_value,int):
            raise TypeError("Invalid type of null_value parameter")

        try:
            spacingX = float(md_dict['RasterInfo.CellInfo.Xdimension'])
            spacingY = float(md_dict['RasterInfo.CellInfo.Ydimension'])
        except:
            raise RuntimeError("Could not determine cell dimensions")

        try:
            if md_dict['CoordinateSpace.CoordinateType']=='EN':
                originX = float(md_dict['RasterInfo.RegistrationCoord.Eastings'])
                originY = float(md_dict['RasterInfo.RegistrationCoord.Northings'])
            elif md_dict['CoordinateSpace.CoordinateType']=='LL':
                originX = float(md_dict['RasterInfo.RegistrationCoord.Longitude'])
                originY = float(md_dict['RasterInfo.RegistrationCoord.Latitude'])
            else:
                raise RuntimeError("Unknown CoordinateType")
        except:
            self.logger.warning("Could not determine coordinate origin. Setting to (0.0, 0.0)")
            originX,originY = 0.0, 0.0

        # data sets have origin in top-left corner so y runs top-down and
        # we need to flip accordingly
        originY-=NY*spacingY

        if 'GEODETIC' in md_dict['CoordinateSpace.Projection']:
            # it appears we have lat/lon coordinates so need to convert
            # origin and spacing. Try using gdal to get the wkt if available:
            try:
                from osgeo import gdal
                ds=gdal.Open(self.__headerfile)
                wkt=str(ds.GetProjection())
            except:
                wkt='GEOGCS["GEOCENTRIC DATUM of AUSTRALIA",DATUM["GDA94",SPHEROID["GRS80",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
                self.logger.warning('GDAL not available or file read error, assuming GDA94 data')
            if self.getReferenceSystem().isCartesian():
                originX_UTM,originY_UTM,zone = LatLonToUTM(originX, originY, wkt)
                op1X,op1Y,_ = LatLonToUTM(originX+spacingX, originY+spacingY, wkt)
                # we are rounding to avoid interpolation issues
                spacingX = np.round(op1X-originX_UTM)
                spacingY = np.round(op1Y-originY_UTM)
                originX = np.round(originX_UTM)
                originY = np.round(originY_UTM)
                self.__utm_zone = zone
            else:
                op1X, op1Y = originX+spacingX, originY+spacingY
                spacingX = np.round(op1X-originX,5)
                spacingY = np.round(op1Y-originY,5)
                originX = np.round(originX,5)
                originY = np.round(originY,5)

        self.__dataorigin=[originX, originY]
        self.__delta = [spacingX, spacingY]
        self.__nPts = [NX, NY]
        self.__origin = [originX, originY]
        ### scale factor
        if self.__scale_factor is None:
            if self.__data_type == self.GRAVITY:
                self.logger.info("Assuming gravity data scale is 1e-6 m/s^2.")
                self.__scale_factor = 1e-6
            else:
                self.logger.info("Assuming magnetic data units are 'nT'.")
                self.__scale_factor = 1e-9

        ### error value
        if self.__error_value is None:
            self.__error_value = 2.
        elif not isinstance(self.__error_value,float) and not isinstance(self.__error_value,int):
            raise TypeError("Invalid type of error parameter")


    def getDataExtents(self):
        """
        returns ( (x0, y0), (nx, ny), (dx, dy) )
        """
        return (list(self.__origin), list(self.__nPts), list(self.__delta))

    def getDataType(self):
        return self.__data_type

    def getSurveyData(self, domain, origin, NE, spacing):
        FS = ReducedFunction(domain)
        nValues=self.__nPts
        # determine base location of this dataset within the domain
        first=[int((self.__origin[i]-origin[i])/spacing[i]) for i in range(len(self.__nPts))]
        # determine the resolution difference between domain and data.
        # If domain has twice the resolution we can double up the data etc.
        multiplier=[int(round(self.__delta[i]/spacing[i])) for i in range(len(self.__nPts))]
        if domain.getDim()==3:
            first.append(int((self.getHeightScale()*self.__altitude-origin[2])/spacing[2]))
            multiplier=multiplier+[1]
            nValues=nValues+[1]

        reverse = [0]*domain.getDim()
        byteorder=BYTEORDER_NATIVE
        self.logger.debug("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s, reverse=%s"%(str(first),str(nValues),str(multiplier),str(reverse)))
        data = readBinaryGrid(self.__datafile, FS, shape=(),
                fill=self.__null_value, byteOrder=byteorder,
                dataType=self.__celltype, first=first, numValues=nValues,
                multiplier=multiplier, reverse=reverse)
        sigma = self.__error_value * esu.whereNonZero(data-self.__null_value)

        data = data * self.__scale_factor
        sigma = sigma * self.__scale_factor
        return data, sigma

    def getUtmZone(self):
        return self.__utm_zone


##############################################################################
class NetCdfData(DataSource):
    """
    Data Source for gridded netCDF data that use CF/COARDS conventions.
    """
    def __init__(self, data_type, filename, altitude=0., data_variable=None,
                       error=None, scale_factor=None, null_value=None, reference_system=None):
        """
        :param filename: file name for survey data in netCDF format
        :type filename: ``str``
        :param data_type: type of data, must be `GRAVITY` or `MAGNETIC`
        :type data_type: ``int``
        :param altitude: altitude of measurements in meters
        :type altitude: ``float``
        :param data_variable: name of the netCDF variable that holds the data.
                              If not provided an attempt is made to determine
                              the variable and an exception thrown on failure.
        :type data_variable: ``str``
        :param error: either the name of the netCDF variable that holds the
                      uncertainties of the measurements or a constant value
                      to use for the uncertainties. If a constant value is
                      supplied, it is scaled by the same factor as the
                      measurements. If not provided the error is assumed to
                      be 2 units for all measurements (i.e. 0.2 mGal and 2 nT
                      for gravity and magnetic, respectively)
        :type error: ``str`` or ``float``
        :param scale_factor: the measurements and error values are scaled by
                             this factor. By default, gravity data is assumed
                             to be given in 1e-6 m/s^2 (0.1 mGal), while
                             magnetic data is assumed to be in 1e-9 T (1 nT).
        :type scale_factor: ``float``
        :param null_value: value that is used in the file to mark undefined
                           areas. This information is usually included in the
                           file.
        :type null_value: ``float``
        :param reference_system: reference coordinate system to be used.
                                 For a Cartesian reference (default) the
                                 appropriate UTM transformation is applied.
        :type reference_system: `ReferenceSystem`
        :note: it is the responsibility of the caller to ensure all data
               sources and the domain builder use the same reference system.
        """
        super(NetCdfData,self).__init__(reference_system, [filename])
        self.__filename=filename
        if not data_type in [self.GRAVITY,self.MAGNETIC]:
            raise ValueError("Invalid value for data_type parameter")
        self.__data_type = data_type
        self.__altitude = altitude
        self.__data_name = data_variable
        self.__scale_factor = scale_factor
        self.__null_value = null_value
        self.__utm_zone = None
        self.__readMetadata(error)

    def __readMetadata(self, error):
        self.logger.debug("Checking Data Source: %s"%self.__filename)
        
        #first check what sort of file we are dealing with
        nct=es.NcFType(self.__filename)
        if nct=='u':
            raise RuntimeError('Unsupported filetype for '+str(self.__filename))            
        elif nct=='C' or nct=='c':
            try:
                import netCDF4 as nc4
                nct=4
                f=nc4.Dataset(self.__filename, 'r')
            except:
                if HAVE_SCIPY_NETCDF:
                    nct=3
                    f=netcdf_file(self.__filename, 'r')
                else:
                    raise RuntimeError('Unable to load file '+str(self.__filename))
        elif nct=='4':
            nct=4
            import netCDF4 as nc4
            f=nc4.Dataset(self.__filename, 'r')
        else:
            raise RuntimeError('Unable to determine type of file '+str(self.__filename)+" "+nct)
        ### longitude- / X-dimension and variable
        NX=0
        for n in ['lon','longitude','x']:
            if n in f.dimensions:
                NX=f.dimensions[n]
                lon_name=n
                break
        if NX==0:
            raise RuntimeError("Could not determine extents of data")

        # CF/COARDS states that coordinate variables have the same name as
        # the dimensions
        if not lon_name in f.variables:
            raise RuntimeError("Could not determine longitude variable")
        longitude=f.variables.pop(lon_name)

        ### latitude- / Y-dimension and variable
        NY=0
        for n in ['lat','latitude','y']:
            if n in f.dimensions:
                NY=f.dimensions[n]
                lat_name=n
                break
        if NY==0:
            raise RuntimeError("Could not determine extents of data")
        if not lat_name in f.variables:
            raise RuntimeError("Could not determine latitude variable")
        latitude=f.variables.pop(lat_name)

        # correct for netCDF4 having a different type
        if nct==4:
            NX=NX.size
            NY=NY.size
        ### data variable
        if self.__data_name is not None:
            try:
                dims = f.variables[self.__data_name].dimensions
                if not ((lat_name in dims) and (lon_name in dims)):
                    raise ValueError("Invalid data variable name supplied")
                
            except KeyError:
                raise ValueError("Invalid data variable name supplied")
        else:
            for n in sorted(f.variables.keys()):
                dims=f.variables[n].dimensions
                if (lat_name in dims) and (lon_name in dims):
                    self.__data_name=n
                    break
        if self.__data_name is None:
            raise RuntimeError("Could not determine data variable")
        self.__data_name=str(self.__data_name)          # try to deal with unicode->str issue in py2
        
        datavar = f.variables[self.__data_name]

        ### error value/variable
        self.__error_name = None
        if isinstance(error,str):
            try:
                dims = f.variables[error].dimensions
                if not ((lat_name in dims) and (lon_name in dims)):
                    raise ValueError("Invalid error variable name supplied")
            except KeyError:
                raise ValueError("Invalid error variable name supplied")
            self.__error_name = error
        elif isinstance(error,float) or isinstance(error,int):
            self.__error_value = float(error)
        elif error is None:
            self.__error_value = 2.
        else:
            raise TypeError("Invalid type of error parameter")

        ### mask/null value
        # note that NaN is always filtered out in ripley
        if self.__null_value is None:
            if hasattr(datavar, 'missing_value'):
                self.__null_value = float(datavar.missing_value)
            elif hasattr(datavar, '_FillValue'):
                self.__null_value = float(datavar._FillValue)
            else:
                self.logger.debug("Could not determine null value, using default.")
                self.__null_value = 99999
        elif not isinstance(self.__null_value,float) and not isinstance(self.__null_value,int):
            raise TypeError("Invalid type of null_value parameter")

        # try to determine units of data - this is disabled until we obtain a
        # file with valid information
        #if hasattr(f.variables[data_name], 'units'):
        #   units=f.variables[data_name].units

        ### scale factor
        if self.__scale_factor is None:
            if self.__data_type == self.GRAVITY:
                self.logger.info("Assuming gravity data scale is 1e-6 m/s^2.")
                self.__scale_factor = 1e-6
            else:
                self.logger.info("Assuming magnetic data units are 'nT'.")
                self.__scale_factor = 1e-9

        # see if there is a WKT string to convert coordinates
        try:
            wkt_string=str(datavar.esri_pe_string)
            self.logger.debug("wkt_string is: %s"%wkt_string)
        except:
            wkt_string=None

        # CF GDAL output: see if there is a grid_mapping attribute which
        # contains the name of a dummy variable that holds information about
        # mapping. GDAL puts the WKT string into spatial_ref:
        if wkt_string is None:
            try:
                mapvar=f.variables[datavar.grid_mapping]
                wkt_string=str(mapvar.spatial_ref)
                self.logger.debug("wkt_string is: %s"%wkt_string)
            except:
                self.logger.debug("no wkt_string found!")
        # actual_range & geospatial_lon_min/max do not always contain correct
        # values so we have to obtain the min/max in a less efficient way:
        if nct==4:
            lon_range=min(longitude),max(longitude)
            lat_range=min(latitude),max(latitude)
        else:
            lon_range=longitude.data.min(),longitude.data.max()
            lat_range=latitude.data.min(),latitude.data.max()
        if self.getReferenceSystem().isCartesian():
            lon_range,lat_range,zone=LatLonToUTM(lon_range, lat_range, wkt_string)
            self.__utm_zone = zone
             
        lengths=[lon_range[1]-lon_range[0], lat_range[1]-lat_range[0]]

        # see if lat or lon is stored in reverse order to domain conventions
        self.__reverse=[False,False]
        d=longitude
        if nct==3:
            d=d.data
        if d[0]>d[-1]:
            self.__reverse[0]=True
        d=latitude
        if nct==3:
            d=d.data
        if d[0]>d[-1]:
            self.__reverse[1]=True
        self.__nPts=[NX, NY]
        self.__origin=[lon_range[0],lat_range[0]]
        # we are rounding to avoid interpolation issues
        if self.getReferenceSystem().isCartesian():
            # rounding will give us about meter-accuracy with UTM coordinates
            r=0
        else:
            # this should give us about meter-accuracy with lat/lon coords
            r=5
        self.__delta=[np.round(lengths[i]/self.__nPts[i],r) for i in range(2)]
        del longitude, latitude, d, datavar
        f.close()

    def getDataExtents(self):
        """
        returns ( (x0, y0), (nx, ny), (dx, dy) )
        """
        return (list(self.__origin), list(self.__nPts), list(self.__delta))

    def getDataType(self):
        return self.__data_type
        
    def getScaleFactor(self):
        """
        returns the data scale factor for adjusting measurement to SI
        """
        return self.__scale_factor
        
    def getSurveyData(self, domain, origin, NE, spacing):
        if not HAS_RIPLEY:
            raise RuntimeError("Ripley module not available for reading")
        FS=ReducedFunction(domain)
        nValues=self.__nPts
        # determine base location of this dataset within the domain
        first=[int((self.__origin[i]-origin[i])/spacing[i]) for i in range(len(self.__nPts))]
        # determine the resolution difference between domain and data.
        # If domain has twice the resolution we can double up the data etc.
        multiplier=[int(round(self.__delta[i]/spacing[i])) for i in range(len(self.__nPts))]
        reverse = [int(self.__reverse[i]) for i in range(len(self.__reverse))]

        if domain.getDim() == 3:
            first.append(int((self.getHeightScale()*self.__altitude-origin[2])/spacing[2]))
            multiplier = multiplier + [1]
            nValues = nValues + [1]
            reverse = reverse + [0]

        self.logger.debug("calling readNcGrid with dataname=%s, first=%s, nValues=%s, multiplier=%s, reverse=%s"%(
            self.__data_name, str(first),str(nValues),str(multiplier),str(reverse)))
        data = ripleycpp._readNcGrid(self.__filename, self.__data_name, FS,
                shape=(), fill=self.__null_value, first=first,
                numValues=nValues, multiplier=multiplier, reverse=reverse)

        if self.__error_name is not None:
            self.logger.debug("calling readNcGrid with dataname=%s, first=%s, nValues=%s, multiplier=%s, reverse=%s"%(
                self.__data_name, str(first),str(nValues),str(multiplier),str(reverse)))
            sigma = ripleycpp._readNcGrid(self.__filename, self.__error_name,
                    FS, shape=(), fill=0., first=first, numValues=nValues,
                    multiplier=multiplier, reverse=reverse)
        else:
            # arithmetics with NaN produces undesired results so we replace
            # NaNs by a large positive number which (hopefully) is not present
            # in the real dataset
            if np.isnan(self.__null_value):
                data.replaceNaN(1.e300)
                self.__null_value = 1.e300
            sigma = self.__error_value * esu.whereNonZero(data-self.__null_value)

        data = data * self.__scale_factor
        sigma = sigma * self.__scale_factor
        return data, sigma

    def getUtmZone(self):
        return self.__utm_zone


##############################################################################
class SourceFeature(object):
    """
    A feature adds a density/susceptibility distribution to (parts of) a
    domain of a synthetic data source, for example a layer of a specific
    rock type or a simulated ore body.
    """
    def getValue(self):
        """
        Returns the value for the area covered by mask. It can be constant
        or a `Data` object with spatial dependency.
        """
        raise NotImplementedError

    def getMask(self, x):
        """
        Returns the mask of the area of interest for this feature. That is,
        mask is non-zero where the value returned by `getValue()` should be
        applied, zero elsewhere.
        """
        raise NotImplementedError

class SmoothAnomaly(SourceFeature):
    """
    A source feature in the form of a blob (roughly gaussian).
    """
    def __init__(self, lx, ly, lz, x, y, depth, v_inner=None, v_outer=None):
        """
        Intializes the smooth anomaly data.

        :param lx: size of blob in x-dimension
        :param ly: size of blob in y-dimension
        :param lz: size of blob in z-dimension
        :param x: location of blob in x-dimension
        :param y: location of blob in y-dimension
        :param depth: depth of blob
        :param v_inner: value in the centre of the blob
        :param v_outer: value in the periphery of the blob
        """
        self.x=x
        self.y=y
        self.lx=lx
        self.ly=ly
        self.lz=lz
        self.depth=depth
        self.v_inner=v_inner
        self.v_outer=v_outer
        self.value=None
        self.mask=None

    def getValue(self,x):
        if self.value is None:
            if self.v_outer is None or self.v_inner is None:
                self.value=0
            else:
                DIM=x.getDomain().getDim()
                alpha=-log(abs(self.v_outer/self.v_inner))*4
                value=exp(-alpha*((x[0]-self.x)/self.lx)**2)
                value=value*exp(-alpha*((x[DIM-1]+self.depth)/self.lz)**2)
                self.value=maximum(abs(self.v_outer), abs(self.v_inner*value))
                if self.v_inner<0: self.value=-self.value

        return self.value

    def getMask(self, x):
        DIM=x.getDomain().getDim()
        m=esu.whereNonNegative(x[DIM-1]+self.depth+self.lz/2) * esu.whereNonPositive(x[DIM-1]+self.depth-self.lz/2) \
            *esu.whereNonNegative(x[0]-(self.x-self.lx/2)) * esu.whereNonPositive(x[0]-(self.x+self.lx/2))
        if DIM>2:
            m*=esu.whereNonNegative(x[1]-(self.y-self.ly/2)) * esu.whereNonPositive(x[1]-(self.y+self.ly/2))
        self.mask = m
        return m

##############################################################################
class SyntheticDataBase(DataSource):
    """
    Base class to define reference data based on a given property distribution
    (density or susceptibility). Data are collected from a square region of
    vertical extent ``length`` on a grid with ``number_of_elements`` cells in
    each direction.

    The synthetic data are constructed by solving the appropriate forward
    problem. Data can be sampled with an offset from the surface at z=0 or
    using the entire subsurface region.
    """
    def __init__(self, data_type,
                 DIM=2,
                 number_of_elements=10,
                 length=1*U.km,
                 B_b=None,
                 data_offset=0,
                 full_knowledge=False):
        """
        :param data_type: data type indicator
        :type data_type: `DataSource.GRAVITY`, `DataSource.MAGNETIC`
        :param DIM: number of spatial dimensions
        :type DIM: ``int`` (2 or 3)
        :param number_of_elements: lateral number of elements in the region
                                   where data are collected
        :type number_of_elements: ``int``
        :param length: lateral extent of the region where data are collected
        :type length: ``float``
        :param B_b: background magnetic flux density [B_r, B_latiude, B_longitude]. Only used for magnetic data.
        :type B_b: ``list`` of ``Scalar``
        :param data_offset: offset of the data collection region from the surface
        :type data_offset: ``float``
        :param full_knowledge: if ``True`` data are collected from the entire
                               subsurface region. This is mainly for testing.
        :type full_knowledge: ``Bool``
        """
        super(SyntheticDataBase, self).__init__()
        if not data_type in [self.GRAVITY, self.MAGNETIC]:
            raise ValueError("Invalid value for data_type parameter")
        self.DIM=DIM
        self.number_of_elements=number_of_elements
        self.length=length
        self.__data_type = data_type
        self.__full_knowledge= full_knowledge
        self.__data_offset=data_offset
        self.__B_b =None
        # this is for Cartesian (FIXME ?)
        if data_type  ==  self.MAGNETIC:
            if self.DIM < 3:
                self.__B_b = np.array([B_b[0], B_b[2]])
            else:
                self.__B_b = ([B_b[0], B_b[1],B_b[2]])
        self.__origin = [0]*(DIM-1)
        self.__delta = [float(length)/number_of_elements]*(DIM-1)
        self.__nPts = [number_of_elements]*(DIM-1)
        self._reference_data=None

    def getUtmZone(self):
        """
        returns a dummy UTM zone since this class does not use real coordinate
        values.
        """
        return 0

    def getDataExtents(self):
        """
        returns the lateral data extend of the data set
        """
        return (list(self.__origin), list(self.__nPts), list(self.__delta))

    def getDataType(self):
        """
        returns the data type
        """
        return self.__data_type

    def getSurveyData(self, domain, origin, number_of_elements, spacing):
        """
        returns the survey data placed on a given domain.

        :param domain: domain on which the data are to be placed
        :type domain: ``Domain``
        :param origin: origin of the domain
        :type origin: ``list`` of ``float``
        :param number_of_elements: number of elements (or cells) in each
                                   spatial direction used to span the domain
        :type number_of_elements: ``list`` of ``int``
        :param spacing: cell size in each spatial direction
        :type spacing: ``list`` of ``float``
        :return: observed gravity field or magnetic flux density for each cell
                 in the domain and for each cell an indicator 1/0 if the data
                 are valid or not.
        :rtype: pair of ``Scalar``
        """
        pde=LinearSinglePDE(domain)
        DIM=domain.getDim()
        x=domain.getX()
        # set the reference data
        k=self.getReferenceProperty(domain)
        # calculate the corresponding potential
        z=x[DIM-1]
        m_psi_ref=esu.whereZero(z-esu.sup(z))
        if self.getDataType()==DataSource.GRAVITY:
            pde.setValue(A=esu.kronecker(domain), Y=-4*np.pi*U.Gravitational_Constant*self._reference_data, q=m_psi_ref)
        else:
            pde.setValue(A=esu.kronecker(domain), X=self._reference_data*self.__B_b, q=m_psi_ref)
        pde.setSymmetryOn()
        #pde.getSolverOptions().setTolerance(1e-13)
        psi_ref=pde.getSolution()
        del pde
        if self.getDataType()==DataSource.GRAVITY:
            data = -esu.grad(psi_ref, ReducedFunction(domain))
        else:
            data = self._reference_data*self.__B_b-esu.grad(psi_ref, ReducedFunction(domain))

        x=ReducedFunction(domain).getX()
        if self.__full_knowledge:
            sigma = esu.whereNegative(x[DIM-1])
        else:
            sigma=1.
            # limit mask to non-padding in horizontal area
            for i in range(DIM-1):
                x_i=x[i]
                sigma=sigma * esu.wherePositive(x_i) * esu.whereNegative(x_i-(esu.sup(x_i)+esu.inf(x_i)))
            # limit mask to one cell thickness at z=0
            z=x[DIM-1]
            oo=int(self.__data_offset/spacing[DIM-1]+0.5)*spacing[DIM-1]
            sigma = sigma * esu.whereNonNegative(z-oo) * esu.whereNonPositive(z-oo-spacing[DIM-1])
        return data,sigma

    def getReferenceProperty(self, domain=None):
        """
        Returns the reference `Data` object that was used to generate
        the gravity/susceptibility anomaly data.

        :return: the density or susceptibility anomaly used to create the
                 survey data
        :note: it can be assumed that in the first call the ``domain``
               argument is present so the actual anomaly data can be created.
               In subsequent calls this may not be true.
        :note: method needs to be overwritten
        """
        raise NotImplementedError()

class SyntheticFeatureData(SyntheticDataBase):
    """
    Uses a list of `SourceFeature` objects to define synthetic anomaly data.
    """
    def __init__(self, data_type,
                       features,
                       DIM=2,
                       number_of_elements=10,
                       length=1*U.km,
                       B_b=None,
                       data_offset=0,
                       full_knowledge=False):
        """
        :param data_type: data type indicator
        :type data_type: `DataSource.GRAVITY`, `DataSource.MAGNETIC`
        :param features: list of features. It is recommended that the features
                         are located entirely below the surface.
        :type features: ``list`` of `SourceFeature`
        :param DIM: spatial dimensionality
        :type DIM: ``int`` (2 or 3)
        :param number_of_elements: lateral number of elements in the region
                                   where data are collected
        :type number_of_elements: ``int``
        :param length: lateral extent of the region where data are collected
        :type length: ``float``
        :param B_b: background magnetic flux density [B_r, B_latiude, B_longitude]. Only used for magnetic data.
        :type B_b: ``list`` of ``Scalar``
        :param data_offset: offset of the data collection region from the surface
        :type data_offset: ``float``
        :param full_knowledge: if ``True`` data are collected from the entire subsurface region. This is mainly for testing.
        :type full_knowledge: ``Bool``
        """
        super(SyntheticFeatureData,self).__init__(
                                 data_type=data_type, DIM=DIM,
                                 number_of_elements=number_of_elements,
                                 length=length, B_b=B_b,
                                 data_offset=data_offset,
                                 full_knowledge=full_knowledge)
        self._features = features


    def getReferenceProperty(self, domain=None):
        """
        Returns the reference `Data` object that was used to generate
        the gravity/susceptibility anomaly data.
        """
        if self._reference_data == None:
            DIM=domain.getDim()
            x=domain.getX()
            k=0.
            for f in self._features:
                m=f.getMask(x)
                k = k * (1-m) + f.getValue(x) * m
            self._reference_data= k
        return self._reference_data

class SyntheticData(SyntheticDataBase):
    """
    Defines synthetic gravity/magnetic data based on harmonic property anomaly

        rho = amplitude * sin(n_depth * pi /depth * (z+depth_offset)) * sin(n_length * pi /length * (x - shift) )

    for all x and z<=0. For z>0 rho = 0.
    """
    def __init__(self, data_type,
                       n_length=1,
                       n_depth=1,
                       depth_offset=0.,
                       depth=None,
                       amplitude=None,
                       DIM=2,
                       number_of_elements=10,
                       length=1*U.km,
                       B_b=None,
                       data_offset=0,
                       full_knowledge=False,
                       s=0.):
        """
        :param data_type: data type indicator
        :type data_type: `DataSource.GRAVITY`, `DataSource.MAGNETIC`
        :param n_length: number of oscillations in the anomaly data within the
                         observation region
        :type n_length: ``int``
        :param n_depth: number of oscillations in the anomaly data below surface
        :type n_depth: ``int``
        :param depth_offset: vertical offset of the data
        :type depth_offset: ``float``
        :param depth: vertical extent in the anomaly data. If not present the
                      depth of the domain is used.
        :type depth: ``float``
        :param amplitude: data amplitude. Default value is 200 U.kg/U.m**3 for
                          gravity and 0.1 for magnetic data.
        :param DIM: spatial dimensionality
        :type DIM: ``int`` (2 or 3)
        :param number_of_elements: lateral number of elements in the region
                                   where data are collected
        :type number_of_elements: ``int``
        :param length: lateral extent of the region where data are collected
        :type length: ``float``
        :param B_b: background magnetic flux density [B_r, B_latiude, B_longitude].
                    Only used for magnetic data.
        :type B_b: ``list`` of ``Scalar``
        :param data_offset: offset of the data collection region from the surface
        :type data_offset: ``float``
        :param full_knowledge: if ``True`` data are collected from the entire
                               subsurface region. This is mainly for testing.
        :type full_knowledge: ``Bool``
        """
        super(SyntheticData,self).__init__(
                                 data_type=data_type, DIM=DIM,
                                 number_of_elements=number_of_elements,
                                 length=length, B_b=B_b,
                                 data_offset=data_offset,
                                 full_knowledge=full_knowledge)
        self.__n_length = n_length
        self.__n_depth = n_depth
        self.depth = depth
        self.depth_offset=depth_offset
        if amplitude == None:
            if data_type == DataSource.GRAVITY:
                amplitude = 200 *U.kg/U.m**3
            else:
                amplitude = 0.1
        self.__amplitude = amplitude
        self.__s=s

    def getReferenceProperty(self, domain=None):
        """
        Returns the reference `Data` object that was used to generate
        the gravity/susceptibility anomaly data.
        """
        if self._reference_data is None:
            DIM=domain.getDim()
            x=domain.getX()
            # set the reference data
            z=x[DIM-1]
            dd=self.depth
            if dd is None: dd=esu.inf(z)
            z2=(z+self.depth_offset)/(self.depth_offset-dd)
            k=esu.sin(self.__n_depth * np.pi  * z2) * esu.whereNonNegative(z2) * esu.whereNonPositive(z2-1.) * self.__amplitude
            for i in range(DIM-1):
               x_i=x[i]
               min_x=esu.inf(x_i)
               max_x=esu.sup(x_i)
               k*= esu.sin(self.__n_length*np.pi*(x_i-min_x-self.__s)/(max_x-min_x))
            self._reference_data= k
        return self._reference_data

##############################################################################
class NumpyData(DataSource):
    """
    """

    def __init__(self, data_type, data, error=1., length=1.*U.km, null_value=-1., tags=[], origin=None):
        """
        A data source that uses survey data from a ``numpy`` object or list
        instead of a file.
        The dimensionality is inferred from the shape of ``data`` (1- and
        2-dimensional data is supported). The data origin is assumed to be
        at the coordinate origin.

        :param data_type: data type indicator
        :type data_type: `DataSource.GRAVITY`, `DataSource.MAGNETIC`
        :param data: the survey data array. Note that for a cartesian coordinate
                     system the shape of the data is considered to be
                     (nz,ny,nx).
        :type data: ``numpy.array`` or ``list``
        :param error: constant value to use for the data uncertainties or a
                      numpy object with uncertainties for every data point.
        :type error: ``float`` or ``list`` or ``ndarray``
        :param length: side length(s) of the data slice/volume. This can be
                       a scalar to indicate constant length in all dimensions
                       or an array/list of values in each coordinate dimension.
        :type length: ``float`` or ``list`` or ``ndarray``
        :param null_value: value that is used in the undefined regions of the
                           survey data object.
        :type null_value: ``float``
        :param tags: a list of tags associated with the data set.
        :type tags: ``list`` of almost any type (typically `str`) 
        :param origin: offset of origin of offset
        :type origin: ``list`` of ``float``
        """
        super(NumpyData, self).__init__(tags=tags)
        if not data_type in [self.GRAVITY, self.MAGNETIC, self.ACOUSTIC, self.MT ]:
            raise ValueError("Invalid value for data_type parameter")
        self.__data_type = data_type
        if not isinstance(data, np.ndarray) or data.dtype not in [ np.float64, np.complex128]:
            self.__data = np.asarray(data, dtype=np.float64)
        else:
            self.__data = data
        DIM = len(self.__data.shape)
        if DIM not in (1,2):
            raise ValueError("NumpyData requires 1- or 2-dimensional data")
        self.__error = np.asarray(error, dtype=np.float64)
        if len(self.__error.shape) > 0 and \
                self.__error.shape != self.__data.shape:
            raise ValueError("error argument must be scalar or match the shape of the data")
        if isinstance(length,float) or isinstance(length,int):
            length = [float(length)] * DIM
        else:
            length=np.asarray(length, dtype=float)
            if len(length.shape) != 1 or length.shape[0] != DIM:
                raise ValueError("length must be scalar or an array with %d values"%DIM)
            length=length.tolist()

        self.__length=length
        self.__null_value = null_value
        self.__nPts = list(reversed(self.__data.shape))
        self.__delta = [length[i]/self.__nPts[i] for i in range(DIM)]
        if origin is None:
            self.__origin = [0.] * DIM 
        else:
            self.__origin = origin[:DIM-1]
            
    def getDataExtents(self):
        return (self.__origin, self.__nPts, self.__delta)

    def getDataType(self):
        return self.__data_type

    def getSurveyData(self, domain, origin, NE, spacing):
        import os
        DIM=domain.getDim()
        if self.getDataType()  == self.ACOUSTIC:
            x=esu.FunctionOnBoundary(domain).getX()
            BBX=boundingBox(domain)
            z=x[DIM-1] 
            mask= esu.whereZero( z - esu.inf(z))  # we don't use BBX[DIM-1][1] due to mountains'
            for i in range(DIM-1):
                x_i=x[i]
                mask+=esu.whereNonPositive( x_i - self.__origin[i] ) + esu.whereNonNegative( x_i - ( self.__origin[i]+self.__length[i] ) ) 
            mask=1-esu.wherePositive(mask)
            
            data=esu.Data(0.,(2,), esu.FunctionOnBoundary(domain))
            step= [ self.__length[i]/self.__data.shape[i] for i in range(DIM-1) ]
            if DIM == 2:
                data[0] = esu.interpolateTable(self.__data.real, x[0],self.__origin[0], step[0])
                data[1] = esu.interpolateTable(self.__data.imag, x[0],self.__origin[0], step[0])
                if len(self.__error.shape) > 0:
                    sigma = esu.interpolateTable(self.__error, x[0],  self.__origin[0], step[0])
                else:
                    sigma = esu.Scalar(self.__error.item(), esu.FunctionOnBoundary(domain))
            else:
                raise ValueError("3D domains are not supported yet.")
            data*=mask
            sigma*=mask
        elif self.getDataType()  == self.MT:
            if DIM == 2:
                step= [ self.__length[i]/self.__data.shape[i] for i in range(DIM-1) ]
                if len(self.__error.shape) > 0:
                    sigma = esu.interpolateTable(self.__error, x[0],  self.__origin[0], step[0])
                else:
                    sigma = esu.Scalar(self.__error.item(), esu.FunctionOnBoundary(domain))
                return self.__data, sigma
            else:
                raise ValueError("3D domains are not supported yet.")

        else:
            FS = ReducedFunction(domain)
            nValues = self.__nPts
            dataDIM = len(nValues)
            # determine base location of this dataset within the domain
            first=[int((self.__origin[i]-origin[i])/spacing[i]) for i in range(dataDIM)]
            # determine the resolution difference between domain and data.
            # If domain has twice the resolution we can double up the data etc.
            multiplier=[int(round(self.__delta[i]/spacing[i])) for i in range(dataDIM)]
            if domain.getDim() > dataDIM:
                first.append(int(-origin[-1]/spacing[-1]))
                multiplier=multiplier+[1]
                nValues=nValues+[1]

            _handle, numpyfile = tempfile.mkstemp()
            os.close(_handle)
            self.__data.tofile(numpyfile)

            reverse=[0]*DIM
            byteorder=BYTEORDER_NATIVE
            datatype=DATATYPE_FLOAT64
            self.logger.debug("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s"%(str(first),str(nValues),str(multiplier)))
            data = readBinaryGrid(numpyfile, FS, shape=(),
                fill=self.__null_value, byteOrder=byteorder, dataType=datatype,
                first=first, numValues=nValues, multiplier=multiplier,
                reverse=reverse)
            if len(self.__error.shape) > 0:
                self.__error.tofile(numpyfile)
                self.logger.debug("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s"%(str(first),str(nValues),str(multiplier)))
                sigma = readBinaryGrid(numpyfile, FS, shape=(),
                    fill=0., byteOrder=byteorder, dataType=datatype,
                    first=first, numValues=nValues, multiplier=multiplier,
                    reverse=reverse)
            else:
                sigma = self.__error.item() * esu.whereNonZero(data-self.__null_value)
            os.unlink(numpyfile)

        return data, sigma

    def getUtmZone(self):
        """
        returns a dummy UTM zone since this class does not use real coordinate
        values.
        """
        return 0

class MT2DTe(object):
    """
    class used to store frequency information accosicated with mt data
    """

    def __init__(self,x, omega=0):
        """
        initiale the MT2DTe tag object

        :param omega: frequency of readings
        :type omega: ``float``
        :param x: coordinates of measurements 
        :type x: ``list`` of ``tuple`` with ``float``
        """
        self.__omega=omega
        self.__x=x

    def getFrequency(self):
        """
        return frequency of measurement
        :rtype: ``float``
        """
        return self.__omega
    
    def getX(self):
        """
        return coordinates of measurement
        :rtype: ``float``
        """
        return self.__x

class SeismicSource(object):
    """
    describes a seimic source by location (x,y), frequency omega, power (if known) and orientation (if known).
    this class is used to tag seismic data sources.
    """ 
    def __init__(self, x, y=0., omega=0., elevation=0.,  power = None, orientation=None):
        """
        initiale the source
        
        :param x: lateral x location
        :param y: lateral y location
        :param omega: frequency of source
        :param elevation: elevation of source above reference level 
        :param power: power of source at frequence
        :param orientation: oriententation of source in 3D or 2D (or None)
        :type x: ``float``
        :type y: ``float``
        :type omega: ``float``
        :type power: ``complex`` or ``None``
        :type orientation: vector of appropriate length or ``None``
        """
        self.__loc=(x,y)
        self.__omega=omega
        self.__power=power
        self.__elevation=elevation
        self.__orientation=orientation
        
    def __eq__(self, other):
        if isinstance(other, SeismicSource):
            return self.__loc == other.getLocation() \
                and self.__omega == other.getFrequency() \
                and self.__elevation == other.getElevation() \
                and self.__power == other.getPower() \
                and self.__orientation == other.getOrientation()
        else:
            return False
    
    def __ne__(self, other):
            return not self.__eq__(other)
            
    def getLocation(self):
        """
        return location of source
        :rtype: ``tuple`` of ``float``
        """
        return self.__loc
    def getFrequency(self):
        """
        return frequency of source
        :rtype: ``float``
        """
        return self.__omega
        
    def getElevation(self):
        """
        return elevation of source
        :rtype: ``float``
        """
        return self.__elevation
        
    def getPower(self):
        """
        return power of source at frequency
        :rtype: ``complex`` or ``None``
        """
        return self.__power
    def getOrientation(self):
        """
        return power of source orientation at frequency
        :rtype: vector type object or ``None``
        """
        return self.__orientation
    
    
