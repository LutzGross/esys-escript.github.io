
##############################################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

"""Data readers/providers for inversions"""

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['simpleGeoMagneticFluxDensity', 'DataSource','ErMapperData', 'SyntheticDataBase' , 'SyntheticFeatureData', 'SyntheticData','SmoothAnomaly']

import logging
import numpy as np
from esys.escript import ReducedFunction
from esys.escript.linearPDEs import LinearSinglePDE
from esys.escript.util import *
import esys.escript.unitsSI as U
from esys.ripley import Brick, Rectangle, ripleycpp

try:
    from scipy.io.netcdf import netcdf_file
    __all__ += ['NetCdfData']
except:
    pass

def LatLonToUTM(lon, lat, wkt_string=None):
    """
    Converts one or more longitude,latitude pairs to the corresponding x,y
    coordinates in the Universal Transverse Mercator projection.

    :note: If the ``pyproj`` module is not installed a warning is printed and
           the input values are scaled by a constant and returned.
    :note: If `wkt_string` is not given or invalid or the ``gdal`` module is
           not available to convert the string, then the input values are
           assumed to be given in the Clarke 1866 projection.

    :param lon: longitude value(s)
    :type lon: `float`, `list`, `tuple`, or ``numpy.array``
    :param lat: latitude value(s)
    :type lat: `float`, `list`, `tuple`, or ``numpy.array``
    :rtype: ``numpy.array``
    """

    # not really optimal: if pyproj is not installed we return the input
    # values scaled by a constant.
    try:
        import pyproj
    except:
        print("Warning, pyproj not available. Domain extents will be wrong")
        return np.array(lon)*1000., np.array(lat)*1000.

    # determine UTM zone from the input data
    zone=int(np.median((np.floor((np.array(lon) + 180)/6) + 1) % 60))
    try:
        import osgeo.osr
        srs = osgeo.osr.SpatialReference()
        srs.ImportFromWkt(wkt_string)
        p_src = pyproj.Proj(srs.ExportToProj4())
    except:
        p_src = pyproj.Proj('+proj=longlat +ellps=clrk66 +no_defs')
    # we assume southern hemisphere here
    p_dest = pyproj.Proj('+proj=utm +zone=%d +south +units=m'%zone)
    x,y=pyproj.transform(p_src, p_dest, lon, lat)
    return x,y

def simpleGeoMagneticFluxDensity(latitude, longitude=0.):
        theta = (90-latitude)/180.*np.pi
        B_0=U.Mu_0  * U.Magnetic_Dipole_Moment_Earth / (4 * np.pi *  U.R_Earth**3)
        B_theta= B_0 * sin(theta)
        B_r= 2 * B_0 * cos(theta)
        return B_r, B_theta, 0.

class DataSource(object):
    """
    A class that provides survey data for the inversion process.
    This is an abstract base class that implements common functionality.
    Methods to be overwritten by subclasses are marked as such.
    This class assumes 2D data which is mapped to a slice of a 3D domain.
    For other setups override the methods as required.
    """

    GRAVITY, MAGNETIC = list(range(2))

    def __init__(self):
        """
        Constructor. Sets some defaults and initializes logger.
        """
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)
        self.__subsampling_factor=1
        self.__background_magnetic_field=None

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
        Subclasses must return `GRAVITY` or `MAGNETIC` as appropriate.
        """
        raise NotImplementedError

    def getSurveyData(self, domain, origin, NE, spacing):
        """
        This method is called by the `DomainBuilder` to retrieve the survey
        data as `Data` objects on the given domain.
        Subclasses should return one or more `Data` objects with survey data
        interpolated on the given ripley domain. The exact return type
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
    def __init__(self, datatype, headerfile, datafile=None, altitude=0.):
        """
        :param datatype: type of data, must be `GRAVITY` or `MAGNETIC`
        :type datatype: ``int``
        :param headerfile: ER Mapper header file (usually ends in .ers)
        :type headerfile: ``str``
        :param datafile: ER Mapper binary data file name. If not supplied the
                         name of the header file without '.ers' is assumed
        :type datafile: ``str``
        :param altitude: altitude of measurements above ground in meters
        :type altitude: ``float``
        """
        super(ErMapperData,self).__init__()
        self.__headerfile=headerfile
        if datafile is None:
            self.__datafile=headerfile[:-4]
        else:
            self.__datafile=datafile
        self.__altitude=altitude
        self.__datatype=datatype
        self.__readHeader()

    def __readHeader(self):
        self.logger.debug("Checking Data Source: %s (header: %s)"%(self.__datafile, self.__headerfile))
        metadata=open(self.__headerfile, 'r').readlines()
        # parse metadata
        start=-1
        for i in range(len(metadata)):
            if metadata[i].strip() == 'DatasetHeader Begin':
                start=i+1
        if start==-1:
            raise RuntimeError('Invalid ERS file (DatasetHeader not found)')

        md_dict={}
        section=[]
        for i in range(start, len(metadata)):
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
            maskval = 99999

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
            self.logger.warn("Could not determine coordinate origin. Setting to (0.0, 0.0)")
            originX,originY = 0.0, 0.0

        if 'GEODETIC' in md_dict['CoordinateSpace.Projection']:
            # it appears we have lat/lon coordinates so need to convert
            # origin and spacing. Try using gdal to get the wkt if available:
            try:
                from osgeo import gdal
                ds=gdal.Open(self.__headerfile)
                wkt=ds.GetProjection()
            except:
                wkt='GEOGCS["GEOCENTRIC DATUM of AUSTRALIA",DATUM["GDA94",SPHEROID["GRS80",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
                self.logger.warn('GDAL not available or file read error, assuming GDA94 data')
            originX_UTM,originY_UTM=LatLonToUTM(originX, originY, wkt)
            op1X,op1Y=LatLonToUTM(originX+spacingX, originY+spacingY, wkt)
            # we are rounding to avoid interpolation issues
            spacingX=np.round(op1X-originX_UTM)
            spacingY=np.round(op1Y-originY_UTM)
            originX=np.round(originX_UTM)
            originY=np.round(originY_UTM)

        # data sets have origin in top-left corner so y runs top-down
        self.__dataorigin=[originX, originY]
        originY-=(NY-1)*spacingY
        self.__delta = [spacingX, spacingY]
        self.__maskval = maskval
        self.__nPts = [NX, NY]
        self.__origin = [originX, originY]
        if self.__datatype == self.GRAVITY:
            self.logger.info("Assuming gravity data scale is 1e-6 m/s^2.")
            self.__scalefactor = 1e-6
        else:
            self.logger.info("Assuming magnetic data units are 'nT'.")
            self.__scalefactor = 1e-9

    def getDataExtents(self):
        """
        returns ( (x0, y0), (nx, ny), (dx, dy) )
        """
        return (list(self.__origin), list(self.__nPts), list(self.__delta))

    def getDataType(self):
        return self.__datatype

    def getSurveyData(self, domain, origin, NE, spacing):
        nValues=self.__nPts
        # determine base location of this dataset within the domain
        first=[int((self.__origin[i]-origin[i])/spacing[i]) for i in range(len(self.__nPts))]
        if domain.getDim()==3:
            first.append(int((self.__altitude-origin[2])/spacing[2]))
            nValues=nValues+[1]

        data=ripleycpp._readBinaryGrid(self.__datafile,
                ReducedFunction(domain),
                first, nValues, (), self.__maskval)
        sigma = whereNonZero(data-self.__maskval)
        data = data*self.__scalefactor
        sigma = sigma * 2. * self.__scalefactor
        return data, sigma


##############################################################################
class NetCdfData(DataSource):
    """
    Data Source for gridded netCDF data that use CF/COARDS conventions.
    """
    def __init__(self, datatype, filename, altitude=0.):
        """
        :param filename: file name for survey data in netCDF format
        :type filename: ``str``
        :param datatype: type of data, must be `GRAVITY` or `MAGNETIC`
        :type datatype: ``int``
        :param altitude: altitude of measurements in meters
        :type altitude: ``float``
        """
        super(NetCdfData,self).__init__()
        self.__filename=filename
        if not datatype in [self.GRAVITY,self.MAGNETIC]:
            raise ValueError("Invalid value for datatype parameter")
        self.__datatype=datatype
        self.__altitude=altitude
        self.__readMetadata()

    def __readMetadata(self):
        self.logger.debug("Checking Data Source: %s"%self.__filename)
        f=netcdf_file(self.__filename, 'r')
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

        # try to figure out data variable name
        data_name=None
        if len(f.variables)==1:
            data_name=f.variables.keys()[0]
        else:
            for n in f.variables.keys():
                dims=f.variables[n].dimensions
                if (lat_name in dims) and (lon_name in dims):
                    data_name=n
                    break
        if data_name is None:
            raise RuntimeError("Could not determine data variable")

        # try to determine value for unused data
        if hasattr(f.variables[data_name], 'missing_value'):
            maskval = float(f.variables[data_name].missing_value)
        elif hasattr(f.variables[data_name], '_FillValue'):
            maskval = float(f.variables[data_name]._FillValue)
        else:
            self.logger.debug("missing_value attribute not found, using default.")
            maskval = 99999

        # try to determine units of data - this is disabled for now
        #if hasattr(f.variables[data_name], 'units'):
        #   units=f.variables[data_name].units
        if self.__datatype == self.GRAVITY:
            self.logger.info("Assuming gravity data scale is 1e-6 m/s^2.")
            self.__scalefactor = 1e-6
        else:
            self.logger.info("Assuming magnetic data units are 'nT'.")
            self.__scalefactor = 1e-9

        # see if there is a wkt string to convert coordinates
        try:
            wkt_string=f.variables[data_name].esri_pe_string
        except:
            wkt_string=None

        # we don't trust actual_range & geospatial_lon_min/max since subset
        # data does not seem to have these fields updated.
        # Getting min/max from the arrays is obviously not very efficient but..
        #lon_range=longitude.actual_range
        #lat_range=latitude.actual_range
        #lon_range=[f.geospatial_lon_min,f.geospatial_lon_max]
        #lat_range=[f.geospatial_lat_min,f.geospatial_lat_max]
        lon_range=longitude.data.min(),longitude.data.max()
        lat_range=latitude.data.min(),latitude.data.max()
        if lon_range[1]<180:
            lon_range,lat_range=LatLonToUTM(lon_range, lat_range, wkt_string)
        lengths=[lon_range[1]-lon_range[0], lat_range[1]-lat_range[0]]
        f.close()

        self.__nPts=[NX, NY]
        self.__origin=[lon_range[0],lat_range[0]]
        # we are rounding to avoid interpolation issues
        self.__delta=[np.round(lengths[i]/self.__nPts[i]) for i in range(2)]
        #self.__wkt_string=wkt_string
        #self.__lon_name=lon_name
        #self.__lat_name=lat_name
        self.__data_name=data_name
        self.__maskval=maskval

    def getDataExtents(self):
        """
        returns ( (x0, y0), (nx, ny), (dx, dy) )
        """
        return (list(self.__origin), list(self.__nPts), list(self.__delta))

    def getDataType(self):
        return self.__datatype

    def getSurveyData(self, domain, origin, NE, spacing):
        nValues=self.__nPts
        # determine base location of this dataset within the domain
        first=[int((self.__origin[i]-origin[i])/spacing[i]) for i in range(len(self.__nPts))]
        if domain.getDim()==3:
            first.append(int((self.__altitude-origin[2])/spacing[2]))
            nValues=nValues+[1]

        data=ripleycpp._readNcGrid(self.__filename, self.__data_name,
                  ReducedFunction(domain), first, nValues, (), self.__maskval)
        sigma=whereNonZero(data-self.__maskval)
        data=data*self.__scalefactor
        sigma=sigma * 2. * self.__scalefactor
        return data, sigma


##############################################################################
class SourceFeature(object):
    """
    A feature adds a density distribution to (parts of) a domain of a synthetic
    data source, for example a layer of a specific rock type or a simulated
    ore body.
    """
    def getValue(self):
        """
        Returns the value for the area covered by mask. It can be constant
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
    def __init__(self, lx, ly, lz, x, y, depth, v_inner=None, v_outer=None):
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
        m=whereNonNegative(x[DIM-1]+self.depth+self.lz/2) * whereNonPositive(x[DIM-1]+self.depth-self.lz/2) \
            *whereNonNegative(x[0]-(self.x-self.lx/2)) * whereNonPositive(x[0]-(self.x+self.lx/2))
        if DIM>2:
            m*=whereNonNegative(x[1]-(self.y-self.ly/2)) * whereNonPositive(x[1]-(self.y+self.ly/2))
        self.mask = m
        return m

##############################################################################
class SyntheticDataBase(DataSource):
  """
  Base class to define reference data based on a given property distribution (density or
  susceptibility). Data are collected from a square region of vertical extend `length` on a 
  grid with ``number_of_elements`` cells in each direction.
  
  The synthetic data are constructed by solving the appropriate forward problem. Data can be sampled
  with an offset from the surface at z=0 or using the entire subsurface region.
  """
  def __init__(self, datatype, 
		DIM=2,
		number_of_elements=10,
		length=1*U.km,
		B_b=None,
		data_offset=0,
		full_knowledge=False,
		spherical=False):
	"""
        :param datatype: data type indicator
        :type datatype: ``DataSource.GRAVITY``, ``DataSource.MAGNETIC``
        :param DIM: spatial dimension
        :type DIM: 2 or 3
        :param number_of_elements: lateral number of elements in the region where data are collected
        :type number_of_elements: ``int``
        :param length: lateral extend of the region where data are collected
        :type length: ``float``
        :param B_b: background magnetic flux density [B_r, B_latiude, B_longitude]. Only used for magnetic data.
        :type B_b: ``list`` of ``Scalar``
        :param data_offset: offset of the data collection region from the surface
        :type data_offset: ``float``
        :param full_knowledge: if ``True`` data are collected from the entire subsurface region. This is mainly for testing.
        :type full_knowledge: ``Bool``
        :param spherical: if ``True`` sperical coordinates are used (ignored)
        :type spherical: ``Bool``
	"""
	super(SyntheticDataBase,self).__init__()
	if not datatype in [self.GRAVITY,self.MAGNETIC]:
	    raise ValueError("Invalid value for datatype parameter")	  
	self.DIM=DIM
	self.number_of_elements=number_of_elements
	self.length=length
	self.__datatype = datatype
	
	self.__spherical = spherical   
	self.__full_knowledge= full_knowledge
	self.__data_offset=data_offset
	self.__B_b =None
	# this is for Cartesian (FIXME ?)
	if datatype  ==  self.MAGNETIC:
	    if self.DIM<3:
	      self.__B_b =  np.array([-B_b[2],  -B_b[0]])
	    else:
	      self.__B_b = ([-B_b[1],  -B_b[2],  -B_b[0]])     
	self.__origin = [0]*(DIM-1)
	self.__delta = [float(length)/number_of_elements]*(DIM-1)
	self.__nPts = [number_of_elements]*(DIM-1)
	self._reference_data=None
	
  def getDataExtents(self):
	"""
	returns the lateral data extend of the data set
	"""
	return (list(self.__origin), list(self.__nPts), list(self.__delta))

  def getDataType(self):
	"""
	returns the data type
	"""
	return self.__datatype
	
  def getSurveyData(self, domain, origin, number_of_elements, spacing):
	"""
	returns the survey data placed on a given domain.
	
	:param domain: domain on which the data are to be placed
	:type param: ``Domain``
	:param origin: origin of the domain
	:type origin: ``list`` of ``float``
	:param number_of_elements: number of elements (or cells) in each spatial direction used 
				  span the domain
	:type number_of_elements: `list`` of ``int``
	:param spacing: cell size in each spatial direction
	:type spacing: ``list`` of ``float``
	:return: observed gravity field or magnetic flux density for each cell in the domain and
	for each cell an indicator 1/0 if the data are valid or not.
	:rtype: pair of ``Scalar``
	"""
	pde=LinearSinglePDE(domain)
	DIM=domain.getDim()
	x=domain.getX()
	# set the reference data
	
	k=self.getReferenceProperty(domain)
	# calculate the corresponding potential
	z=x[DIM-1]
	m_psi_ref=whereZero(z-sup(z))
	if self.getDataType()==DataSource.GRAVITY:
	    pde.setValue(A=kronecker(domain), Y=-4*np.pi*U.Gravitational_Constant*self._reference_data, q=m_psi_ref)
	else:
	    pde.setValue(A=kronecker(domain), X=self._reference_data*self.__B_b, q=m_psi_ref)
	pde.setSymmetryOn()
	psi_ref=pde.getSolution()
	del pde
	if self.getDataType()==DataSource.GRAVITY:
	    data = -grad(psi_ref, ReducedFunction(domain))
	else:
	    data = self._reference_data*self.__B_b-grad(psi_ref, ReducedFunction(domain))
	  
	x=ReducedFunction(domain).getX()    
	if self.__full_knowledge:
	    sigma = whereNegative(x[DIM-1])
	else:
	  
	    sigma=1.
	    # limit mask to non-padding in horizontal area        
	    for i in range(DIM-1):
		x_i=x[i]
		sigma=sigma * wherePositive(x_i) * whereNegative(x_i-(sup(x_i)+inf(x_i)))
	    # limit mask to one cell thickness at z=0
	    z=x[DIM-1]
	    oo=int(self.__data_offset/spacing[DIM-1]+0.5)*spacing[DIM-1]
	    sigma = sigma * whereNonNegative(z-oo) * whereNonPositive(z-oo-spacing[DIM-1])
	return data,sigma
	
  def getReferenceProperty(self, domain=None):
	"""
	Returns the reference density Data object that was used to generate
	the gravity/susceptibility anomaly data.
	
	:return: the density or susceptibility anomaly used to create the survey data.
	:note: it can be assumed that in the first call the ``domain`` argument is present so the 
	actual anomaly data can be created. In subsequent calls this may not be true.
	:note: method needs to be overwritten
	"""
	raise NotImplementedError()       
      
class SyntheticFeatureData(SyntheticDataBase):
    """
    uses a list of ``SourceFeature`` to define synthetic anomaly data
    """
    def __init__(self, datatype, 
                       features,
                       DIM=2,
                       number_of_elements=10,
                       length=1*U.km,
                       B_b=None,
                       data_offset=0,
                       full_knowledge=False,
                       spherical=False):
        """
        :param datatype: data type indicator
        :type datatype: ``DataSource.GRAVITY``, ``DataSource.MAGNETIC``
        :param features: list of features. It is recommended that the features do entirely lay below surface.
        :type features: ``list`` of ``SourceFeature``
        :param DIM: spatial dimension
        :type DIM: 2 or 3
        :param number_of_elements: lateral number of elements in the region where data are collected
        :type number_of_elements: ``int``
        :param length: lateral extend of the region where data are collected
        :type length: ``float``
        :param B_b: background magnetic flux density [B_r, B_latiude, B_longitude]. Only used for magnetic data.
        :type B_b: ``list`` of ``Scalar``
        :param data_offset: offset of the data collection region from the surface
        :type data_offset: ``float``
        :param full_knowledge: if ``True`` data are collected from the entire subsurface region. This is mainly for testing.
        :type full_knowledge: ``Bool``
        :param spherical: if ``True`` sperical coordinates are used (ignored)
        :type spherical: ``Bool``
        """
        super(SyntheticFeatureData,self).__init__(
                                 datatype=datatype, DIM=DIM, number_of_elements=number_of_elements, 
                                 length=length, B_b=B_b, 
                                 data_offset=data_offset,
                                 full_knowledge=full_knowledge,
                                 spherical=spherical)
        self._features = features


    def getReferenceProperty(self, domain=None):
        """
	Returns the reference density Data object that was used to generate
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
    defines synthetic  gravity/magnetic data based on harmonic property anomaly
    
        rho = mean + amplitude * sin(n_depth * pi /depth * z) * sin(n_length * pi /length * x - shift )
        
    for all x and z<=0. for z>0 rho = 0.        
    """
    def __init__(self, datatype, 
                       n_length=1,
                       n_depth=1, 
                       shift=0.,
                       amplitude=None, 
                       DIM=2,
                       number_of_elements=10,
                       length=1*U.km,
                       B_b=None,
                       data_offset=0,
                       full_knowledge=False,
                       spherical=False):
        """
        :param datatype: data type indicator
        :type datatype: ``DataSource.GRAVITY``, ``DataSource.MAGNETIC``
        :param n_length: number of oscillation in the anomaly data within the observation region.
        :type n_length: ``int`` 
        :param n_depth: number of oscillation in the anomaly data below surface
        :param shift: lateral phase shift in the  anomaly data
        :param amplitude: data amplitude. Default value is 200 *U.kg/U.m**3 for gravity and 0.1 for magnetic data.
        :param features: list of features. It is recommended that the features do entirely lay below surface.
        :type features: ``list`` of ``SourceFeature``
        :param DIM: spatial dimension
        :type DIM: 2 or 3
        :param number_of_elements: lateral number of elements in the region where data are collected
        :type number_of_elements: ``int``
        :param length: lateral extend of the region where data are collected
        :type length: ``float``
        :param B_b: background magnetic flux density [B_r, B_latiude, B_longitude]. Only used for magnetic data.
        :type B_b: ``list`` of ``Scalar``
        :param data_offset: offset of the data collection region from the surface
        :type data_offset: ``float``
        :param full_knowledge: if ``True`` data are collected from the entire subsurface region. This is mainly for testing.
        :type full_knowledge: ``Bool``
        :param spherical: if ``True`` sperical coordinates are used (ignored)
        :type spherical: ``Bool``
        """      
        super(SyntheticData,self).__init__(
                                 datatype=datatype, DIM=DIM, number_of_elements=number_of_elements, 
                                 length=length, B_b=B_b, 
                                 data_offset=data_offset,
                                 full_knowledge=full_knowledge,
                                 spherical=spherical)
        self.__n_length = n_length
        self.__n_depth = n_depth 
        self.__shift = shift
        if amplitude == None:
	    if datatype == DataSource.GRAVITY:
	        amplitude = 200 *U.kg/U.m**3
	    else:
	         amplitude =0.1
        self.__amplitude = amplitude



    def getReferenceProperty(self, domain=None):
        """
        Returns the reference density Data object that was used to generate
        the gravity anomaly data.
        """
        if self._reference_data == None:
            DIM=domain.getDim()
            x=domain.getX()
            # set the reference data
            z=x[DIM-1]
            k=sin(self.__n_depth * np.pi/inf(z) * z) * whereNegative(z) * self.__amplitude 
            for i in xrange(DIM-1):
	       k*= sin(self.__n_length * np.pi /self.__length * x[i] - shift )
            self._reference_data= k
        return self._reference_data


