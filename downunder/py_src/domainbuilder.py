
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

"""Domain construction from survey data for inversions"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['DomainBuilder']

import logging
import numpy as np
import esys.escript.util as esu
from esys.escript import unitsSI as U
from .datasources import DataSource
from .coordinates import ReferenceSystem, CartesianReferenceSystem

try:
    from esys.ripley import Rectangle, Brick
    HAVE_RIPLEY = True
except ImportError:
    HAVE_RIPLEY = False

class DomainBuilder(object):
    """
    This class is responsible for constructing an escript Domain object with
    suitable extents and resolution for survey data (`DataSource` objects)
    that are added to it.

    The domain covers a region above and below the Earth surface. The
    East-West direction is used as the x- or longitudinal or x[0] direction,
    the North-South direction is used as the y- or latitudinal or x[1]
    direction, the vertical direction is denoted by z or radial or x[2]
    direction. The corresponding terms are used synonymously. 
    """

    def __init__(self, dim=3, reference_system=None):
        """
        Constructor.

        :param dim: Dimensionality (2 or 3) of the target domain.
                    This has implications for the survey data than can be
                    added. By default a 3D domain is created.
        :type dim: ``int``
        :param reference_system: reference coordinate system. By default the 
                                 Cartesian coordinate system is used.
        :type reference_system: `ReferenceSystem`
        """

        if not HAVE_RIPLEY:
            raise ImportError("Ripley module not available")
        self.logger = logging.getLogger('inv.%s'%self.__class__.__name__)
        if dim not in (2,3):
            raise ValueError("Number of dimensions must be 2 or 3")
        if not reference_system:
            self.__reference_system=CartesianReferenceSystem()
        else:
            self.__reference_system=reference_system 

        if self.__reference_system.isCartesian():
            self.__v_scale=1.
        else:
            self.__v_scale=1./self.getReferenceSystem().getHeightUnit()

        self.__domain=None
        self.__dim=dim
        self.__sources=[]
        self.__background_magnetic_field=None
        # list of all tags used by all data sources being attached
        self.__tags=[]
        self.setElementPadding()
        self.setVerticalExtents()
        self.fixDensityBelow()
        self.fixSusceptibilityBelow()
        self.fixVelocityBelow()
        
    def getReferenceSystem(self):
        """
        returns the reference coordinate system
        
        :rtype: `ReferenceSystem`
        """
        return self.__reference_system

    def getTags(self):
        """
        returns a list of all tags in use by the attached data sources.
        The list may be empty.
        """
        return self.__tags
            
    def addSource(self, source):
        """
        Adds a survey data provider to the domain builder.
        An exception is raised if the domain has already been built.
        An exception is also reported if the reference system used is
        cartesian and the UTM zone of `source` does not match the UTM zone of
        sources already added to the domain builder (see Inversion Cookbook
        for more information).
        The dimensionality of the data source must be compatible with this
        domain builder. That is, the dimensionality of the data must be one
        less than the dimensionality of the domain (specified in the
        constructor).

        :param source: The data source to be added. Its reference system needs
                       to match the reference system of the DomainBuilder.
        :type source: `DataSource`
        """
        if self.__domain is not None:
            raise RuntimeError("Invalid call to addSource(). Domain is already built.")
        if not isinstance(source, DataSource):
            raise TypeError("source is not a DataSource")
        if not source.getReferenceSystem() == self.getReferenceSystem():
           raise ValueError("source reference system does not match.")

        DATA_DIM = len(source.getDataExtents()[0])
        if DATA_DIM != self.__dim-1:
            raise ValueError("Data must be %d-dimensional."%(self.__dim-1))
        if len(self.__sources)>0 and self.getReferenceSystem().isCartesian():
            if self.__sources[0].getUtmZone() != source.getUtmZone():
                raise ValueError("It is not possible to combine data sources located in different UTM zones at the moment.")

        self.__sources.append(source)
        if source.getTags(): self.__tags=list(set(self.__tags + source.getTags()))

    def setFractionalPadding(self, pad_x=None, pad_y=None, pad_lat=None, pad_lon=None):
        """
        Sets the amount of padding around the dataset as a fraction of the
        dataset side lengths.

        For example, calling ``setFractionalPadding(0.2, 0.1)`` with a data
        source of size 10x20 will result in the padded data set size
        14x24 (10*(1+2*0.2), 20*(1+2*0.1))

        :param pad_x: Padding per side in x direction (default: no padding)
        :type pad_x: ``float``
        :param pad_y: Padding per side in y direction (default: no padding)
        :type pad_y: ``float``
        :param pad_lat: Padding per side in latitudinal direction (default: no padding)
        :type pad_lat: ``float``
        :param pad_lon: Padding per side in longitudinal direction (default: no padding)
        :type pad_lon: ``float``        
        :note: `pad_y` is ignored for 2-dimensional domains. 
        """
        if not pad_lat == None:
            if not pad_x == None:
               raise ValueError("Either pad_lat or pad_x can be set.")
            else:
              pad_x = pad_lat
        if not pad_lon == None:
            if not pad_y == None:
              raise ValueError("Either pad_lon or pad_y can be set.")
            else:
              pad_y = pad_lan         
        if self.__domain is not None:
            raise RuntimeError("Invalid call to setFractionalPadding(). Domain is already built.")
        if pad_x is not None:
            if pad_x < 0:
                raise ValueError("setFractionalPadding: Arguments must be non-negative")
            if pad_x > 10:
                raise ValueError("setFractionalPadding: Argument too large")
        if pad_y is not None:
            if pad_y < 0:
                raise ValueError("setFractionalPadding: Arguments must be non-negative")
            if pad_y > 10:
                raise ValueError("setFractionalPadding: Argument too large")
        self._padding = [pad_x,pad_y], 'f'

    def setPadding(self, pad_x=None, pad_y=None,  pad_lat=None, pad_lon=None):
        """
        Sets the amount of padding around the dataset in absolute length units.

        The final domain size will be the length in x (in y) of the dataset
        plus twice the value of `pad_x` (`pad_y`). The arguments must be
        non-negative.

        :param pad_x: Padding per side in x direction (default: no padding)
        :type pad_x: ``float`` in units of length (meter)
        :param pad_y: Padding per side in y direction (default: no padding)
        :type pad_y: ``float`` in units of length (meter)
        :note: `pad_y` is ignored for 2-dimensional domains.
        :note: this function can only be used if the reference system is Cartesian
        """
        if not self.getReferenceSystem().isCartesian():
            raise RuntimeError("setPadding can be called for the Cartesian reference system only.")
        if self.__domain is not None:
            raise RuntimeError("Invalid call to setPadding(). Domain is already built.")
        if pad_x is not None:
            if pad_x < 0:
                raise ValueError("setPadding: Arguments must be non-negative")
        if pad_y is not None:
            if pad_y < 0:
                raise ValueError("setPadding: Arguments must be non-negative")
        self._padding = [pad_x,pad_y], 'l'
        
    def setGeoPadding(self, pad_lat=None, pad_lon=None):
        """
        Sets the amount of padding around the dataset in longitude and latitude.

        The final domain size will be the extent in the latitudinal (in
        longitudinal) direction of the dataset plus twice the value of
        `pad_lat` (`pad_lon`). The arguments must be non-negative.

        :param pad_lat: Padding per side in latitudinal direction (default: 0)
        :type pad_lat: ``float`` in units of degree 
        :param pad_lon: Padding per side in longitudinal direction (default: 0)
        :type pad_lon: ``float``  in units of degree  
        :note: `pad_lon` is ignored for 2-dimensional domains.
        :note: this function can only be used if the reference system is not Cartesian
        """
        if self.getReferenceSystem().isCartesian():
            raise RuntimeError("setGeoPadding can be called for non-Cartesian reference systems only.")
        if self.__domain is not None:
            raise RuntimeError("Invalid call to setPadding(). Domain is already built.")
        if pad_lat is not None:
            if pad_lat < 0:
                raise ValueError("setPadding: Arguments must be non-negative")
        if pad_lon is not None:
            if pad_lon < 0:
                raise ValueError("setPadding: Arguments must be non-negative")
        self._padding = [pad_lat,pad_lon], 'd'
        
    def setElementPadding(self, pad_x=None, pad_y=None, pad_lat=None, pad_lon=None):
        """
        Sets the amount of padding around the dataset in number of elements
        (cells).

        When the domain is constructed `pad_x` (`pad_y`) elements are added
        on each side of the x- (y-) dimension. The arguments must be
        non-negative.

        :param pad_x: Padding per side in x direction (default: no padding)
        :type pad_x: ``int``
        :param pad_y: Padding per side in y direction (default: no padding)
        :type pad_y: ``int``
        :note: `pad_y` is ignored for 2-dimensional datasets.
        """
        if not pad_lat == None:
            if not pad_x == None:
              raise ValueError("Either pad_lat or pad_x can be set.")
            else:
              pad_x = pad_lat
        if not pad_lon == None:
            if not pad_y == None:
              raise ValueError("Either pad_lon or pad_y can be set.")
            else:
              pad_y = pad_lan
              
        if self.__domain is not None:
            raise RuntimeError("Invalid call to setElementPadding(). Domain is already built.")
        if pad_x is not None:
            if type(pad_x) is not int:
                raise TypeError("setElementPadding expects integer arguments")
            if pad_x < 0:
                raise ValueError("setElementPadding: Arguments must be non-negative")
        if pad_y is not None:
            if type(pad_y) is not int:
                raise TypeError("setElementPadding expects integer arguments")
            if pad_y < 0:
                raise ValueError("setElementPadding: Arguments must be non-negative")
        self._padding = [pad_x,pad_y], 'e'

    def getGravitySurveys(self):
        """
        Returns a list of gravity surveys, see `getSurveys` for details.
        """
        return self.getSurveys(DataSource.GRAVITY)

    def getMagneticSurveys(self):
        """
        Returns a list of magnetic surveys, see `getSurveys` for details.
        """
        return self.getSurveys(DataSource.MAGNETIC)

        
    def fixDensityBelow(self, depth=None):
        """
        Defines the depth below which the density anomaly is set to a given
        value. If no value is given zero is assumed. 
        
        :param depth: depth below which the density is fixed. If not set, no
                      constraint at depth is applied.
        :type depth: ``float``
        """
        self.__fix_density_below=depth

    def fixSusceptibilityBelow(self, depth=None):
        """
        Defines the depth below which the susceptibility anomaly is set to a
        given value. If no value is given zero is assumed. 
        
        :param depth: depth below which the susceptibility is fixed. If not
                      set, no constraint at depth is applied.
        :type depth: ``float``
        """
        self.__fix_susceptibility_below=depth

    def fixVelocityBelow(self, depth=None):
        """
        Defines the depth below which the velocity and Q index is set to a
        given value. If no value is given zero is assumed. 
        
        :param depth: depth below which the velocity is fixed. If not
                      set, no constraint at depth is applied.
        :type depth: ``float``
        """
        self.__fix_velocity_below=depth


    def getSurveys(self, datatype, tags=None):
        """
        Returns a list of `Data` objects for all surveys of type `datatype`
        available to this domain builder. If a list of `tags` is given 
        only data sources whose tag matches the tag list are returned.

        :return: List of surveys which are tuples (anomaly,error).
        :rtype: ``list``
        """
        surveys=[]
        for src in self.__sources:
            if src.getDataType()==datatype:
                if tags is None or ( src.getTags() is not None and all( [ t in tags for t in src.getTags() ] )  ) :
                    surveys.append(src.getSurveyData(self.getDomain(), self._dom_origin, self._dom_NE, self._spacing))
        return surveys

    def setBackgroundMagneticFluxDensity(self, B):
        """
        Sets the background magnetic flux density B=(B_East, B_North, B_Vertical)
        """
        self.__background_magnetic_field=B

    def getBackgroundMagneticFluxDensity(self):
        """
        Returns the background magnetic flux density.
        """
        B = self.__background_magnetic_field
        if B is None:
            raise ValueError("No background magnetic flux density set!")

        if self.__dim < 3 :
            return np.array([B[0], B[2]])
        else:
            return np.array(B)

    def getSetDensityMask(self):
        """
        Returns the density mask data object which is non-zero for cells
        whose density value is fixed, zero otherwise.
        """
        z=self.getDomain().getX()[self.__dim-1]
        m = esu.whereNonNegative(z)
        if self.__fix_density_below:
            m += esu.whereNonPositive(z+self.__v_scale*self.__fix_density_below)
        return m

    def getSetSusceptibilityMask(self):
        """
        Returns the susceptibility mask data object which is non-zero for
        cells whose susceptibility value is fixed, zero otherwise.
        """
        z=self.getDomain().getX()[self.__dim-1]
        m = esu.whereNonNegative(z)
        if self.__fix_susceptibility_below:
            m += esu.whereNonPositive(z+self.__v_scale*self.__fix_susceptibility_below)
        return m

    def getDomain(self):
        """
        Returns a domain that spans the data area plus padding.

        The domain is created the first time this method is called,
        subsequent calls return the same domain so anything that affects
        the domain (such as padding) needs to be set beforehand.

        :return: The escript domain for this data source
        :rtype: `esys.escript.Domain`
        """
        if self.__domain is None:
            self.__domain=self.__createDomain()
        return self.__domain

    def setVerticalExtents(self, depth=40000., air_layer=10000., num_cells=25):
        """
        This method sets the target domain parameters for the vertical
        dimension.

        :param depth: Depth of the domain (in meters)
        :type depth: ``float``
        :param air_layer: Depth of the layer above sea level (in meters)
        :type air_layer: ``float``
        :param num_cells: Number of domain elements for the entire vertical
                          dimension
        :type num_cells: ``int``
        """
        if self.__domain is not None:
            raise RuntimeError("Invalid call to setVerticalExtents(). Domain is already built.")
        self._v_depth=depth
        self._v_air_layer=air_layer
        self._v_num_cells=num_cells

    def __getTotalExtentsWithPadding(self):
        """
        Helper method that computes origin and number of data elements
        after adding padding to the bounding box of all available survey data.
        """
        X0, NX, DX = self.__getTotalExtents()
        DATA_DIM=len(X0)
        frac=[]
        # padding is applied to each side so multiply by 2 to get the total
        # amount of padding per dimension
        pad, pt = self._padding
        for i in range(DATA_DIM):
            if pad[i] is None:
                frac.append(0.)
                continue
            if pt == 'f' : # fraction of side length
                frac.append(2.*pad[i])
            elif pt == 'e': # number of elements
                frac.append(2.*pad[i]/float(NX[i]))
            else: # absolute length
                f=pad[i]/DX[i]
                frac.append(2.*f/float(NX[i]))

        # calculate new number of elements
        NX_padded=[int(round(NX[i]*(1+frac[i]))) for i in range(DATA_DIM)]
        NXdiff=[NX_padded[i]-NX[i] for i in range(DATA_DIM)]
        X0_padded=[X0[i]-NXdiff[i]/2.*DX[i] for i in range(DATA_DIM)]
        return X0_padded, NX_padded, DX

    def __getTotalExtents(self):
        """
        Helper method that computes the origin, number of elements and
        minimal element spacing taking into account all available survey data.
        """
        if len(self.__sources)==0:
            raise ValueError("No data")
        X0, NE, DX = self.__sources[0].getDataExtents()
        # do not mess with the values if only one source used
        if len(self.__sources)>1:
            XN=[X0[i]+NE[i]*DX[i] for i in range(len(NE))]

            for src in self.__sources[1:]:
                d_x0, d_ne, d_dx = src.getDataExtents()
                for i in range(len(d_x0)):
                    X0[i]=min(X0[i], d_x0[i])
                for i in range(len(d_dx)):
                    DX[i]=min(DX[i], d_dx[i])
                for i in range(len(d_ne)):
                    XN[i]=max(XN[i], d_x0[i]+d_ne[i]*d_dx[i])
            # FIXME: should this be rounded up instead?
            NE=[int((XN[i]-X0[i])/DX[i]) for i in range(len(XN))]
        return X0, NE, DX

    def __createDomain(self):
        """
        Creates and returns an escript domain that spans the entire area of
        available data plus a padding zone. This method is called only once
        the first time `getDomain()` is invoked.

        :return: The escript domain
        :rtype: `esys.escript.Domain`
        """
        X0, NX, DX = self.__getTotalExtentsWithPadding()

        # number of domain elements
        NE = NX + [self._v_num_cells]

        # origin of domain
        origin = X0 + [-self._v_depth*self.__v_scale]

        if self.getReferenceSystem().isCartesian():
            # rounding will give us about meter-accuracy with UTM coordinates
            self._dom_origin = [np.floor(oi) for oi in origin]
        else:
            # this should give us about meter-accuracy with lat/lon coords
            self._dom_origin = [1e-5*np.floor(oi*1e5) for oi in origin]

        # cell size / point spacing
        spacing = DX + [self.__v_scale*np.floor((self._v_depth+self._v_air_layer)/self._v_num_cells)]
        #self._spacing = [float(np.floor(si)) for si in spacing]
        self._spacing = spacing

        lo=[(self._dom_origin[i], self._dom_origin[i]+NE[i]*self._spacing[i]) for i in range(self.__dim)]

        if self.__dim==3:
            dom=Brick(*NE, l0=lo[0], l1=lo[1], l2=lo[2])
        else:
            dom=Rectangle(*NE, l0=lo[0], l1=lo[1])

        # ripley may internally adjust NE and length, so recompute
        self._dom_len=[esu.sup(dom.getX()[i])-esu.inf(dom.getX()[i]) for i in range(self.__dim)]
        self._dom_NE=[int(self._dom_len[i]/self._spacing[i]) for i in range(self.__dim)]

        self.logger.debug("Domain size: "+str(self._dom_NE))
        self.logger.debug("     length: "+str(self._dom_len))
        self.logger.debug("     origin: "+str(self._dom_origin))
        return dom

