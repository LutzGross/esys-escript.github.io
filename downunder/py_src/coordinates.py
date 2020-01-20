
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

"""Functions to deal with coordinate systems"""

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['ReferenceSystem', 'CartesianReferenceSystem',
    'GeodeticReferenceSystem', 'SphericalReferenceSystem',
    'WGS84ReferenceSystem', 'GRS80ReferenceSystem',
    'SpatialCoordinateTransformation', 'GeodeticCoordinateTransformation',
    'CartesianCoordinateTransformation', 'makeTransformation']

from esys.escript import unitsSI as U
import esys.escript as esc

class ReferenceSystem(object):
    """
    Generic identifier for coordinate systems.
    """
    def __init__(self,name='none'):
        """
        initialization of reference system

        :param name: name of the reference system
        :type name: ``str``
        """
        self.__name=name

    def getName(self):
        """
        returns the name of the reference system
        """
        return self.__name

    def __str__(self):
        return ("%s (id %s)"%(self.getName(),id(self)))

    def __eq__(self, other):
        return self.isTheSame(other)

    def __ne__(self, other):
        return not self.isTheSame(other)

    def isTheSame(self, other):
        """
        test if argument ``other`` defines the same reference system

        :param other: a second reference system
        :type other: `ReferenceSystem`
        :returns: ``True`` if other defines the same reference system
        :rtype: ``bool``
        
        .. note:: needs to be overwritten by a particular reference system
        """
        raise NotImplementedError()

    def isCartesian(self):
        """
        returns if the reference system is Cartesian
        
        .. note:: needs to be overwritten by a particular reference system
        
        :rtype: ``bool``
        """
        raise NotImplementedError()
      
    def createTransformation(self, domain):
        """
        creates an appropriate coordinate transformation on a given domain

        .. note:: needs to be overwritten by a particular reference system
                
        :param domain: domain of transformation
        :type domain: `esys.escript.AbstractDomain`
        :rtype: `SpatialCoordinateTransformation`
        """
        raise NotImplementedError()

class CartesianReferenceSystem(ReferenceSystem):
    """
    Identifies the Cartesian coordinate system
    """
    def __init__(self, name="CARTESIAN"):
        """
        set up Cartesian coordinate system
        """
        super(CartesianReferenceSystem, self).__init__(name)

    def isTheSame(self, other):
        """
        test if argument ``other`` defines the same reference system

        :param other: a second reference system
        :type other: `ReferenceSystem`
        :returns: ``True`` if ``other`` is a `CartesianReferenceSystem` instance.
        :rtype: ``bool``
        :note: every two `CartesianReferenceSystem` instances are considered
               as being the same.
        """
        return isinstance(other, CartesianReferenceSystem)

    def createTransformation(self, domain):
        """
        creates an appropriate coordinate transformation on a given domain

        :param domain: domain of transformation
        :type domain: `esys.escript.AbstractDomain`
        :rtype: `SpatialCoordinateTransformation`
        """
        return SpatialCoordinateTransformation(domain, reference=self)
      
    def isCartesian(self):
        """
        returns if the reference system is Cartesian
        
        :rtype: ``bool``
        """
        return True
      
class GeodeticReferenceSystem(ReferenceSystem):
    """
    Identifies a Geodetic coordinate system
    """
    def __init__(self, a=6378137.0 *U.m, f=1/298.257223563, angular_unit=1.*U.DEG, height_unit=1.*U.km, name="WGS84"):
        """
        initializes a geodetic reference system

        :param a: semi-major axis in meter
        :type a: positive ``double``
        :param f: flattening
        :type f: non-negative ``double``, less than one
        :param name: name of the reference system
        :type name: ``str``
        :param angular_unit: factor to scale the unit of latitude and
                             longitude to radians.
        :type  angular_unit: positive ``double``
        :param height_unit: factor to scale the unit of latitude and
                             longitude to radians.
        :type  height_unit: positive ``double``
        """
        if not a>0:
            raise ValueError("length of semi-major axis a must be positive.")
        if not ( f>=0 and f<1 ):
            raise ValueError("flattening f must be non-negative and less than one.")
        if not angular_unit > 0:
            raise ValueError("angular_unit must be positive.")
        if not height_unit > 0:
            raise ValueError("height_unit must be positive.")  
        super(GeodeticReferenceSystem, self).__init__(name)
      
        self.__a=a
        self.__f=f
        self.__angular_unit=angular_unit
        self.__height_unit=height_unit
        

    def isCartesian(self):
        """
        returns if the reference system is Cartesian
        
        :rtype: ``bool``
        """
        return False
      
    def getAngularUnit(self):
        """
        returns the angular unit
        """
        return self.__angular_unit
        
    def getHeightUnit(self):
        """
        returns the height unit
        """
        return self.__height_unit
        
    def getSemiMajorAxis(self):
        """
        returns the length of semi major axis
        """
        return self.__a

    def getSemiMinorAxis(self):
        """
        returns the length of semi minor axis
        """
        a=self.getSemiMajorAxis()
        f=self.getFlattening()
        return a*(1-f)

    def getFlattening(self):
        """
        returns the flattening
        """
        return self.__f

    def isTheSame(self, other):
        """
        test if ``other`` defines the same reference system

        :param other: a second reference system
        :type other: `ReferenceSystem`
        :returns: ``True`` if other defines then same reference system
        :rtype: ``bool``
        
        .. note:: two `GeodeticReferenceSystem` are considered to be the same
               if the use the same semi major axis, the same flattening
               and the same angular unit.
        """
        if isinstance(other, GeodeticReferenceSystem):
            if self.getSemiMajorAxis() == other.getSemiMajorAxis() \
                      and self.getFlattening() == other.getFlattening() \
                      and self.getAngularUnit() == other.getAngularUnit():
                return True
            else:
                return False
        else:
            return False

    def createTransformation(self, domain):
        """
        creates an appropriate coordinate transformation on a given domain

        :param domain: domain of transformation
        :type domain: `esys.escript.AbstractDomain`
        :rtype: `SpatialCoordinateTransformation`
        """
        return GeodeticCoordinateTransformation(domain, reference=self)

def SphericalReferenceSystem(R=6378137.0*U.m):
    """
    returns the `GeodeticReferenceSystem` of a sphere
    :param R: sphere radius
    :type R: positive ``double``
    """
    return GeodeticReferenceSystem(a=R, f=0, angular_unit=1*U.DEG,  height_unit=1.*U.km, name="SPHERE")

def WGS84ReferenceSystem():
    """
    returns the `GeodeticReferenceSystem` for the WGS84 Ellipsoid
    """
    return GeodeticReferenceSystem(a=6378137.0 *U.m, f=1/298.257223563, angular_unit=1*U.DEG,  height_unit=100.*U.km, name="WGS84")

def GRS80ReferenceSystem():
    """
    returns the `GeodeticReferenceSystem` for the GRS80 Ellipsoid eg. used by Geocentric Datum of Australia GDA94 
    """
    return GeodeticReferenceSystem(a=6378137.0 *U.m, f=1/298.257222101, angular_unit=1*U.DEG,  height_unit=1.*U.km, name="GRS80")


class SpatialCoordinateTransformation(object):
    """
    Defines an orthogonal coordinate transformation from a domain into the
    Cartesian domain using a coordinate transformation.

    The default implementation is the identity transformation (i.e.
    no changes are applied to the domain). Overwrite the appropriate
    methods to define other coordinate systems.
    """
    def __init__(self, domain, reference=CartesianReferenceSystem()):
        """
        set up the orthogonal coordinate transformation.

        :param domain: domain in the domain of the coordinate transformation
        :type domain: `esys.escript.AbstractDomain`
        :param reference: the reference system
        :type reference: `ReferenceSystem`
        """
        self.__domain = domain
        self.__reference_system=reference
        self._volumefactor=esc.Scalar(1., esc.Function(domain))
        self._scaling_factors = esc.Vector(1., esc.Function(domain))

    def __eq__(self, other):
         return self.isTheSame(other)

    def __ne__(self, other):
         return not self.isTheSame(other)

    def isTheSame(self, other):
        """
        test if argument ``other`` defines the same coordinate transformation

        :param other: a second coordinate transformation
        :type other: `SpatialCoordinateTransformation`
        :returns: ``True`` if other defines then same coordinate transformation
        :rtype: ``bool``
        """
        if isinstance(other, SpatialCoordinateTransformation):
            if self.getDomain() == other.getDomain() \
                    and self.getReferenceSystem() == other.getReferenceSystem():
                return True

        return False

    def isCartesian(self):
        """
        returns ``True`` if the scaling factors (and the volume factor) are equal to 1
        
        :rtype: ``bool``
        """
        return self.__reference_system.isCartesian()

    def getDomain(self):
        """
        returns the domain of the coordinate transformation.
        
        :rtype: `esys.escript.AbstractDomain`
        """
        return self.__domain

    def getReferenceSystem(self):
        """
        returns the reference system used to to define the coordinate transformation
        
        :rtype: `ReferenceSystem`
        """
        return self.__reference_system

    def getVolumeFactor(self):
        """
        returns the volume factor for the coordinate transformation
        
        :rtype:  `esys.escript.Scalar`
        """
        return self._volumefactor

    def getScalingFactors(self):
        """
        returns the scaling factors
        
        :rtype: `esys.escript.Vector`
        """
        return self._scaling_factors

    def getGradient(self, u):
        """
        returns the gradient of a scalar function in direction of the
        coordinate axis.
        
        :rtype: `esys.escript.Vector`
        """
        g=esc.grad(u)
        if not self.isCartesian():
                d=self.getScalingFactors()
                g*=d
        return g



def CartesianCoordinateTransformation(domain, reference=CartesianReferenceSystem() ):

    return SpatialCoordinateTransformation(domain, reference)

class GeodeticCoordinateTransformation(SpatialCoordinateTransformation):
    """
    A geodetic coordinate transformation
    """
    def __init__(self, domain, reference=WGS84ReferenceSystem() ):
        """
        set up the orthogonal coordinate transformation.

        :param domain: domain in the domain of the coordinate transformation
        :type domain: `esys.escript.AbstractDomain`
        :param reference: the reference system
        :type reference: `ReferenceSystem`

        """
        DIM=domain.getDim()
        super(GeodeticCoordinateTransformation, self).__init__(domain, reference )

        a=reference.getSemiMajorAxis()
        f=reference.getFlattening()
        f_a=reference.getAngularUnit()
        f_h=reference.getHeightUnit()

        x=esc.Function(domain).getX()
        if DIM == 2:
            phi=0.
        else:
            phi=x[1] * f_a
        h=x[DIM-1] * f_h

        e = esc.sqrt(2*f-f**2)
        N = a/esc.sqrt(1 - e**2 * esc.sin(phi)**2 )
        M = ( a*(1-e**2) ) /esc.sqrt(1 - e**2 * esc.sin(phi)**2 )**3
        v_phi = f_a * (M + h) 
        v_lam = f_a * (N + h) * esc.cos(phi)
        v_h = f_h
        s= esc.Vector(1., esc.Function(domain)) 
        if DIM == 2:
            v= v_phi * v_h 
            s[0]=1/v_lam
            s[1]=1/v_h
        else:
            v= v_phi * v_lam * v_h
            s[0]=1/v_lam
            s[1]=1/v_phi
            s[2]=1/v_h

        self._volumefactor=v
        self._scaling_factors = s

def makeTransformation(domain, coordinates=None):
    """
    returns a `SpatialCoordinateTransformation` for the given domain
    
    :param domain: domain in the domain of the coordinate transformation
    :type domain: `esys.escript.AbstractDomain`
    :param coordinates: the reference system or spatial coordinate system.
    :type coordinates: `ReferenceSystem` or `SpatialCoordinateTransformation`
    :return: the spatial coordinate system for the given domain of the specified 
             reference system ``coordinates``. If ``coordinates`` is already spatial coordinate system based on the 
             riven domain ``coordinates`` is returned. Otherwise an appropriate spatial coordinate system 
             is created.
    :rtype: `SpatialCoordinateTransformation` 
    """
    if coordinates == None:
        return CartesianCoordinateTransformation(domain)
    elif isinstance(coordinates, ReferenceSystem):
        return coordinates.createTransformation(domain)
    else:
        if not coordinates.getDomain() == domain:
            raise ValueError("Domain of spatial coordinate system and given domain don't match.")
        else:
            return coordinates

