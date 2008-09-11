#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
Geometrical Primitives

the concept is inspired by gmsh and very much focused on the fact that
the classes are used to wrk with gmsh.

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""


__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"

try:
   import numpy
   numpyImported=True
except:
   numpyImported=False   

import numarray
from transformations import _TYPE, Translation, Dilation, Transformation
from math import sqrt


def resetGlobalPrimitiveIdCounter():
   """
   initializes the global primitive ID counter
   """
   global global_primitive_id_counter
   global_primitive_id_counter=1

def setToleranceForColocation(tol=1.e-11):
   """
   set the global tolerance for colocation checks to tol
   """
   global global_tolerance_for_colocation
   global_tolerance_for_colocation=tol

def getToleranceForColocation():
   """
   returns the global tolerance for colocation checks
   """
   return global_tolerance_for_colocation

resetGlobalPrimitiveIdCounter()
setToleranceForColocation()


class PrimitiveBase(object):
    """
    template for a set of primitives 
    """
    def __init__(self): 
       """
       initializes PrimitiveBase instance object with id
       """ 
       pass

    def __cmp__(self,other):
       """
       compares object with other by comparing the absolute value of the ID
       """
       if isinstance(other, PrimitiveBase):
           return cmp(self.getID(),other.getID())
       else:
           return False
    def getConstructionPoints(self):
        """
        returns the points used to construct the primitive
        """
        out=[]
        for i in self.getPrimitives(): 
           if isinstance(i,Point): out.append(i)
        return out

    def getPrimitives(self):
        """
        returns a list of primitives used to construct the primitive with no double entries
        """
        out=[]
        for p in self.collectPrimitiveBases():
            if not p  in out: out.append(p)
        return out

    def copy(self):
       """
       returns a deep copy of the object
       """
       return self.substitute({})

    def modifyBy(self,transformation):
       """
       modifies the coordinates by applying a transformation 
       """
       for p in self.getConstructionPoints(): p.modifyBy(transformation)

    def __add__(self,other):
        """
        returns a new object shifted by other
        """
        return self.apply(Translation(numarray.array(other,_TYPE)))

    def __sub__(self,other):
        """
        returns a new object shifted by other
        """
        return self.apply(Translation(-numarray.array(other,_TYPE)))

    def __iadd__(self,other):
        """
        shifts the point by other
        """
        self.modifyBy(Translation(numarray.array(other,_TYPE)))
        return self

    def __isub__(self,other):
        """
        shifts the point by -other
        """
        self.modifyBy(Translation(-numarray.array(other,_TYPE)))
        return self

    def __imul__(self,other):
        """
        modifies object by applying L{Transformation} other. If other is not a L{Transformation} it will try convert it.
        """
        if isinstance(other,int) or isinstance(other,float):
            trafo=Dilation(other)
        elif isinstance(other,numarray.NumArray):
            trafo=Translation(other)
        elif isinstance(other,Transformation):
            trafo=other
        else:
            raise TypeError, "cannot convert argument to Trnsformation class object."
        self.modifyBy(trafo)
        return self

    def __rmul__(self,other):
        """
        applies L{Transformation} other to object. If other is not a L{Transformation} it will try convert it.
        """
        if isinstance(other,int) or isinstance(other,float):
            trafo=Dilation(other)
        elif isinstance(other,numarray.NumArray):
            trafo=Translation(other)
        elif isinstance(other,Transformation):
            trafo=other
        else:
            raise TypeError, "cannot convert argument to Transformation class object."
        return self.apply(trafo)


    def setLocalScale(self,factor=1.):
       """
       sets the local refinement factor
       """
       for p in self.getConstructionPoints(): p.setLocalScale(factor)

    def apply(self,transformation):
        """
        returns a new object by applying the transformation
        """
        out=self.copy()
        out.modifyBy(transformation)
        return out

class Primitive(object):
    """
    A general primitive
    """
    def __init__(self): 
       """
       initializes PrimitiveBase instance object with id
       """ 
       global global_primitive_id_counter
       self.__ID=global_primitive_id_counter
       global_primitive_id_counter+=1

    def getID(self):
       """
       returns the primitive ID
       """
       return self.__ID

    def getDirectedID(self):
        """
        returns the primitive ID where a negative signs means that the reversed ordring is used.
        """
        return self.getID()

    def __repr__(self):
       return "%s(%s)"%(self.__class__.__name__,self.getID())

    def getUnderlyingPrimitive(self):
        """
        returns the underlying primitive
        """
        return self
    def hasSameOrientation(self,other):
        """
        returns True if other is the same primitive and has the same orientation
        """
        return self == other and isinstance(other,Primitive)

    def __neg__(self):
        """
        returns a view onto the curve with reversed ordering

        @note: this class is overwritten by subclass
        """
        raise NotImplementedError("__neg__ is not implemented.")

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.

        @note: this class is overwritten by subclass
        """
        raise NotImplementedError("substitute is not implemented.")

    def collectPrimitiveBases(self):
        """
        returns a list of primitives used to construct the primitive. It may contain primitives twice
        
        @note: this class is overwritten by subclass
        """
        raise NotImplementedError("collectPrimitiveBases is not implemented.")

    def isColocated(self,primitive):
       """
       returns True is the two primitives are located at the smae position

       @note: this class is overwritten by subclass
       """
       raise NotImplementedError("isColocated is not implemented.")


class ReversePrimitive(object):
    """
    A view onto a primitive creating an reverse orientation
    """
    def __init__(self,primitive): 
       """
       instantiate a view onto primitve
       """ 
       if not isinstance(primitive, Primitive):
           raise ValueError("argument needs to be a Primitive class object.")
       self.__primitive=primitive

    def getID(self):
       """
       returns the primitive ID
       """
       return self.__primitive.getID()

    def getUnderlyingPrimitive(self):
        """
        returns the underlying primitive
        """
        return self.__primitive

    def hasSameOrientation(self,other):
        """
        returns True if other is the same primitive and has the same orientation
        """
        return self == other and isinstance(other,ReversePrimitive)

    def __repr__(self):
       return "-%s(%s)"%(self.__primitive.__class__.__name__,self.getID())

    def getDirectedID(self):
        """
        returns the primitive ID where a negative signs means that the reversed ordring is used.
        """
        return -self.__primitive.getID()

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            sub_dict[self]=-self.getUnderlyingPrimitive().substitute(sub_dict)
        return sub_dict[self]
            
    def __neg__(self):
          """
          returns a view onto the curve with reversed ordering
          """
          return self.__primitive

    def collectPrimitiveBases(self):
        """
        returns a list of primitives used to construct the primitive. It may contain primitives twice
        """
        return self.__primitive.collectPrimitiveBases()

    def isColocated(self,primitive):
       """
       returns True is the two primitives are located at the smae position

       @note: this class is overwritten by subclass
       """
       return self.__primitive.isColocated(primitive)

class Point(Primitive, PrimitiveBase):
    """
    a three dimensional point
    """
    def __init__(self,x=0.,y=0.,z=0.,local_scale=1.): 
       """
       creates a point with coorinates x,y,z with the local refinement factor local_scale
       """ 
       PrimitiveBase.__init__(self)
       Primitive.__init__(self)
       self.setCoordinates(numarray.array([x,y,z],_TYPE))
       self.setLocalScale(local_scale)

    def setLocalScale(self,factor=1.):
       """
       sets the local refinement factor
       """
       if factor<=0.:
          raise ValueError("scaling factor must be positive.")
       self.__local_scale=factor

    def getLocalScale(self):
       """
       returns the local refinement factor
       """
       return self.__local_scale
    def getCoordinates(self):
       """
       returns the coodinates of the point as L{numarray.NumArray} object
       """
       return self._x
    def setCoordinates(self,x):
       """
       returns the coodinates of the point as L{numarray.NumArray} object
       """
       if not isinstance(x, numarray.NumArray):
          self._x=numarray.array(x,_TYPE)
       else:
          self._x=x

    def collectPrimitiveBases(self):
       """
       returns primitives used to construct the primitive
       """
       return [self]
 
    def isColocated(self,primitive):
       """
       returns True if L{Point} primitive is colocation (same coordinates) 
       that means if |self-primitive| <= tol * max(|self|,|primitive|)
       """
       if isinstance(primitive,Point):
          primitive=primitive.getCoordinates()
          c=self.getCoordinates()
          d=c-primitive
          if numpyImported:
            return numpy.dot(d,d)<=getToleranceForColocation()**2*max(numpy.dot(c,c),numpy.dot(primitive,primitive))
          else:
            return numarray.dot(d,d)<=getToleranceForColocation()**2*max(numarray.dot(c,c),numarray.dot(primitive,primitive))
       else:
          return False

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
           c=self.getCoordinates()
           sub_dict[self]=Point(c[0],c[1],c[2],local_scale=self.getLocalScale())
        return sub_dict[self]

    def modifyBy(self,transformation):
        """
        modifies the coordinates by applying a transformation 
        """
        self.setCoordinates(transformation(self.getCoordinates()))


    def __neg__(self):
        """
        returns a view of the object with reverse orientiention. As a point has no direction the object itself is returned.
        """
        return self
       
class Manifold1D(PrimitiveBase):
    """
    general one-dimensional minifold in 3D defined by a start and end point.
    """
    def __init__(self):
        """
        create a one-dimensional manifold
        """
        PrimitiveBase.__init__(self)

    def getStartPoint(self):
         """
         returns start point
         """
         raise NotImplementedError()

    def getEndPoint(self):
         """
         returns end point
         """
         raise NotImplementedError()
    def getBoundary(self):
        """
        returns a list of the zero-dimensional manifolds forming the boundary of the curve 
        """
        return [ self.getStartPoint(), self.getEndPoint()]

class CurveBase(Manifold1D):
    """
    A Curve is defined by a set of control points 
    """
    def __init__(self):
          """
          create curve 
          """
          Manifold1D.__init__(self)

    def __len__(self):
          """
          returns the number of control points
          """
          return len(self.getControlPoints())

    def getStartPoint(self):
         """
         returns start point
         """
         return self.getControlPoints()[0]

    def getEndPoint(self):
         """
         returns end point
         """
         return self.getControlPoints()[-1]

    def getControlPoints(self):
         """
         returns a list of the points
         """
         raise NotImplementedError()

class Curve(CurveBase, Primitive):
    """
    a curve defined through a list of control points. 
    """
    def __init__(self,*points):
       """
       defines a curve form control points 
       """
       if len(points)<2:
           raise ValueError("Curve needs at least two points")
       i=0
       for p in points:
              i+=1
              if not isinstance(p,Point): raise TypeError("%s-th argument is not a Point object."%i)
       self.__points=points
       CurveBase.__init__(self)
       Primitive.__init__(self)

    def getControlPoints(self):
         """
         returns a list of the points
         """
         return self.__points
      
    def __neg__(self):
          """
          returns a view onto the curve with reversed ordering
          """
          return ReverseCurve(self)

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            new_p=[]
            for p in self.getControlPoints(): new_p.append(p.substitute(sub_dict))
            sub_dict[self]=self.__class__(*tuple(new_p))
        return sub_dict[self]

    def collectPrimitiveBases(self):
       """
       returns primitives used to construct the Curve
       """
       out=[self]
       for p in self.getControlPoints(): out+=p.collectPrimitiveBases()
       return out

    def isColocated(self,primitive):
       """
       returns True curves are on the same position
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
         if isinstance(primitive.getUnderlyingPrimitive(),self.__class__): 
           if len(primitive) == len(self):
             cp0=self.getControlPoints()
             cp1=primitive.getControlPoints()
             match=True
             for i in range(len(cp0)):
                if not cp0[i].isColocated(cp1[i]):
                   match=False
                   break
             if not match:
                for i in range(len(cp0)):
                   if not cp0[i].isColocated(cp1[len(cp0)-1-i]):
                      return False
             return True
       return False

class ReverseCurve(CurveBase, ReversePrimitive):
    """
    a curve defined through a list of control points. 
    """
    def __init__(self,curve):
       """
       defines a curve form control points 
       """
       if not isinstance(curve, Curve):
           raise TypeError("ReverseCurve needs to be an instance of Curve")
       CurveBase.__init__(self)
       ReversePrimitive.__init__(self,curve)

    def getControlPoints(self):
         """
         returns a list of the points
         """
         out=[p for p in self.getUnderlyingPrimitive().getControlPoints()]
         out.reverse()
         return out

class Spline(Curve):
    """
    a spline curve defined through a list of control points. 
    """
    pass

class BezierCurve(Curve):
    """
    a Bezier curve
    """
    pass

class BSpline(Curve):
    """
    a BSpline curve. Control points may be repeated.
    """
    pass

class Line(Curve):
    """
    a line is defined by two pointDirecteds
    """
    def __init__(self,*points):
        """
        defines a line with start and end point
        """
        if len(points)!=2:
           raise TypeError("Line needs two points")
        Curve.__init__(self,*points)

class ArcBase(Manifold1D):
    def __init__(self):
          """
          create curve 
          """
          Manifold1D.__init__(self)
    def collectPrimitiveBases(self):
       """
       returns the primitives used to construct the Curve
       """
       out=[self]
       out+=self.getStartPoint().collectPrimitiveBases()
       out+=self.getEndPoint().collectPrimitiveBases()
       out+=self.getCenterPoint().collectPrimitiveBases()
       return out


    def getCenterPoint(self):
         """
         returns center
         """
         raise NotImplementedError()

class Arc(ArcBase, Primitive):
    """
    defines an arc which is strictly, smaller than Pi
    """
    def __init__(self,center,start,end):
       """
       creates an arc by the start point, end point and center
       """
       if not isinstance(center,Point): raise TypeError("center needs to be a Point object.")
       if not isinstance(end,Point): raise TypeError("end needs to be a Point object.")
       if not isinstance(start,Point): raise TypeError("start needs to be a Point object.")
       if center.isColocated(end): raise TypeError("center and start point are colocated.")
       if center.isColocated(start): raise TypeError("center end end point are colocated.")
       if start.isColocated(end): raise TypeError("start and end are colocated.")
       # TODO: check length of circle.
       ArcBase.__init__(self)
       Primitive.__init__(self)
       self.__center=center
       self.__start=start
       self.__end=end
    def __neg__(self):
          """
          returns a view onto the curve with reversed ordering
          """
          return ReverseArc(self)

    def getStartPoint(self):
       """
       returns start point
       """
       return self.__start

    def getEndPoint(self):
       """
       returns end point
       """
       return self.__end

    def getCenterPoint(self):
       """
       returns center
       """
       return self.__center

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            sub_dict[self]=Arc(self.getCenterPoint().substitute(sub_dict),self.getStartPoint().substitute(sub_dict),self.getEndPoint().substitute(sub_dict))
        return sub_dict[self]


    def isColocated(self,primitive):
       """
       returns True curves are on the same position
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
          if isinstance(primitive.getUnderlyingPrimitive(),Arc):
            return (self.getCenterPoint().isColocated(primitive.getCenterPoint())) and ( \
                   (self.getEndPoint().isColocated(primitive.getEndPoint()) and self.getStartPoint().isColocated(primitive.getStartPoint()) ) \
                or (self.getEndPoint().isColocated(primitive.getStartPoint()) and self.getStartPoint().isColocated(primitive.getEndPoint()) ) )
       return False

class ReverseArc(ArcBase, ReversePrimitive):
    """
    defines an arc which is strictly, smaller than Pi
    """
    def __init__(self,arc):
       """
       creates an arc by the start point, end point and center
       """
       if not isinstance(arc, Arc):
           raise TypeError("ReverseCurve needs to be an instance of Arc")
       ArcBase.__init__(self)
       ReversePrimitive.__init__(self,arc)

    def getStartPoint(self):
       """
       returns start point
       """
       return self.getUnderlyingPrimitive().getEndPoint()

    def getEndPoint(self):
       """
       returns end point
       """
       return self.getUnderlyingPrimitive().getStartPoint()

    def getCenterPoint(self):
       """
       returns center
       """
       return self.getUnderlyingPrimitive().getCenterPoint()

class EllipseBase(Manifold1D):
    def __init__(self):
          """
          create ellipse
          """
          Manifold1D.__init__(self)
    def collectPrimitiveBases(self):
       """
       returns the primitives used to construct the Curve
       """
       out=[self]
       out+=self.getStartPoint().collectPrimitiveBases()
       out+=self.getEndPoint().collectPrimitiveBases()
       out+=self.getCenterPoint().collectPrimitiveBases()
       out+=self.getPointOnMainAxis().collectPrimitiveBases()
       return out


class Ellipse(EllipseBase, Primitive):
    """
    defines an ellipse which is strictly, smaller than Pi
    """
    def __init__(self,center,point_on_main_axis,start,end):
       """
       creates an arc by the start point, end point, the center and a point on a main axis.
       """
       if not isinstance(center,Point): raise TypeError("center needs to be a Point object.")
       if not isinstance(end,Point): raise TypeError("end needs to be a Point object.")
       if not isinstance(start,Point): raise TypeError("start needs to be a Point object.")
       if not isinstance(point_on_main_axis,Point): raise TypeError("point on main axis needs to be a Point object.")
       if center.isColocated(end): raise TypeError("center and start point are colocated.")
       if center.isColocated(start): raise TypeError("center end end point are colocated.")
       if center.isColocated(point_on_main_axis): raise TypeError("center and point on main axis are colocated.")
       if start.isColocated(end): raise TypeError("start and end point are colocated.")
       # TODO: check length of circle.
       EllipseBase.__init__(self)
       Primitive.__init__(self)
       self.__center=center
       self.__start=start
       self.__end=end
       self.__point_on_main_axis=point_on_main_axis

    def __neg__(self):
          """
          returns a view onto the curve with reversed ordering
          """
          return ReverseEllipse(self)

    def getStartPoint(self):
       """
       returns start point
       """
       return self.__start

    def getEndPoint(self):
       """
       returns end point
       """
       return self.__end

    def getCenterPoint(self):
       """
       returns center
       """
       return self.__center

    def getPointOnMainAxis(self):
       """
       returns a point on a main axis
       """
       return self.__point_on_main_axis

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            sub_dict[self]=Ellipse(self.getCenterPoint().substitute(sub_dict),
                                   self.getPointOnMainAxis().substitute(sub_dict),
                                   self.getStartPoint().substitute(sub_dict),
                                   self.getEndPoint().substitute(sub_dict))
        return sub_dict[self]


    def isColocated(self,primitive):
       """
       returns True curves are on the same position
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
          if isinstance(primitive.getUnderlyingPrimitive(),Ellipse):
            self_c=self.getCenterPoint().getCoordinates()
            p=self.getPointOnMainAxis().getCoordinates()-self_c
            q=primitive.getPointOnMainAxis().getCoordinates()-self_c
            # are p and q orthogonal or collinear?
            len_p=sqrt(p[0]**2+p[1]**2+p[2]**2)
            len_q=sqrt(q[0]**2+q[1]**2+q[2]**2)
            p_q= abs(p[0]*q[0]+p[1]*q[1]+p[2]*q[2])
            return ((p_q <= getToleranceForColocation() * len_q * p_q) or \
                                   (abs(p_q - len_q * p_q) <= getToleranceForColocation())) and \
                   self.getCenterPoint().isColocated(primitive.getCenterPoint()) and \
                   (                                                                  \
                        (self.getEndPoint().isColocated(primitive.getEndPoint()) and \
                         self.getStartPoint().isColocated(primitive.getStartPoint()) ) \
                                       or                                              \
                         (self.getEndPoint().isColocated(primitive.getStartPoint()) and \
                          self.getStartPoint().isColocated(primitive.getEndPoint()) ) \
                   )
       return False

class ReverseEllipse(EllipseBase, ReversePrimitive):
    """
    defines an arc which is strictly, smaller than Pi
    """
    def __init__(self,arc):
       """
       creates an instance of a reverse view to an ellipse
       """
       if not isinstance(arc, Ellipse):
           raise TypeError("ReverseCurve needs to be an instance of Ellipse")
       EllipseBase.__init__(self)
       ReversePrimitive.__init__(self,arc)

    def getStartPoint(self):
       """
       returns start point
       """
       return self.getUnderlyingPrimitive().getEndPoint()

    def getEndPoint(self):
       """
       returns end point
       """
       return self.getUnderlyingPrimitive().getStartPoint()

    def getCenterPoint(self):
       """
       returns center
       """
       return self.getUnderlyingPrimitive().getCenterPoint()

    def getPointOnMainAxis(self):
       """
       returns a point on a main axis
       """
       return self.getUnderlyingPrimitive().getPointOnMainAxis()


class CurveLoop(Primitive, PrimitiveBase):
    """
    An oriented loop of one-dimensional manifolds (= curves and arcs)

    The loop must be closed and the L{Manifold1D}s should be oriented consistently.
    """
    def __init__(self,*curves):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       if len(curves)<2:
            raise ValueError("at least two curves have to be given.")
       for i in range(len(curves)):
           if not isinstance(curves[i],Manifold1D):
              raise TypeError("%s-th argument is not a Manifold1D object."%i)
       # for the curves a loop:
       used=[ False for i in curves]
       self.__curves=list(curves)
       Primitive.__init__(self)
       PrimitiveBase.__init__(self)

    def getCurves(self):
       """
       returns the curves defining the CurveLoop
       """
       return self.__curves

    def __neg__(self):
       """
       returns a view onto the curve with reversed ordering
       """
       return ReverseCurveLoop(self)

    def __len__(self):
       """
       return the number of curves in the CurveLoop
       """
       return len(self.getCurves())


    def collectPrimitiveBases(self):
       """
       returns primitives used to construct the CurveLoop
       """
       out=[self]
       for c in self.getCurves(): out+=c.collectPrimitiveBases()
       return out

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            new_c=[]
            for c in self.getCurves(): new_c.append(c.substitute(sub_dict))
            sub_dict[self]=CurveLoop(*tuple(new_c))
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True if each curve is colocted with a curve in primitive
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
          if isinstance(primitive.getUnderlyingPrimitive(),CurveLoop):
             if len(primitive) == len(self):
                cp0=self.getCurves()
                cp1=primitive.getCurves()
                for c0 in cp0: 
                    colocated = False
                    for c1 in cp1:
                         colocated = colocated or c0.isColocated(c1)
                    if not colocated: return False
                return True
       return False

class ReverseCurveLoop(ReversePrimitive, PrimitiveBase):
    """
    An oriented loop of one-dimensional manifolds (= curves and arcs)

    The loop must be closed and the one-dimensional manifolds should be oriented consistently.
    """
    def __init__(self,curve_loop):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       if not isinstance(curve_loop, CurveLoop):
           raise TypeError("arguments need to be an instance of CurveLoop.")
       ReversePrimitive.__init__(self, curve_loop)
       PrimitiveBase.__init__(self)

    def getCurves(self):
       """
       returns the curves defining the CurveLoop
       """
       return [ -c for c in  self.getUnderlyingPrimitive().getCurves() ]

    def __len__(self):
        return len(self.getUnderlyingPrimitive())

#=
class Manifold2D(PrimitiveBase):
    """
    general two-dimensional manifold
    """
    def __init__(self):
       """
       create a two-dimensional manifold
       """
       PrimitiveBase.__init__(self)

    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        raise NotImplementedError()

class RuledSurface(Primitive, Manifold2D):
    """
    A ruled surface, i.e., a surface that can be interpolated using transfinite interpolation
    """
    def __init__(self,loop):
       """
       creates a ruled surface with boundary loop

       @param loop: L{CurveLoop} defining the boundary of the surface. 
       """
       if not isinstance(loop.getUnderlyingPrimitive(),CurveLoop):
           raise TypeError("argument loop needs to be a CurveLoop object.")
       if len(loop)<2:
           raise ValueError("the loop must contain at least two Curves.")
       if len(loop)>4:
           raise ValueError("the loop must contain at least three Curves.")
       Primitive.__init__(self)
       Manifold2D.__init__(self)
       self.__loop=loop

    def __neg__(self):
          """
          returns a view onto the suface with reversed ordering
          """
          return ReverseRuledSurface(self)

    def getBoundaryLoop(self):
        """
        returns the loop defining the outer boundary
        """
        return self.__loop

    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        return self.getBoundaryLoop().getCurves()

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            sub_dict[self]=RuledSurface(self.getBoundaryLoop().substitute(sub_dict))
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True if each curve is colocted with a curve in primitive
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
          if isinstance(primitive.getUnderlyingPrimitive(),RuledSurface):
             return self.getBoundaryLoop().isColocated(primitive.getBoundaryLoop())
       return False

    def collectPrimitiveBases(self):
        """
        returns primitives used to construct the Surface
        """
        return [self] + self.getBoundaryLoop().collectPrimitiveBases()

def createRuledSurface(*curves):
      """
      an easier way to create a L{RuledSurface} from given curves.
      """
      return RuledSurface(CurveLoop(*curves))


class ReverseRuledSurface(ReversePrimitive, Manifold2D):
    """
    creates a view onto a L{RuledSurface} but with the reverse orientation
    """
    def __init__(self,surface):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       if not isinstance(surface, RuledSurface):
           raise TypeError("arguments need to be an instance of CurveLoop.")
       ReversePrimitive.__init__(self, surface)
       Manifold2D.__init__(self)

    def getBoundaryLoop(self):
       """
       returns the CurveLoop defining the RuledSurface
       """
       return -self.getUnderlyingPrimitive().getBoundaryLoop()

    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        return self.getBoundaryLoop().getCurves()
#==============================
class PlaneSurface(Primitive, Manifold2D):
    """
    a plane surface with holes
    """
    def __init__(self,loop,holes=[]):
       """
       creates a  plane surface with a hole

       @param loop: L{CurveLoop} defining the boundary of the surface
       @param holes: list of L{CurveLoop} defining holes in the surface. 
       @note: A CurveLoop defining a hole should not have any lines in common with the exterior CurveLoop.  
              A CurveLoop defining a hole should not have any lines in common with another CurveLoop defining a hole in the same surface.
       """
       if not isinstance(loop.getUnderlyingPrimitive(),CurveLoop):
           raise TypeError("argument loop needs to be a CurveLoop object.")
       for i in range(len(holes)):
            if not isinstance(holes[i].getUnderlyingPrimitive(), CurveLoop):
                 raise TypeError("%i-th hole needs to be a CurveLoop object.")
       #TODO: check if lines and holes are in a plane
       #TODO: are holes really holes?
       Primitive.__init__(self)
       Manifold2D.__init__(self)
       self.__loop=loop
       self.__holes=holes
    def getHoles(self):
       """
       returns the holes
       """
       return self.__holes

    def getBoundaryLoop(self):
        """
        returns the loop defining the boundary
        """
        return self.__loop

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            sub_dict[self]=PlaneSurface(self.getBoundaryLoop().substitute(sub_dict),[ h.substitute(sub_dict) for h in self.getHoles()])
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True if each curve is colocted with a curve in primitive
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
          if isinstance(primitive.getUnderlyingPrimitive(),PlaneSurface):
             if self.getBoundaryLoop().isColocated(primitive.getBoundaryLoop()):
                hs0=self.getHoles()
                hs1=primitive.getHoles()
                if len(hs0) == len(hs1):
                    for h0 in hs0:
                       colocated = False
                       for h1 in hs1: 
                         colocated = colocated or h0.isColocated(h1)
                       if not colocated: return False
                    return True
       return False
    def collectPrimitiveBases(self):
        """
        returns primitives used to construct the Surface
        """
        out=[self] + self.getBoundaryLoop().collectPrimitiveBases()
        for i in self.getHoles(): out+=i.collectPrimitiveBases()
        return out
    def __neg__(self):
          """
          returns a view onto the curve with reversed ordering
          """
          return ReversePlaneSurface(self)
    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        out = []+ self.getBoundaryLoop().getCurves()
        for h in self.getHoles(): out+=h.getCurves()
        return out

class ReversePlaneSurface(ReversePrimitive, Manifold2D):
    """
    creates a view onto a L{PlaneSurface} but with the reverse orientation
    """
    def __init__(self,surface):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       if not isinstance(surface, PlaneSurface):
           raise TypeError("arguments need to be an instance of PlaneSurface.")
       ReversePrimitive.__init__(self, surface)
       Manifold2D.__init__(self)

    def getBoundaryLoop(self):
       """
       returns the CurveLoop defining the RuledSurface
       """
       return -self.getUnderlyingPrimitive().getBoundaryLoop()

    def getHoles(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        return [ -h for h in self.getUnderlyingPrimitive().getHoles() ]

    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        out = [] + self.getBoundaryLoop().getCurves()
        for h in self.getHoles(): out+=h.getCurves()
        return out


#=========================================================================
class SurfaceLoop(Primitive, PrimitiveBase):
    """
    a loop of 2D primitives. It defines the shell of a volume. 

    The loop must represent a closed shell, and the primitives should be oriented consistently.
    """
    def __init__(self,*surfaces):
       """
       creates a surface loop
       """
       if len(surfaces)<2:
            raise ValueError("at least two surfaces have to be given.")
       for i in range(len(surfaces)):
           if not isinstance(surfaces[i].getUnderlyingPrimitive(),Manifold2D):
              raise TypeError("%s-th argument is not a Manifold2D object."%i)
       self.__surfaces=list(surfaces)
       Primitive.__init__(self)
       PrimitiveBase.__init__(self)
    def __len__(self):
       """
       return the number of curves in the SurfaceLoop
       """
       return len(self.__surfaces)

    def __neg__(self):
       """
       returns a view onto the curve with reversed ordering
       """
       return ReverseSurfaceLoop(self)

    def getSurfaces(self):
       """
       returns the surfaces defining the SurfaceLoop
       """
       return self.__surfaces

    def collectPrimitiveBases(self):
       """
       returns primitives used to construct the SurfaceLoop
       """
       out=[self]
       for c in self.getSurfaces(): out+=c.collectPrimitiveBases()
       return out

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            new_s=[]
            for s in self.getSurfaces(): new_s.append(s.substitute(sub_dict))
            sub_dict[self]=SurfaceLoop(*tuple(new_s))
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True if each surface is colocted with a curve in primitive and vice versa.
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
         if isinstance(primitive.getUnderlyingPrimitive(),SurfaceLoop):
            if len(primitive) == len(self):
                sp0=self.getSurfaces()
                sp1=primitive.getSurfaces()
                for s0 in sp0: 
                    colocated = False
                    for s1 in sp1:
                         colocated = colocated or s0.isColocated(s1)
                    if not colocated: return False
                return True
       return False

class ReverseSurfaceLoop(ReversePrimitive, PrimitiveBase):
    """
    a view to SurfaceLoop with reverse orientaion

    The loop must represent a closed shell, and the primitives should be oriented consistently.
    An oriented loop of 2-dimensional manifolds (= RuledSurface, PlaneSurface)

    The loop must be closed and the one-dimensional manifolds should be oriented consistently.
    """
    def __init__(self,surface_loop):
       """
       creates a polygon from a list of line surfaces. The curves must form a closed loop.
       """
       if not isinstance(surface_loop, SurfaceLoop):
           raise TypeError("arguments need to be an instance of SurfaceLoop.")
       ReversePrimitive.__init__(self, surface_loop)
       PrimitiveBase.__init__(self)

    def getSurfaces(self):
       """
       returns the surfaces defining the SurfaceLoop
       """
       return [ -s for s in  self.getUnderlyingPrimitive().getSurfaces() ]

    def __len__(self):
        return len(self.getUnderlyingPrimitive())

#==============================
class Manifold3D(PrimitiveBase):
    """
    general three-dimensional manifold
    """
    def __init__(self):
       """
       create a three-dimensional manifold
       """
       PrimitiveBase.__init__(self)

    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the volume (including holes)
        """
        raise NotImplementedError()

class Volume(Manifold3D, Primitive):
    """
    a volume with holes.
    """
    def __init__(self,loop,holes=[]):
       """
       creates a volume

       @param loop: L{SurfaceLoop} defining the boundary of the surface
       @param holes: list of L{SurfaceLoop} defining holes in the surface. 
       @note: A SurfaceLoop defining a hole should not have any surfaces in common with the exterior SurfaceLoop.  
              A SurfaceLoop defining a hole should not have any surfaces in common with another SurfaceLoop defining a hole in the same volume.
       """
       if not isinstance(loop.getUnderlyingPrimitive(), SurfaceLoop):
           raise TypeError("argument loop needs to be a SurfaceLoop object.")
       for i in range(len(holes)):
            if not isinstance(holes[i].getUnderlyingPrimitive(), SurfaceLoop):
                 raise TypeError("%i th hole needs to be a SurfaceLoop object.")
       Primitive.__init__(self)
       Manifold3D.__init__(self)
       self.__loop=loop
       self.__holes=holes
    def getHoles(self):
       """
       returns the hole in the volume
       """
       return self.__holes
    def getSurfaceLoop(self):
       """
       returns the loop forming the surface
       """
       return self.__loop

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            sub_dict[self]=Volume(self.getSurfaceLoop().substitute(sub_dict),[ h.substitute(sub_dict) for h in self.getHoles()])
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True if each curve is colocted with a curve in primitive
       """
       if hasattr(primitive,"getUnderlyingPrimitive"): 
          if isinstance(primitive.getUnderlyingPrimitive(),Volume):
             if self.getSurfaceLoop().isColocated(primitive.getSurfaceLoop()):
                hs0=self.getHoles()
                hs1=primitive.getHoles()
                if len(hs0) == len(hs1):
                    for h0 in hs0:
                       colocated = False
                       for h1 in hs1: 
                         colocated = colocated or h0.isColocated(h1)
                       if not colocated: return False
                    return True
       return False
    def collectPrimitiveBases(self):
        """
        returns primitives used to construct the Surface
        """
        out=[self] + self.getSurfaceLoop().collectPrimitiveBases()
        for i in self.getHoles(): out+=i.collectPrimitiveBases()
        return out
    def getBoundary(self):
        """
        returns a list of the one-dimensional manifolds forming the boundary of the Surface (including holes)
        """
        out = []+ self.getSurfaceLoop().getSurfaces()
        for h in self.getHoles(): out+=h.getSurfaces()
        return out

class PropertySet(Primitive, PrimitiveBase):
    """
    defines a group of L{Primitive} which can be accessed through a name
    """
    def __init__(self,name,*items):
       Primitive.__init__(self)
       self.__dim=None
       self.clearItems()
       self.addItem(*items)
       self.setName(name)

    def getDim(self):
       """
       returns the dimension of the items
       """ 
       if self.__dim == None:
           items=self.getItems()
           if len(items)>0:
                if isinstance(items[0] ,Manifold1D): 
                     self.__dim=1
                elif isinstance(items[0] ,Manifold2D): 
                     self.__dim=2
                elif isinstance(items[0] ,Manifold3D): 
                    self.__dim=3
                else:
                    self.__dim=0
       return self.__dim
    def __repr__(self):
       """
       returns a string representation
       """
       return "%s(%s)"%(self.getName(),self.getID())
    def getManifoldClass(self):
        """
        returns the manifold class expected from items
        """
        d=self.getDim()
        if d == None:
           raise ValueError("undefined spatial diemnsion.")
        else:
           if d==0:
              return Point
           elif d==1:
              return Manifold1D
           elif d==2:
              return Manifold2D
           else:
              return Manifold3D

    def getName(self):
        """
        returns the name of the set
        """
        return self.__name
    def setName(self,name):
        """
        sets the name.
        """
        self.__name=str(name)

    def addItems(self,*items):
        """
        adds items. An item my be any L{Primitive} but no L{PropertySet}
        """
        self.addItem(*items)

    def addItem(self,*items): 
        """
        adds items. An item my be any L{Primitive} but no L{PropertySet}
        """
        for i in items: 
            if not i in self.__items: 
               if len(self.__items)>0:
                  m=self.getManifoldClass()
                  if not isinstance(i, m):
                     raise TypeError("argument %s is not a %s class object."%(i, m.__name__))
               self.__items.append(i)
    def getNumItems(self):
        """
        returns the number of items in the property set
        """ 
        return len(self.__items)

    def getItems(self):
        """
        returns the list of items
        """
        return self.__items

    def clearItems(self):
        """
        clears the list of items 
        """
        self.__items=[]
    def collectPrimitiveBases(self):
        """
        returns primitives used to construct the PropertySet
        """
        out=[self] 
        for i in self.getItems(): out+=i.collectPrimitiveBases()
        return out

    def getTag(self):
         """
         returns the tag used for this property set
         """
         return self.getID()
