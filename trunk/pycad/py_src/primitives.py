# $Id:$

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

import numarray
from transformations import _TYPE, Translation, Dilation, Transformation


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
    template for elementary geometrical object
    """
    def __init__(self,id): 
       """
       initializes PrimitiveBase instance object with id
       """ 
       self.__ID=id

    def getID(self):
       """
       returns the primitive ID
       """
       return self.__ID

    def __repr__(self):
       return "%s(%s)"%(self.getClass().__name__,self.getID())

    def __cmp__(self,other):
       """
       compares object with other by comparing the absolute value of the ID
       """
       if isinstance(self, PrimitiveBase):
           return cmp(abs(self.getID()),abs(other.getID()))
       else:
           return False

    def getConstructionPoints(self):
        """
        returns the points used to construct the primitive
        """
        out=set()
        for i in self.getPrimitives(): 
           if isinstance(i,Point): out.add(i)
        return list(out)

    def getPrimitives(self):
        """
        returns a list of primitives used to construct the primitive with no double entries
        """
        out=set()
        for p in self.collectPrimitives():
           if isintance(p,ReversePrimitive):
              out.add(-p)
           elif isintance(p,Primitive):
              out.add(p)
        return list(set(self.collectPrimitives()))

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

    def __neg__(self):
        """
        returns a view of the object with reverse orientiention
        """
        if isinstance(self,ReversePrimitive):
            return self.getUnderlyingPrimitive()
        else:
            return ReversedPrimitive(self)

    def apply(self,transformation):
        """
        returns a new object by applying the transformation
        """
        out=self.copy()
        out.modifyBy(transformation)
        return out

    def collectPrimitives(self):
        """
        returns a list of primitives used to construct the primitive. It may contain primitives twice
        
        @note: this class is overwritten by subclass
        """
        return [self]

    def getGmshCommand(self, local_scaling_factor=1.):
        """
        returns the Gmsh command(s) to create the primitive

        @note: this class is overwritten by subclass
        """
        raise NotImplementedError("getGmshCommand is not implemented for this class %s."%self.__class__.__name__)

    def isColocated(self,primitive):
       """
       returns True is the two primitives are located at the smae position

       @note: this class is overwritten by subclass
       """
       raise NotImplementedError("isColocated is not implemented for this class %s."%self.__class__.__name__)

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.

        @note: this class is overwritten by subclass
        """
        if not sub_dict.has_key(self):
           sub_dict[self]=self.getClass()
        return sub_dict[self]

class Primitive(PrimitiveBase):
    """
    A general primitive
    """
    def __init__(self): 
       """
       instantiate a primitve
       """ 
       global global_primitive_id_counter
       self.__ID=global_primitive_id_counter
       PrimitiveBase.__init__(self,global_primitive_id_counter)

    def getClass(self):
        """
        returns the class of the object
        """
        return self.__class__

class ReversePrimitive(PrimitiveBase):
    """
    A view onto a primitive creating an reverse orientation
    """
    def __init__(self,primitive): 
       """
       instantiate a view onto primitve
       """ 
       PrimitiveBase.__init__(self,-primitive.getID())
       self.__primitive=primitive

    def getUnderlyingPrimitive(self):
        """
        returns the underlying primitive
        """
        return self.__primitive

    def getClass(self):
        """
        returns the class of the underlying object
        """
        return self.getUnderlyingPrimitive().__class__()

class Point(Primitive):
    """
    a three dimensional point
    """
    def __init__(self,x=0.,y=0.,z=0.,local_scale=1.): 
       """
       creates a point with coorinates x,y,z with the local refinement factor local_scale
       """ 
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

    def collectPrimitives(self):
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


    def getGmshCommand(self, local_scaling_factor=1.):
        """
        returns the Gmsh command(s) to create the primitive
        """
        c=self.getCoordinates()
        return "Point(%s) = {%s , %s, %s , %s };"%(self.getID(),c[0],c[1],c[2], self.getLocalScale()*local_scaling_factor)

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
        PrimitiveBass.__init__(self)

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

class CurveBase(Manifold1D):
    def __init__(self):
          """
          create a one-dimensional primitive
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
         return self.getControlPoints[0]

    def getEndPoint(self):
         """
         returns end point
         """
         return self.getControlPoints[-1]

    def collectPrimitives(self):
       """
       returns primitives used to construct the Curve
       """
       out=[self]
       for p in self.getControlPoints(): out+=p.collectPrimitives()
       return out

    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            new_p=[]
            for p in self.getControlPoints(): new_p.append(p.substitute(sub_dict))
            sub_dict[self]=self.getGeneratingClass()(*tuple(new_p))
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True curves are on the same position
       """
       if isinstance(primitive,self.__class__):
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
          else:
             return False
       else:
          return False

      
class Curve(CurveBase):
    """
    a curve defined through a list of control points. 
    """
    def __init__(self,*points):
       """
       defines a curve form control points 
       """
       CurveBase.__init__(self)
       if len(points)<2:
           raise TypeError("Curve needs at least two points")
       super(Curve, self).__init__()
       i=0
       for p in points:
              i+=1
              if not isinstance(p,Point): raise TypeError("%s-th argument is not a Point object."%i)
       self.__points=points

    def __neg__(self):
          """
          returns a view onto the curce with reversed ordering
          """
          return ReverseCurve(self)
    def getControlPoints(self):
         """
         returns a list of the points
         """
         return self.__points

    def getGeneratingClass(self):
         return self.__class__

class ReverseCurve(CurveBase):
    """
    a curve defined through a list of control points. 
    """
    def __init__(self,curve):
       """
       defines a curve form control points 
       """
       CurveBase.__init__(self)
       if not isinstance(curve, CurveBase):
           raise TypeError("curve needs to be an instance of Curve")
       self.__curve=curve

    def __neg__(self):
          """
          returns a view onto the curce with reversed ordering
          """
          return self.__curve
    def getControlPoints(self):
         """
         returns a list of the points
         """
         out=[p for p in self.__curve.getControlPoints()]
         out.reverse()
         return out
    def getGeneratingClass(self):
         return self.__curve.getGeneratingClass()

class Spline(Curve):
    """
    a spline curve defined through a list of control points. 
    """
    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the Curve
        """
        out=""
        for i in self.getControlPoints():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Spline(%s) = {%s};"%(self.getID(),out)
    

class BezierCurve(Curve):
    """
    a Bezier curve
    """
    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the Curve
        """
        out=""
        for i in self.getControlPoints():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Bezier(%s) = {%s};"%(self.getID(),out)

class BSpline(Curve):
    """
    a BSpline curve. Control points may be repeated.
    """
    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the Curve
        """
        out=""
        for i in self.getControlPoints():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "BSpline(%s) = {%s};"%(self.getID(),out)

class Line(Curve):
    """
    a line is defined by two L{Point}s
    """
    def __init__(self,*points):
        """
        defines a line with start and end point
        """
        if len(points)!=2:
           raise TypeError("Line needs two points")
        super(Line, self).__init__(*points)
    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the Curve
        """
        return "Line(%s) = {%s, %s};"%(self.getID(),self.getStartPoint().getID(),self.getEndPoint().getID())


class ArcBase(Manifold1D):
    def collectPrimitives(self):
       """
       returns the primitives used to construct the Curve
       """
       out=[self]
       out+=self.getStartPoint().collectPrimitives()
       out+=self.getEndPoint().collectPrimitives()
       out+=self.getCenterPoint().collectPrimitives()
       return out

    def getGmshCommand(self):
       """
       returns the Gmsh command(s) to create the primitive
       """
       return "Circle(%s) = {%s, %s, %s};"%(self.getID(),self.getStartPoint().getID(),self.getCenterPoint().getID(),self.getEndPoint().getID())

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
       if isinstance(primitive,Arc):
            return (self.getCenterPoint().isColocated(primitive.getCenterPoint())) and ( \
                   (self.getEndPoint().isColocated(primitive.getEndPoint()) and self.getStartPoint().isColocated(primitive.getStartPoint()) ) \
                or (self.getEndPoint().isColocated(primitive.getStartPoint()) and self.getStartPoint().isColocated(primitive.getEndPoint()) ) )
       else:
          return False

class Arc(ArcBase):
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
       # TODO: check length of circle.
       ArcBase.__init__(self)
       self.__center=center
       self.__start=start
       self.__end=end
    def __neg__(self):
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

class ReverseArc(_ArcBase):
    """
    defines an arc which is strictly, smaller than Pi
    """
    def __init__(self,arc):
       """
       creates an arc by the start point, end point and center
       """
       if not isinstance(arc,Arc): raise TypeError("arc needs to be an instance of Arc.")
       ArcBase.__init__(self)
       self.__arc=arc

    def getStartPoint(self):
       """
       returns start point
       """
       return self.__arc.getEndPoint()

    def getEndPoint(self):
       """
       returns end point
       """
       return self.__arc.getStartPoint()

    def getCenterPoint(self):
       """
       returns center
       """
       return self.__arc.getCenterPoint()

    def __neg__(self):
       return self.__arc

#====================================
class CurveLoop(Primitive):
    """
    An oriented loop of 1D primitives (= curves and arcs)

    The loop must be closed and the L{Manifold1D}s should be oriented consistently.
    """
    def __init__(self,*curves):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       super(CurveLoop, self).__init__()
       if len(curves)<2:
            raise TypeError("at least two curves have to be given.")
       for i in range(len(curves)):
           if not isinstance(curves[i],Manifold1D):
              raise TypeError("%s-th argument is not a Manifold1D object."%i)
       # for the curves a loop:
       used=[ False for i in curves]
       self.__curves=[curves[0]]
       used[0]=True
       while not min(used):
          found=False
          for i in xrange(len(curves)):
             if not used[i]:
                if self.__curves[-1].getEndPoint() == curves[i].getStartPoint():
                   self.__curves.append(curves[i])
                   used[i]=True
                   found=True
                   break
          if not found:
             raise ValueError("loop is not closed.")
       if not self.__curves[0].getStartPoint() == self.__curves[-1].getEndPoint():
          raise ValueError("loop is not closed.")

    def getCurves(self):
       """
       returns the curves defining the CurveLoop
       """
       return self.__curves

    def __len__(self):
       """
       return the number of curves in the CurveLoop
       """
       return len(self.getCurves())

    def collectPrimitives(self):
       """
       returns primitives used to construct the CurveLoop
       """
       out=[self]
       for c in self.getCurves(): out+=c.collectPrimitives()
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
       returns True if each curve is collocted with a curve in primitive
       """
       if isinstance(primitive,CurveLoop):
          if len(primitive) == len(self):
             cp0=self.getCurves()
             cp1=primitive.getCurves()
             for c0 in cp0: 
                 collocated = False
                 for c1 in cp1: 
                      collocated = collocated or c0.isColocated(c1)
                 if not collocated: return False
             return True
          else:
             return False
       else:
          return False

    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the primitive
        """
        out=""
        for i in self.getCurves():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Line Loop(%s) = {%s};"%(self.getID(),out)

class Primitive2D(Primitive):
    """
    general two-dimensional primitive
    """
    def __init__(self):
          """
          create a two-dimensional primitive
          """
          super(Primitive2D, self).__init__()

    def getBoundary(self):
        """
        returns a list of the 1D primitives forming the boundary of the Surface (including holes)
        """
        out=[]
        for i in self.getPrimitives():
             if isinstance(i, Manifold1D): out.append(i)
        return out

    def getBoundary(self):
        """
        returns a list of the 1D primitives forming the boundary of the Surface (including holes)
        """
        out=[]
        for i in self.getPrimitives():
             if isinstance(i, Manifold1D): out.append(i)
        return out

    def getBoundaryLoop(self):
        """
        returns the loop defining the outer boundary
        """
        raise NotImplementedError("getBoundaryLoop is not implemented for this class %s."%self.__class__.__name__)

    def getHoles(self):
       """
       returns the holes
       """
       raise NotImplementedError("getHoles is not implemented for this class %s."%self.__class__.__name__)

    def collectPrimitives(self):
        """
        returns primitives used to construct the Surface
        """
        out=[self] + self.getBoundaryLoop().collectPrimitives()
        for i in self.getHoles(): out+=i.collectPrimitives()
        return out

class RuledSurface(Primitive2D):
    """
    A ruled surface, i.e., a surface that can be interpolated using transfinite interpolation
    """
    def __init__(self,loop):
       """
       creates a ruled surface with boundary loop

       @param loop: L{CurveLoop} defining the boundary of the surface. 
       """
       super(RuledSurface, self).__init__()
       if not isinstance(loop,CurveLoop):
           raise TypeError("argument loop needs to be a CurveLoop object.")
       if len(loop)<2:
           raise TypeError("the loop must contain at least two Curves.")
       if len(loop)>4:
           raise TypeError("the loop must contain at least three Curves.")

       self.__loop=loop

    def getBoundaryLoop(self):
        """
        returns the loop defining the outer boundary
        """
        return self.__loop

    def getHoles(self):
        """
        returns the holes
        """
        return []

    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the primitive
        """
        return "Ruled Surface(%s) = {%s};"%(self.getID(),self.getBoundaryLoop().getID())

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
       returns True if each curve is collocted with a curve in primitive
       """
       if isinstance(primitive,RuledSurface):
          return self.getBoundaryLoop().isColocated(primitive.getBoundaryLoop())
       else:
          return False

def createRuledSurface(*curves):
      """
      an easier way to create a L{RuledSurface} from given curves.
      """
      return RuledSurface(CurveLoop(*curves))

class PlaneSurface(Primitive2D):
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
       super(PlaneSurface, self).__init__()
       if not isinstance(loop,CurveLoop):
           raise TypeError("argument loop needs to be a CurveLoop object.")
       for l in loop.getCurves():
           if not isinstance(l,Line):
             raise TypeError("loop may be formed by Lines only.")
       for i in range(len(holes)):
            if not isinstance(holes[i], CurveLoop):
                 raise TypeError("%i-th hole needs to be a CurveLoop object.")
            for l in holes[i].getCurves():
               if not isinstance(l,Line):
                  raise TypeError("holes may be formed by Lines only.")
       #TODO: check if lines and holes are in a plane
       #TODO: are holes really holes?
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

    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the primitive
        """
        out=""
        for i in self.getHoles():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        if len(out)>0:
          return "Plane Surface(%s) = {%s, %s};"%(self.getID(),self.getBoundaryLoop().getID(), out)
        else:
          return "Plane Surface(%s) = {%s};"%(self.getID(),self.getBoundaryLoop().getID())

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
       returns True if each curve is collocted with a curve in primitive
       """
       if isinstance(primitive,PlaneSurface):
          if self.getBoundaryLoop().isColocated(primitive.getBoundaryLoop()):
             hs0=self.getHoles()
             hs1=primitive.getHoles()
             if len(hs0) == len(hs1):
                 for h0 in hs0:
                    collocated = False
                    for h1 in hs1: 
                      collocated = collocated or h0.isColocated(h1)
                    if not collocated: return False
                 return True
             return False
          else:
             return False
       else:
          return False

class SurfaceLoop(Primitive):
    """
    a loop of 2D primitives. It defines the shell of a volume. 

    The loop must represent a closed shell, and the primitives should be oriented consistently.
    """
    def __init__(self,*surfaces):
       """
       creates a surface loop
       """
       super(SurfaceLoop, self).__init__()
       if len(surfaces)<2:
            raise TypeError("at least two surfaces have to be given.")
       for i in range(len(surfaces)):
           if not isinstance(surfaces[i],Primitive2D):
              raise TypeError("%s-th argument is not a Primitive2D object."%i)
       # for the curves a loop:
       used=[ True for s in surfaces]
       self.__surfaces=[surfaces[0]]
       edges=[ e in surfaces[0].getBoundary() ]
       used_edges=[ False in surfaces[0].getBoundary() ]
       while min(used):
          found=False
          for i in xrange(len(surfaces)):
             if not used[i]:
                i_boundary=surfaces[i].getBoundary()
                for ib in xrange(i_boundary):  
                    if i_boundary[ib] in edges:
                         if used_edges[edges.index(i_boundary[ib])]:
                            raise TypeError("boundary segment %s is shared by more than one surface."%str(i_boundary[ib]))
                         used_edges[edges.index(i_boundary[ib])]=True
                         self.__surfaces.append(surfaces[i])
                         for b in i_boundary:
                            if not b in edges:
                                edges.append(b)
                                used_edges.append(False)
                         found=True
                         used[i]=True
                         break
             if found: break
          if not found:
             raise ValueError("loop is not closed.")
       if min(used_edges): 
          raise ValueError("loop is not closed. Surface is missing.")
    def __len__(self):
       """
       return the number of curves in the SurfaceLoop
       """
       return len(self.__surfaces)

    def getSurfaces(self):
       """
       returns the surfaces defining the SurfaceLoop
       """
       return self.__curves

    def collectPrimitives(self):
       """
       returns primitives used to construct the SurfaceLoop
       """
       out=[self]
       for c in self.getSurfaces(): out+=c.collectPrimitives()
       return out
    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the primitive
        """
        out=""
        for i in self.getSurfaces():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Surface Loop(%s) = {%s};"%(self.getID(),out)
    def substitute(self,sub_dict):
        """
        returns a copy of self with substitutes for the primitives used to construct it given by the dictionary C{sub_dict}.
        If a substitute for the object is given by C{sub_dict} the value is returned, otherwise a new instance 
        with substituted arguments is returned.
        """
        if not sub_dict.has_key(self):
            new_s=[]
            for s in self.getCurves(): new_s.append(s.substitute(sub_dict))
            sub_dict[self]=SurfaceLoop(*tuple(new_s))
        return sub_dict[self]

    def isColocated(self,primitive):
       """
       returns True if each surface is collocted with a curve in primitive and vice versa.
       """
       if isinstance(primitive,SurfaceLoop):
         if len(primitive) == len(self):
             sp0=self.getSurfaces()
             sp1=primitive.getCurves()
             for s0 in sp0: 
                 collocated = False
                 for s1 in sp1: 
                      collocated = collocated or s0.isColocated(s1)
                 if not collocated: return False
             return True
         else:
             return False
       else:
         return False

#==========================
class Volume(Primitive):
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
       super(Volume, self).__init__()
       if not loop.isSurfaceLoop():
           raise TypeError("argument loop needs to be a SurfaceLoop object.")
       for i in range(len(holes)):
            if not holes[i].isSurfaceLoop():
                 raise TypeError("%i th hole needs to be a SurfaceLoop object.")
       self.__loop=loop
       self.__holes=holes
    def getHoles(self):
       return self.__holes
    def getSurfaceLoop(self):
       return self.__loop
    def __add__(self,other):
       return Volume(self.getSurfaceLoop()+other, holes=[h+other for h in self.getHoles()])
    def collectPrimitives(self):
        out=[self] + self.getSurfaceLoop().collectPrimitives()
        for i in self.getHoles(): out+=i.collectPrimitives()
        return out
    def getConstructionPoints(self):
        out=self.getSurfaceLoop().getConstructionPoints()
        for i in self.getHoles(): out|=i.Points()
        return out
    def getGmshCommand(self):
        """
        returns the Gmsh command(s) to create the primitive
        """
        out=""
        for i in self.getHoles():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        if len(out)>0:
          return "Volume(%s) = {%s, %s};"%(self.getID(),self.getSurfaceLoop().getID(), out)
        else:
          return "Volume(%s) = {%s};"%(self.getID(),self.getSurfaceLoop().getID())

class ReversedPrimitive(object):
    def __init__(self,prim):
       self.__prim=prim
    def __getattr__(self,name):
       if name == "getID":
          return self.getReverseID
       else:
          return getattr(self.__prim,name)
    def getReverseID(self):
        return -self.__prim.getID()

class PropertySet(Primitive):
    """
    defines a group L{Primitive} objects. 
    """
    def __init__(self,tag=None,*items):
       super(PropertySet, self).__init__()
       self.__items=items
       self.__tag=tag
    def collectPrimitives(self):
        out=[self]+self.getBoundaryLoop().collectPrimitives()
        for i in self.getHoles(): out+=i.collectPrimitives()
        return out

class PrimitiveStack(object):
      def __init__(self,*items):
        self.__prims=set()
        for i in items:
            self.__prims|=i.getPrimitives()
        self.__prims=list(self.__prims)
        self.__prims.sort()

      def getGmshCommands(self):
        out=""
        for i in self.__prims:
           out+=i.getGmshCommand()+"\n"
        return out
