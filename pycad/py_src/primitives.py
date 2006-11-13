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

global global_primitive_id_counter
global_primitive_id_counter=0

class Primitive(object):
    """
    template for elementary geometrical object
    """
    def __init__(self): 
       """
       """ 
       global global_primitive_id_counter
       self.__ID=global_primitive_id_counter
       global_primitive_id_counter+=1
    def getID(self):
       return self.__ID
    def __repr__(self):
       return "%s(%s)"%(self.__class__.__name__,self.getID())
    def __cmp__(self,other):
       return cmp(self.getID(),other.getID())
    def getGmshCommand(self):
        raise NotImplementedError("getGmshCommand is not implemented for this class %s."%self.__class__.__name__)
    def isPoint(self):
        return False
    def isCurve(self):
        return False
    def isSurface(self):
        return False
    def isCurveLoop(self):
        return False
    def isSurfaceLoop(self):
        return False
    def __neg__(self):
        return ReversedPrimitive(self)

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

class Point(Primitive):
    """
    a three dimensional point
    """
    def __init__(self,x=0.,y=0.,z=0.,local_scale=1.): 
       """
       creates a point with coorinates x,y,z with a relative refinement factor 
       """ 
       super(Point, self).__init__()
       if not local_scale > 0.:
           raise ValueError("local_scale needs to be positive.")
       self._x=numarray.array([x,y,z],numarray.Float64)
       self.setLocalScale(local_scale)
    def setLocalScale(self,factor=1.):
       self.__local_scale=factor
    def isPoint(self):
        return True
    def getLocalScale(self):
       return self.__local_scale
    def getCoordinates(self):
       return self._x
    def __add__(self,other):
       c=self.getCoordinates()+numarray.array(other,numarray.Float64)
       return Point(x=c[0],y=c[1],z=c[2],local_scale=self.getLocalScale())
    def getGmshCommand(self):
        c=self.getCoordinates()
        return "Point(%s) = {%e , %e, %e , %e * scale};"%(self.getID(),c[0],c[1],c[2], self.getLocalScale())
    def getHistory(self):
        return set([self])

class Curve(Primitive):
      """
      a curve
      """
      def __init__(self,*args):
          """
          defines a curve form a set of control points
          """
          super(Curve, self).__init__()
          l=len(args)
          for i in range(l):
              if not args[i].isPoint():
                 raise TypeError("%s-th argument is not a Point object."%i)
          self.__nodes=args
      def __len__(self):
          return len(self.__nodes)
      def isCurve(self):
        return True
      def getHistory(self):
          out=set([self])
          for i in self.getNodes(): out|=i.getHistory()
          return out

      def getStart(self):
         """
         returns start point
         """
         return self.__nodes[0]

      def getEnd(self):
         """
         returns end point
         """
         return self.__nodes[-1]

      def getNodes(self):
         """
         returns a list of the nodes
         """
         return self.__nodes
      def __add__(self,other):
         return Curve([p+other for p in self.getNodes()])
      def getGmshCommand(self):
        out=""
        for i in self.getNodes():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Spline(%s) = {%s};"%(self.getID(),out)

class BezierCurve(Curve):
    """
    a Bezier curve
    """
    def __neg__(self):
         """
         returns the line segment with swapped start and end points
         """
         return BezierCurve(self.getNodes()[::-1])
    def __add__(self,other):
         return BezierCurve([p+other for p in self.getNodes()])
    def getGmshCommand(self):
        out=""
        for i in self.getNodes():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Bezier(%s) = {%s};"%(self.getID(),out)

class BSplineCurve(Curve):
    """
    a BSpline curve. Control points may be repeated.
    """
    def __neg__(self):
         """
         returns the line segment with swapped start and end points
         """
         return BSplineCurve(self.getNodes()[::-1])
    def __add__(self,other):
         return BSplineCurve([p+other for p in self.getNodes()])
    def getGmshCommand(self):
        out=""
        for i in self.getNodes():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "BSpline(%s) = {%s};"%(self.getID(),out)

class Line(Curve):
    """
    a line is defined by two L{Point}s
    """
    def __init__(self,start,end):
        """
        defines a curve form a set of control points
        """
        super(Line, self).__init__(start,end)
    def __neg__(self):
       return ReversedPrimitive(self)
    def __add__(self,other):
       return Line(self.getEnd()+other,self.getStart()+other)
    def getGmshCommand(self):
        return "Line(%s) = {%s, %s};"%(self.getID(),self.getStart().getID(),self.getEnd().getID())

class Arc(Curve):
    """
    defines an arc
    """
    def __init__(self,center,start,end):
       """
       creates an arc by the start point, end point and center
       """
       if center.isPoint():
           raise TypeError("center needs to be a Point object.")
       super(Arc, self).__init__(start,end)
       self.__center=center

    def getCenter(self):
       """
       returns center
       """
       return self.__center
    def __add__(self,other):
       return Arc(self.getCenter()+other,self.getStart()+other,self.getEnd()+other)

    def getHistory(self):
          out=set([self])
          out|=self.getCenter().getHistory()
          for i in self.getNodes(): out|=i.getHistory()
          return out
    def getGmshCommand(self):
        return "Circle(%s) = {%s, %s, %s};"%(self.getID(),self.getStart().getID(),self.getCenter().getID(),self.getEnd().getID())

class CurveLoop(Primitive):
    """
    An oriented loop of curves. 

    The loop must be closed and the L{Curves}s should be oriented consistently.
    """
    def __init__(self,*curves):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       super(CurveLoop, self).__init__()
       self.__curves=[]
       self.addCurve(*curves)
    def addCurve(self,*curves):
       for i in range(len(curves)):
           if not curves[i].isCurve():
              raise TypeError("%s-th argument is not a Curve object."%i)
       self.__curves+=curves

    def isCurveLoop(self):
        return True
    def getCurves(self):
       return self.__curves
    def __add__(self,other):
       return CurveLoop(*tuple([c+other for c in self.getCurves()[::-1]]))
    def __len__(self):
       return len(self.__curves)
    def getHistory(self):
          out=set([self])
          for i in self.getCurves(): out|=i.getHistory()
          return out
    def getGmshCommand(self):
        out=""
        for i in self.getCurves():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Line Loop(%s) = {%s};"%(self.getID(),out)

class Surface(Primitive):
    """
    a surface
    """
    def __init__(self,loop):
       """
       creates a  surface with boundary loop

       @param loop: L{CurveLoop} defining the boundary of the surface
       """
       super(Surface, self).__init__()
       if not loop.isCurveLoop():
           raise TypeError("argument loop needs to be a CurveLoop object.")
       self.__loop=loop
    def isSurface(self):
        return True
    def getBoundaryLoop(self):
       return self.__loop
    def __add__(self,other):
       return Surface(self.getBoundaryLoop()+other)
    def getHistory(self):
        out=set([self]) | self.getBoundaryLoop().getHistory()
        return out
    def getGmshCommand(self):
        return "Ruled Surface(%s) = {%s};"%(self.getID(),self.getBoundaryLoop().getID())

class PlaneSurface(Surface):
    """
    a plane surface with holes
    """
    def __init__(self,loop,holes=[]):
       """
       creates a  plane surface. 

       @param loop: L{CurveLoop} defining the boundary of the surface
       @param holes: list of L{CurveLoop} defining holes in the surface. 
       @note: A CurveLoop defining a hole should not have any lines in common with the exterior CurveLoop.  
              A CurveLoop defining a hole should not have any lines in common with another CurveLoop defining a hole in the same surface.
       """
       super(PlaneSurface, self).__init__(loop)
       for i in range(len(holes)):
            if not holes[i].inCurveLoop():
                 raise TypeError("%i th hole needs to be a CurveLoop object.")
       self.__holes=holes
    def getHoles(self):
       return self.__holes
    def __add__(self,other):
       return PlaneSurface(self.getBoundaryLoop()+other, holes=[h+other for h in self.getHoles()])
    def getHistory(self):
        out=set([self]) | self.getBoundaryLoop().getHistory()
        for i in self.getHoles(): out|=i.getHistory()
        return out
    def getGmshCommand(self):
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

class RuledSurface(Surface):
    """
   A ruled surface, i.e., a surface that can be interpolated using transfinite interpolation
    """
    def __init__(self,loop):
       """
       creates a ruled surface from a 

       @param loop: L{CurveLoop} defining the boundary of the surface. There is a restriction of composed of either three or four L{Curve} objects.
       """
       if not loop.isCurveLoop():
           raise TypeError("argument loop needs to be a CurveLoop object.")
       if len(loop)<3:
           raise TypeError("the loop must contain at least three Curves.")
       super(RuledSurface, self).__init__(loop)
    def __add__(self,other):
       return RuledSurface(self.getBoundaryLoop()+other)
    def getGmshCommand(self):
        return "Ruled Surface(%s) = {%s};"%(self.getID(),self.getBoundaryLoop().getID())

class SurfaceLoop(Primitive):
    """
    a surface loop. It defines the shell of a volume. 

    The loop must represent a closed shell, and the L{Surface}s should be oriented consistently.
    """
    def __init__(self,*surfaces):
       """
       creates a surface loop
       """
       super(SurfaceLoop, self).__init__()
       self.__surfaces=[]
       self.addSurface(*surfaces)
    def addSurface(self,*surfaces):
       for i in range(len(surfaces)):
           if not surfaces[i].isSurface():
              raise TypeError("%s-th argument is not a Surface object."%i)
       self.__surfaces+=surfaces

    def isSurfaceLoop(self):
        return True
    def getSurfaces(self):
       return self.__surfaces
    def __add__(self,other):
       return SurfaceLoop([c+other for c in self.getSurfaces])
    def __len__(self):
       return len(self.__surfaces)
    def getHistory(self):
          out=set([self])
          for i in self.getSurfaces(): out|=i.getHistory()
          return out
    def getGmshCommand(self):
        out=""
        for i in self.getSurfaces():
            if len(out)>0: 
                out+=", %s"%i.getID()
            else:
                out="%s"%i.getID()
        return "Surface Loop(%s) = {%s};"%(self.getID(),out)

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
    def getHistory(self):
        out=set([self]) | self.getSurfaceLoop().getHistory()
        for i in self.getHoles(): out|=i.getHistory()
        return out
    def getGmshCommand(self):
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

class PropertySet(Primitive):
    """
    defines a group L{Primitive} objects. 
    """
    def __init__(self,tag=None,*items):
       super(PropertySet, self).__init__()
       self.__items=items
       self.__tag=tag
    def getHistory(self):
        out=set([self, self.getBoundaryLoop().getHistory()])
        for i in self.getHoles(): out|=i.getHistory()
        return out

class PrimitiveStack(object):
      def __init__(self,*items):
        self.__prims=set()
        for i in items:
            self.__prims|=i.getHistory()
        self.__prims=list(self.__prims)
        self.__prims.sort()

      def getGmshCommands(self):
        out=""
        for i in self.__prims:
           out+=i.getGmshCommand()+"\n"
        return out

      def getPycadCommands(self,name="TMP_"):
        out=""
        for i in self.__prims:
           out+=i.getPycadCommand(name)+"\n"
