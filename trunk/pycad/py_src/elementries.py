# $Id:$

"""
Geometrical Elementries

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

class Elementary(object):
    """
    template for elementary geometrical object
    """
    def __init__(self): 
       """
       """ 
       pass

class Point(Elementary):
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
       self.setLocalScale(0.,local_scale)
    def __str__(self):
       return str(self._x)
    def setLocalScale(self,factor=1.):
       self.__local_scale=factor
 

class Curve(Elementary):
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
              if not isinstance(args[i],Point):
                 raise TypeError("%s-th argument is not a Point object.")
          super(Curve, self).__init__(args[0],args[l-1])
          self.__nodes=list(args)
      def __len__(self):
          return len(self.__nodes)

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

      def __neg__(self):
         """
         returns the line segment with swapped start and end points
         """
         return Curve(self.__nodes[::-1])

class BezierCurve(Curve):
    """
    a Bezier curve
    """
    def __neg__(self):
       """
       returns the line segment with swapped start and end points
       """
       return BezierCurve(self.getNodesInReversedOrder())

class BSplineCurve(Curve):
    """
    a BSpline curve. Control points may be repeated.
    """
    def __neg__(self):
       """
       returns the line segment with swapped start and end points
       """
       return BSplineCurve(self.getNodesInReversedOrder())


class Line(Curve):
    """
    a line is defined by two L{Point}s
    """
    def __neg__(self):
       """
       returns the line segment with swapped start and end points
       """
       return Line(self.getEnd(),self.getStart())

class Arc(Curve):
    """
    defines an arc
    """
    def __init__(self,center,start,end):
       """
       creates an arc by the start point, end point and center
       """
       if not isinstance(center,Point):
           raise TypeError("center needs to be a Point object.")
       super(Arc, self).__init__(start,end)
       self.__center=center

    def getCenter(self):
       """
       returns center
       """
       return self.__center
    def __neg__(self):
       """
       returns the line segment with swapped start and end points
       """
       return Arc(self.getCenter(),self.getEnd(),self.getStart())

class CurveLoop(Elementary):
    """
    An oriented loop of curves. 

    The loop must be closed and the L{Curves}s should be oriented consistently.
    """
    def __init__(self,*curves):
       """
       creates a polygon from a list of line curves. The curves must form a closed loop.
       """
       super(CurveLoop, self).__init__()
       for i in curves:
           if not isinstance(curves[i],Curve):
              raise TypeError("%s-th argument is not a Curve object.")
       if length(curves)>0:
          for i in range(len(curves)-1):
             if not curves[i].getEnd() == curves[i+1].getStart():
                 raise ValueError("start point of %s-th curve does not match end point of %s-th curve."%(i,i+1))
          if curves[0].getStart() == curves[len(curves)-1].getEnd():
                 raise ValueError("loop is not closed.")
       self.__curves=curves

    def getCurves(self):
       return self.__curves

    def __neg__(self):
       return CurveLoop([-c for c in self.__curves[::-1]])

    def __len__(self):
       return len(self.__curves)

class Surface(Elementary):
    """
    a surface
    """
    def __init__(self,loop):
       """
       creates a  surface with boundary loop

       @param loop: L{CurveLoop} defining the boundary of the surface
       """
       super(Surface, self).__init__()
       if not isinstance(loop,CurveLoop):
           raise TypeError("argument loop needs to be a CurveLoop object.")
       self.__loop=loop
    def getBoundaryLoop(self):
       return self.__loop
    def __neg__(self):
       return Surface(-self.getBoundaryLoop())
       
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
            if not isinstance(holes[i],CurveLoop):
                 raise TypeError("%i th hole needs to be a CurveLoop object.")
       self.__holes=holes
    def getHoles(self):
       return self.__holes
    def __neg__(self):
       return PlaneSurface(-self.getBoundaryLoop(),holes=[-h for h in self.getHoles()] )

class RuledSurface(Surface):
    """
   A ruled surface, i.e., a surface that can be interpolated using transfinite interpolation
    """
    def __init__(self,loop):
       """
       creates a ruled surface from a 

       @param loop: L{CurveLoop} defining the boundary of the surface. There is a restriction of composed of either three or four L{Curve} objects.
       """
       if not isinstance(loop,CurveLoop):
           raise TypeError("argument loop needs to be a CurveLoop object.")
       if len(loop)<3:
           raise TypeError("the loop must contain at least three Curves.")
       super(RuledSurface, self).__init__(loop)
    def __neg__(self):
       return RuledSurface(-self.getBoundaryLoop())

class SurfaceLoop(Elementary):
    """
    a surface loop. It defines the shell of a volume. 

    The loop must represent a closed shell, and the L{Surface}s should be oriented consistently.
    """
    def __init__(self,*surfaces):
       """
       creates a surface loop
       """
       super(SurfaceLoop, self).__init__()
       for i in surfaces:
           if not isinstance(surfaces[i],Curve):
              raise TypeError("%s-th argument is not a Curve object.")
       self.__surfaces=surfaces

    def getSurfaces(self):
       return self.__surfaces

    def __neg__(self):
       return SurfaceLoop(-c for c in self.__surfaces[::-1])

    def __len__(self):
       return len(self.__surfaces)

class Volume(Elementary):
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
       if not isinstance(loop,SurfaceLoop):
           raise TypeError("argument loop needs to be a SurfaceLoop object.")
       for i in range(len(holes)):
            if not isinstance(holes[i],SurfaceLoop):
                 raise TypeError("%i th hole needs to be a SurfaceLoop object.")
       self.__loop=loop
       self.__holes=holes
    def getHoles(self):
       return self.__holes
    def getSurfaceLoop(self):
       return self.__loop
    def __neg__(self):
       return PlaneSurface(-self.getSurfaceLoop(),holes=[-h for h in self.getHoles()] )

class PropertySet(Elementary):
    """
    defines a group L{Elementary} objects. 
    """
    def __init__(self,tag=None,*items):
       super(PropertySet, self).__init__()
       self.__items=items
       self.__tag=tag
