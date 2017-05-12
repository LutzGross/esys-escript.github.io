
##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
mesh generation using Triangle

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Brett Tully, b.tully@uq.edu.au"

import tempfile
import os
import glob
import esys.pycad.design as design
import math 
from esys.pycad.primitives import Point, Spline, BezierCurve, BSpline, Line, Arc, CurveLoop, RuledSurface, PlaneSurface, SurfaceLoop, Volume, PropertySet
from esys.escript import getMPIWorldMax, getMPIRankWorld

class Design(design.AbstractDesign):
    """
    Design for Triangle.
    """
    def __init__(self,dim=2,keep_files=False):
       """
       Initializes the Triangle design.

       :param dim: spatial dimension
       :param keep_files: flag to keep work files
       """
       if dim != 2:
           raise ValueError("only dimension 2 is supported by Triangle.")
       design.AbstractDesign.__init__(self,dim=dim,keep_files=keep_files)
       self.__scriptname=""
       self.setScriptFileName()
       self.setMeshFileName()
       self.setOptions()

    def setScriptFileName(self,name=None):
       """
       Sets the filename for the Triangle input script. If no name is given
       a name with extension `poly` is generated.
       """
       if self.__scriptname:
           os.unlink(self.__scriptname)
       if name == None:
           self.__scriptname_set=False
           self.__scriptname=tempfile.mkstemp(suffix=".poly")[1]
       else:
           self.__scriptname_set=True
           self.__scriptname=name
           self.setMeshFileName(name)

    def getScriptFileName(self):
       """
       Returns the name of the gmsh script file.
       """
       return self.__scriptname

    def setMeshFileName(self,name=None):
       """
       Sets the name of the Triangle mesh file.
       """
       if self.__mshname:
           os.unlink(self.__mshname)
       if name == None:
           self.__mshname=tempfile.mkstemp(suffix="")[1]
       else:
           # Triangle creates a default filename so all we want to pass is
           # the basename
           if (".poly" in name) or (".ele" in name) or (".node" in name):
               s=name.split(".")[:-1]
               name=s[0]
               if len(s) > 1: name+=".%s"*(len(s)-1)%tuple(s[1:])
           self.__mshname=name
           self.setKeepFilesOn()

    def getMeshFileName(self):
       """
       Returns the name of the Triangle mesh file.
       """
       return self.__mshname

    def setOptions(self,cmdLineArgs=""):
        """
        Sets command line options for the mesh generator::

            triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file

            see http://www.cs.cmu.edu/~quake/triangle.switch.html

        :param cmdLineArgs: the switches you would ordinarily use at the
                            command line (e.g. cmdLineArgs="pq25a7.5")
        """
        self.__cmdLineArgs=cmdLineArgs

    def __del__(self):
        """
        Cleans up.
        """
        if not self.keepFiles():
            if not self.__scriptname_set:
                os.unlink(self.getScriptFileName())
            if not self.__mshname_set:
                os.unlink(self.getMeshFileName())

    def getCommandString(self):
        """
        Returns the Triangle command line::

            triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file

            see http://www.cs.cmu.edu/~quake/triangle.switch.html
        """
        if self.__cmdLineArgs == "":
            print("warning: using default command line arguments for Triangle")
        exe="triangle %s %%s"%self.__cmdLineArgs
        return exe

    def getMeshHandler(self):
        """
        Returns a handle to a mesh meshing the design. In the current
        implementation a mesh file name in Triangle format is returned.
        """
        args = self.getCommandString().split()
        args[-1]=args[-1]%self.getScriptFileName()
        if getMPIRankWorld() == 0:
            import subprocess
            open(self.getScriptFileName(),"w").write(self.getScriptString())
            ret = subprocess.call(args) / 256
        else:
            ret=0
        ret=getMPIWorldMax(ret)
        if ret > 0:
          raise RuntimeError("Could not build mesh: %s"%" ".join(args))
        else:
            # <hack> so that users can set the mesh filename they want.
            name=self.getScriptFileName()
            if (".poly" in name) or (".ele" in name) or (".node" in name):
                s=name.split(".")[:-1]
                name=s[0]
                if len(s) > 1: name+=".%s"*(len(s)-1)%tuple(s[1:])
            files=glob.glob("%s.1.*"%name)
            for f in files:
                sufx=f.split(".")[-1]
                mshName=self.getMeshFileName()+"."+sufx
                os.rename(f,mshName)
            # </hack>
            return self.getMeshFileName()

    def getScriptString(self):
        """
        Returns the Triangle script to generate the mesh.
        """
        # comment lines are prefixed by a '#' in triangle *.poly files
        out="# generated by esys.pycad\n"
        vertices="";vertCnt=0
        segments="";segCnt=0
        holestr="";holeCnt=0
        ctrlPts={}
        for prim in self.getItems():

            p=prim.getUnderlyingPrimitive()

            if isinstance(p, PropertySet):
               PSprims=p.getItems()
               for PSp in PSprims:
                   if isinstance(PSp, Point):
                       c=PSp.getCoordinates()
                       vertCnt+=1
                       vertices+="%s %s %s %s\n"%(vertCnt,c[0],c[1],p.getID())

                   elif isinstance(PSp, Line):
                       # get end points and add them as vertices
                       # add line segments - bnd_mkr's are p.getID()
                       # and p.getID() should be mapped to the FID's from the GIS
                       pts=list(PSp.getControlPoints())
                       for pt in pts:
                           c=pt.getCoordinates()
                           if pt not in list(ctrlPts.keys()):
                               vertCnt+=1
                               vertices+="%s %s %s %s\n"%(vertCnt,c[0],c[1],p.getID())
                               ctrlPts[pt]=vertCnt
                           ptIndx=pts.index(pt)
                           if ptIndx != 0:
                               segCnt+=1
                               segments+="%s %s %s %s\n"%(segCnt,
                                                          ctrlPts[pts[ptIndx-1]],
                                                          ctrlPts[pts[ptIndx]],
                                                          p.getID())

                   elif isinstance(PSp, Spline) or isinstance(PSp, BezierCurve) or \
                        isinstance(PSp, BSpline) or isinstance(PSp, Arc):
                       TypeError("Triangle can only handle linear curves: not %s objects."%str(type(p)))

                   elif isinstance(PSp, PlaneSurface):

                       outerBnd=PSp.getBoundaryLoop()
                       holes=PSp.getHoles()

                       for curve in outerBnd.getCurves():
                           pts=list(curve.getControlPoints())
                           for pt in pts:
                               c=pt.getCoordinates()
                               if pt not in list(ctrlPts.keys()):
                                   vertCnt+=1
                                   vertices+="%s %s %s %s\n"%(vertCnt,c[0],c[1],p.getID())
                                   ctrlPts[pt]=vertCnt
                               ptIndx=pts.index(pt)
                               if ptIndx != 0:
                                   segCnt+=1
                                   segments+="%s %s %s %s\n"%(segCnt,
                                                              ctrlPts[pts[ptIndx-1]],
                                                              ctrlPts[pts[ptIndx]],
                                                              p.getID())
# in order to deal with holes in Triangle, you must place a hole node inside
# the boundary of the polygon that describes the hole. For a convex polygon
# (all internal angles < 180 degrees) this is easy, however, for concave
# polygons (one or more internal angles > 180 degrees) this is much
# more difficult. Easiest method is to find the smallest internal angle, and
# place the hole node inside the triangle formed by the three vertices
# associated with that internal angle.
                       for hole in holes:
                           holePts=[]
                           for curve in hole.getCurves():
                               pts=list(curve.getControlPoints())
                               for pt in pts:
                                   c=pt.getCoordinates()
                                   if pt not in list(ctrlPts.keys()):
                                       vertCnt+=1
                                       vertices+="%s %s %s %s\n"%(vertCnt,c[0],c[1],p.getID())
                                       ctrlPts[pt]=vertCnt
                                   ptIndx=pts.index(pt)
                                   if ptIndx != 0:
                                       segCnt+=1
                                       segments+="%s %s %s %s\n"%(segCnt,
                                                                  ctrlPts[pts[ptIndx-1]],
                                                                  ctrlPts[pts[ptIndx]],
                                                                  p.getID())
                                   if pt not in holePts:
                                       holePts.append(pt)
                           vectors=[] # the key corresponds to the ctrlPts index
                           # create vectors
                           for i in range(len(holePts)):
                               A=holePts[i]
                               vectors.append([])
                               if i == 0:
                                   B=holePts[1]
                                   C=holePts[-1]
                               elif i == len(holePts)-1:
                                   B=holePts[0]
                                   C=holePts[-2]
                               else:
                                   B=holePts[i+1]
                                   C=holePts[i-1]
                               vectors[i].append(self.__getVector(A,B))
                               vectors[i].append(self.__getVector(A,C))
                           # get angle between vectors at each vertex
                           for i in range(len(vectors)):
                               angle=self.__getAngle(vectors[i][0],vectors[i][1])
                               vectors[i].append(angle)
                           # find the vertex with the smallest angle
                           minAngle=360.
                           indx=0
                           for i in range(len(vectors)):
                               if vectors[i][2] < minAngle:
                                   indx=i
                                   minAngle=vectors[i][2]
                           # find a node inside the triangle stemming from the
                           # vertex with the smallest internal angle. Do this by
                           # adding 5% of each of the vectorsesys.pycad. to the current point
                           A=holePts[indx]
                           cA=A.getCoordinates()
                           x=cA[0]+(vectors[indx][0][0]+vectors[indx][1][0])/20.
                           y=cA[1]+(vectors[indx][0][1]+vectors[indx][1][1])/20.
                           # use this node to define the hole.
                           holeCnt+=1
                           holestr+="%s %s %s\n"%(holeCnt,x,y)

                   else:
                       raise TypeError("please pass correct PropertySet objects to Triangle design")
            else:
                raise TypeError("please pass only PropertySet objects to Triangle design")
        out+="# vertices #\n"
        out+="%i 2 0 1\n"%(vertCnt)
        out+=vertices
        out+="# segments #\n"
        out+="%i 1\n"%(segCnt)
        out+=segments
        out+="# holes #\n"
        out+="%i\n"%(holeCnt)
        out+=holestr
        return out

    def __getVector(self,A,B):
        # get the vector between two points.
        cA=A.getCoordinates()
        cB=B.getCoordinates()
        x=cB[0]-cA[0]
        y=cB[1]-cA[1]
        return [x,y]

    def __getAngle(self,v,w):
        # get internal (directional) angle between two vectors (in degrees).
        theta=atan2(v[1],v[0]) - atan2(w[1],w[0])
        theta=theta*180./pi
        if theta < 0.:
            theta+=360.
        return theta

