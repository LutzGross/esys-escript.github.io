
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
mesh generation using gmsh

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

import design
import tempfile
import os
from primitives import Point, Spline, BezierCurve, BSpline, Line, Arc, CurveLoop, RuledSurface, PlaneSurface, SurfaceLoop, Volume, PropertySet, Ellipse

class Design(design.Design):
    """
    design fo gmsh
    """
    DELAUNAY="iso"
    NETGEN="netgen"
    TETGEN="tetgen"

    __script_file = None

    def __init__(self,dim=3,element_size=1.,order=1,keep_files=False):
       """
       initializes the gmsh design

       @param dim: patial dimension
       @param element_size: global element size
       @param order: element order
       @param keep_files: flag to keep work files.
       """ 
       design.Design.__init__(self,dim=dim,element_size=element_size,order=order,keep_files=keep_files)
       self.setScriptFileName()
       self.setMeshFileName()
       self.setOptions()
    def setScriptFileName(self,name=None):
       """
       set the filename for the gmsh input script. if no name is given a name with extension geo is generated.
       """
       if name == None:
           self.__scriptname=tempfile.mkstemp(suffix=".geo")[1]
       else:
           self.__scriptname=name
           self.setKeepFilesOn()
    def getScriptFileName(self):
       """
       returns the name of the file for the gmsh script
       """
       return self.__scriptname
    def setMeshFileName(self, name=None):
       """
       sets the name for the gmsh mesh file. if no name is given a name with extension msh is generated.
       """
       if name == None:
           self.__mshname=tempfile.mkstemp(suffix=".msh")[1]
       else:
           self.__mshname=name
           self.setKeepFilesOn()
    def getMeshFileName(self):
       """
       returns the name of the file for the gmsh msh
       """
       return self.__mshname
    def setOptions(self,algorithm=None,optimize_quality=True,smoothing=1, curvature_based_element_size=False):
        """
        sets options for the mesh generator
        """ 
        if algorithm==None: algorithm=self.DELAUNAY
        self.__curvature_based_element_size=curvature_based_element_size
        self.__algo=algorithm
        self.__optimize_quality=optimize_quality
        self.__smoothing=smoothing
    def __del__(self):
        """
        clean up
        """
        if not self.keepFiles() and (self.scriptFile() != None) :
            self.scriptFile().close()
            os.unlink(self.getScriptFileName())
            os.unlink(self.getMeshFileName())

    def getCommandString(self):
        """
        returns the gmsh comand
        """
        if self.__optimize_quality:
              opt="-optimize "
        else:
              opt=""
        if self.__curvature_based_element_size:
              clcurv="-clcurv "
        else: 
              clcurv=""
 
        exe="gmsh -%s -algo %s %s-smooth %s %s-v 0 -order %s -o %s %s"%(self.getDim(),
                                                                       self.__algo,
                                                                       clcurv,
                                                                       self.__smoothing,
                                                                       opt,
                                                                       self.getElementOrder(),
                                                                       self.getMeshFileName(),
                                                                       self.getScriptFileName())
        return exe
    def getMeshHandler(self):
        """
        returns a handle to a mesh meshing the design. In the current implementation 
        a mesh file name in gmsh format is returned.
        """
        f = open(self.getScriptFileName(),"w")
        self.setScriptFile(f)
        f.write(self.getScriptString())
        cmd = self.getCommandString()
        ret = os.system(cmd) / 256
	if ret > 0:
	  raise RuntimeError, "Could not build mesh: %s"%cmd
	else:
          return self.getMeshFileName()

    def setScriptFile(self,f=None):
        self.__script_file=f
        return

    def scriptFile(self):
        return self.__script_file



    def getScriptString(self):
        """
        returns the gmsh script to generate the mesh
        """
        h=self.getElementSize()
        out="// generated by esys.pycad\n"
        for prim in self.getAllPrimitives():
           p=prim.getUnderlyingPrimitive()
           if isinstance(p, Point):
               c=p.getCoordinates()
               out+="Point(%s) = {%s , %s, %s , %s };\n"%(p.getID(),c[0],c[1],c[2], p.getLocalScale()*h)
         
           elif isinstance(p, Spline):
               out+="Spline(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getControlPoints()))
    
           elif isinstance(p, BezierCurve):
               out+="Bezier(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getControlPoints()))

           elif isinstance(p, BSpline):
               out+="BSpline(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getControlPoints()))

           elif isinstance(p, Line):
               out+="Line(%s) = {%s, %s};\n"%(p.getID(),p.getStartPoint().getDirectedID(),p.getEndPoint().getDirectedID())

           elif isinstance(p, Arc):
              out+="Circle(%s) = {%s, %s, %s};\n"%(p.getID(),p.getStartPoint().getDirectedID(),p.getCenterPoint().getDirectedID(),p.getEndPoint().getDirectedID())
       
           elif isinstance(p, Ellipse):
              out+="Ellipse(%s) = {%s, %s, %s, %s};\n"%(p.getID(),p.getStartPoint().getDirectedID(),p.getCenterPoint().getDirectedID(),p.getPointOnMainAxis().getDirectedID(), p.getEndPoint().getDirectedID())

           elif isinstance(p, CurveLoop):
               out+="Line Loop(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getCurves()))
       
           elif isinstance(p, RuledSurface):
               out+="Ruled Surface(%s) = {%s};\n"%(p.getID(),p.getBoundaryLoop().getDirectedID())
       
           elif isinstance(p, PlaneSurface):
               line=self.__mkArgs(p.getHoles())
               if len(line)>0:
                 out+="Plane Surface(%s) = {%s, %s};\n"%(p.getID(),p.getBoundaryLoop().getDirectedID(), line)
               else:
                 out+="Plane Surface(%s) = {%s};\n"%(p.getID(),p.getBoundaryLoop().getDirectedID())
       
           elif isinstance(p, SurfaceLoop):
               out+="Surface Loop(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getSurfaces()))
       
           elif isinstance(p, Volume):
               line=self.__mkArgs(p.getHoles())
               if len(line)>0:
                 out+="Volume(%s) = {%s, %s};\n"%(p.getID(),p.getSurfaceLoop().getDirectedID(), line)
               else:
                 out+="Volume(%s) = {%s};\n"%(p.getID(),p.getSurfaceLoop().getDirectedID())

           elif isinstance(p, PropertySet):
               if p.getNumItems()>0: 
	          dim=p.getDim()
                  line="Physical "
                  if dim==0: 
                      line+="Point"
                  elif dim==1: 
                      line+="Line"
                  elif dim==2: 
                      line+="Surface"
                  else:
                      line+="Volume"
                  out+=line+"(" + str(p.getID()) + ") = {"+self.__mkArgs(p.getItems())+"};\n"

           else:
               raise TypeError("unable to pass %s object to gmsh."%str(type(p)))
        return out


    def __mkArgs(self,args):
        line=""
        for i in args:
            if len(line)>0: 
                line+=", %s"%i.getDirectedID()
            else:
                line="%s"%i.getDirectedID()
        return line 
