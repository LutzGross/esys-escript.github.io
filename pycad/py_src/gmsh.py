# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
mesh generation using gmsh

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from . import design
import tempfile
import os
from .primitives import Point, Spline, BezierCurve, BSpline, Line, Arc, CurveLoop, RuledSurface, PlaneSurface, SurfaceLoop, Volume, PropertySet, Ellipse
from esys.escript import getMPIWorldMax, getMPIRankWorld, gmshGeo2Msh
from .transformations import DEG

class Design(design.AbstractDesign):
    """
    Design for gmsh.
    """
    DELAUNAY="Delauny"
    MESHADAPT="MeshAdapt"
    FRONTAL="Frontal"
    NETGEN="Frontal"
    TETGEN="Delauny"

    def __init__(self, dim=3, element_size=1., order=1, keep_files=False):
        """
        Initializes the gmsh design.

        :param dim: spatial dimension
        :param element_size: global element size
        :param order: element order
        :param keep_files: flag to keep work files
        """
        design.AbstractDesign.__init__(self,dim=dim,element_size=element_size,order=order,keep_files=keep_files)
        self.__mshname_set = False
        self.__scriptname=""
        self.setScriptFileName()
        self.setOptions()
        self.setFileFormat(self.GMSH)

    def setScriptFileName(self, name=None):
        """
        Sets the filename for the gmsh input script. If no name is given a name
        with extension `geo` is generated.
        """
        if self.__scriptname:
            os.unlink(self.__scriptname)
        if name == None:
            self.__scriptname_set=False
            tmp_f_id=tempfile.mkstemp(suffix=".geo")
            self.__scriptname=tmp_f_id[1]
            os.close(tmp_f_id[0])
        else:
            self.__scriptname=name
            self.__scriptname_set=True

    def getScriptFileName(self):
        """
        Returns the name of the gmsh script file.
        """
        return self.__scriptname

    def setOptions(self, algorithm=None, optimize_quality=True, smoothing=1,
                   curvature_based_element_size=False, algorithm2D=None,
                   algorithm3D=None, generate_hexahedra=False,
                   random_factor=None):
        """
        Sets options for the mesh generator.

        :param algorithm: selects 2D meshing algorithm
        :type algorithm: in self.DELAUNAY, self.MESHADAPT, self.FRONTAL
        :param algorithm2D: must be equal to algorithm
        :type algorithm2D: in self.DELAUNAY, self.MESHADAPT, self.FRONTAL
        :param algorithm3D: selects 3D meshing algorithm
        :type algorithm3D: in self.DELAUNAY, self.FRONTAL
        :param curvature_based_element_size: switch for curvature based element size adaption
        :type curvature_based_element_size: ```bool```
        :param smoothing: number of smoothing steps
        :type smoothing: non-negative ```int```
        :param optimize_quality: switch for mesh quality optimization
        :type optimize_quality: ```bool```
        :param generate_hexahedra: switch for using quadrangles/hexahedra elements everywhere.
        :type generate_hexahedra: ```bool```
        :param random_factor: used in the 2D meshing algorithm (should be increased if RandomFactor * size(triangle)/size(model) approaches machine accuracy)
        :type random_factor: positive ```float```
        """
        if random_factor==None: random_factor=1.e-9
        if not random_factor > 0:
            raise ValueError("random_factor must be positive.")
        smoothing=int(smoothing)
        if not smoothing > 0:
            raise ValueError("smoothing must be positive.")

        if algorithm3D is None:
            algorithm3D=self.FRONTAL
        if algorithm is None:
            if algorithm2D is None:
                algorithm2D=self.MESHADAPT
        else:
            if not algorithm2D is None:
                if not algorithm == algorithm2D:
                    raise ValueError("argument algorithm (=%s) and algorithm2D (=%s) must have the same value if set."%(algorithm, algorithm2D))
            algorithm2D = algorithm
        if not algorithm2D in [ self.DELAUNAY, self.MESHADAPT, self.FRONTAL ]:
            raise ValueError("illegal 2D meshing algorithm %s."%algorithm2D)
        if not algorithm3D in [ self.DELAUNAY, self.FRONTAL ]:
            raise ValueError("illegal 3D meshing algorithm %s."%algorithm3D)

        self.__curvature_based_element_size=curvature_based_element_size
        self.__algo2D=algorithm2D
        self.__algo3D=algorithm3D
        self.__optimize_quality=optimize_quality
        self.__smoothing=smoothing
        self.__generate_hexahedra=generate_hexahedra
        self.__random_factor=random_factor

    def getOptions(self, name=None):
        """
        Returns the current options for the mesh generator.
        """
        if name is None:
            return {"optimize_quality" : self.__optimize_quality ,
                    "smoothing" : self.__smoothing,
                    "curvature_based_element_size" : self.__curvature_based_element_size,
                    "generate_hexahedra" : self.__generate_hexahedra,
                    "algorithm2D" : self.__algo2D,
                    "algorithm3D" : self.__algo3D ,
                    "random_factor" : self.__random_factor }
        else:
            return self.getOption()[name]

    def __del__(self):
        """
        Cleans up.
        """
        try:
            if not self.keepFiles():
                if not self.__scriptname_set: #i.e. it's a tempfile
                    os.unlink(self.getScriptFileName())
                if not self.__mshname_set: #i.e. it's a tempfile
                    os.unlink(self.getMeshFileName())
        except OSError:
            pass # The file might not have been created and there is nothing
                 # to do about a "failure" here anyway

    def getScriptHandler(self):
        """
        Returns a handler to the script file to generate the geometry.
        In the current implementation a script file name is returned.
        """
        if getMPIRankWorld() == 0:
            open(self.getScriptFileName(),"w").write(self.getScriptString())
        return self.getScriptFileName()

    def getMeshHandler(self):
        """
        Returns a handle to a mesh meshing the design. In the current
        implementation a mesh file name in gmsh format is returned.
        """

        verbosity = 3
        ret = gmshGeo2Msh(self.getScriptHandler(), self.getMeshFileName(),
                          self.getDim(), self.getElementOrder(), verbosity)
        if ret > 0:
            self.setKeepFilesOn() #no files to delete, so don't try to
            raise RuntimeError("Could not build mesh using gmsh.\n" + \
                    "Is gmsh available?")
        return self.getMeshFileName()

    def getScriptString(self):
        """
        Returns the gmsh script to generate the mesh.
        """
        h=self.getElementSize()
        out='// generated by esys.pycad\nGeneral.Terminal = 1;\nGeneral.ExpertMode = 1;\n'
        options=self.getOptions()
        if options["optimize_quality"]:
            out += "Mesh.Optimize = 1;\n"
        else:
            out += "Mesh.Optimize = 0;\n"

        if options["curvature_based_element_size"]:
            out += "Mesh.CharacteristicLengthFromCurvature = 1;\n"
        else:
            out += "Mesh.CharacteristicLengthFromCurvature = 0;\n"

        if options["generate_hexahedra"]:
            if self.getDim() == 2:
                out += "Mesh.SubdivisionAlgorithm = 1;\n"
            else:
                out += "Mesh.SubdivisionAlgorithm = 2;\n"
        else:
            out += "Mesh.SubdivisionAlgorithm = 0;\n"

        out += "Mesh.Smoothing = %d;\n"%options["smoothing"]
        out += "Mesh.RandomFactor = %.14e;\n"%options["random_factor"]
        if options["algorithm2D"] == self.MESHADAPT:
            out += "Mesh.Algorithm = 1; // = MeshAdapt\n"
        elif options["algorithm2D"] == self.DELAUNAY:
            out += "Mesh.Algorithm = 5; // = Delaunay\n"
        elif options["algorithm2D"] == self.FRONTAL:
            out += "Mesh.Algorithm = 6; // = Frontal\n"

        if options["algorithm3D"] == self.DELAUNAY:
            out += "Mesh.Algorithm3D = 1; // = Delaunay\n"
        elif options["algorithm3D"] == self.FRONTAL:
            out += "Mesh.Algorithm3D = 4; // = Frontal\n"

        for prim in self.getAllPrimitives():
            p=prim.getUnderlyingPrimitive()
            if isinstance(p, Point):
                c=p.getCoordinates()
                #out+="Point(%s) = {%f , %f, %f , %f };\n"%(p.getID(),c[0],c[1],c[2], p.getLocalScale()*h)
                out += "Point(%s) = {%.14e, %.14e, %.14e, %.14e};\n"%(p.getID(),c[0],c[1],c[2], p.getLocalScale()*h)

            elif isinstance(p, Spline):
                out += "Spline(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getControlPoints()))+self.__mkTransfiniteLine(p)

            elif isinstance(p, BezierCurve):
                out += "Bezier(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getControlPoints()))+self.__mkTransfiniteLine(p)

            elif isinstance(p, BSpline):
                out += "BSpline(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getControlPoints()))+self.__mkTransfiniteLine(p)

            elif isinstance(p, Line):
                out += "Line(%s) = {%s, %s};\n"%(p.getID(),p.getStartPoint().getDirectedID(),p.getEndPoint().getDirectedID())+self.__mkTransfiniteLine(p)

            elif isinstance(p, Arc):
                out += "Circle(%s) = {%s, %s, %s};\n"%(p.getID(),p.getStartPoint().getDirectedID(),p.getCenterPoint().getDirectedID(),p.getEndPoint().getDirectedID())+self.__mkTransfiniteLine(p)

            elif isinstance(p, Ellipse):
                out += "Ellipse(%s) = {%s, %s, %s, %s};\n"%(p.getID(),p.getStartPoint().getDirectedID(),p.getCenterPoint().getDirectedID(),p.getPointOnMainAxis().getDirectedID(), p.getEndPoint().getDirectedID())+self.__mkTransfiniteLine(p)

            elif isinstance(p, CurveLoop):
                out += "Line Loop(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getCurves()))

            elif isinstance(p, RuledSurface):
                out += "Ruled Surface(%s) = {%s};\n"%(p.getID(),p.getBoundaryLoop().getDirectedID())+self.__mkTransfiniteSurface(p)

            elif isinstance(p, PlaneSurface):
                line = self.__mkArgs(p.getHoles())
                if len(line) > 0:
                    out += "Plane Surface(%s) = {%s, %s};\n"%(p.getID(),p.getBoundaryLoop().getDirectedID(), line)+self.__mkTransfiniteSurface(p)
                else:
                    out += "Plane Surface(%s) = {%s};\n"%(p.getID(),p.getBoundaryLoop().getDirectedID())+self.__mkTransfiniteSurface(p)

            elif isinstance(p, SurfaceLoop):
                out += "Surface Loop(%s) = {%s};\n"%(p.getID(),self.__mkArgs(p.getSurfaces()))

            elif isinstance(p, Volume):
                line = self.__mkArgs(p.getHoles())
                if len(line)>0:
                    out += "Volume(%s) = {%s, %s};\n"%(p.getID(),p.getSurfaceLoop().getDirectedID(), line)+self.__mkTransfiniteVolume(p)
                else:
                    out += "Volume(%s) = {%s};\n"%(p.getID(),p.getSurfaceLoop().getDirectedID())+self.__mkTransfiniteVolume(p)

            elif isinstance(p, PropertySet):
                if p.getNumItems() > 0:
                    dim=p.getDim()
                    line = "Physical "
                    if dim==0:
                        line += "Point"
                    elif dim==1:
                        line += "Line"
                    elif dim==2:
                        line += "Surface"
                    else:
                        line += "Volume"
                    out += line+"(" + str(p.getID()) + ") = {"+self.__mkArgs(p.getItems(),useAbs=True)+"};\n"

            else:
                raise TypeError("unable to pass %s object to gmsh."%str(type(p)))
        return out

    def __mkArgs(self, args, useAbs=False):
        line = ""
        for i in args:
            id = i.getDirectedID()
            if useAbs: id=abs(id)
            if len(line) > 0:
                line += ", %s"%id
            else:
                line = "%s"%id
        return line

    def __mkTransfiniteLine(self, p):
        s = p.getElementDistribution()
        if not s == None:
            if s[2]:
                out="Transfinite Line{%d} = %d Using Bump %s;\n"%(p.getID(),s[0],s[1])
            else:
                out="Transfinite Line{%d} = %d Using Progression %s;\n"%(p.getID(),s[0],s[1])
        else:
            out=""
        return out

    def __mkTransfiniteSurface(self, p):
        out = ""
        o = p.getRecombination()
        s = p.getTransfiniteMeshing()
        if not s == None:
            out2 = ""
            if not s[0] is None:
                for q in s[0]:
                    if len(out2)==0:
                        out2 = "%s"%q.getID()
                    else:
                        out2 = "%s,%s"%(out2, q.getID())
            if s[1] is None:
                out += "Transfinite Surface{%s} = {%s};\n"%(p.getID(),out2)
            else:
                out += "Transfinite Surface{%s} = {%s} %s;\n"%(p.getID(),out2,s[1])
        if not o is None:
            out += "Recombine Surface {%s} = %f;\n"%(p.getID(), o/DEG)
        return out

    def __mkTransfiniteVolume(self, p):
        out=""
        s=p.getTransfiniteMeshing()
        if not s == None:
            if len(s)>0:
                out2=""
                for q in s[0]:
                    if len(out2)==0:
                        out2="%s"%q.getID()
                    else:
                        out2="%s,%s"%(out2,q.getID())
                out+="Transfinite Volume{%s} = {%s};\n"%(p.getID(),out2)
            else:
                out+="Transfinite Volume{%s};\n"%(p.getID(),)
        return out

