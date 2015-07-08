from __future__ import print_function
##############################################################################
#
# Copyright (c) 2003-2015 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from esys.escript import *
from math import pi
import tempfile, os
import logging
logger=logging.getLogger('inv.DCResDomGenerator')

HAS_FINLEY = True
try:
    from esys.finley import ReadGmsh, ReadMesh
except ImportError as e:
    HAS_FINLEY = False

HAVE_GMSH = getEscriptParamInt("GMSH_SUPPORT")

class DCResDomGenerator(object):
    
    """
    This class is used to generate escript domains which are suitable for
    solving dc resistivity problems. The mesh will be refined close to the 
    points provided in electrodeLst. The tags will be included in the domain
    as diracTags so that set tagged value can be used. If prism is set to None
    the mesh will be a homogeneous halfspace. interfaces can be set when calling
    getDom to get a layerd earth model.

    note that the coordinate system is such that -z is the direction of the 
    subsurface and the domain run from -x extents/2,-y extents/2 to 
    x extents/2,y extents/2

    usage: setup dom generator object then call getDom to get the domain.
    """

    def __init__(self, extents, electrodeLst, cl=0.1, tmpDir=None, prism=None, bufferThickness=None, interfaces=None):
        """
        :param extents: x,y,z extents of the domain
        :type extents: list or tuple, len should = 3
        :param electrodeLst: A list of tuples of the form (tag,coords) for each electrode
        :type electrodeLst: list of tuples
        :param cl: Characteristic Length, Contrains the size of elements
        :type cl:float
        :param tmpDir: directory in which to save temporary files:
        :type tmpDir: string
        :param prism: provide start point,extents and a extrude depth for a cubic prism
            note that the prism will be labeled prism for the purpose of set tagged value
        :type prims: [(x,y,z)_start,(x,y,z)_extent]
        :param bufferThickness: a buffer can be added to reduce
        boundary effects this parameter defines this thickness:
        :type bufferThickness:
        :param interfaces: Specify a list of interfaces for a layered model.
                           Doing this will ignore the z-extent. The layers
                           will be tagged iteratively from volume-0 to
                           volume-(n-1).
        :type interfaces: list
        """
        if not HAS_FINLEY:
            raise RuntimeError("Finley module not available")
        if(len(extents)==3 or len(extents)==4):
            self.__extents=extents
        else:
            raise ValueError("extents should be of length 3 or 4")
        for electrode in electrodeLst:
            if len(electrode[1]) != 4:
                raise ValueError("currently only 3d domains are supported electrodeLst elements must be of length 4)")
        self.__extentLen = len(self.__extents)
        self.__electrodeLst=electrodeLst
        self.__lc=cl
        self.__scriptString=""
        self.__scriptArray=[]
        self.__pntList=""
        self.__prism=prism
        self.__tags=[]
        self.__points=[]
        self.__tmpDir=tmpDir
        self.__bufferThickness=bufferThickness
        self.__pointCount=5 # the point from which to start added embeded points
        self.nodeFieldSize=0
        self.nodeMeshSize=0
        self.filename=""
        self.interfaces=interfaces
        for i in electrodeLst:
            self.__tags.append(i[0])
            self.__points.append(i[1][:-1])

    def generateScriptFile(self, nodeFieldSize):
        """
        create a tempfile in which to store the gmsh code. 
        generate the scriptstring and then write it to the file.
        """
        fd, filename = tempfile.mkstemp(suffix=".geo", dir=self.__tmpDir)
        os.close(fd)
        
        if self.interfaces is None:
            self.generateScriptString()

        else:
            raise ValueError("layered script not implemented")
            # self.generateLayedScriptString(self.interfaces, nodeFieldSize)
        
        open(filename, "w").write(self.__scriptString)
        return filename

    def genPrism(self):
        pc = self.__pointCount
        pris=self.__prism
        self.__scriptArray.append("Point(%d)={%f, %f, %f, lc};\n"%(pc ,pris[0][0],pris[0][1],pris[0][2]))
        self.__scriptArray.append("Point(%d)={%f, %f, %f, lc};\n"%(pc + 1,(pris[0][0] + pris[1][0]), pris[0][1],pris[0][2]))
        self.__scriptArray.append("Point(%d)={%f, %f, %f, lc};\n"%(pc + 2,(pris[0][0] + pris[1][0]), (pris[0][1] + pris[1][1]),pris[0][2]))
        self.__scriptArray.append("Point(%d)={%f, %f, %f, lc};\n"%(pc + 3,pris[0][0], (pris[0][1] + pris[1][1]),pris[0][2]))
        self.__scriptArray.append("l1=newreg; Line(l1) = {%d,%d};\n"%(pc,pc + 1))
        self.__scriptArray.append("l2=newreg; Line(l2) = {%d,%d};\n"%(pc + 2,pc + 1))
        self.__scriptArray.append("l3=newreg; Line(l3) = {%d,%d};\n"%(pc + 2,pc + 3))
        self.__scriptArray.append("l4=newreg; Line(l4) = {%d,%d};\n"%(pc + 3,pc))
        self.__scriptArray.append("ll2=newreg; Line Loop(ll2) = {l4,l1,-l2,l3} ;\n")
        self.__scriptArray.append("ps2=newreg; Plane Surface(ps2) = {ll2}; \n")
        self.__scriptArray.append("prism[]=Extrude {0, 0, %f} { Surface {ps2};};\n"%(pris[1][2]))
        self.__scriptArray.append("Physical Volume(\"prism\") = {%d} ;\n"%(2))

    def genScriptWithBuffer(self):
        """
        generate the scipt with buffer added
        """
        self.__scriptArray.append("lc=%f;\n"%self.__lc)
        self.__scriptArray.append("Point(1)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.-self.__bufferThickness) , (-self.__extents[1]/2.-self.__bufferThickness)))
        self.__scriptArray.append("Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.+self.__bufferThickness)  , (-self.__extents[1]/2.-self.__bufferThickness)))
        self.__scriptArray.append("Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.+self.__bufferThickness)  , (self.__extents[1]/2.+self.__bufferThickness)))
        self.__scriptArray.append("Point(4)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.-self.__bufferThickness) , (self.__extents[1]/2.+self.__bufferThickness)))
        self.__scriptArray.append("Line(1) = {1,2} ;\n")
        self.__scriptArray.append("Line(2) = {3,2} ;\n")
        self.__scriptArray.append("Line(3) = {3,4} ;\n")
        self.__scriptArray.append("Line(4) = {4,1} ;\n")
        self.__scriptArray.append("Line Loop(5) = {4,1,-2,3} ; \n")
        self.__scriptArray.append("Plane Surface(6) = {5} ; \n")
        pctmp=self.__pointCount
        for i in self.__electrodeLst:
            pntInfo=i[1]
            self.__scriptArray.append("Point(%d)={%f,%f,%f,%f};\n"%(self.__pointCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
            self.__pointCount+=1
        self.__scriptArray.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%(self.__extents[2]+self.__bufferThickness))
        self.__pntList=str([i for i in range(pctmp,self.__pointCount)])[1:-1]
        self.__scriptArray.append("Point{%s} In Surface{6};\n"%self.__pntList)
        self.__scriptArray.append("Physical Volume(\"hs-%d\") = {%d} ;\n"%(1,1))
        
        

    def genScriptWithoutBuffer(self):
        """
        generate the scipt with no buffer added
        """

        self.__scriptArray.append("lc=%f;\n"%self.__lc)
        self.__scriptArray.append("Point(1)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.) , (-self.__extents[1]/2.)))
        self.__scriptArray.append("Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.)  , (-self.__extents[1]/2.)))
        self.__scriptArray.append("Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.)  , (self.__extents[1]/2.)))
        self.__scriptArray.append("Point(4)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.) , (self.__extents[1]/2.)))
        self.__scriptArray.append("Line(1) = {1,2} ;\n")
        self.__scriptArray.append("Line(2) = {3,2} ;\n")
        self.__scriptArray.append("Line(3) = {3,4} ;\n")
        self.__scriptArray.append("Line(4) = {4,1} ;\n")
        self.__scriptArray.append("Line Loop(5) = {4,1,-2,3} ; \n")
        self.__scriptArray.append("Plane Surface(6) = {5} ; \n")

        pctmp=self.__pointCount
        for i in self.__electrodeLst:
            pntInfo=i[1]
            pc = self.__pointCount
            self.__scriptArray.append("p%d=newp; Point(p%d)={%f,%f,%f,%f};\n"%(pc,pc,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
            self.__pointCount+=1
        self.__scriptArray.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%(self.__extents[2]))
        
        self.__pntList=[i for i in range(pctmp,self.__pointCount)]
        tmplst = []
        for i in self.__pntList:
            tmplst.append("p%s"%i)
        self.__scriptArray.append("Point{%s} In Surface{6};\n"%','.join(tmplst))
        self.__scriptArray.append("Physical Volume(\"hs-%d\") = {%d} ;\n"%(1,1))

    def generateScriptString(self):
        #the str's that follow are used for specifying the boundary faces of the domain
        leftStr ="-out0[2],"
        rightStr="-out0[4],"
        frontStr="-out0[5],"
        backStr ="-out0[3],"
        out=[]

        if not self.__bufferThickness is None:
            self.genScriptWithBuffer()
        else:
            self.genScriptWithoutBuffer()

        if self.__prism != None:
            self.genPrism()

        self.__scriptArray.append("Physical Surface(\"Top\") = { -6 };\n")
        self.__scriptArray.append("Physical Surface(\"Bottom\") = { -out%d[0] };\n"%0)
        self.__scriptArray.append("Physical Surface(\"Left\") = { %s };\n"%leftStr[:-1])
        self.__scriptArray.append("Physical Surface(\"Right\") = { %s };\n"%rightStr[:-1])
        self.__scriptArray.append("Physical Surface(\"Front\") = { %s };\n"%frontStr[:-1])
        self.__scriptArray.append("Physical Surface(\"Back\") = { %s };\n"%backStr[:-1])
        
        if not self.__bufferThickness is None:
            self.__scriptArray.append("Field[1] = Box;\n")
            self.__scriptArray.append("Field[1].VIn=lc;\n")    
            self.__scriptArray.append("Field[1].VOut=5*lc;\n")
            # define field one bounding box
            self.__scriptArray.append("Field[1].XMax=%f;\n"%(self.__extents[0]/2.))
            self.__scriptArray.append("Field[1].XMin=%f;\n"%(-self.__extents[0]/2.))
            self.__scriptArray.append("Field[1].YMax=%f;\n"%(self.__extents[1]/2.))
            self.__scriptArray.append("Field[1].YMin=%f;\n"%(-self.__extents[1]/2.))
            self.__scriptArray.append("Field[1].ZMax=0;\n")
            self.__scriptArray.append("Field[1].ZMin=-%f;\n"%self.__extents[2])

            self.__scriptArray.append("Field[2] = Attractor;\n")
            self.__scriptArray.append("Field[2].NodesList = {%s};\n"%self.__pntList)
            self.__scriptArray.append("Field[3] = Threshold;\n")
            self.__scriptArray.append("Field[3].IField = 2;\n")
            if self.nodeMeshSize is None:
                self.__scriptArray.append("Field[3].LcMin = lc / 5;\n")
            else:
                self.__scriptArray.append("Field[3].LcMin = %g;\n"%self.nodeMeshSize)
            self.__scriptArray.append("Field[3].LcMax = 100*lc;\n") # this value is so high because It should not play a role in field 4
            self.__scriptArray.append("Field[3].DistMin = %g;\n"%self.nodeFieldSize[0])
            self.__scriptArray.append("Field[3].DistMax = %g;\n"%self.nodeFieldSize[1])
            self.__scriptArray.append("Field[4] = Min;\n")
            self.__scriptArray.append("Field[4].FieldsList = {1, 3};\n")
            self.__scriptArray.append("Background Field = 4;\n")
        self.__scriptArray.append("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
        self.__scriptString = "".join(self.__scriptArray)

    # def generateLayedScriptString(self, interfaces):
    #     self.__pointCount=5
    #     leftStr ="-out0[2],"
    #     rightStr="-out0[4],"
    #     frontStr="-out0[5],"
    #     backStr ="-out0[3],"
    #     out = []
    #     extentCount=float(interfaces[0])

    #     if not self.__bufferThickness == None:
    #         out.append("lc=%f;\n"%self.__lc)
    #         out.append("Point(1)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.-self.__bufferThickness) , (-self.__extents[1]/2.-self.__bufferThickness)))
    #         out.append("Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.+self.__bufferThickness)  , (-self.__extents[1]/2.-self.__bufferThickness)))
    #         out.append("Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.+self.__bufferThickness)  , (self.__extents[1]/2.+self.__bufferThickness)))
    #         out.append("Point(4)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.-self.__bufferThickness) , (self.__extents[1]/2.+self.__bufferThickness)))
    #         out.append("Line(1) = {1,2} ;\n")
    #         out.append("Line(2) = {3,2} ;\n")
    #         out.append("Line(3) = {3,4} ;\n")
    #         out.append("Line(4) = {4,1} ;\n")
    #         out.append("Line Loop(5) = {4,1,-2,3} ; \n")
    #         out.append("Plane Surface(6) = {5} ; \n")
    #         for i in self.__electrodeLst:
    #             pntInfo=i[1]
    #             out.append("Point(%d)={%f,%f,%f,%f};\n"%(self.__pointCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
    #             self.__pointCount+=1
    #         #out.append("out[]=Extrude {0, 0, -%f} { Surface {6};};\n"%self.__bufferThickness)
    #         #out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(0,1))
    #         #leftStr+= "-out[2],"
    #         #rightStr+="-out[4],"
    #         #frontStr+="-out[5],"
    #         #backStr+= "-out[3],"
    #         out.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%interfaces[0])
    #         out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1))
    #         #out.append("Point{%s} In Surface{out[0]};\n"%str(range(5,self.__pointCount))[1:-1])
    #         self.__pntList=str([i for i in range(5,self.__pointCount)])[1:-1]
    #         out.append("Point{%s} In Surface{6};\n"%self.__pntList)
    #         for i in range(1,len(interfaces)):
    #             extentCount+=float(interfaces[i])
    #             out.append("out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,interfaces[i],i-1))
    #             out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(i+1,i+1))
    #             leftStr+= "-out%d[2],"%i
    #             rightStr+="-out%d[4],"%i
    #             frontStr+="-out%d[5],"%i
    #             backStr+= "-out%d[3],"%i
    #         i+=1
    #         out.append("out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,self.__bufferThickness,i-1))
    #         out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(i+1,i+1))
    #         leftStr+= "-out%d[2],"%i
    #         rightStr+="-out%d[4],"%i
    #         frontStr+="-out%d[5],"%i
    #         backStr+= "-out%d[3],"%i
    #     else:
    #         out.append("lc=%f;\n"%self.__lc)
    #         out.append("Point(1)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.) , (-self.__extents[1]/2.)))
    #         out.append("Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.)  , (-self.__extents[1]/2.)))
    #         out.append("Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]/2.)  , (self.__extents[1]/2.)))
    #         out.append("Point(4)={%f, %f, 0, lc};\n"%( (-self.__extents[0]/2.) , (self.__extents[1]/2.)))
    #         out.append("Line(1) = {1,2} ;\n")
    #         out.append("Line(2) = {3,2} ;\n")
    #         out.append("Line(3) = {3,4} ;\n")
    #         out.append("Line(4) = {4,1} ;\n")
    #         out.append("Line Loop(5) = {4,1,-2,3} ; \n")
    #         out.append("Plane Surface(6) = {5} ; \n")
    #         for i in self.__electrodeLst:
    #             pntInfo=i[1]
    #             out.append("Point(%d)={%f,%f,%f,%f};\n"%(self.__pointCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
    #             self.__pointCount+=1
    #         self.__pntList=str([i for i in range(5,self.__pointCount)])[1:-1]
    #         out.append("Point{%s} In Surface{6};\n"%self.__pntList)
    #         out.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%interfaces[0])
    #         out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(0,1))

    #         for i in range(1,len(interfaces)):
    #             out.append("out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,interfaces[i],i-1))
    #             out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(i,i+1))
    #             leftStr+= "-out%d[2],"%i
    #             rightStr+="-out%d[4],"%i
    #             frontStr+="-out%d[5],"%i
    #             backStr+= "-out%d[3],"%i


    #     # out.append("Physical Surface(\"Top\") = { -6 };\n")
    #     # out.append("Physical Surface(\"Bottom\") = { -out%d[0] };\n"%i)
    #     # out.append("Physical Surface(\"Left\") = { %s };\n"%leftStr[:-1])
    #     # out.append("Physical Surface(\"Right\") = { %s };\n"%rightStr[:-1])
    #     # out.append("Physical Surface(\"Front\") = { %s };\n"%frontStr[:-1])
    #     # out.append("Physical Surface(\"Back\") = { %s };\n"%backStr[:-1])
    #     if not self.__bufferThickness is None:
    #         out.append("Field[1] = Box;\n")
    #         out.append("Field[1].VIn=lc;\n")
    #         out.append("Field[1].VOut=5*lc;\n")
    #         out.append("Field[1].XMax=%f;\n"%(self.__extents[0]/2.))
    #         out.append("Field[1].XMin=%f;\n"%(-self.__extents[0]/2.))
    #         out.append("Field[1].YMax=%f;\n"%(self.__extents[1]/2.))
    #         out.append("Field[1].YMin=%f;\n"%(-self.__extents[1]/2.))
    #         out.append("Field[1].ZMax=0;\n")
    #         out.append("Field[1].ZMin=-%f;\n"%extentCount)
    #         out.append("Field[2] = Attractor;\n")
    #         out.append("Field[2].NodesList = {%s};\n"%self.__pntList)
    #         out.append("Field[3] = Threshold;\n")
    #         out.append("Field[3].IField = 2;\n")
    #         if self.nodeMeshSize is None:
    #             out.append("Field[3].LcMin = lc / 5;\n")
    #         else:
    #             out.append("Field[3].LcMin = %g;\n"%self.nodeMeshSize)
    #         out.append("Field[3].LcMax = 100*lc;\n")
    #         out.append("Field[3].DistMin = %g;\n"%self.nodeFieldSize[0])
    #         out.append("Field[3].DistMax = %g;\n"%self.nodeFieldSize[1])

    #         out.append("Field[4] = Min;\n")
    #         out.append("Field[4].FieldsList = {1, 3};\n")
    #         out.append("Background Field = 4;\n")
    #         out.append("Mesh.CharacteristicLengthExtendFromBoundary = 0;\n")

    #     self.__scriptString = "".join(out)

    def getDom(self, nodeFieldSize , nodeMeshSize = None, mshName = None , reUse = False):
        """
        Generates and returns a valid escript domain.

        :param nodeFieldSize: size of field around the electrodes for which the
        mesh size is reduced to nodeMeshSize.
        :type nodeFieldSize: [float,float]
        :param nodeMeshSize: size to which the mesh is reduced around electrodes if specified 
        otherwise size will be cl/5 where cl is the Characteristiclength
        :type nodeMeshSize: float
        :param mshName: Specify a name for output mesh this is usefull if you would like
         to reuse the mesh rather than generating a new mesh every time.

        :type mshName: string

        :param reUse: should the msh be reused or should a new file be generated
        :type reUse: bool
        """
        self.nodeFieldSize = nodeFieldSize
        self.nodeMeshSize = nodeMeshSize

        filename = "" # won't be used by non-0 ranks
        self.filename=filename
        if (mshName is not None and os.path.isfile(mshName) and reUse==True):
            if mshName[-4:]=='.msh':
                dom=ReadGmsh(mshName, 3, diracTags=self.__tags, diracPoints=self.__points)
            else:
                raise ValueError("mshName must end in .msh")
            return dom
        elif mshName is None and reUse == True:
            raise ValueError("Can't reuse msh if mshname is not specified")

        # early exit so we don't even create files if we don't have to
        if not HAVE_GMSH:
            raise RuntimeError("gmsh is not available to build meshfiles")
        
        if getMPIRankWorld() == 0:
            filename = self.generateScriptFile(nodeFieldSize)
            self.filename = filename
        
        verbosity = 3
        if mshName is None:
            mshName = filename[:-4]+".msh"

        if gmshGeo2Msh(filename, mshName, 3, 1, verbosity)!=0:
            raise RuntimeError("Call out to gmsh failed")
        dom=ReadGmsh(mshName, 3, diracTags=self.__tags, diracPoints=self.__points)
        return dom

    def getFileName(self):
        return self.filename
