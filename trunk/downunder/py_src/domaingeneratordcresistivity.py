from __future__ import print_function
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
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
    This class is used to generate an escript domain which is suitable for
    solving dc resistivity problems
    """
    def __init__(self, extents, electrodeDict, lc=0.1, tmpDir=None, prism=None, bufferThickness=None):
        """
        :param extents: x,y,z extents of the domain
        :type extents: list or tuple, len should=3
        :param electrodeDict: dictionary to hold coords and characteristic length of
        points to be used as electrodes.
        :type electrodeDict: dictionary
        :param lc:
        :type float
        :param prism: provide start point,extents and a extrude depth for a cubic prism
        :type [(x,y,z)_start,(x,y,z)_extent]
        :
        """
        if not HAS_FINLEY:
            raise RuntimeError("Finley module not available")
        if(len(extents)==3 or len(extents)==4):
            self.__extents=extents
        else:
            raise ValueError("extents should be of length 3 or 4")
        self.__extentLen = len(self.__extents)
        self.__electrodeDict=electrodeDict
        self.__lc=lc
        self.__scriptString=""
        self.__pntList=""
        self.__prism=prism
        self.__tags=[]
        self.__points=[]
        self.__tmpDir=tmpDir
        self.__bufferThickness=bufferThickness
        # logger.debug(electrodeDict)

        for i in electrodeDict:
            self.__tags.append(i)
            self.__points.append(electrodeDict[i][:-1])

    def generateScriptFile(self, interfaces=None):
        fd, filename = tempfile.mkstemp(suffix=".geo", dir=self.__tmpDir)
        os.close(fd)
        
        if interfaces is None:
            self.generateScriptString()
        else:
            self.generateLayedScriptString(interfaces)
        
        open(filename, "w").write(self.__scriptString)
        return filename

    def generateScriptString(self):
        pntCount=5
        leftStr ="-out0[2],"
        rightStr="-out0[4],"
        frontStr="-out0[5],"
        backStr ="-out0[3],"
        out=[]

        if not self.__bufferThickness == None:
            out.append("lc=%f;\n"%self.__lc)
            out.append("Point(1)={%f, %f, 0, lc};\n"%(-self.__bufferThickness , -self.__bufferThickness))
            out.append("Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), -self.__bufferThickness))
            out.append("Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), (self.__extents[1]+self.__bufferThickness)))
            out.append("Point(4)={%f, %f, 0, lc};\n"%( -self.__bufferThickness, (self.__extents[1]+self.__bufferThickness)))
            out.append("Line(1) = {1,2} ;\n")
            out.append("Line(2) = {3,2} ;\n")
            out.append("Line(3) = {3,4} ;\n")
            out.append("Line(4) = {4,1} ;\n")
            out.append("Line Loop(5) = {4,1,-2,3} ; \n")
            out.append("Plane Surface(6) = {5} ; \n")

            if self.__prism != None:
                pntCount+=4
                out.append("Point(5)={%f, %f, -%f, lc};\n"%self.__prism[0])
                out.append("Point(6)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]), self.__prism[0][1],self.__prism[0][2]))
                out.append("Point(7)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]), (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2]))
                out.append("Point(8)={%f, %f, -%f, lc};\n"%(self.__prism[0][0]                     , (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2]))
                out.append("Line(5) = {5,6};\n")
                out.append("Line(6) = {7,6};\n")
                out.append("Line(7) = {7,8};\n")
                out.append("Line(8) = {8,5};\n")
                out.append("Line Loop(6) = {8,5,-6,7} ;\n")
                out.append("Plane Surface(7) = {6} ; \n")

            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out.append("Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
                pntCount+=1
            out.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%(self.__extents[2]+self.__bufferThickness))
            self.__pntList=str([i for i in range(5,pntCount)])[1:-1]
            out.append("Point{%s} In Surface{6};\n"%self.__pntList)
            out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1))
            if self.__prism != None:
                out.append("out[]=Extrude {0, 0, -10.000000} { Surface {7};};\n")
                out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(2,2))
                out.append("s=newreg;\n")
                out.append("Compound Volume(s) = {1,2};\n")
                out.append("Physical Volume(\"volume-%d\") = {s} ;\n"%(3))

        else:
            out.append("lc=%f;\n"%self.__lc)
            out.append("Point(1)={0,-1,0,lc};\n")
            out.append("Point(2)={%f,-1,0,lc};\n"%self.__extents[0])
            out.append("Point(3)={%f,%f,0,lc};\n"%(self.__extents[0],(self.__extents[1]+1)))
            out.append("Point(4)={0,%f,0,lc};\n"%(self.__extents[1]+1))
            out.append("Line(1) = {1,2} ;\n")
            out.append("Line(2) = {3,2} ;\n")
            out.append("Line(3) = {3,4} ;\n")
            out.append("Line(4) = {4,1} ;\n")
            out.append("Line Loop(5) = {4,1,-2,3} ; \n")
            out.append("Plane Surface(6) = {5} ; \n")

            if self.__prism != None:
                pntCount+=4
                out.append("Point(5)={%f, %f, -%f, lc};\n"%self.__prism[0])
                out.append("Point(6)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]),
                        self.__prism[0][1],self.__prism[0][2]))
                out.append("Point(7)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]),
                        (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2]))
                out.append("Point(8)={%f, %f, -%f, lc};\n"%(self.__prism[0][0],
                        (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2]))
                out.append("Line(5) = {5,6};\n")
                out.append("Line(6) = {7,6};\n")
                out.append("Line(7) = {7,8};\n")
                out.append("Line(8) = {8,5};\n")
                out.append("Line Loop(6) = {8,5,-6,7} ;\n")
                out.append("Plane Surface(7) = {6} ; \n")

            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out.append("Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
                pntCount+=1

            out.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%(self.__extents[2]+bufferThickness))
            out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1))
            if self.__prism != None:
                out.append("out[]=Extrude {0, 0, -10.000000} { Surface {7};};\n")
                out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(2,2))
                out.append("s=newreg;\n")
                out.append("Compound Volume(s) = {1,2};\n")
                out.append("Physical Volume(\"volume-%d\") = {s} ;\n"%(3))
            self.__pntList=str([i for i in range(5,pntCount)])[1:-1]
            out.append("Point{%s} In Surface{6};\n"%self.__pntList)

        out.append("Physical Surface(\"Top\") = { -6 };\n")
        out.append("Physical Surface(\"Bottom\") = { -out%d[0] };\n"%0)
        out.append("Physical Surface(\"Left\") = { %s };\n"%leftStr[:-1])
        out.append("Physical Surface(\"Right\") = { %s };\n"%rightStr[:-1])
        out.append("Physical Surface(\"Front\") = { %s };\n"%frontStr[:-1])
        out.append("Physical Surface(\"Back\") = { %s };\n"%backStr[:-1])
        if not self.__bufferThickness is None:
            out.append("Field[1] = Box;\n")
            out.append("Field[1].VIn=lc;\n")
            out.append("Field[1].VOut=5*lc;\n")
            out.append("Field[1].XMax=%f;\n"%self.__extents[0])
            out.append("Field[1].XMin=0;\n")
            out.append("Field[1].YMax=%f;\n"%self.__extents[1])
            out.append("Field[1].YMin=0;\n")
            out.append("Field[1].ZMax=0;\n")
            out.append("Field[1].ZMin=-%f;\n"%self.__extents[2])
            out.append("Field[2] = Attractor;\n")
            out.append("Field[2].NodesList = {%s};\n"%self.__pntList)
            out.append("Field[3] = Threshold;\n")
            out.append("Field[3].IField = 2;\n")
            out.append("Field[3].LcMin = lc / 5;\n")
            out.append("Field[3].LcMax = 100*lc;\n") # this value is so high because It should not play a role in field 4
            out.append("Field[3].DistMin = 50;\n")
            out.append("Field[3].DistMax = 100;\n")
            out.append("Field[4] = Min;\n")
            out.append("Field[4].FieldsList = {1, 3};\n")
            out.append("Background Field = 4;\n")
            out.append("Mesh.CharacteristicLengthExtendFromBoundary = 0;\n")
        self.__scriptString = "".join(out)

    def generateLayedScriptString(self, interfaces):
        pntCount=5
        leftStr ="-out0[2],"
        rightStr="-out0[4],"
        frontStr="-out0[5],"
        backStr ="-out0[3],"
        out = []
        extentCount=float(interfaces[0])

        if not self.__bufferThickness == None:
            out.append("lc=%f;\n"%self.__lc)
            out.append("Point(1)={%f, %f, 0, lc};\n"%(-self.__bufferThickness , -self.__bufferThickness))
            out.append("Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), -self.__bufferThickness))
            out.append("Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), (self.__extents[1]+ self.__bufferThickness)))
            out.append("Point(4)={%f, %f, 0, lc};\n"%(-self.__bufferThickness, (self.__extents[1]+self.__bufferThickness)))
            out.append("Line(1) = {1,2} ;\n")
            out.append("Line(2) = {3,2} ;\n")
            out.append("Line(3) = {3,4} ;\n")
            out.append("Line(4) = {4,1} ;\n")
            out.append("Line Loop(5) = {4,1,-2,3} ; \n")
            out.append("Plane Surface(6) = {5} ; \n")
            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out.append("Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
                pntCount+=1
            #out.append("out[]=Extrude {0, 0, -%f} { Surface {6};};\n"%self.__bufferThickness)
            #out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(0,1))
            #leftStr+= "-out[2],"
            #rightStr+="-out[4],"
            #frontStr+="-out[5],"
            #backStr+= "-out[3],"
            out.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%interfaces[0])
            out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1))
            #out.append("Point{%s} In Surface{out[0]};\n"%str(range(5,pntCount))[1:-1])
            self.__pntList=str([i for i in range(5,pntCount)])[1:-1]
            out.append("Point{%s} In Surface{6};\n"%self.__pntList)
            for i in range(1,len(interfaces)):
                extentCount+=float(interfaces[i])
                out.append("out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,interfaces[i],i-1))
                out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(i+1,i+1))
                leftStr+= "-out%d[2],"%i
                rightStr+="-out%d[4],"%i
                frontStr+="-out%d[5],"%i
                backStr+= "-out%d[3],"%i
            i+=1
            out.append("out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,self.__bufferThickness,i-1))
            out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(i+1,i+2))
            leftStr+= "-out%d[2],"%i
            rightStr+="-out%d[4],"%i
            frontStr+="-out%d[5],"%i
            backStr+= "-out%d[3],"%i
        else:
            out.append("lc=%f;\n"%self.__lc)
            out.append("Point(1)={0,-1,0,lc};\n")
            out.append("Point(2)={%f,-1,0,lc};\n"%self.__extents[0])
            out.append("Point(3)={%f,%f,0,lc};\n"%(self.__extents[0],(self.__extents[1]+1)))
            out.append("Point(4)={0,%f,0,lc};\n"%(self.__extents[1]+1))
            out.append("Line(1) = {1,2} ;\n")
            out.append("Line(2) = {3,2} ;\n")
            out.append("Line(3) = {3,4} ;\n")
            out.append("Line(4) = {4,1} ;\n")
            out.append("Line Loop(5) = {4,1,-2,3} ; \n")
            out.append("Plane Surface(6) = {5} ; \n")
            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out.append("Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3]))
                pntCount+=1
            self.__pntList=str([i for i in range(5,pntCount)])[1:-1]
            out.append("Point{%s} In Surface{6};\n"%self.__pntList)
            out.append("out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%interfaces[0])
            out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(0,1))

            for i in range(1,len(interfaces)):
                out.append("out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,interfaces[i],i-1))
                out.append("Physical Volume(\"volume-%d\") = {%d} ;\n"%(i,i+1))
                leftStr+= "-out%d[2],"%i
                rightStr+="-out%d[4],"%i
                frontStr+="-out%d[5],"%i
                backStr+= "-out%d[3],"%i


        # out.append("Physical Surface(\"Top\") = { -6 };\n")
        # out.append("Physical Surface(\"Bottom\") = { -out%d[0] };\n"%i)
        # out.append("Physical Surface(\"Left\") = { %s };\n"%leftStr[:-1])
        # out.append("Physical Surface(\"Right\") = { %s };\n"%rightStr[:-1])
        # out.append("Physical Surface(\"Front\") = { %s };\n"%frontStr[:-1])
        # out.append("Physical Surface(\"Back\") = { %s };\n"%backStr[:-1])
        if not self.__bufferThickness == None:
            out.append("Field[1] = Box;\n")
            out.append("Field[1].VIn=lc;\n")
            out.append("Field[1].VOut=5*lc;\n")
            out.append("Field[1].XMax=%f;\n"%self.__extents[0])
            out.append("Field[1].XMin=0;\n")
            out.append("Field[1].YMax=%f;\n"%self.__extents[1])
            out.append("Field[1].YMin=0;\n")
            out.append("Field[1].ZMax=0;\n")
            out.append("Field[1].ZMin=-%f;\n"%extentCount)
            out.append("Field[2] = Attractor;\n")
            out.append("Field[2].NodesList = {%s};\n"%self.__pntList)
            out.append("Field[3] = Threshold;\n")
            out.append("Field[3].IField = 2;\n")
            out.append("Field[3].LcMin = lc / 5;\n")
            out.append("Field[3].LcMax = 100*lc;\n")
            out.append("Field[3].DistMin = 50;\n")
            out.append("Field[3].DistMax = 100;\n")
            out.append("Field[4] = Min;\n")
            out.append("Field[4].FieldsList = {1, 3};\n")
            out.append("Background Field = 4;\n")
            out.append("Mesh.CharacteristicLengthExtendFromBoundary = 0;\n")

        self.__scriptString = "".join(out)

    def getDom(self, mshName=None, interfaces=None, reUse=False):
        """
        Generates the domain.
        :param interfaces: Specify a list of interfaces for a layered model.
                           Doing this will ignore the z-extent. The layers
                           will be tagged iteratively from volume-0 to
                           volume-(n-1).
        :type interfaces: list
        :param reUse: should the msh be reused or should a new file be generated
        :type reUse: bool
        """

        if (mshName is not None and os.path.isfile(mshName) and reUse==True):
            if mshName[-4:]=='.msh':
                dom=ReadGmsh(mshName, 3, diracTags=self.__tags, diracPoints=self.__points)
            elif mshName[-4:]=='.fly':
                dom=ReadMesh(mshName, diracTags=self.__tags, diracPoints=self.__points)
            return dom

        # early exit so we don't even create files if we don't have to
        if not HAVE_GMSH:
            raise RuntimeError("gmsh is not available to build meshfiles")

        filename = "" # won't be used by non-0 ranks
        if getMPIRankWorld() == 0:
            filename = self.generateScriptFile(interfaces)
        verbosity = 3
        if mshName is None:
            mshName = filename[:-4]+".msh"
        else:
            mshName = mshName[:-4]+".msh"

        if gmshGeo2Msh(filename, mshName, 3, 1, verbosity)!=0:
            raise RuntimeError("Call out to gmsh failed")
        dom=ReadGmsh(mshName, 3, diracTags=self.__tags, diracPoints=self.__points)
        return dom

