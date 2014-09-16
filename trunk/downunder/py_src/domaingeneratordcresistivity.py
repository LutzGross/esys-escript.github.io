# from esys.escript import getMPIWorldMax, getMPIRankWorld
from esys.escript import *
from esys.escript.linearPDEs import LinearPDE, SolverOptions
from esys.finley import ReadGmsh, ReadMesh
from esys.weipa import saveSilo
from math import pi
import esys.escript.pdetools           as pdetools
import tempfile, os
import shlex


import logging
logger=logging.getLogger('inv.DCResDomGenerator')


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
        if(len(extents)==3 or len(extents)==4):
            self.__extents=extents
        else:
            raise ValueError("extents should be of length 3 or 4")
        self.__extentLen = len(self.__extents)
        self.__electrodeDict=electrodeDict
        self.__lc=lc
        self.__scriptname=""
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
    
    def generateScriptFile(self):
       if self.__scriptname:
           os.unlink(self.__scriptname)
       self.__scriptname_set=False
       if(self.__tmpDir==None):
           tmp_f_id=tempfile.mkstemp(suffix=".geo")
       else:
           tmp_f_id=tempfile.mkstemp(suffix=".geo", dir=self.__tmpDir)
       self.__scriptname=tmp_f_id[1]
       os.close(tmp_f_id[0])

    def generateScriptString(self):
        print "setting up contant half space"
        pntCount=5
        leftStr ="-out0[2],"
        rightStr="-out0[4],"
        frontStr="-out0[5],"
        backStr ="-out0[3],"
        out=""
        
        if not self.__bufferThickness == None:
            out+="lc=%f;\n"%self.__lc
            out+="Point(1)={%f, %f, 0, lc};\n"%(-self.__bufferThickness , -self.__bufferThickness   )
            out+="Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), -self.__bufferThickness  )
            out+="Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), (self.__extents[1]+self.__bufferThickness)  )
            out+="Point(4)={%f, %f, 0, lc};\n"%( -self.__bufferThickness, (self.__extents[1]+self.__bufferThickness) )
            out+="Line(1) = {1,2} ;\n"
            out+="Line(2) = {3,2} ;\n"
            out+="Line(3) = {3,4} ;\n"
            out+="Line(4) = {4,1} ;\n"
            out+="Line Loop(5) = {4,1,-2,3} ; \n"
            out+="Plane Surface(6) = {5} ; \n"

            if self.__prism != None:
                pntCount+=4
                out+="Point(5)={%f, %f, -%f, lc};\n"%self.__prism[0]
                out+="Point(6)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]), self.__prism[0][1],self.__prism[0][2]                   )
                out+="Point(7)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]), (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2])
                out+="Point(8)={%f, %f, -%f, lc};\n"%(self.__prism[0][0]                     , (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2])
                out+="Line(5) = {5,6};\n"
                out+="Line(6) = {7,6};\n"
                out+="Line(7) = {7,8};\n"
                out+="Line(8) = {8,5};\n"
                out+="Line Loop(6) = {8,5,-6,7} ;\n"
                out+="Plane Surface(7) = {6} ; \n"

            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out+="Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3])
                pntCount+=1
            out+="out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%(self.__extents[2]+self.__bufferThickness)
            self.__pntList=str(range(5,pntCount))[1:-1]
            out+="Point{%s} In Surface{6};\n"%self.__pntList
            out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1)
            if self.__prism != None:    
                out+="out[]=Extrude {0, 0, -10.000000} { Surface {7};};\n"
                out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(2,2)
                out+="s=newreg;\n"
                out+="Compound Volume(s) = {1,2};\n"
                out+="Physical Volume(\"volume-%d\") = {s} ;\n"%(3)
            
        else:
            out+="lc=%f;\n"%self.__lc
            out+="Point(1)={0,-1,0,lc};\n"
            out+="Point(2)={%f,-1,0,lc};\n"%self.__extents[0]
            out+="Point(3)={%f,%f,0,lc};\n"%(self.__extents[0],(self.__extents[1]+1))
            out+="Point(4)={0,%f,0,lc};\n"%(self.__extents[1]+1)
            out+="Line(1) = {1,2} ;\n"
            out+="Line(2) = {3,2} ;\n"
            out+="Line(3) = {3,4} ;\n"
            out+="Line(4) = {4,1} ;\n"
            out+="Line Loop(5) = {4,1,-2,3} ; \n"
            out+="Plane Surface(6) = {5} ; \n"

            if self.__prism != None:
                pntCount+=4
                out+="Point(5)={%f, %f, -%f, lc};\n"%self.__prism[0]
                out+="Point(6)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]), self.__prism[0][1],self.__prism[0][2]                   )
                out+="Point(7)={%f, %f, -%f, lc};\n"%((self.__prism[0][0]+self.__prism[1][0]), (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2])
                out+="Point(8)={%f, %f, -%f, lc};\n"%(self.__prism[0][0]                     , (self.__prism[0][1]+self.__prism[1][1]),self.__prism[0][2])
                out+="Line(5) = {5,6};\n"
                out+="Line(6) = {7,6};\n"
                out+="Line(7) = {7,8};\n"
                out+="Line(8) = {8,5};\n"
                out+="Line Loop(6) = {8,5,-6,7} ;\n"
                out+="Plane Surface(7) = {6} ; \n"

            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out+="Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3])
                pntCount+=1

            out+="out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%(self.__extents[2]+bufferThickness)
            out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1)
            if self.__prism != None:    
                out+="out[]=Extrude {0, 0, -10.000000} { Surface {7};};\n"
                out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(2,2)
                out+="s=newreg;\n"
                out+="Compound Volume(s) = {1,2};\n"
                out+="Physical Volume(\"volume-%d\") = {s} ;\n"%(3)
            self.__pntList=str(range(5,pntCount))[1:-1]
            out+="Point{%s} In Surface{6};\n"%self.__pntList
        
        out+="Physical Surface(\"Top\") = { -6 };\n"
        out+="Physical Surface(\"Bottom\") = { -out%d[0] };\n"%0
        out+="Physical Surface(\"Left\") = { %s };\n"%leftStr[:-1]
        out+="Physical Surface(\"Right\") = { %s };\n"%rightStr[:-1]
        out+="Physical Surface(\"Front\") = { %s };\n"%frontStr[:-1]
        out+="Physical Surface(\"Back\") = { %s };\n"%backStr[:-1]
        if not self.__bufferThickness == None:
            out+="Field[1] = Box;\n"
            out+="Field[1].VIn=lc;\n"
            out+="Field[1].VOut=5*lc;\n"
            out+="Field[1].XMax=%f;\n"%self.__extents[0]
            out+="Field[1].XMin=0;\n"
            out+="Field[1].YMax=%f;\n"%self.__extents[1]
            out+="Field[1].YMin=0;\n"
            out+="Field[1].ZMax=0;\n"
            out+="Field[1].ZMin=-%f;\n"%self.__extents[2]
            out+="Field[2] = Attractor;\n"
            out+="Field[2].NodesList = {%s};\n"%self.__pntList
            out+="Field[3] = Threshold;\n"
            out+="Field[3].IField = 2;\n"
            out+="Field[3].LcMin = lc / 5;\n"
            out+="Field[3].LcMax = 100*lc;\n" # this value is so high because It should not play a role in field 4
            out+="Field[3].DistMin = 50;\n"
            out+="Field[3].DistMax = 100;\n"
            out+="Field[4] = Min;\n"
            out+="Field[4].FieldsList = {1, 3};\n"
            out+="Background Field = 4;\n"
            out+="Mesh.CharacteristicLengthExtendFromBoundary = 0;\n"
        self.__scriptString=out   

    def generateLayedScriptString(self, interfaces):
        pntCount=5
        leftStr ="-out0[2],"
        rightStr="-out0[4],"
        frontStr="-out0[5],"
        backStr ="-out0[3],"
        out=""
        extentCount=float(interfaces[0])

        if not self.__bufferThickness == None:
            out+="lc=%f;\n"%self.__lc
            out+="Point(1)={%f, %f, 0, lc};\n"%(-self.__bufferThickness , -self.__bufferThickness   )
            out+="Point(2)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), -self.__bufferThickness  )
            out+="Point(3)={%f, %f, 0, lc};\n"%( (self.__extents[0]+self.__bufferThickness), (self.__extents[1]+self.__bufferThickness)  )
            out+="Point(4)={%f, %f, 0, lc};\n"%( -self.__bufferThickness, (self.__extents[1]+self.__bufferThickness) )
            out+="Line(1) = {1,2} ;\n"
            out+="Line(2) = {3,2} ;\n"
            out+="Line(3) = {3,4} ;\n"
            out+="Line(4) = {4,1} ;\n"
            out+="Line Loop(5) = {4,1,-2,3} ; \n"
            out+="Plane Surface(6) = {5} ; \n"
            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out+="Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3])
                pntCount+=1
            #out+="out[]=Extrude {0, 0, -%f} { Surface {6};};\n"%self.__bufferThickness
            #out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(0,1)
            #leftStr+= "-out[2],"
            #rightStr+="-out[4],"
            #frontStr+="-out[5],"
            #backStr+= "-out[3],"
            out+="out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%interfaces[0]
            out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(1,1)
            #out+="Point{%s} In Surface{out[0]};\n"%str(range(5,pntCount))[1:-1]
            self.__pntList=str(range(5,pntCount))[1:-1]
            out+="Point{%s} In Surface{6};\n"%self.__pntList
            for i in range(1,len(interfaces)):
                extentCount+=float(interfaces[i])
                out+="out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,interfaces[i],i-1)
                out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(i+1,i+1)
                leftStr+= "-out%d[2],"%i
                rightStr+="-out%d[4],"%i
                frontStr+="-out%d[5],"%i
                backStr+= "-out%d[3],"%i
            i+=1
            out+="out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,self.__bufferThickness,i-1)
            out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(i+1,i+2)
            leftStr+= "-out%d[2],"%i
            rightStr+="-out%d[4],"%i
            frontStr+="-out%d[5],"%i
            backStr+= "-out%d[3],"%i
        else:
            out+="lc=%f;\n"%self.__lc
            out+="Point(1)={0,-1,0,lc};\n"
            out+="Point(2)={%f,-1,0,lc};\n"%self.__extents[0]
            out+="Point(3)={%f,%f,0,lc};\n"%(self.__extents[0],(self.__extents[1]+1))
            out+="Point(4)={0,%f,0,lc};\n"%(self.__extents[1]+1)
            out+="Line(1) = {1,2} ;\n"
            out+="Line(2) = {3,2} ;\n"
            out+="Line(3) = {3,4} ;\n"
            out+="Line(4) = {4,1} ;\n"
            out+="Line Loop(5) = {4,1,-2,3} ; \n"
            out+="Plane Surface(6) = {5} ; \n"
            for i in self.__electrodeDict:
                pntInfo=self.__electrodeDict[i]
                out+="Point(%d)={%f,%f,%f,%f};\n"%(pntCount,pntInfo[0],pntInfo[1],pntInfo[2],pntInfo[3])
                pntCount+=1
            self.__pntList=str(range(5,pntCount))[1:-1]
            out+="Point{%s} In Surface{6};\n"%self.__pntList
            out+="out0[]=Extrude {0, 0, -%f} { Surface {6};};\n"%interfaces[0]
            out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(0,1)
            
            for i in range(1,len(interfaces)):
                out+="out%d[]=Extrude {0, 0, -%f} { Surface {out%d[0]};};\n"%(i,interfaces[i],i-1)
                out+="Physical Volume(\"volume-%d\") = {%d} ;\n"%(i,i+1)
                leftStr+= "-out%d[2],"%i
                rightStr+="-out%d[4],"%i
                frontStr+="-out%d[5],"%i
                backStr+= "-out%d[3],"%i
        

        out+="Physical Surface(\"Top\") = { -6 };\n"
        out+="Physical Surface(\"Bottom\") = { -out%d[0] };\n"%i
        out+="Physical Surface(\"Left\") = { %s };\n"%leftStr[:-1]
        out+="Physical Surface(\"Right\") = { %s };\n"%rightStr[:-1]
        out+="Physical Surface(\"Front\") = { %s };\n"%frontStr[:-1]
        out+="Physical Surface(\"Back\") = { %s };\n"%backStr[:-1]
        if not self.__bufferThickness == None:
            out+="Field[1] = Box;\n"
            out+="Field[1].VIn=lc;\n"
            out+="Field[1].VOut=5*lc;\n"
            out+="Field[1].XMax=%f;\n"%self.__extents[0]
            out+="Field[1].XMin=0;\n"
            out+="Field[1].YMax=%f;\n"%self.__extents[1]
            out+="Field[1].YMin=0;\n"
            out+="Field[1].ZMax=0;\n"
            out+="Field[1].ZMin=-%f;\n"%extentCount
            out+="Field[2] = Attractor;\n"
            out+="Field[2].NodesList = {%s};\n"%self.__pntList
            out+="Field[3] = Threshold;\n"
            out+="Field[3].IField = 2;\n"
            out+="Field[3].LcMin = lc / 5;\n"
            out+="Field[3].LcMax = 100*lc;\n"
            out+="Field[3].DistMin = 50;\n"
            out+="Field[3].DistMax = 100;\n"
            out+="Field[4] = Min;\n"
            out+="Field[4].FieldsList = {1, 3};\n"
            out+="Background Field = 4;\n"
            out+="Mesh.CharacteristicLengthExtendFromBoundary = 0;\n"
        self.__scriptString=out        

    def runGmsh(self, args):
        if getMPIRankWorld() == 0:
            import subprocess
            try:
                ret = subprocess.call(args) // 256
            except:
                ret = 1
        else:
            ret = 0
        ret=getMPIWorldMax(ret)
        return ret    

    def getDom(self, mshName=None, interfaces=None):
        """
        generate the domain
        :param interfaces, specify a list of interfaces for a layered model. Doing this will
            ignore the z-extent. the layers will be tagged iterative from volume-0 to volume-(n-1)
        :type interfaces, list    
        """
        if mshName!=None:
            try:
                print "atempting to read mesh %s"%mshName
                if mshName[-4:]=='.msh':
                    dom=ReadGmsh(mshName ,3,diracTags=self.__tags, diracPoints=self.__points)     
                elif mshName[-4:]=='.fly':
                    dom=ReadMesh(mshName ,3,diracTags=self.__tags, diracPoints=self.__points)  
                print "mesh was read directly"
                return dom
            except RuntimeError:
                # this might happen even if the cause was not the non existence of the
                # file this should get checked at a later stage
                print "msh file does not exist will try generate it now"

        self.generateScriptFile()
        if getMPIRankWorld() == 0:
            if interfaces == None:
                self.generateScriptString()
            else:
                self.generateLayedScriptString(interfaces)
            print self.__scriptString
            open(self.__scriptname,"w").write(self.__scriptString)
        if mshName == None:    
            exe="gmsh -3 -order 1 -v 3 %s"%self.__scriptname
        else:
            #exe="gmsh -3 -order 1 -v 3 -o %s %s"%(mshName[:-4]+".msh", self.__scriptname)
            exe="gmsh -3 -order 1 -o %s %s"%(mshName[:-4]+".msh", self.__scriptname)
        args=shlex.split(exe)
        print args
        ret=self.runGmsh(args)
        if ret > 0:
            raise RuntimeError("Could not build mesh using: " + \
                    "%s"%" ".join(args) + "\nCheck gmsh is available")
        if mshName == None:
            dom=ReadGmsh(self.__scriptname[:-4]+".msh" ,3,diracTags=self.__tags, diracPoints=self.__points)                  
            # os.unlink(self.__scriptname[:-4]+".msh")
            print "unlinking;%s"%self.__scriptname[:-4]+".msh"
        else:
            if mshName[-4:]=='.msh':
                dom=ReadGmsh(mshName ,3,diracTags=self.__tags, diracPoints=self.__points)     
            elif mshName[-4:]=='.fly':
                #"dom=ReadMesh(mshName ,3,diracTags=self.__tags, diracPoints=self.__points)"
                raise RuntimeError("please convert msh to fly")
            # os.unlink(mshName+".msh")
            # print "unlinking;%s"%mshName
        #os.unlink(self.__scriptname)
        return dom


