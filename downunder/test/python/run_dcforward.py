from __future__ import print_function
##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
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

from esys.downunder import *
from esys.escript import *
from esys.escriptcore.testing import *
import esys.escriptcore.utestselect as unittest
import numpy as np

try:
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT='ref_data'

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

try:
    from esys.finley import Rectangle, Brick
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

HAVE_GMSH = getEscriptParamInt("GMSH_SUPPORT")

mpisize = getMPISizeWorld()
mpirank = getMPIRankWorld()

@unittest.skipIf(not HAVE_FINLEY, "Finley module not available")
class TestDCResistivityForward(unittest.TestCase):
    def test_getPotential2d(self):
        dom=Rectangle(20,20,l0=1000,l1=-1000,d1=mpisize)
        extents=[1000,1000]
        primaryConductivity=Data(1/100., ContinuousFunction(dom))
        secondaryConductivity=Data(1/130., ContinuousFunction(dom))
        current=1000.
        a=0.05*extents[0]
        start = [0.25*extents[0]]
        directionVector=[1]
        numElectrodes = 10

        self.assertRaises(NotImplementedError,
                lambda:PolePoleSurvey(dom, primaryConductivity,
                        secondaryConductivity, current, a, start,
                        directionVector, numElectrodes))

    def test_getpotential3dPolePole(self):
        structured=False
        if structured:
            extents=[1000,1000,1000]
            dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=30.0
            bufferThickness=3000
            extents=[1000,2000,2000]
            electrodeDict={}
            lcDiv=10
            electrodeDict["e1"]=[-0.4*extents[0], 0, 0,lc/lcDiv]
            electrodeDict["e2"]=[-0.2*extents[0], 0, 0,lc/lcDiv]
            electrodeDict["e3"]=[0.2*extents[0], 0, 0,lc/lcDiv]
            electrodeDict["e4"]=[0.4*extents[0], 0, 0,lc/lcDiv]
            runName=os.path.join(WORKDIR, "dcResPolePole%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=WORKDIR,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh")
            if mpirank==0: 
                os.unlink(runName+".msh")
        totalApparentRes = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        a=4*0.05*extents[0]
        midPoint = [0,0]
        directionVector=[1.,0.]
        numElectrodes = 4

        pps=PolePoleSurvey(dom, primaryConductivity, secondaryConductivity, current, a, midPoint, directionVector, numElectrodes)
        delPhi=pps.getPotential()
        totalApparentResList = pps.getApparentResistivityTotal()
        for i in totalApparentResList:
            res_a = abs(i-totalApparentRes)
            res_b = 0.075 * totalApparentRes
            self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

    def test_getPotential3dSchlumberger(self):
        structured=False
        numElectrodes = 12
        directionVector=[1.,0.]
        midPoint=[]
        totalApparentResVal = 130.
        current=0.5
        interval_a = 5
        interval_n = 5


        if structured:
            #does not work because finley does not allow the specification of domain origin
            extents=[200,200,200]
            dom=Brick(25,25,25,l0=(-extents[0]/2,extents[0]/2),l1=(-extents[1]/2,extents[1]/2),l2=-extents[2])
            midPoint = [0,0]
        else:
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=10.0
            bufferThickness=100
            extents=[200,200,200]
            midPoint = [0,0]
            lcDiv=10.0
            electrodes=[]
            start=[]
            start.append(midPoint[0] - (((numElectrodes-1)*interval_a)/2. * directionVector[0]))
            start.append(midPoint[1] - (((numElectrodes-1)*interval_a)/2. * directionVector[1]))
            electrodeTags=[]
            electrodeDict={}
            for i in range(numElectrodes):
                electrodes.append([start[0]+(directionVector[0]*i*interval_a), start[1]+(directionVector[1]*i*interval_a),0])
                electrodeTags.append("e%d"%i)
                electrodeDict[electrodeTags[i]]=[electrodes[i][0], electrodes[i][1], electrodes[i][2], lc/lcDiv]
            runName=os.path.join(WORKDIR, "dcResSchlum%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=WORKDIR,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh",fieldSize=[70,100])
            fn = domGen.getFileName()
            if mpirank==0: 
                os.unlink(runName+".msh")
            
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))    
        schs=SchlumbergerSurvey(dom, primaryConductivity, secondaryConductivity,
                current, interval_a, interval_n, midPoint, directionVector,
                numElectrodes)
        potentials=schs.getPotentialAnalytic()
        totalApparentRes=schs.getApparentResistivity(potentials[2])
        for i in totalApparentRes:
            for j in i:
                res_a = abs(j-totalApparentResVal)
                res_b = 0.05 * totalApparentResVal
                self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

    def test_getPotentialDipDip(self):
        structured=False
        totalApparentRes = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(25,25,25,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=10.0
            bufferThickness=300
            extents=[100,100,100]
            electrodeDict={}
            lcDiv=10
            electrodeDict["e0" ] = [-22.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e1" ] = [-18.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e2" ] = [-14.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e3" ] = [-10.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e4" ] = [-6.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e5" ] = [-2.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e6" ] = [ 2.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e7" ] = [ 6.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e8" ] = [ 10.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e9" ] = [ 14.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e10"] = [ 18.0, 0.0, 0, lc/lcDiv]
            electrodeDict["e11"] = [ 22.0, 0.0, 0, lc/lcDiv]
            runName=os.path.join(WORKDIR, "dcResdipdip%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=WORKDIR,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh",reUse=False)
            if mpirank==0: 
                os.unlink(runName+".msh")
        n=5
        totalApparentResVal = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 12
        # a=(.8*extents[0])/numElectrodes
        a=4
        midPoint = [0,0]
        directionVector=[1.,0.]
        dipdips=DipoleDipoleSurvey(dom, primaryConductivity, secondaryConductivity, current, a,n, midPoint, directionVector, numElectrodes)
        dipdips.getPotential()
        primaryApparentRes=dipdips.getApparentResistivityPrimary()
        SecondaryApparentRes=dipdips.getApparentResistivitySecondary()
        totalApparentRes=dipdips.getApparentResistivityTotal()
        for i in totalApparentRes:
            for j in i:
                res_a = abs(j-totalApparentResVal)
                res_b = 0.075 * totalApparentResVal
                self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

    def test_getPotentialWenner(self):
        structured=True
        totalApparentResVal = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=50.0
            bufferThickness=3000
            extents=[1000,2000,2000]
            electrodeDict={}
            lcDiv=10

            electrodeDict["e0" ] = [28.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e1" ] = [32.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e2" ] = [36.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e3" ] = [40.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e4" ] = [44.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e5" ] = [48.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e6" ] = [52.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e7" ] = [56.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e8" ] = [60.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e9" ] = [64.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e10"] = [68.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e11"] = [72.0, 48.0, 0, lc/lcDiv]

            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=WORKDIR,bufferThickness=bufferThickness,prism=None)
            runName=os.path.join(WORKDIR, "wenner%d-%d"%(lc,lc/lcDiv))
            dom = domGen.getDom(mshName=runName+".msh")
            if mpirank==0: 
                os.unlink(runName+".msh")
        totalApparentRes = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 8
        # a=(.8*extents[0])/numElectrodes
        a=2
        midPoint = [0.5*extents[0]+1,0.5*extents[1]]
        directionVector=[1.,0.]
        wenSurv=WennerSurvey(dom, primaryConductivity, secondaryConductivity,
                current, a, midPoint, directionVector, numElectrodes)
        wenSurv.getPotential()
        primaryApparentRes=wenSurv.getApparentResistivityPrimary()
        SecondaryApparentRes=wenSurv.getApparentResistivitySecondary()
        totalApparentRes=wenSurv.getApparentResistivityTotal()
        for i in totalApparentRes:
            res_a = abs(i-totalApparentResVal)
            res_b = 0.05 * totalApparentResVal
            self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

    def test_getPotentialPolDip(self):
        structured=False
        totalApparentRes = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(25,25,25,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=10.0
            bufferThickness=300
            extents=[100,100,100]
            electrodeDict={}
            lcDiv=10
            electrodeDict["e0" ] = [28.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e1" ] = [32.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e2" ] = [36.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e3" ] = [40.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e4" ] = [44.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e5" ] = [48.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e6" ] = [52.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e7" ] = [56.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e8" ] = [60.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e9" ] = [64.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e10"] = [68.0, 48.0, 0, lc/lcDiv]
            electrodeDict["e11"] = [72.0, 48.0, 0, lc/lcDiv]
            runName=os.path.join(WORKDIR, "dcRespoldip%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=WORKDIR,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh")
            if mpirank==0: 
                os.unlink(runName+".msh")
        n=5
        totalApparentResVal = 130.
        primaryConductivity   =  Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity =  Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 12
        # a=(.8*extents[0])/numElectrodes
        a=4
        midPoint = [0.5*extents[0],0.5*extents[1] - 2]
        directionVector=[1.,0.]
        poldips=PoleDipoleSurvey(dom, primaryConductivity,
                secondaryConductivity, current, a,n, midPoint,
                directionVector, numElectrodes)
        poldips.getPotential()
        primaryApparentRes=poldips.getApparentResistivityPrimary()
        SecondaryApparentRes=poldips.getApparentResistivitySecondary()
        totalApparentRes=poldips.getApparentResistivityTotal()
        for i in totalApparentRes:
            for j in i:
                res_a = abs(j-totalApparentResVal)
                res_b = 0.075 * totalApparentResVal
                self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

################################
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

    
