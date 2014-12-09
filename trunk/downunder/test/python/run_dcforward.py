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

from esys.downunder import *
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
import numpy as np
from esys.escript import *
from esys.weipa import saveSilo
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE
from esys.escript.pdetools import Locator

try:
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT='ref_data'

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

tmpDir=os.path.join(TEST_DATA_ROOT, "dc_forward")

try:
    from esys.finley import Rectangle, Brick
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

HAVE_GMSH = getEscriptParamInt("GMSH_SUPPORT")
GMSH_MPI = HAVE_GMSH and getEscriptParamInt("GMSH_MPI")

mpisize = getMPISizeWorld()

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
                lambda:polepoleSurvey(dom, primaryConductivity,
                        secondaryConductivity, current, a, start,
                        directionVector, numElectrodes))

    def test_getpotential3dPolePole(self):
        structured=False
        if structured:
            extents=[1000,1000,1000]
            dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if GMSH_MPI:
                raise unittest.SkipTest("MPI gmsh not currently supported")
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=50.0
            bufferThickness=3000
            extents=[1000,2000,2000]
            electrodeDict={}
            lcDiv=10
            electrodeDict["e1"]=[0.2*extents[0], 0.5*extents[1], 0,lc/lcDiv]
            electrodeDict["e2"]=[0.4*extents[0], 0.5*extents[1], 0,lc/lcDiv]
            electrodeDict["e3"]=[0.6*extents[0], 0.5*extents[1], 0,lc/lcDiv]
            electrodeDict["e4"]=[0.8*extents[0], 0.5*extents[1], 0,lc/lcDiv]
            runName=os.path.join(TEST_DATA_ROOT, "dc_forward/dcResPolePole%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=tmpDir,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh")
            os.unlink(runName+".msh")
        totalApparentRes = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        a=4*0.05*extents[0]
        midPoint = [0.5*extents[0],0.5*extents[1]]
        directionVector=[1.,0.]
        numElectrodes = 4

        pps=polepoleSurvey(dom, primaryConductivity, secondaryConductivity, current, a, midPoint, directionVector, numElectrodes)
        delPhi=pps.getPotential()
        totalApparentResList = pps.getApparentResistivityTotal()
        for i in totalApparentResList:
            self.assertTrue(abs(i-totalApparentRes) < 0.05 * totalApparentRes)

    def test_getPotential3dSchlumberger(self):
        structured=True
        totalApparentRes = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if GMSH_MPI:
                raise unittest.SkipTest("MPI gmsh not currently supported")
            if not HAVE_GMSH:
                raise unittest.SkipTest("gmsh required for test")
            lc=50.0
            bufferThickness=3000
            extents=[1000,2000,2000]
            electrodeDict={}
            lcDiv=10
            electrodeDict["e1"]=[440., 0.5*extents[1], 0,lc/lcDiv]
            electrodeDict["e2"]=[480., 0.5*extents[1], 0,lc/lcDiv]
            electrodeDict["e3"]=[520., 0.5*extents[1], 0,lc/lcDiv]
            electrodeDict["e4"]=[560., 0.5*extents[1], 0,lc/lcDiv]
            runName=os.path.join(TEST_DATA_ROOT, "dc_forward/dcResSchlum%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=tmpDir,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh")
            os.unlink(runName+".msh")
        totalApparentResVal = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 12
        midPoint = [0.5*extents[0]+1,0.5*extents[1]]
        directionVector=[1.,0.]
        interval_n = 5
        interval_a = 2
        schs=SchlumbergerSurvey(dom, primaryConductivity, secondaryConductivity,
                current, interval_a, interval_n, midPoint, directionVector,
                numElectrodes)
        schs.getPotential()
        primaryApparentRes=schs.getApparentResistivityPrimary()
        SecondaryApparentRes=schs.getApparentResistivitySecondary()
        totalApparentRes=schs.getApparentResistivityTotal()
        for i in totalApparentRes:
            for j in i:
                self.assertTrue(abs(j-totalApparentResVal) < 0.05 * totalApparentResVal)

    def test_getPotentialDipDip(self):
        structured=False
        totalApparentRes = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(25,25,25,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if GMSH_MPI:
                raise unittest.SkipTest("MPI gmsh not currently supported")
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
            runName=os.path.join(TEST_DATA_ROOT, "dc_forward/dcResdipdip%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=tmpDir,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh")
            os.unlink(runName+".msh")
        n=5
        totalApparentResVal = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 12
        # a=(.8*extents[0])/numElectrodes
        a=4
        midPoint = [0.5*extents[0],0.5*extents[1] - 2]
        directionVector=[1.,0.]
        dipdips=DipoleDipoleSurvey(dom, primaryConductivity, secondaryConductivity, current, a,n, midPoint, directionVector, numElectrodes)
        dipdips.getPotential()
        primaryApparentRes=dipdips.getApparentResistivityPrimary()
        SecondaryApparentRes=dipdips.getApparentResistivitySecondary()
        totalApparentRes=dipdips.getApparentResistivityTotal()
        for i in totalApparentRes:
            for j in i:
                self.assertTrue(abs(j-totalApparentResVal) < 0.05 * totalApparentResVal)

    def test_getPotentialWenner(self):
        structured=True
        totalApparentResVal = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if GMSH_MPI:
                raise unittest.SkipTest("MPI gmsh not currently supported")
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

            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=tmpDir,bufferThickness=bufferThickness,prism=None)
            runName="dc_forward/wenner%d-%d"%(lc,lc/lcDiv)
            dom = domGen.getDom(mshName=runName+".msh")
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
                self.assertTrue(abs(i-totalApparentResVal ) < 0.05 * totalApparentResVal)


    def test_getPotentialPolDip(self):
        structured=False
        totalApparentRes = 130.
        if structured:
            extents=[100,100,100]
            dom=Brick(25,25,25,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            if GMSH_MPI:
                raise unittest.SkipTest("MPI gmsh not currently supported")
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
            runName=os.path.join(TEST_DATA_ROOT, "dc_forward/dcRespoldip%d-%d"%(lc,lc/lcDiv))
            domGen=DCResDomGenerator(extents, electrodeDict,lc=lc,tmpDir=tmpDir,bufferThickness=bufferThickness,prism=None)
            dom = domGen.getDom(mshName=runName+".msh")
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
        poldips=poledipoleSurvey(dom, primaryConductivity,
                secondaryConductivity, current, a,n, midPoint,
                directionVector, numElectrodes)
        poldips.getPotential()
        primaryApparentRes=poldips.getApparentResistivityPrimary()
        SecondaryApparentRes=poldips.getApparentResistivitySecondary()
        totalApparentRes=poldips.getApparentResistivityTotal()
        for i in totalApparentRes:
            for j in i:
                self.assertTrue(abs(j-totalApparentResVal) < 0.05 * totalApparentResVal)

################################
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

    
