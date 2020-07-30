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

from esys.downunder import *
from esys.escript import *
from esys.escriptcore.testing import *
import esys.escriptcore.utestselect as unittest
import numpy as np

try:
    TEST_DATA_ROOT=os.environ['DOWNUNDER_TEST_DATA_ROOT']
except KeyError:
    TEST_DATA_ROOT=os.path.join(os.getcwd(), "downunder/test/python/ref_data")

try:
    WORKDIR=os.environ['DOWNUNDER_WORKDIR']
except KeyError:
    WORKDIR='.'

try:
    from esys.finley import Rectangle, Brick, ReadGmsh
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

HAVE_GMSH = hasFeature("gmsh")

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

    @unittest.skipIf(not HAVE_GMSH, "gmsh not available")
    def test_getpotential3dPolePole(self):
        structured=False
        if structured:
            extents=[1000,1000,1000]
            dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        else:
            extents=[1000,2000,2000]
            tags=[]
            points=[]
            tags.append("e1")
            tags.append("e2")
            tags.append("e3")
            tags.append("e4")
            points.append([-0.4*extents[0], 0, 0])
            points.append([-0.2*extents[0], 0, 0])
            points.append([0.2*extents[0], 0,  0])
            points.append([0.4*extents[0], 0,  0])
            verbosity = 3 
            filename  = os.path.join(TEST_DATA_ROOT, "pole.geo")
            meshname  = os.path.join(TEST_DATA_ROOT, "dcResPolePole50-5.msh")
            gmshGeo2Msh(filename, meshname, 3, 1, verbosity)
            dom = ReadGmsh(meshname, 3, diracTags=tags, diracPoints=points)
        totalApparentRes = 130.
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        a=4*0.05*extents[0]
        midPoint = [0.5*extents[0],0.5*extents[1]]
        directionVector=[1.,0.]
        numElectrodes = 4

        pps=PolePoleSurvey(dom, primaryConductivity, secondaryConductivity, current, a, midPoint, directionVector, numElectrodes)
        delPhi=pps.getPotential()
        totalApparentResList = pps.getApparentResistivityTotal()
        for i in totalApparentResList:
            res_a = abs(i-totalApparentRes)
            res_b = 0.05 * totalApparentRes
            self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

    def test_getPotential3dSchlumberger(self):
        numElectrodes = 12
        directionVector=[1.,0.]
        midPoint=[]
        totalApparentResVal = 130.
        current=0.5
        interval_a = 5
        interval_n = 5
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
        electrodeLst=[]
        for i in range(numElectrodes):
            electrodes.append([start[0]+(directionVector[0]*i*interval_a), start[1]+(directionVector[1]*i*interval_a),0])
            electrodeTags.append("e%d"%i)
        runName=os.path.join(WORKDIR, "dcResSchlum%d-%d"%(lc,lc/lcDiv))
        filename  = os.path.join(TEST_DATA_ROOT, "schlum.geo")
        meshname  = os.path.join(TEST_DATA_ROOT, "dcResSchlum10-1.msh")
        verbosity=3
        gmshGeo2Msh(filename, meshname, 3, 1, verbosity)
        dom = ReadGmsh(meshname, 3, diracTags=electrodeTags, diracPoints=electrodes)
        if mpirank==0: 
            os.unlink(meshname)
            
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
                res_b = 0.1 * totalApparentResVal
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
            electrodes=[]
            tags=[]
            lcDiv=10
            tags.append("e0" )
            tags.append("e1" )
            tags.append("e2" )
            tags.append("e3" )
            tags.append("e4" )
            tags.append("e5" )
            tags.append("e6" )
            tags.append("e7" )
            tags.append("e8" )
            tags.append("e9" )
            tags.append("e10")
            tags.append("e11")
            electrodes.append([-22.0, 0.0, 0])
            electrodes.append([-18.0, 0.0, 0])
            electrodes.append([-14.0, 0.0, 0])
            electrodes.append([-10.0, 0.0, 0])
            electrodes.append([-6.0, 0.0, 0])
            electrodes.append([-2.0, 0.0, 0])
            electrodes.append([ 2.0, 0.0, 0])
            electrodes.append([ 6.0, 0.0, 0])
            electrodes.append([ 10.0, 0.0, 0])
            electrodes.append([ 14.0, 0.0, 0])
            electrodes.append([ 18.0, 0.0, 0])
            electrodes.append([ 22.0, 0.0, 0])
            filename  = os.path.join(TEST_DATA_ROOT, "dip.geo")
            meshname  = os.path.join(TEST_DATA_ROOT, "dcResdipdip-1.msh")
            verbosity=3
            gmshGeo2Msh(filename, meshname, 3, 1, verbosity)
            dom = ReadGmsh(meshname, 3, diracTags=tags, diracPoints=electrodes)
            if mpirank==0: 
                os.unlink(meshname)
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
                res_b = 0.1 * totalApparentResVal
                self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

    def test_getPotentialWenner(self):
        totalApparentResVal = 130.
        extents=[100,100,100]
        dom=Brick(50,50,50,l0=extents[0],l1=extents[1],l2=-extents[2])
        primaryConductivity=Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity=Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 8
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
            res_b = 0.1 * totalApparentResVal
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

            extents=[100,100,100]
            electrodes=[]
            tags=[]
            tags.append("e0" )
            tags.append("e1" )
            tags.append("e2" )
            tags.append("e3" )
            tags.append("e4" )
            tags.append("e5" )
            tags.append("e6" )
            tags.append("e7" )
            tags.append("e8" )
            tags.append("e9" )
            tags.append("e10")
            tags.append("e11")
            electrodes.append([-22.0, 0.0, 0])
            electrodes.append([-18.0, 0.0, 0])
            electrodes.append([-14.0, 0.0, 0])
            electrodes.append([-10.0, 0.0, 0])
            electrodes.append([-6.0, 0.0, 0])
            electrodes.append([-2.0, 0.0, 0])
            electrodes.append([ 2.0, 0.0, 0])
            electrodes.append([ 6.0, 0.0, 0])
            electrodes.append([ 10.0, 0.0, 0])
            electrodes.append([ 14.0, 0.0, 0])
            electrodes.append([ 18.0, 0.0, 0])
            electrodes.append([ 22.0, 0.0, 0])
            filename  = os.path.join(TEST_DATA_ROOT, "dip.geo")
            meshname  = os.path.join(TEST_DATA_ROOT, "dcRespoldip10-1.msh")
            verbosity=3
            gmshGeo2Msh(filename, meshname, 3, 1, verbosity)
            dom = ReadGmsh(meshname, 3, diracTags=tags, diracPoints=electrodes)

            if mpirank==0: 
                os.unlink(meshname)
        n=5
        totalApparentResVal = 130.
        primaryConductivity   =  Scalar(1/100., ContinuousFunction(dom))
        secondaryConductivity =  Scalar(1/130., ContinuousFunction(dom))
        current=1000.
        numElectrodes = 12
        # a=(.8*extents[0])/numElectrodes
        a=4
        midPoint = [0,0]
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
                res_b = 0.1 * totalApparentResVal
                self.assertLess(res_a, res_b, "result of %g greater than tolerance of %g"%(res_a, res_b))

################################
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

    
