
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

from run_camera import TestCamera2D, TestCamera3D 
from run_carpet import TestCarpetScalarWarp, TestCarpetVectorWarp 
from run_carpet_with_lazy_evaluation import TestCarpet
from run_contour import TestContour
from run_contour_with_lazy_evaluation import TestContourGenerate, \
		TestContourOnPlaneCut, TestContourOnPlaneClip
from run_datacollector import TestSourceXml2DPointData, \
		TestSourceXml2DCellData,TestSourceXml3DPointData, \
		TestSourceXml3DCellData 
from run_ellipsoid import TestEllipsoid
from run_ellipsoid_with_lazy_evaluation import TestEllipsoidScaleResolution, \
		TestEllipsoidOnPlaneCut, TestEllipsoidOnPlaneClip
from run_escript_with_lazy_evaluation import TestEscriptMap, \
		TestEscriptVelocity, TestEscriptEllipsoid 
from run_image import TestImageOnAMap
from run_image_with_lazy_evaluation import TestImage
from run_light import TestLight2D, TestLight3D
from run_logo import TestAccessLogo
from run_logo_with_lazy_evaluation import TestLogo
from run_map import TestMapOneViewport, TestMapFourViewports, \
		TestMap2DCellDataWithCellToPointConversion, \
		TestMap2DCellDataWithoutCellToPointConversion, TestMap3DPointData, \
		TestMap3DCellDataWithCellToPointConversion, \
		TestMap3DCellDataWithoutCellToPointConversion, TestMap3DSecondOrder, \
		TestMapGreyScaleLut, TestMapOnPlaneCut, TestMapOnPlaneClip, \
		TestMapOnScalarClip, \
		TestMapOnScalarClipWithRotation
from run_map_with_lazy_evaluation import TestMapLazy, TestMapOnPlaneCutLazy, \
		TestMapOnPlaneClipLazy, TestMapOnScalarClipLazy, \
		TestMapOnScalarClipWithRotationLazy
from run_scene import TestSceneOneViewport, TestSceneFourViewports, TestSceneInit
from run_streamline import TestStreamLinePointSource, TestStreamLineModule, \
		TestStreamLineTube
from run_streamline_with_lazy_evaluation import TestStreamLine
from run_text import TestText2D
from run_velocity import TestVelocity2DArrowVectorColor, \
		TestVelocity2DArrowScalarColor, TestVelocity3DSecondOrder
from run_velocity_with_lazy_evaluation import TestVelocity, TestVelocityOnPlaneCut, TestVelocityOnPlaneClip
from run_exporter import TestVRMLExporter, TestIVExporter
from run_legend import TestLegend
from run_legend_with_lazy_evaluation import TestLegendWithLazyEvaluation
from run_movie_with_lazy_evaluation import TestGenerateMovie
from run_rectangle import TestRectangleOnAMap
from esys.escript import getMPISizeWorld, MPIBarrierWorld

import unittest



if __name__ == '__main__':
    suite = unittest.TestSuite()
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera2D))
    else:
	print "test TestCamera2D is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera3D))
    else:
	print "test TestCamera3D is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpetScalarWarp))
    else:
	print "test TestCarpetScalarWarp is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpetVectorWarp))
    else:
	print "test TestCarpetVectorWarp is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContour))
    else:
	print "test TestContour is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpet))
    else:
	print "test TestCarpet is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContourGenerate))
    else:
	print "test TestContourGenerate is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContourOnPlaneCut))
    else:
	print "test TestContourOnPlaneCut is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContourOnPlaneClip))
    else:
	print "test TestContourOnPlaneClip is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml2DPointData))
    else:
	print "test TestSourceXml2DPointData is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml2DCellData))
    else:
	print "test TestSourceXml2DCellData is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml3DPointData))
    else:
	print "test TestSourceXml3DPointData is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml3DCellData))
    else:
	print "test TestSourceXml3DCellData is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoid))
    else:
	print "test TestEllipsoid is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidScaleResolution))
    else:
	print "test TestEllipsoidScaleResolution is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidOnPlaneCut))
    else:
	print "test TestEllipsoidOnPlaneCut is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidOnPlaneClip))
    else:
	print "test TestEllipsoidOnPlaneClip is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImageOnAMap))
    else:
	print "test TestImageOnAMap is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImage))
    else:
	print "test TestImage is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLight2D))
    else:
	print "test TestLight2D is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLight3D))
    else:
	print "test TestLight3D is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAccessLogo))
    else:
	print "test TestAccessLogo is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLogo))
    else:
	print "test TestLogo is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOneViewport))
    else:
	print "test TestMapOneViewport is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapFourViewports))
    else:
	print "test TestMapFourViewports is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap2DCellDataWithCellToPointConversion))
    else:
	print "test TestMap2DCellDataWithCellToPointConversion is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap2DCellDataWithoutCellToPointConversion))
    else:
	print "test TestMap2DCellDataWithoutCellToPointConversion is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DPointData))
    else:
	print "test TestMap3DPointData is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DCellDataWithCellToPointConversion))
    else:
	print "test TestMap3DCellDataWithCellToPointConversion is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DCellDataWithoutCellToPointConversion))
    else:
	print "test TestMap3DCellDataWithoutCellToPointConversion is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapGreyScaleLut))
    else:
	print "test TestMapGreyScaleLut is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneCut))
    else:
	print "test TestMapOnPlaneCut is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneClip))
    else:
	print "test TestMapOnPlaneClip is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClip))
    else:
	print "test TestMapOnScalarClip is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClipWithRotation))
    else:
	print "test TestMapOnScalarClipWithRotation is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DSecondOrder))
    else:
	print "test TestMap3DSecondOrder is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapLazy))
    else:
	print "test TestMapLazy is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneCutLazy))
    else:
	print "test TestMapOnPlaneCutLazy is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneClipLazy))
    else:
	print "test TestMapOnPlaneClipLazy is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClipLazy))
    else:
	print "test TestMapOnScalarClipLazy is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClipWithRotationLazy))
    else:
	print "test TestMapOnScalarClipWithRotationLazy is dropped as MPI size > 1"
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneInit))
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneOneViewport))
    else:
	print "test TestSceneOneViewport is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneFourViewports))
    else:
	print "test TestSceneFourViewports is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestText2D))
    else:
	print "test TestText2D is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity2DArrowScalarColor))
    else:
	print "test TestVelocity2DArrowScalarColor is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity2DArrowVectorColor))
    else:
	print "test TestVelocity2DArrowVectorColor is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity3DSecondOrder))
    else:
	print "test TestVelocity3DSecondOrder is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity))
    else:
	print "test TestVelocity is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocityOnPlaneCut))
    else:
	print "test TestVelocityOnPlaneCut is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocityOnPlaneClip))
    else:
	print "test TestStreamLinePointSource is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLinePointSource))
    else:
	print "test TestStreamLineModule is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineModule))
    else:
	print "test TestStreamLineTube is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineTube))
    else:
	print "test TestStreamLine is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLine))
    else:
	print "test TestEscriptMap is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptMap))
    else:
	print "test TestEscriptVelocity is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptVelocity))
    else:
	print "test TestEscriptEllipsoid is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptEllipsoid))
    else:
	print "test TestVRMLExporter is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVRMLExporter))
    else:
	print "test TestIVExporter is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIVExporter))
    else:
	print "test TestLegend is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLegend))
    else:
	print "test TestLegendWithLazyEvaluation is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLegendWithLazyEvaluation))
    else:
	print "test TestGenerateMovie is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	if os.system('ppmtompeg')/256==127:
	    print "test TestGenerateMovie is dropped as ppmtompeg is not available"
	else:
	    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestGenerateMovie))
    else:
	print "test TestRectangleOnAMap is dropped as MPI size > 1"
    if getMPISizeWorld() == 1:
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestRectangleOnAMap))
    else:
	print "test TestRectangleOnAMap is dropped as MPI size > 1"
    unittest.TextTestRunner(verbosity=2).run(suite)
