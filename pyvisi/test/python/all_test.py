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
		TestMapOnScalarClip
from run_map_with_lazy_evaluation import TestMapLazy, TestMapOnPlaneCutLazy, \
		TestMapOnPlaneClipLazy, TestMapOnScalarClipLazy
from run_scene import TestSceneOneViewport, TestSceneFourViewports
from run_streamline import TestStreamLinePointSource, TestStreamLineModule, \
		TestStreamLineTube
from run_streamline_with_lazy_evaluation import TestStreamLine
from run_text import TestText2D
from run_velocity import TestVelocity2DArrowVectorColor, \
		TestVelocity2DArrowScalarColor, TestVelocity3DSecondOrder
from run_velocity_with_lazy_evaluation import TestVelocity, TestVelocityOnPlaneCut, TestVelocityOnPlaneClip

import unittest

if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera2D))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCamera3D))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpetScalarWarp))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpetVectorWarp))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCarpet))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContour))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContourGenerate))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContourOnPlaneCut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestContourOnPlaneClip))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml2DPointData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml2DCellData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml3DPointData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSourceXml3DCellData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoid))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidScaleResolution))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidOnPlaneCut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEllipsoidOnPlaneClip))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImageOnAMap))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImage))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLight2D))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLight3D))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestAccessLogo))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestLogo))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOneViewport))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapFourViewports))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap2DCellDataWithCellToPointConversion))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap2DCellDataWithoutCellToPointConversion))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DPointData))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DCellDataWithCellToPointConversion))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DCellDataWithoutCellToPointConversion))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapGreyScaleLut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneCut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneClip))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClip))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMap3DSecondOrder))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneCutLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnPlaneClipLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestMapOnScalarClipLazy))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneOneViewport))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestSceneFourViewports))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestText2D))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity2DArrowVectorColor))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity2DArrowScalarColor))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity3DSecondOrder))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocity))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocityOnPlaneCut))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestVelocityOnPlaneClip))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLinePointSource))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineModule))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLineTube))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestStreamLine))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptMap))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptVelocity))
	suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestEscriptEllipsoid))


	unittest.TextTestRunner(verbosity=2).run(suite)

