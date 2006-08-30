"""
class that shows a vector field by arrows

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane, L. Gross"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision:$"
__date__="$Date:$"

import vtk

class Arrows:

	def __init__(self, open_scene, data_collector):
		self.open_scene = open_scene
		self.data_collector = data_collector
		self.vtk_glyph = None
		self.vtk_arrows_mapper = None
		self.vtk_arrows_actor = None		

		self.setArrows()
		self.setMapper()
		self.setActor()
	
	# set up the glyph and use arrows as the source	
	def setArrows(self):
		vtk_arrows = vtk.vtkArrowSource()
		
		self.vtk_glyph = vtk.vtkGlyph3D()
		self.vtk_glyph.SetInput(self.data_collector.getReader().GetOutput())
		self.vtk_glyph.SetSource(vtk_arrows.GetOutput())
		self.vtk_glyph.SetVectorModeToUseVector()
		self.vtk_glyph.SetScaleModeToScaleByVector()
		self.vtk_glyph.SetColorModeToColorByScalar()
		self.vtk_glyph.SetScaleFactor(0.2)
	
	# set up the mapper and data	
	def setMapper(self):
		self.vtk_arrows_mapper = vtk.vtkPolyDataMapper()
		self.vtk_arrows_mapper.SetInput(
			self.vtk_glyph.GetOutput())
	
	# set up the actor and add the actor to the scene	
	def setActor(self):
		self.vtk_arrows_actor = vtk.vtkActor()
		self.vtk_arrows_actor.SetMapper(self.vtk_arrows_mapper)

		self.open_scene.getRenderer().AddActor(self.vtk_arrows_actor)		




#class ArrowsOnPlane:
"""
shows a vector field by arrows on a plane
"""
pass
