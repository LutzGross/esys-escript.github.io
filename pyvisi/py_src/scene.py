"""
defines a scene in which items are shown

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

class OpenScene:
	def __init__(self, renderer = "vtk-online", x_size = 800, y_size = 600):
		self.renderer = renderer
		self.x_size = x_size
		self.y_size = y_size
		self.vtk_renderer = None
		self.vtk_render_window = None

		self.setupRenderingWindow()

	def setupRenderingWindow(self):
		self.vtk_renderer = vtk.vtkRenderer()
		self.vtk_render_window = vtk.vtkRenderWindow()
		self.vtk_render_window.AddRenderer(self.vtk_renderer)
		self.vtk_render_window.SetSize(self.x_size, self.y_size)
		self.vtk_renderer.SetBackground(1, 1, 1)	

	def render(self):
		vtk_render_window_interactor = vtk.vtkRenderWindowInteractor()
		vtk_render_window_interactor.SetRenderWindow(self.vtk_render_window)
		vtk_render_window_interactor.Initialize()
		self.vtk_render_window.Render()
		vtk_render_window_interactor.Start()

	def getRenderer(self):
		return self.vtk_renderer
		
		


	



