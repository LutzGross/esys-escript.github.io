"""
Class and functions associated with lights

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

class Light:

	def __init__(self, open_scene):
		self.open_scene = open_scene
		self.vtk_light = None

		self.setLight()

	def setLight(self):
		self.vtk_light = vtk.vtkLight()
		self.open_scene.getRenderer().AddLight(self.vtk_light)

	def setColor(self, red, green, blue):
		self.vtk_light.SetColor(red, green, blue)

	def setFocalPoint(self, x_coor, y_coor, z_coor):
		self.vtk_light.SetFocalPoint(x_coor, y_coor, z_coor)

	def setPosition(self, x_coor, y_coor, z_coor):
		self.vtk_light.SetPosition(x_coor, y_coor, z_coor)

	def setIntensity(self, intensity):
		self.vtk_light.SetIntensity(intensity)
