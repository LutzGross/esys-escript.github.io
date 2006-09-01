"""
Class and functions associated with cameras and views

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


class Camera:
	def __init__(self, open_scene):
		self.open_scene = open_scene
		self.vtk_camera = None 

		self.setCamera()

	def setCamera(self):
		self.vtk_camera = vtk.vtkCamera()		
		self.open_scene.getRenderer().SetActiveCamera(self.vtk_camera)

	def setClippingRange(self, near_clipping, far_clipping):
		self.vtk_camera.SetClippingRange(near_clipping, far_clipping)

	def setFocalPoint(self, x_coor, y_coor, z_coor):
		self.vtk_camera.SetFocalPoint(x_coor, y_coor, z_coor)	

	def setPosition(self, x_coor, y_coor, z_coor):
		self.vtk_camera.SetPosition(x_coor, y_coor, z_coor)

	def setViewUp(self, x_view, y_view, z_view):
		self.vtk_camera.SetViewUp(x_view, y_view, z_view)






class FrontView(Camera):
    pass

class BackView(Camera):
    pass

class TopView(Camera):
    pass

class BottomView(Camera):
    pass

class LeftView(Camera):
    pass

class RightView(Camera):
    pass

class IsometricView(Camera):
    pass
