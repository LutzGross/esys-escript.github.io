
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

"""
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class Sphere:
	"""
	Class that defines a sphere.
	"""

	def __init__(self):
		"""
		Initialise the sphere.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.Sphere is not running on more than one processor."
		self.__vtk_sphere = vtk.vtkSphereSource()

	def setThetaResolution(self, resolution):
		"""
		Set the theta resolution of the sphere.

		:type resolution: Number
		:param resolution: Theta resolution
		"""

		self.__vtk_sphere.SetThetaResolution(resolution)

	def setPhiResolution(self, resolution):
		"""
		Set the phi resolution of the sphere.

		:type resolution: Number
		:param resolution: Phi resolution
		"""

		self.__vtk_sphere.SetPhiResolution(resolution)

	def _getSphereOutput(self):
		"""
		Return the output of the sphere.

		:rtype: vtkPolyData
		:return: Polygonal data
		"""

		return self.__vtk_sphere.GetOutput()

