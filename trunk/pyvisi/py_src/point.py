
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
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from position import GlobalPosition
from esys.escript import getMPISizeWorld
if getMPISizeWorld()==1: import vtk

class PointSource:
	"""
	Class that defines the source (location) to generate points. The points are
	generated within the radius of a sphere.
	"""

	def __init__(self):
		"""
		Initialise the point source.
		"""
                if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.PointSource is not running on more than one processor."
		self.__vtk_point_source = vtk.vtkPointSource()
		self.__center = None

	def _setupPointSource(self, object):
		"""
		Setup the point source.

		@type object: vtkPolyData, etc
		@param object: Used to define the location of the sphere 
		"""

		self.__object = object

		# Default number of points to generate within the sphere is 10.
		self.setPointSourceNumberOfPoints(10)
		# Default radius of the sphere is 0.5.
		self.setPointSourceRadius(0.5)

	def setPointSourceRadius(self, radius):
		"""
		Set the radius of the sphere.

		@type radius: Number
		@param radius: Radius of the sphere
		"""

		self.__vtk_point_source.SetRadius(radius)

	# This method is used to delay the execution of setting the point source
	# center.
	def setPointSourceCenter(self, center):
		"""
		Specity the sphere's center.
		
		@type center: L{GLobalPosition <position.GlobalPosition>} object
		@param center: Center of the sphere
		"""
		
		self.__center = center

	def _setPointSourceCenter(self):
		"""
		Set the sphere's center.
		"""

		self.__vtk_point_source.SetCenter(self.__center._getGlobalPosition())

	def setPointSourceNumberOfPoints(self, points):
		"""
		Set the number of points to generate within the sphere (the larger the
		number of points, the more streamlines are generated)

		@type points: Number
		@param points: Number of points to generate
		"""

		self.__vtk_point_source.SetNumberOfPoints(points)	

	def _getPointSourceOutput(self):
		"""
		Return the output of the point source.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_point_source.GetOutput()

	def _isPointSourceCenterSet(self):
		"""
		Return whether the center has been specified.

		@rtype: Boolean
		@return: True or False
		"""

		if(self.__center != None):
			return True
		else:
			return False


##############################################################################


class MaskPoints:
	"""
	Class that defines the masking of points 
	every n'th point.  This is useful to prevent the rendered object 
	from being cluttered with arrows or ellipsoids.
	"""

	def __init__(self):
		"""
		Initialise the mask points.
		"""
		if getMPISizeWorld()>1:
                     raise ValueError,"pyvisi.MaskPoints is not running on more than one processor."
                self.__vtk_mask_points = vtk.vtkMaskPoints()

	def _setupMaskPoints(self, object):
		"""
		Setup the mask points.

		@type object: vtkDataSet (i.e. vtkUnstructuredGrid, etc)
		@param object: Data source from which to mask points 
		"""

		self.__object = object
		self.__setInput()

	def __setInput(self):
		"""
		Set the input for the mask points.
		"""

		self.__vtk_mask_points.SetInput(self.__object)

	def setRatio(self, ratio):
		"""
		Mask every n'th point.

		@type ratio: Number
		@param ratio: Masking ratio
		"""

		self.__vtk_mask_points.SetOnRatio(ratio)

	def randomOn(self):
		"""
		Enables the randomization of the points selected for masking.
		"""

		self.__vtk_mask_points.RandomModeOn()

	def randomOff(self):
		"""
		Disables the randomization of the points selected for masking.
		"""

		self.__vtk_mask_points.RandomModeOff()
	
	def _getMaskPointsOutput(self):
		"""
		Return the output of the masked points.

		@rtype: vtkPolyData
		@return: Polygonal datda
		"""

		return self.__vtk_mask_points.GetOutput()


