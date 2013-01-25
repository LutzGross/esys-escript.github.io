#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


import vtk

class Glyph3D:
	"""
	Class that defines 3D glyphs.
	"""

	def __init__(self):
		"""
		Initialise the 3D glyph.
		"""

		self.__vtk_glyph3D = vtk.vtkGlyph3D()

	def _setupGlyph3D(self, object, source):
		"""
		Setup the 3D glyph.

		@type object: vtkDataSet, etc
		@param object: Input for the 3D glyph 
		@type source: vtkPolyData 	
		@param source: Source for the 3D glyph (i.e. Arrow2D, Arrow3D, etc)
		"""

		self.__object = object
		self.__source = source

		self.__setInput()
		self.__setSource()
		self.setScaleModeByVector()	
		self.__setVectorModeByVector()
		self.__setClampingOn()

		self.__setScalingOn()
		self.__setOrientOn()

	def __setInput(self):
		"""
		Set the input for the 3D glyph.	
		"""

		self.__vtk_glyph3D.SetInput(self.__object)		

	def __setSource(self):
		"""
		Set the source for the 3D glyph.	
		"""

		self.__vtk_glyph3D.SetSource(self.__source)

	def setScaleModeByVector(self):
		"""
		Set the 3D glyph to scale according to the vector.
		"""

		self.__vtk_glyph3D.SetScaleModeToScaleByVector()

	def setScaleModeByScalar(self):
		"""
		Set the 3D glyph to scale according to the scalar.
		"""

		self.__vtk_glyph3D.SetScaleModeToScaleByScalar()

	def _setColorModeByVector(self):
		"""
		Set the 3D glyph to color according to the vector.
		"""

		self.__vtk_glyph3D.SetColorModeToColorByVector()

	def _setColorModeByScalar(self):
		"""
		Set the 3D glyph to color according to the scalar.
		"""

		self.__vtk_glyph3D.SetColorModeToColorByScalar()

	def __setVectorModeByVector(self):
		"""
		Set the 3D glyph vector mode according to the vector.
		"""

		self.__vtk_glyph3D.SetVectorModeToUseVector()

	def setScaleFactor(self, scale_factor):
		"""
		Set the 3D glyph scale factor.
		
		@type scale_factor: Number
		@param scale_factor: Scale factor
		"""

		self.__vtk_glyph3D.SetScaleFactor(scale_factor)

	def __setClampingOn(self):
		"""
		Enable clamping of "scalar" values to range.
		"""

		self.__vtk_glyph3D.SetClamping(1)

	def __setScalingOn(self):
		"""
		Enable the scaling of the rendered object.	
		"""

		self.__vtk_glyph3D.ScalingOn()

	def __setOrientOn(self):
		"""
		Enable the orientation of the rendered object along the vector/normal.	
		"""

		self.__vtk_glyph3D.OrientOn()
		
	def _setRange(self, range):
		"""
		Set the range to map scalar values.

		@type range: Two column tuple containing numbers
		@param range: Range to map scalar values 
		"""

		self.__vtk_glyph3D.SetRange(range)

	def _getGlyph3DOutput(self):
		"""
		Return the output of the 3D glyph.

		@rtype: vtkPolyData
		@return Polygonal data
		"""

		return self.__vtk_glyph3D.GetOutput()


###############################################################################


class TensorGlyph:
	"""
	Class that defines tensor glyphs.
	"""

	def __init__(self):
		"""
		Initialise the tensor glyph.
		"""

		self.__vtk_tensor_glyph = vtk.vtkTensorGlyph()

	def _setupTensorGlyph(self, object, source):
		"""
		Setup the tensor glyph.

		@type object: vtkDataSet, etc
		@param object: Input for the 3D glyph 
		@type source: vtkPolyData 	
		@param source: Source for the 3D glyph (i.e. Sphere, etc)
		"""

		self.__object = object
		self.__source = source

		self.__setInput()
		self.__setSource()
		self.__vtk_tensor_glyph.ClampScalingOn()

	def __setInput(self):
		"""
		Set the input for the tensor glyph.
		"""

		self.__vtk_tensor_glyph.SetInput(self.__object)

	def __setSource(self):
		"""
		Set the source for the tensor glyph.
		"""
													
		self.__vtk_tensor_glyph.SetSource(self.__source)

	def setScaleFactor(self, scale_factor):
		"""
		Set the scale factor for the tensor glyph.
		
		@type scale_factor: Number
		@param scale_factor: Scale factor
		"""

		self.__vtk_tensor_glyph.SetScaleFactor(scale_factor)

	def setMaxScaleFactor(self, max_scale_factor):
		"""
		Set the maximum allowable scale factor for the tensor glyph.	

		@type max_scale_factor: Number
		@param max_scale_factor: Maximum allowable scale factor.
		"""

		self.__vtk_tensor_glyph.SetMaxScaleFactor(max_scale_factor)

	def _getTensorGlyphOutput(self):
		"""
		Return the output of the tensor glyph.

		@rtype: vtkPolyData
		@return: Polygonal data
		"""

		return self.__vtk_tensor_glyph.GetOutput()

