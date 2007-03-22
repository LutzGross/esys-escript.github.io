"""
@author: John NGUI
"""

import vtk

class MaskPoints:

	def __init__(self, object):
		self.__object = object
		self.__vtk_mask_points = vtk.vtkMaskPoints

		self.__setupMaskPoints()

	def __setupMaskPoints(self):
		self.__setInput()
		self.randomOn()

	def __setInput(self):
		self.__vtk_mask_points.SetInput(self.__object)

	def setRatio(self, ratio):
		self.__vtk_mask_points.SetOnRation(ratio)

	def randomOn(self):
		self.__vtk_mask_points.RandomModeOn()

	def randomOff(self):
		self.__vtk_mask_points.RandomModeOff()
	
	def _getOutput(self):
		return self.__vtk_mask_points.GetOutput()
