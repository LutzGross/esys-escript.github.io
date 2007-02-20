"""
@author: John NGUI
"""

import vtk
import tempfile, os, sys
from constant import Source 
try:
	import esys.escript 
except ImportError:
	print "Warning: importing esys.escript failed."

class DataCollector:
	"""
	Class that defines a data collector which dealrs with the source 
	of data for the visualisation.
	"""

	def __init__(self, source = Source.XML):
		"""
		Initialise the data collector.

		@type source: L{Source <constant.Source>} constant
		@param source: Source data type
		"""

		self.__source=source

		if(source == Source.XML): # Source is an XML file.
			self.__vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
		# This is a temporary solution using a file: 
		elif (self.__source == Source.ESCRIPT):
			self.__vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
			self.__tmp_file=tempfile.mkstemp(suffix=".xml")[1]
			self.__vtk_xml_reader.SetFileName(self.__tmp_file)
			self.__output = self.__vtk_xml_reader.GetOutput()

	def __del__(self):
		"""
		Perform some clean up of the temporary file.
		"""

		if (self.__source == Source.ESCRIPT):
			if os.access(self.__tmp_file,os.F_OK): os.unlink(self.__tmp_file)

	def setFileName(self, file_name):
		"""
		Set the source file name to read. 

		@type file_name: String
		@param file_name: Name of the file to read
		"""

		if (self.__source == Source.XML):
			self.__vtk_xml_reader.SetFileName(file_name)

			# NOTE: Update must be called after SetFileName to make the reader 
			# up to date. Otherwise, some output values may be incorrect.
			self.__vtk_xml_reader.Update() 
			self.__output = self.__vtk_xml_reader.GetOutput()
			self.__get_attribute_lists()
		else:
			raise ValueError("source type %s does not support setFileName"%self.__source)

	def setData(self,**args):
		"""
		sets the data using. <name>=<data> sets the data tagged by <name> 
		to the object <data>. It is expected that the data are given in 
		an appropriate source type.
		"""

		if self.__source == Source.ESCRIPT:
			esys.escript.saveVTK(self.__tmp_file,**args)
			self.__vtk_xml_reader.Update()
		else:
			raise ValueError("source type %s does not support setData"%self.__source)

	def _setActiveScalar(self, scalar):
		"""
		Specify the scalar field to load from the source file.

		@type scalar: String
		@param scalar: Scalar field to load from the file. 
		"""
		
		if scalar in self.__point_attribute['scalars']:
			self._getOutput().GetPointData().SetActiveScalars(scalar)
		elif scalar in self.__cell_attribute['scalars']:
			self._getOutput().GetCellData().SetActiveScalars(scalar)
		else:
			print "\nSorry, no scalar called '%s' is available.\n" % scalar
			sys.exit(0)	
			

	def _setActiveVector(self, vector):
		"""
		Specify the vector field to load from the source file. 

		@type vector: String
		@param vector: Vector field to load from the file. 
		"""
		
		if vector in self.__point_attribute['vectors']:
			self._getOutput().GetPointData().SetActiveVectors(vector)
		elif vector in self.__cell_attribute['vectors']:
			self._getOutput().GetCellData().SetActiveVectors(vector)
		else:
			print "\nSorry, no vector called '%s' is available.\n" % vector
			sys.exit(0)	
			

	def _setActiveTensor(self, tensor):
		"""
		Specify the tensor field to load from the source file.

		@type tensor: String
		@param tensor: Tensor field to load from the file. 
		"""

		if tensor in self.__point_attribute['tensors']:
			self._getOutput().GetPointData().SetActiveTensors(tensor)
		elif tensor in self.__cell_attribute['tensors']:
			self._getOutput().GetCellData().SetActiveTensors(tensor)
		else:
			print "\nSorry, no tensor called '%s' is available.\n" % tensor
			sys.exit(0)	


	def __get_array_type(self, arr):
		num_components = arr.GetNumberOfComponents()

		if num_components == 1:
			return 'scalars'
		elif num_components == 3:
			return 'vectors'
		elif num_components == 9:
			return 'tensors'

	def __get_attribute_list(self, data): 
		attribute = {'scalars':[], 'vectors':[], 'tensors':[]}
		if data:
			num_arrays = data.GetNumberOfArrays()
			for i in range(num_arrays):
				name = data.GetArrayName(i)
				type = self.__get_array_type(data.GetArray(i))
				attribute[type].extend([name])

		return attribute

	def __get_attribute_lists(self):
		self.__point_attribute = \
				self.__get_attribute_list(self._getOutput().GetPointData())
		
		self.__cell_attribute = \
				self.__get_attribute_list(self._getOutput().GetCellData())

	def _getScalarRange(self):
		"""
		Return the scalar range.

		@rtype: Two column tuple containing numbers
		@return: Scalar range
		"""
	
		if(len(self.__point_attribute['scalars']) != 0):
			return self._getOutput().GetPointData().GetScalars().GetRange(-1)
		elif(len(self.__cell_attribute['scalars']) != 0):
			return self._getOutput().GetCellData().GetScalars().GetRange(-1)
		else:
			print "\nSorry, no scalar is available.\n"	
			sys.exit(0)

	def _getVectorRange(self):
		"""
		Return the vector range.
	
		@rtype: Two column tuple containing numbers
		@return: Vector range
		"""
		

		# NOTE: Generally GetRange(-1) returns the correct vector range. 
		# However, there are certain data sets where GetRange(-1) seems 
		# to return incorrect mimimum vector although the maximum vector is 
		# correct. As a result, the mimimum vector has been hard coded to 0.0
		# to accommodate those incorrect cases.
		if(len(self.__point_attribute['vectors']) != 0):
			vector_range = \
					self._getOutput().GetPointData().GetVectors().GetRange(-1)
			return (0.0, vector_range[1])
		elif(len(self.__cell_attribute['vectors']) != 0): 
			vector_range = \
					self._getOutput().GetCellData().GetVectors().GetRange(-1)
			return (0.0, vector_range[1])
		else:
			print "\nSorry, no vector is available.\n"	
			sys.exit(0)

	def _getTensorRange(self):
		"""
		Return the tensor range.

		@rtype: Two column table
		@return: Tensor range
		"""

		if(len(self.__point_attribute['tensors']) != 0):
			return self._getOutput().GetPointData().GetTensors().GetRange(-1)
		elif(len(self.__cell_attribute['tensors']) != 0): 
			return self._getOutput().GetCellData().GetTensors().GetRange(-1)
		else:
			print "\nSorry, no tensor is available.\n"	
			sys.exit(0)

	def _getOutput(self):
		"""
		Return the output of the data collector.

		@rtype: vtkUnstructuredGrid
		@return: Unstructured grid
		"""
		#cell_to_point_data = vtk.vtkCellDataToPointData()
		#cell_to_point_data.SetInput(self.__output)
		#cell_to_point_data.Update()
		#return (cell_to_point_data.GetOutput())
		return self.__output

	
###############################################################################


from constant import ImageFormat

class ImageReader:
	"""
	Class that defines an image reader.
	"""

	def __init__(self, format):
		"""	
		Initialise the image reader.

		@type format: String
		@param format: Format of the image 
		"""

		self.__format = format
		self.__vtk_image_reader = self.getImageReader()

	def getImageReader(self):
		"""
		Return the appropriate image reader based on the supplied image 
		format.

		@rtype: vtkImageReader2 (i.e. vtkJPEGReader, etc)
		@return: Image reader 
		"""

		if(self.__format == ImageFormat.JPG):
			return vtk.vtkJPEGReader()	
		elif(self.__format == ImageFormat.BMP):
			return vtk.vtkBMPReader()
		elif(self.__format == ImageFormat.PNM):
			return vtk.vtkPNMReader()
		elif(self.__format == ImageFormat.PNG):
			return vtk.vtkPNGReader()
		elif(self.__format == ImageFormat.TIF):
			return vtk.vtkTIFFReader()

	def setFileName(self, file_name):
		"""
		Set the image file name.

		@type file_name: String
		@param file_name: Image file name to be read 
		"""

		self.__vtk_image_reader.SetFileName(file_name)

	def _getOutput(self):
		"""
		Return the output of the image reader.

		@rtype: vtkImageData
		@return Image data
		"""

		return self.__vtk_image_reader.GetOutput()

