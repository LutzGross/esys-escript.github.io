"""
@author: John NGUI
"""

import vtk
import tempfile, os, sys
from constant import Source, ColorMode
try:
	import esys.escript 
except ImportError:
	print "Warning: importing esys.escript failed."

class DataCollector:
	"""
	Class that defines a data collector. A data collector is used to read 
	data from an XML file or from an escript object directly. Writing XML 
	files are expensive, but this approach has the advantage given that the 
	results can be analyzed easily after the simulation has completed.   

	@attention: A DataCollector instance can only be used to specify one 
	scalar, vector and tensor attribute from a source at any one time. If a 
	second scalar, vector or tensor attribute needs to be specified from the 
	same source, a second DataCollector instance must be created. 
	"""

	def __init__(self, source = Source.XML):
		"""
		Initialise the data collector.

		@type source: L{Source <constant.Source>} constant
		@param source: Source type
		"""

		self.__source = source
		# Keeps track on whether DataCollector have been modified.
		self.__modified = True
		# Keeps track on whether any specific scalar, vector or tensor 
		# field have been specified.
		self.__set_scalar = False
		self.__set_vector= False
		self.__set_tensor= False

		if(source == Source.XML): # Source is an XML file.
			self.__vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
		# Source is a escript data object using a temp file in the background. 
		elif (self.__source == Source.ESCRIPT):
			self.__vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
			# Create a temporary .xml file and retrieve its path.
			self.__tmp_file = tempfile.mkstemp(suffix=".xml")[1]

	def __del__(self):
		"""
		Perform some clean up of the temporary file.
		"""

		if (self.__source == Source.ESCRIPT):
			if os.access(self.__tmp_file,os.F_OK): os.unlink(self.__tmp_file)

	def setFileName(self, file_name):
		"""
		Set the XML file name to read. 

		@type file_name: String
		@param file_name: Name of the file to read
		"""

		self.__modified = True

		if(self.__source == Source.XML):
			# Check whether the specified file exists, otherwise exit.
			if not(os.access(file_name, os.F_OK)):
				raise IOError("\nERROR: '%s' file does NOT exists.\n" % \
						file_name)

			self.__vtk_xml_reader.SetFileName(file_name)
			# Update must be called after SetFileName to make the reader 
			# up-to-date. Otherwise, some output values may be incorrect.
			self.__vtk_xml_reader.Update() 
			self.__get_attribute_lists()
		else:
			raise ValueError("Source type %s does not support \
			'setFileName'\n" % self.__source)

	def setData(self,**args):
		"""
		Create data using the <name>=<data> pairing. Assumption is made
		that the data will be given in the appropriate format.
		"""

		self.__modified = True

		if self.__source == Source.ESCRIPT:
			esys.escript.saveVTK(self.__tmp_file,**args)
			self.__vtk_xml_reader.SetFileName(self.__tmp_file)

			# Modified must be called for setData but NOT for
			# setFileName. If Modified is not called, only the first file 
			# will always be displayed. The reason Modified is used is 
			# because the same temporary file name is always used 
			# (previous file is overwritten). Modified MUST NOT be used in 
			# setFileName, it can cause incorrect output such as map.
			self.__vtk_xml_reader.Modified()

			# Update must be called after Modified. If Update is called before
			# Modified, then the first/second image(s) may not be updated
			# correctly.
			self.__vtk_xml_reader.Update()
			self.__get_attribute_lists()
		else:
			raise ValueError("Source type %s does not support 'setData'\n" \
					% self.__source)

	# This method is used to delay the execution of setting the active scalar
	# until 'setFileName' or 'setData' have been executed.
	def setActiveScalar(self, scalar):
		"""
		Specify the scalar field to load. 

		@type scalar: String
		@param scalar: Scalar field to load from the file. 
		"""

		self.__set_scalar = True
		self.__active_scalar = scalar

	def _setActiveScalar(self):
		"""
		Load the specified scalar field. 
		"""

		# Check whether the specified scalar is available in either point
		# or cell data. If not available, program exits.

		# NOTE: This check is similar to the check used in _getScalarRange 
		# but this is used only when a scalar attribute has been specified.
		if self.__active_scalar in self.__point_attribute['scalars']:
			self._getDataCollectorOutput().GetPointData().SetActiveScalars(
					self.__active_scalar)
		elif self.__active_scalar in self.__cell_attribute['scalars']:
			self._getDataCollectorOutput().GetCellData().SetActiveScalars(
					self.__active_scalar)
		else:
			raise IOError("ERROR: No scalar called '%s' is available." % \
					self.__active_scalar)

	# This method is used to delay the execution of setting the active vector
	# until 'setFileName' or 'setData' have been executed.
	def setActiveVector(self, vector):
		"""
		Specify the vector field to load.

		@type vector: String
		@param vector: Vector field to load from the file. 
		"""

		self.__set_vector = True
		self.__active_vector = vector

	def _setActiveVector(self):
		"""
		Load the specified vector field.
		"""

		# Check whether the specified vector is available in either point
		# or cell data. If not available, program exits.

		# NOTE: This check is similar to the check used in _getVectorRange 
		# but this is used only when a vector attribute has been specified.
		if self.__active_vector in self.__point_attribute['vectors']:
			self._getDataCollectorOutput().GetPointData().SetActiveVectors(
					self.__active_vector)
		elif self.__active_vector in self.__cell_attribute['vectors']:
			self._getDataCollectorOutput().GetCellData().SetActiveVectors(
					self.__active_vector)
		else:
			raise IOError("ERROR: No vector called '%s' is available." % \
					self.__active_vector)

	# This method is used to delay the execution of setting the active tensor
	# until 'setFileName' or 'setData' have been executed.
	def setActiveTensor(self, tensor):
		"""
		Specify the tensor field to load.

		@type tensor: String
		@param tensor: Tensor field to load from the file. 
		"""

		self.__set_tensor = True
		self.__active_tensor = tensor

	def _setActiveTensor(self):
		"""
		Load the the specified tensor field.
		"""

		# Check whether the specified tensor is available in either point
		# or cell data. If not available, program exits.

		# NOTE: This check is similar to the check used in _getTensorRange 
		# but this is used only when a tensor attribute has been specified.
		if self.__active_tensor in self.__point_attribute['tensors']:
			self._getDataCollectorOutput().GetPointData().SetActiveTensors(
					self.__active_tensor)
		elif self.__active_tensor in self.__cell_attribute['tensors']:
			self._getDataCollectorOutput().GetCellData().SetActiveTensors(
					self.__active_tensor)
		else:
			raise IOError("ERROR: No tensor called '%s' is available." % \
					self.__active_tensor)

	def __get_array_type(self, arr):
		"""
		Return whether an array type is scalar, vector or tensor by looking 
		at the number of components in the array.

		@type arr: vtkDataArray
		@param arr: An array from the source.
		@rtype: String
		@return: Array type ('scalar', vector' or 'tensor')
		"""

		# Number of components in an array.
		num_components = arr.GetNumberOfComponents() 

		if num_components == 1:
			return 'scalars'
		elif num_components == 3:
			return 'vectors'
		elif num_components == 9:
			return 'tensors'

	def __get_attribute_list(self, data): 
		"""
		Return the available scalar, vector and tensor attributes 
		(either point or cell data).

		@type data: vtkPointData or vtkCellData
		@param data: Available point data or cell data from the source
		@rtype: Dictionary
		@return: Dictionary containing the available scalar, vector and \
				tensor attributes
		"""

		attribute = {'scalars':[], 'vectors':[], 'tensors':[]} 
		if data:
			num_arrays = data.GetNumberOfArrays() # Number of arrays.
			for i in range(num_arrays):
				name = data.GetArrayName(i) # Get an array name.
				type = self.__get_array_type(data.GetArray(i)) # Get array type.
				attribute[type].extend([name]) # Add array name to dictionary.

		return attribute

	def __get_attribute_lists(self):
		"""
		Get all the available point and cell data attributes from the source.
		"""

		# Get all the available point data attributes into a list.
		self.__point_attribute = \
				self.__get_attribute_list(
				self._getDataCollectorOutput().GetPointData())

		# Get all the available cell data attribute into another list.	
		self.__cell_attribute = \
				self.__get_attribute_list(
				self._getDataCollectorOutput().GetCellData())

	def _getScalarRange(self):
		"""
		Return the scalar range.

		@rtype: Two column tuple containing numbers
		@return: Scalar range
		"""

		# Check whether any scalar is available in either point or cell data. 
		# If not available, program exits.

		# NOTE: This check is similar to the check used in _setActiveScalar 
		# but this is used only when no scalar attribute has been specified.
		if(len(self.__point_attribute['scalars']) != 0):
			return self._getDataCollectorOutput().GetPointData().\
					GetScalars().GetRange(-1)
		elif(len(self.__cell_attribute['scalars']) != 0):
			return self._getDataCollectorOutput().GetCellData().\
					GetScalars().GetRange(-1)
		else:
			raise IOError("\nERROR: No scalar is available.\n")	

	def _getVectorRange(self):
		"""
		Return the vector range.
	
		@rtype: Two column tuple containing numbers
		@return: Vector range
		"""
		
		# Check whether any vector is available in either point or cell data. 
		# If not available, program exits.

		# NOTE: This check is similar to the check used in _setActiveVector 
		# but this is used only when no vector attribute has been specified.

		# NOTE: Generally GetRange(-1) returns the correct vector range. 
		# However, there are certain data sets where GetRange(-1) seems 
		# to return incorrect mimimum vector although the maximum vector is 
		# correct. As a result, the mimimum vector has been hard coded to 0.0
		# to accommodate for the incorrect cases.
		if(len(self.__point_attribute['vectors']) != 0):
			vector_range = \
					self._getDataCollectorOutput().GetPointData().\
					GetVectors().GetRange(-1)
			return (0.0, vector_range[1])
		elif(len(self.__cell_attribute['vectors']) != 0): 
			vector_range = \
					self._getDataCollectorOutput().GetCellData().\
					GetVectors().GetRange(-1)
			return (0.0, vector_range[1])
		else:
			raise IOError("\nERROR: No vector is available.\n")	

	def _getTensorRange(self):
		"""
		Return the tensor range.

		@rtype: Two column tuple containing numbers
		@return: Tensor range
		"""

		# Check whether any tensor is available in either point or cell data. 
		# If not available, program exits.

		# NOTE: This check is similar to the check used in _setActiveTensor 
		# but this is used only when no tensor attribute has been specified.
		if(len(self.__point_attribute['tensors']) != 0):
			return self._getDataCollectorOutput().GetPointData().\
					GetTensors().GetRange(-1)
		elif(len(self.__cell_attribute['tensors']) != 0): 
			return self._getDataCollectorOutput().GetCellData().\
					GetTensors().GetRange(-1)
		else:
			raise IOError("\nERROR: No tensor is available.\n")	

	def _getDataCollectorOutput(self):
		"""
		Return the output of the data collector.

		@rtype: vtkUnstructuredGrid
		@return: Unstructured grid
		"""
		return self.__vtk_xml_reader.GetOutput()

	def _isModified(self):
		"""
		Return whether the DataCollector has been modified.

		@rtype: Boolean
		@return: True or False
		"""

		if(self.__modified == True):
			self.__modified = False
			return True
		else:
			return False
	
	def _isScalarSet(self):
		"""
		Return whether a specific scalar field has been specified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__set_scalar

	def _isVectorSet(self):
		"""
		Return whether a specific vector field has been specified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__set_vector 

	def _isTensorSet(self):
		"""
		Return whether a specific tensor field has been specified.

		@rtype: Boolean
		@return: True or False
		"""

		return self.__set_tensor 

	def _getCenter(self):
		"""
		Return the center of the rendered object.

		@rtype: Three column tuple containing numbers
		@return: Center of the rendered object
		"""

		return self._getDataCollectorOutput().GetCenter()
		
