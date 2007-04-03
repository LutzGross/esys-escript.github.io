"""
@author: John NGUI
"""

import vtk
import tempfile, os, sys
from constant import Source, VizType, ColorMode
try:
	import esys.escript 
except ImportError:
	print "Warning: importing esys.escript failed."

class DataCollector:
	"""
	Class that defines a data collector which deals with the source 
	of data for the visualisation.

	@attention: One DataCollector instance can only be used to specify one 
	scalar, vector and tensor attribute from a source at any one time. If a 
	second scalar, vector or tensor attribute needs to be specified from the 
	same source, a second DataCollector instance must be created. 

	@attention: When a series of XML files / ESCRIPT objects are read 
	(using 'setFileName' in a for-loop), the 'setActiveScalar' / 
	'setActiveVector' / 'setActiveTensor' have to be called after loading each 
	new file (if a specific field needs to be loaded) as all active fields 
	specified for the previous file goes back to the default once a new file 
	is read.
	"""

	def __init__(self, source = Source.XML):
		"""
		Initialise the data collector.

		@type source: L{Source <constant.Source>} constant
		@param source: Source type
		"""

		self.__source = source
		self.__count = 0 # Keeps track of the number of files/sources read.

		if(source == Source.XML): # Source is an XML file.
			self.__vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
		# Source is a escript data object using a temp file in the background. 
		elif (self.__source == Source.ESCRIPT):
			self.__vtk_xml_reader = vtk.vtkXMLUnstructuredGridReader()
			# Create a temporary .xml file and retrieves its path.
			self.__tmp_file = tempfile.mkstemp(suffix=".xml")[1]
			#self.__vtk_xml_reader.SetFileName(self.__tmp_file)

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

		if(self.__source == Source.XML):
			self.__vtk_xml_reader.SetFileName(file_name)

			# Update must be called after SetFileName to make the reader 
			# up to date. Otherwise, some output values may be incorrect.
			self.__vtk_xml_reader.Update() 
			self.__output = self.__vtk_xml_reader.GetOutput()
			self.__get_attribute_lists()
			
			# Count has to be larger than zero because when setFileName is 
			# called for the first time, the data set mapper has not yet been
			# instantiated. Therefore, the range of the mapper can only be
			# updated after the first file/source has been read.
			if(self.__count > 0):
				self._updateRange(None)

			self.__count+=1

		else:
			raise ValueError("Source type %s does not support 'setFileName'\n" \
					% self.__source)

	def setData(self,**args):
		"""
		Create data using the <name>=<data> pairing. Assumption is made
		that the data will be given in the appropriate format.
		"""

		if self.__source == Source.ESCRIPT:
			esys.escript.saveVTK(self.__tmp_file,**args)
			self.__vtk_xml_reader.SetFileName(self.__tmp_file)
			# Modified must be called for setData but NOT for
			# setFileName. If Modified is not called, only the first file 
			# will always be displayed. The reason Modified is used is 
			# because the same temporary file name is always used 
			# (old file is overwritten). Modified MUST NOT be used in 
			# setFileName, it can cause incorrect output such as map.
			self.__vtk_xml_reader.Modified()
			# Update must be called after Modified. If Update is called before
			# Modified, then the first/second image(s) may not be updated
			# correctly.
			self.__vtk_xml_reader.Update()
			self.__output = self.__vtk_xml_reader.GetOutput()
			self.__get_attribute_lists()
		else:
			raise ValueError("Source type %s does not support 'setData'\n" \
					% self.__source)

	def setActiveScalar(self, scalar):
		"""
		Specify the scalar field to load. 

		@type scalar: String
		@param scalar: Scalar field to load from the file. 
		"""
		
		# Check whether the specified scalar is available in either point
		# or cell data. If not available, program exits.

		# NOTE: This check is similar to the check used in _getScalarRange 
		# but is used only when a scalar attribute has been specified.
		if scalar in self.__point_attribute['scalars']:
			self._getOutput().GetPointData().SetActiveScalars(scalar)
		elif scalar in self.__cell_attribute['scalars']:
			self._getOutput().GetCellData().SetActiveScalars(scalar)
		else:
			print "\nERROR: No scalar called '%s' is available.\n" % scalar
			sys.exit(1)	

	def setActiveVector(self, vector):
		"""
		Specify the vector field to load.

		@type vector: String
		@param vector: Vector field to load from the file. 
		"""
		
		# Check whether the specified vector is available in either point
		# or cell data. If not available, program exits.

		# NOTE: This check is similar to the check used in _getVectorRange 
		# but is used only when a vector attribute has been specified.
		if vector in self.__point_attribute['vectors']:
			self._getOutput().GetPointData().SetActiveVectors(vector)
		elif vector in self.__cell_attribute['vectors']:
			self._getOutput().GetCellData().SetActiveVectors(vector)
		else:
			print "\nERROR: No vector called '%s' is available.\n" % vector
			sys.exit(1)	
			
	def setActiveTensor(self, tensor):
		"""
		Specify the tensor field to load.

		@type tensor: String
		@param tensor: Tensor field to load from the file. 
		"""

		# Check whether the specified tensor is available in either point
		# or cell data. If not available, program exits.

		# NOTE: This check is similar to the check used in _getTensorRange 
		# but is used only when a tensor attribute has been specified.
		if tensor in self.__point_attribute['tensors']:
			self._getOutput().GetPointData().SetActiveTensors(tensor)
		elif tensor in self.__cell_attribute['tensors']:
			self._getOutput().GetCellData().SetActiveTensors(tensor)
		else:
			print "\nERROR: No tensor called '%s' is available.\n" % tensor
			sys.exit(0)	

	# 'object' is set to 'None' because some types of visualization have
	# two different ranges that need to be updated while others only have one.
	def _paramForUpdatingMultipleSources(self, viz_type, color_mode, mapper,
			object = None):
		"""
		Parameters required to update the necessary range when two or more 
		files are read.

		@type viz_type: : L{VizType <constant.VizType>} constant
		@param viz_type: Type if visualization (i.e. Map, Velocity, etc)
		@type color_mode: L{ColorMode <constant.ColorMode>} constant
		@param color_mode: Type of color mode
		@type mapper: vtkDataSetMapper
		@param mapper: Mapped data
		@type object: vtkPolyDataAlgorith (i.e. vtkContourFilter, vtkGlyph3D, \
				etc)
		@param object: Poly data
		"""

		self.__viz_type = viz_type
		self.__color_mode = color_mode
		self.__mapper = mapper
		self.__object = object

	def _updateRange(self, range):
		"""
		Update the necessary range when two or more sources are read in.
		"""

		if self.__viz_type == VizType.MAP or \
				self.__viz_type == VizType.ELLIPSOID or \
				self.__viz_type == VizType.CARPET:
			self.__mapper.SetScalarRange(self._getScalarRange())
		elif self.__viz_type == VizType.VELOCITY:
			if self.__color_mode == ColorMode.VECTOR:
				self.__object.SetRange(self._getVectorRange())
				self.__mapper.SetScalarRange(self._getVectorRange())
			elif self.__color_mode == ColorMode.SCALAR:
				self.__object.SetRange(self._getScalarRange())
				self.__mapper.SetScalarRange(self._getScalarRange())
		elif self.__viz_type == VizType.CONTOUR:
			self.__object.GenerateValues(
					self.__object.GetNumberOfContours(),
					self._getScalarRange()[0],
					self._getScalarRange()[1])
			self.__mapper.SetScalarRange(self._getScalarRange())
		elif self.__viz_type == VizType.STREAMLINE:
			if self.__color_mode == ColorMode.VECTOR:
				self.__mapper.SetScalarRange(self._getVectorRange())
			elif self.__color_mode == ColorMode.SCALAR:
				self.__mapper.SetScalarRange(self._getScalarRange())

	def __get_array_type(self, arr):
		"""
		Return if the array type is a scalar, vector or tensor by looking 
		at the number of components in the array.

		@type arr: vtkDataArray
		@param arr: An array from the source.
		@rtype: String
		@return: Array type ('scalar', vector', 'tensor')
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
				self.__get_attribute_list(self._getOutput().GetPointData())
		# Get all the available cell data attribute into another list.	
		self.__cell_attribute = \
				self.__get_attribute_list(self._getOutput().GetCellData())

	def _getScalarRange(self):
		"""
		Return the scalar range.

		@rtype: Two column tuple containing numbers
		@return: Scalar range
		"""

		# Check whether any tensor is available in either point or cell data. 
		# If not available, program exits.

		# NOTE: This check is similar to the check used in _setActiveScalar 
		# but is used only when no scalar attribute has been specified.
		if(len(self.__point_attribute['scalars']) != 0):
			return self._getOutput().GetPointData().GetScalars().GetRange(-1)
		elif(len(self.__cell_attribute['scalars']) != 0):
			return self._getOutput().GetCellData().GetScalars().GetRange(-1)
		else:
			print "\nERROR: No scalar is available.\n"	
			sys.exit(1)

	def _getVectorRange(self):
		"""
		Return the vector range.
	
		@rtype: Two column tuple containing numbers
		@return: Vector range
		"""
		
		# Check whether any vector is available in either point or cell data. 
		# If not available, program exits.

		# NOTE: This check is similar to the check used in _setActiveVector 
		# but is used only when no vector attribute has been specified.

		# NOTE: Generally GetRange(-1) returns the correct vector range. 
		# However, there are certain data sets where GetRange(-1) seems 
		# to return incorrect mimimum vector although the maximum vector is 
		# correct. As a result, the mimimum vector has been hard coded to 0.0
		# to accommodate the incorrect cases.
		if(len(self.__point_attribute['vectors']) != 0):
			vector_range = \
					self._getOutput().GetPointData().GetVectors().GetRange(-1)
			return (0.0, vector_range[1])
		elif(len(self.__cell_attribute['vectors']) != 0): 
			vector_range = \
					self._getOutput().GetCellData().GetVectors().GetRange(-1)
			return (0.0, vector_range[1])
		else:
			print "\nERROR: No vector is available.\n"	
			sys.exit(0)

	def _getTensorRange(self):
		"""
		Return the tensor range.

		@rtype: Two column table
		@return: Tensor range
		"""

		# Check whether any tensor is available in either point or cell data. 
		# If not available, program exits.

		# NOTE: This check is similar to the check used in _setActiveTensor 
		# but is used only when no tensor attribute has been specified.
		if(len(self.__point_attribute['tensors']) != 0):
			return self._getOutput().GetPointData().GetTensors().GetRange(-1)
		elif(len(self.__cell_attribute['tensors']) != 0): 
			return self._getOutput().GetCellData().GetTensors().GetRange(-1)
		else:
			print "\nERROR: No tensor is available.\n"	
			sys.exit(1)

	def _getOutput(self):
		"""
		Return the output of the data collector.

		@rtype: vtkUnstructuredGrid
		@return: Unstructured grid
		"""

		return self.__output

