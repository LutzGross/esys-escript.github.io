"""
@author: John NGUI
"""

from esys.pyvisi.constant import ImageFormat
from os import system

class Movie:
	"""
	Class that creates a file called 'make_movie' by default (if a parameter 
	file name is not speficied) which contains a list of parameters required 
	by the 'ppmtompeg' command to generate a movie from a series of images.
	"""
	
	def __init__(self, parameter_file = "make_movie"):	
		"""
		Initialise the generation of the movie.

		@type parameter_file: String
		@param parameter_file: Name of the file containing the list of 
				parameters required by the 'ppmtompeg' command.
		"""

		self.__parameter_file = parameter_file

	def makeMovie(self, input_directory, first_image, last_image, movie):
		"""
		Coordinate the generation of the movie.

		@type input_directory: String
		@param input_directory: Directory in which the series of images can 
		        be found
		@type first_image: String
		@param first_image: First image name (including the extension)
		@type last_image: String
		@param last_image: Last image name (including the extension)
		@type movie : String
		@param movie: Movie name (including the extension)
		"""

		self.__input_directory = input_directory
		self.__first_image = first_image
		self.__last_image = last_image
		self.__movie = movie

		self.__splitInput()
		self.__retrieveFirstImageDetails()
		self.__retrieveLastImageDetails()
		self.__retrieveConversionCommand()
		self.__generateParameterFile()

		# If a parameter file name was not specified, then the default file 
		# will be deleted automatically once the movie has been generated. 
		# However, if a paramter file name was specified, the file will be 
		# maintained.
		system('ppmtompeg ' + self.__parameter_file)
		if(self.__parameter_file == "make_movie"):
			system('rm ' + self.__parameter_file)

	def __splitInput(self):
		"""
		Split the image label (i.e. temp-0001.jpg) to separate the 
		image name (i.e. temp-0001) and image format (i.e. jpg).
		"""

		# Image format (i.e. jpg)
		self.__image_format = self.__first_image.split('.')[1]
		# First image name.
		self.__first_image = self.__first_image.split('.')[0]
		# Last image name.
		self.__last_image = self.__last_image.split('.')[0]

	def __retrieveFirstImageDetails(self):
		"""
		Retrieve the first image name details.
		"""

		# For a image called 'temp-0001.jpg', self.__first_image_number 
		# will be '0001', while self.__image_prefix will be 'temp-'.
		self.__first_image_number = ""
		self.__image_prefix = ""

		for i in range(len(self.__first_image)):
			if(self.__first_image[i].isdigit()): # Retrieve first image number.
				self.__first_image_number = \
						self.__first_image_number + self.__first_image[i]
			else: # Retrieve image prefix.
				self.__image_prefix = \
						self.__image_prefix + self.__first_image[i]
			
	def __retrieveLastImageDetails(self):
		"""
		Retrieve the last image name details.
		"""

		self.__last_image_number = ""

		for i in range(len(self.__last_image)):
			if(self.__last_image[i].isdigit()): # Retrieve last image number.
				self.__last_image_number = \
						self.__last_image_number + self.__last_image[i]

		# The last image number is deducted by one. This allows the user
		# to use the same image range when generating a movie as the 
		# ones used for the for-loop.
		#self.__last_image_number = unicode(int(self.__last_image_number) - 1)
		
	def __retrieveConversionCommand(self):
		"""
		Retrieve the conversion command (depending on the image format)
		required by the 'ppmtompeg' command.
		"""
		print "-----", self.__image_format
		if(self.__image_format.endswith(ImageFormat.JPG)):
			self.__command = 'jpeg'
		elif(self.__image_format.endswith(ImageFormat.BMP)):
			self.__command = ImageFormat.BMP
		elif(self.__image_format.endswith(ImageFormat.PNM)):
			self.__command = ImageFormat.PNM
		elif(self.__image_format.endswith(ImageFormat.PNG)):
			self.__command = ImageFormat.PNG
		elif(self.__image_format.endswith(ImageFormat.TIF)):
			self.__command = ImageFormat.TIF
		else:
			raise IOError("ERROR: Invalid image format.")


	def __generateParameterFile(self):
		"""
		Write the list of parameters into the file. 
		"""

		parameter_file = open(self.__parameter_file, 'w')

		parameter_file.write('PATTERN IBBPBBPBBPBBPBBP\n' +
			'OUTPUT ' + self.__movie + '\n' 
			'BASE_FILE_FORMAT PNM\n' + 
			'INPUT_CONVERT ' +  self.__command + 'topnm *\n' +
			'GOP_SIZE 16\n' +
			'SLICES_PER_FRAME 10\n' +
			'INPUT_DIR ' + self.__input_directory + '\n' +
			'INPUT\n' +
			self.__image_prefix + '*.' + self.__image_format + ' [' + \
					self.__first_image_number + '-' + \
					self.__last_image_number + ']\n'
			'END_INPUT\n' +
			'PIXEL HALF\n' +
			'RANGE 10\n' +
			'PSEARCH_ALG LOGARITHMIC\n' + 
			'BSEARCH_ALG CROSS2\n' +
			'IQSCALE 8\n' + 
			'PQSCALE 10\n' +
			'BQSCALE 25\n' + 
			'REFERENCE_FRAME DECODED\n' +
			'FORCE_ENCODE_LAST_FRAME\n' +
			'ASPECT_RATIO 1\n' +
			'FRAME_RATE 24')

		parameter_file.close()

