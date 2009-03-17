
########################################################
#
# Copyright (c) 2003-2008 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2008 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.uq.edu.au/esscc/escript-finley"

"""
@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="John Ngui, john.ngui@uq.edu.au"


from esys.pyvisi.constant import ImageFormat
from esys.escript import getMPIWorldMax, getMPIRankWorld
import os


class Movie:
    """
    Class that creates a file called 'make_movie' by default (if a parameter 
    file name is not speficied) which contains a list of parameters required 
    by the 'ppmtompeg' command to generate a movie from a series of images.

    @attention: A movie cannot be generated from postscript (.ps) images.
    """

    def __init__(self, parameter_file = "make_movie"):      
        """
        Initialise the generation of the movie. If a parameter file name
        is supplied, the file will not be deleted once the movie has been 
        generated. Otherwise, the temporary parameter file created will 
        automatically be removed at the end.

        @type parameter_file: String
        @param parameter_file: Name of the file containing the list of 
            parameters required by the 'ppmtompeg' command.
        """

        self.__parameter_file = parameter_file
        return


    def imageRange(self, input_directory, first_image, last_image):
        """
        The image range from which the movie is to be generated from.
        
        @type input_directory: String
        @param input_directory: Directory in which the series of images can 
        be found
        @type first_image: String
        @param first_image: First image name (including the extension)
        @type last_image: String
        @param last_image: Last image name (including the extension)
        """

        # Keeps track whether an image range or image list was provided as
        # the source for generating the movie.
        self.__image_range_used = True
        self.__input_directory = input_directory
        self.__first_image = first_image
        self.__last_image = last_image

        self.__splitInput()
        self.__retrieveFirstImageDetails()
        self.__retrieveLastImageDetails()
        self.__retrieveConversionCommand()
        return

    def imageList(self, input_directory, image_list):
        """
        The image list from which the movie is to be generated from.

        @type input_directory: String
        @param input_directory: Directory in which the series of images can 
        be found
        @type image_list: List
        @param image_list: List of images name (including the extension)
        """

        self.__image_range_used = False 
        self.__input_directory = input_directory
        self.__images = ""
        
        self.__first_image = image_list[0] # Get first image in the list.               
        self.__splitInput()
        self.__retrieveConversionCommand()
        
        for i in image_list:
            self.__images += i  + '\n'
        return

    def makeMovie(self, movie):
        """
        Generate the movie.

        @type movie : String
        @param movie: Movie name (including the .mpg extension)
        """

        self.__movie = movie
        self.__generateParameterFile()

        # If a parameter file name was not specified, then the default file 
        # will be deleted automatically once the movie has been generated. 
        # However, if a paramter file name was specified, the file will be 
        # maintained.
        if getMPIRankWorld():
            ret = os.system('ppmtompeg ' + self.__parameter_file) / 256
        else:
            ret=0
        ret=getMPIWorldMax(ret)
        if ret > 0:
          raise RuntimeError, "Could not  generate movie %s"%'ppmtompeg ' + self.__parameter_file
        
        if(self.__parameter_file == "make_movie"):
            if getMPIRankWorld():
            	ret = os.system('rm ' + self.__parameter_file) / 256
            else:
            	ret=0
            ret=getMPIWorldMax(ret)
            if ret > 0:
          	raise RuntimeError, 'Could not remove file '+ self.__parameter_file
        return

    def __splitInput(self):
        """
        Split the image label (i.e. temp-0001.jpg) to separate the 
        image name (i.e. temp-0001) and image format (i.e. jpg).
        """

        first_image_split = self.__first_image.split('.')
        # Image format (i.e. jpg)
        self.__image_format = first_image_split[len(first_image_split) - 1]
        # First image name.
        self.__first_image = \
                           self.__first_image.rstrip('.' + self.__image_format)

        if(self.__image_range_used == True):
            # Last image name.
            self.__last_image = \
                              self.__last_image.rstrip('.' + self.__image_format)
        return

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
        return

    def __retrieveLastImageDetails(self):
        """
        Retrieve the last image name details.
        """

        self.__last_image_number = ""

        for i in range(len(self.__last_image)):
            if(self.__last_image[i].isdigit()): # Retrieve last image number.
                self.__last_image_number = \
                                         self.__last_image_number + self.__last_image[i]

        return

    def __retrieveConversionCommand(self):
        """
        Retrieve the conversion command (depending on the image format)
        required by the 'ppmtompeg' command.
        """

        if(self.__image_format.endswith(ImageFormat.JPG)):
            self.__command = 'jpeg'
        elif(self.__image_format.endswith(ImageFormat.BMP)):
            self.__command = ImageFormat.BMP
        elif(self.__image_format.endswith(ImageFormat.PNM)):
            self.__command = ImageFormat.PNM
        elif(self.__image_format.endswith(ImageFormat.PNG)):
            self.__command = ImageFormat.PNG
        elif(self.__image_format.endswith(ImageFormat.TIF)):
            self.__command = 'tiff'
        else:
            raise IOError("ERROR: '" + self.__image_format + \
                          "' is an invalid image format.")

        return


    def __generateParameterFile(self):
        """
        Write the list of parameters into the file. 
        """

        if(self.__image_range_used == True): # Image range was provided.
            input = self.__image_prefix + '*.' + self.__image_format + ' [' + \
                    self.__first_image_number + '-' + \
                    self.__last_image_number + ']\n'
        else: # Image list was provided
            input = self.__images

        
        ret=0
        if getMPIRankWorld() == 0:
             try: 
                  parameter_file = open(self.__parameter_file, 'w')

                  if os.name == 'nt' :
                      tmp_input_name = self.__input_directory.replace('\\','/')
                      tmp_movie_name = self.__movie.replace('\\','/')
                  else:
                      tmp_input_name = self.__input_directory
                      tmp_movie_name = self.__movie
      
                  parameter_file.write('PATTERN IBBPBBPBBPBBPBBP\n' +
                                       'OUTPUT ' + tmp_movie_name + '\n' 
                                       'BASE_FILE_FORMAT PNM\n' + 
                                       'INPUT_CONVERT ' +  self.__command + 'topnm *\n' +
                                       'GOP_SIZE 16\n' +
                                       'SLICES_PER_FRAME 10\n' +
                                       'INPUT_DIR ' + tmp_input_name + '\n' +
                                       'INPUT\n' +
                                       input +
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
                  if os.name == 'nt' :
                      tmp_name = self.__parameter_file.replace('\\','/')
                      os.system('dos2unix ' + tmp_name) 
             except Exception, e:
                 ret=1
        ret=getMPIWorldMax(ret)
        if ret > 0:
            if getMPIRankWorld() == 0:
                raise e
            else:
                raise RuntimeError, "Movie generation on processor 0 failed."
        return
