
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


class Color:
	"""
	Constants that define the colors using RGB values.

	:cvar RED: Constant representing red color
	:cvar GREEN: Constant representing green color
	:cvar BLUE: Constant representing blue color
	:cvar BLACK: Constant representing black color
	:cvar WHITE: Constant representing white color
	:cvar YELLOW: Constant representing yellow color
	:cvar PINK: Constant represnting pink color
	:cvar ORANGE: Constant representing orange color
	:cvar PURPLE: Constant representing purple color
	:cvar GREY: Constant representing grey color
	:cvar BROWN: Constant representing brown color
	"""

	RED    = [1, 0, 0]
	GREEN  = [0, 1, 0]
	BLUE   = [0, 0, 1]
	BLACK  = [0, 0, 0]
	WHITE  = [1, 1, 1]
	YELLOW = [1,1,0]
	PINK   = [1, 0.0784, 0.4588]
	ORANGE = [1, 0.2706, 0]
	PURPLE = [0.5412, 0.1680, 0.8828]
	GREY   = [0.6602, 0.6602, 0.6602]
	BROWN  = [0.5430, 0.2700, 0.0742]

class Viewport:
	"""
	Constants that define the four viewports in a window.

	:cvar SOUTH_WEST: Constant representing the bottom-left viewport of a 
			window
	:cvar NORTH_WEST: Constant representing the upper-left viewport of a 
			window
	:cvar NORTH_EAST: Constatnt representing the upper-right viewport of a 
			window
	:cvar SOUTH_EAST: Constant representing the bottom-right viewport of a 				window
	"""

	SOUTH_WEST = 0
	NORTH_WEST = 1
	NORTH_EAST = 2
	SOUTH_EAST = 3

class Source:
	"""
	Constants that define the data source type.

	:cvar XML: Constant representing xml as the source type
	:cvar ESCRIPT: Constant representing escript data objects the source type
	"""
	
	XML     = "xml"
	ESCRIPT = "escript"

class Renderer:
	"""
	Constants that define the renderer type.

	:cvar ONLINE: Constant representing the online renderer
	:cvar ONLINE_JPG: Constant representing the JPG online renderer
	:cvar ONLINE_BMP: Constant representing the BMP online renderer
	:cvar ONLINE_PNM: Constant representing the PNM online renderer
	:cvar ONLINE_PNG: Constant representing the PNG online renderer
	:cvar ONLINE_TIF: Constant representing the TIF online renderer
	:cvar ONLINE_PS: Constant representing the PS online renderer
	:cvar ONLINE_VRML: Constant representing the VRML online renderer
	:cvar ONLINE_IV: Constant representing the OpenInventor online renderer

	:cvar OFFLINE_JPG: Constant representing the JPG offline renderer
	:cvar OFFLINE_BMP: Constant representing the BMP offline renderer
	:cvar OFFLINE_PNM: Constant representing the PNM offline renderer
	:cvar OFFLINE_PNG: Constant representing the PNG offline renderer
	:cvar OFFLINE_TIF: Constant representing the TIF offline renderer
	:cvar OFFLINE_PS: Constant representing the PS offline renderer
	:cvar OFFLINE_VRML: Constant representing the VRML offline renderer
	:cvar OFFLINE_IV: Constant representing the OpenInventor offline renderer

	:cvar DISPLAY: Constant representing the display renderer
	:cvar DISPLAY_JPG: Constant representing the JPG display renderer
	:cvar DISPLAY_BMP: Constant representing the BMP display renderer
	:cvar DISPLAY_PNM: Constant representing the PNM display renderer
	:cvar DISPLAY_PNG: Constant representing the PNG display renderer
	:cvar DISPLAY_TIF: Constant representing the TIF display renderer
	:cvar DISPLAY_PS: Constant representing the PS display renderer
	:cvar DISPLAY_VRML: Constant representing the VRML display renderer
	:cvar DISPLAY_IV: Constant representing the OpenInventor display renderer
	"""

	ONLINE       = "online"
	ONLINE_JPG   = "online_jpg"
	ONLINE_BMP   = "online_bmp"
	ONLINE_PNM   = "online_pnm"
	ONLINE_PNG   = "online_png"
	ONLINE_TIF   = "online_tif"
	ONLINE_PS    = "online_ps"
	ONLINE_VRML  = "online_vrml"
	ONLINE_IV    = "online_iv"

	OFFLINE_JPG   = "offline_jpg"
	OFFLINE_BMP   = "offline_bmp"
	OFFLINE_PNM   = "offline_pnm"
	OFFLINE_PNG   = "offline_png"
	OFFLINE_TIF   = "offline_tif"
	OFFLINE_PS    = "offline_ps"
	OFFLINE_VRML  = "offline_vrml"
	OFFLINE_IV    = "offline_iv"

	DISPLAY       = "display"
	DISPLAY_JPG   = "display_jpg"
	DISPLAY_BMP   = "display_bmp"
	DISPLAY_PNM   = "display_pnm"
	DISPLAY_PNG   = "display_png"
	DISPLAY_TIF   = "display_tif"
	DISPLAY_PS    = "display_ps"
	DISPLAY_VRML  = "display_vrml"
	DISPLAY_IV    = "display_iv"

class Arrow:
	"""
	Constants that define the arrow type.

	:cvar TWO_D: Constant representing the two dimensional arrow type
	:cvar THREE_D: Constant representing the three dimensional arrow type
	"""

	TWO_D   = "2d"
	THREE_D = "3d"

class ColorMode:
	"""
	Constants that define the color mode used to color the data. 

	:cvar VECTOR: Constant representing the vector color mode 
	:cvar SCALAR: Constant representing the scalar color mode
	"""

	VECTOR = "vector"
	SCALAR = "scalar"

class WarpMode:
	"""
	Constants that define the warp mode used to deform the scalar data.

	:cvar VECTOR: Constant representing the vector deformation mode 
	:cvar SCALAR: Constant representing the scalar deformation mode
	"""

	VECTOR = "vector"
	SCALAR = "scalar"

class ImageFormat:
	"""
	Constants that define the image formats.

	:cvar JPG: Constant representing the JPG image format (.jpg)
	:cvar BMP: Constant representing the BMP image format (.bmp)
	:cvar PNM: Constant representing the PNM image format (.pnm)
	:cvar PNG: Constant representing the PNG image format (.png)
	:cvar TIF: Constant representing the TIF image format (.tif)
	:cvar PS:  Constant representing the PS image format  (.ps)
	"""

	JPG = "jpg"
	BMP = "bmp"
	PNM = "pnm"
	PNG = "png"
	TIF = "tif"
	PS  = "ps"

class Lut:
	"""
	Constants that define the type of mapping scheme for the lookup table.

	:cvar COLOR: Constant representing the color scheme
	:cvar GREY_SCALE: Constant representing the grey scale scheme
	"""
	
	COLOR      = "color"
	GREY_SCALE = "grey_scale"

class LegendType:
	"""
	Constants that define the type of legend for the scalar bar.

	:cvar SCALAR: Constant representing the legend using scalar data
	:cvar VECTOR: Constant representing the legend using  vector data
	"""

	SCALAR = "scalar"
	VECTOR = "vector"
