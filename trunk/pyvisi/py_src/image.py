# Copyright (C) 2004-2005 Paul Cochrane
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# $Id: image.py,v 1.22 2005/03/07 07:45:48 paultcochrane Exp $

## @file image.py

"""
Class and functions associated with a pyvisi Image objects

Base implementation.
"""

# generic import
from common import debugMsg, overrideWarning, unsupportedError, fileCheck

from item import Item

__revision__ = '$Revision: 1.22 $'

class Image(Item):
    """
    Image class.  Generic class to handle image data.
    """
    def __init__(self, scene=None):
        """
        Initialises the Image class object
        
        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Item.__init__(self)
        debugMsg("Called Image.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

    def load(self, fname):
        """
        Loads image data from file.

        @param fname: The filename from which to load image data
        @type fname: string
        """
        debugMsg("Called Image.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("Image.load")

        return

class JpegImage(Image):
    """
    Subclass of Image class to explicitly handle jpeg images
    """
    def __init__(self, scene=None):
        """
        Initialises the JpegImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called JpegImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "jpeg"

    def load(self, fname):
        """
        Loads jpeg image data from file.

        @param fname: The filename from which to load jpeg image data
        @type fname: string
        """
        debugMsg("Called JpegImage.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("JpegImage.load")

        return

    def render(self):
        """
        Does JpegImage object specific (pre)rendering stuff
        """
        debugMsg("Called JpegImage.render()")

        # print a warning message if get to here
        overrideWarning("JpegImage.render")

        return

class PngImage(Image):
    """
    Subclass of Image class to explicitly handle png images
    """
    def __init__(self, scene=None):
        """
        Initialises the PngImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called PngImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "png"

    def load(self, fname):
        """
        Loads png image data from file.

        @param fname: The filename from which to load png image data
        @type fname: string
        """
        debugMsg("Called PngImage.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("PngImage.load")

        return

    def render(self):
        """
        Does PngImage object specific (pre)rendering stuff
        """
        debugMsg("Called PngImage.render()")

        # print a warning message if get to here
        overrideWarning("PngImage.render")

        return

class BmpImage(Image):
    """
    Subclass of Image class to explicitly handle bmp images
    """
    def __init__(self, scene=None):
        """
        Initialises the BmpImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called BmpImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "bmp"

    def load(self, fname):
        """
        Loads bmp image data from file.

        @param fname: The filename from which to load bmp image data
        @type fname: string
        """
        debugMsg("Called BmpImage.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("BmpImage.load")

        return

    def render(self):
        """
        Does BmpImage object specific (pre)rendering stuff
        """
        debugMsg("Called BmpImage.render()")

        # print a warning message if get to here
        overrideWarning("BmpImage.render")

        return

class TiffImage(Image):
    """
    Subclass of Image class to explicitly handle tiff images
    """
    def __init__(self, scene=None):
        """
        Initialises the TiffImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called TiffImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "tiff"

    def load(self, fname):
        """
        Loads tiff image data from file.

        @param fname: The filename from which to load tiff image data
        @type fname: string
        """
        debugMsg("Called TiffImage.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("TiffImage.load")

        return

    def render(self):
        """
        Does TiffImage object specific (pre)rendering stuff
        """
        debugMsg("Called TiffImage.render()")

        # print a warning message if get to here
        overrideWarning("TiffImage.render")

        return

class PnmImage(Image):
    """
    Subclass of Image class to explicitly handle pnm images
    """
    def __init__(self, scene=None):
        """
        Initialises the PnmImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called PnmImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "pnm"

    def load(self, fname):
        """
        Loads pnm image data from file.

        @param fname: The filename from which to load pnm image data
        @type fname: string
        """
        debugMsg("Called PnmImage.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("PnmImage.load")

        return

    def render(self):
        """
        Does PnmImage object specific (pre)rendering stuff
        """
        debugMsg("Called PnmImage.render()")

        # print a warning message if get to here
        overrideWarning("PnmImage.render")

        return

class PbmImage(Image):
    """
    Subclass of Image class to explicitly handle pbm images
    """
    def __init__(self, scene=None):
        """
        Initialises the PbmImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called PbmImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "pbm"

    def load(self, fname):
        """
        Loads pbm image data from file.

        @param fname: The filename from which to load pbm image data
        @type fname: string
        """
        debugMsg("Called PbmImage.load()")

        # do a check to see if the file exists
        fileCheck(fname)

        # print a warning message if get to here
        overrideWarning("PbmImage.load")

        return

    def render(self):
        """
        Does PbmImage object specific (pre)rendering stuff
        """
        debugMsg("Called PbmImage.render()")

        # print a warning message if get to here
        overrideWarning("PbmImage.render")

        return

class PsImage(Image):
    """
    Subclass of Image class to explicitly handle ps images
    """
    def __init__(self, scene=None):
        """
        Initialises the PsImage class object

        This object is B{only} used for generating postscript output

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called PsImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "ps"

    def load(self, fname):
        """
        Loads ps image data from file.

        B{NOT} supported by this renderer module

        @param fname: The filename from which to load ps image data
        @type fname: string
        """
        debugMsg("Called PsImage.load()")

        # this ability not handled by this renderer module
        unsupportedError("BASE")
        
        # do a check to see if the file exists
        fileCheck(fname)

        return

    def render(self):
        """
        Does PsImage object specific (pre)rendering stuff
        """
        debugMsg("Called PsImage.render()")

        return

class PdfImage(Image):
    """
    Subclass of Image class to explicitly handle pdf images
    """
    def __init__(self, scene=None):
        """
        Initialises the PdfImage class object

        This object is B{only} used for generating pdf output

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        Image.__init__(self)
        debugMsg("Called PdfImage.__init__()")

        if scene is not None:
            self.renderer = scene.renderer

        self.format = "pdf"

    def load(self, fname):
        """
        Loads pdf image data from file.

        B{NOT} supported by this renderer module

        @param fname: The filename from which to load pdf image data
        @type fname: string
        """
        debugMsg("Called PdfImage.load()")

        # need to check that the file exists and is readable etc here
        # *before* we add to the evalString, better to get the feedback
        # now rather than at the end of the script

        # this ability not handled by this renderer module
        unsupportedError("BASE")
        
        # do a check to see if the file exists
        fileCheck(fname)

        return

    def render(self):
        """
        Does PdfImage object specific (pre)rendering stuff
        """
        debugMsg("Called PdfImage.render()")

        return

# vim: expandtab shiftwidth=4:
