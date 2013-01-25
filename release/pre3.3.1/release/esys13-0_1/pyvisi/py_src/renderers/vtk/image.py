"""
Class and functions associated with a pyvisi Image objects

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


# generic imports
from common import debugMsg, unsupportedError

from esys.pyvisi.common import fileCheck

# module specific imports
from item import Item

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
        debugMsg("Called Image.__init__()")
        Item.__init__(self)

        if scene is not None:
            self.renderer = scene.renderer
        
    def load(self, fname):
        """
        Loads image data from file.

        @param fname: The filename from which to load image data
        @type fname: string
        """
        debugMsg("Called Image.load()")

        fileCheck(fname)

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
        debugMsg("Called JpegImage.__init__()")
        Image.__init__(self)

        if scene is not None:
            self.renderer = scene.renderer
            self.renderer.runString("# JpegImage.__init__()")
            self.renderer.runString("_jpegReader = vtk.vtkJPEGReader()")
            self.readerName = "_jpegReader"
            # add the image to the scene
            scene.add(self)

        self.format = "jpeg"

    def load(self, fname):
        """
        Loads jpeg image data from file.

        @param fname: The filename from which to load jpeg image data
        @type fname: string
        """
        debugMsg("Called JpegImage.load()")

        # need to check that the file exists and is readable etc here
        # *before* we add to the evalString, better to get the feedback
        # now rather than at the end of the script
        
        self.renderer.runString("# JpegImage.load()")
        evalString = "_jpegReader.SetFileName(\"%s\")" % fname
        self.renderer.runString(evalString)
        return

    def render(self):
        """
        Does JpegImage object specific (pre)rendering stuff
        """
        debugMsg("Called JpegImage.render()")

        self.renderer.runString("# JpegImage.render()")
        self.renderer.runString("_imgActor = vtk.vtkImageActor()")
        self.renderer.runString(\
                "_imgActor.SetInput(_jpegReader.GetOutput())")
        self.renderer.runString("_renderer.AddActor(_imgActor)")
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
        debugMsg("Called PngImage.__init__()")
        Image.__init__(self)

        if scene is not None:
            self.renderer = scene.renderer
            self.renderer.runString("# PngImage.__init__()")
            self.renderer.runString("_pngReader = vtk.vtkPNGReader()")
            self.readerName = "_pngReader"
            # add the image to the scene
            scene.add(self)

        self.format = "png"

    def load(self, fname):
        """
        Loads png image data from file.

        @param fname: The filename from which to load png image data
        @type fname: string
        """
        debugMsg("Called PngImage.load()")

        # check to see if the file exists
        fileCheck(fname)
        
        self.renderer.runString("# PngImage.load()")
        evalString = "_pngReader.SetFileName(\"%s\")" % fname
        self.renderer.runString(evalString)
        return

    def render(self):
        """
        Does PngImage object specific (pre)rendering stuff
        """
        debugMsg("Called PngImage.render()")

        self.renderer.runString("# PngImage.render()")
        self.renderer.runString("_imgActor = vtk.vtkImageActor()")
        self.renderer.runString(\
                "_imgActor.SetInput(_pngReader.GetOutput())")
        self.renderer.runString("_renderer.AddActor(_imgActor)")
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
        debugMsg("Called BmpImage.__init__()")
        Image.__init__(self)

        if scene is not None:
            self.renderer = scene.renderer
            self.renderer.runString("# BmpImage.__init__()")
            self.renderer.runString("_bmpReader = vtk.vtkBMPReader()")
            self.readerName = "_bmpReader"
            # add the image to the scene
            scene.add(self)

        self.format = "bmp"

    def load(self, fname):
        """
        Loads bmp image data from file.

        @param fname: The filename from which to load bmp image data
        @type fname: string
        """
        debugMsg("Called BmpImage.load()")

        # check to see if the file exists
        fileCheck(fname)
        
        self.renderer.runString("# BmpImage.load()")
        evalString = "_bmpReader.SetFileName(\"%s\")" % fname
        self.renderer.runString(evalString)
        return

    def render(self):
        """
        Does BmpImage object specific (pre)rendering stuff
        """
        debugMsg("Called BmpImage.render()")

        self.renderer.runString("# BmpImage.render()")
        self.renderer.runString("_imgActor = vtk.vtkImageActor()")
        self.renderer.runString(\
                "_imgActor.SetInput(_bmpReader.GetOutput())")
        self.renderer.runString("_renderer.AddActor(_imgActor)")
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
        debugMsg("Called TiffImage.__init__()")
        Image.__init__(self)

        if scene is not None:
            self.renderer = scene.renderer
            self.renderer.runString("# TiffImage.__init__()")
            self.renderer.runString("_tiffReader = vtk.vtkTIFFReader()")
            self.readerName = "_tiffReader"
            # add the image to the scene
            scene.add(self)

        self.format = "tiff"

    def load(self, fname):
        """
        Loads tiff image data from file.

        @param fname: The filename from which to load tiff image data
        @type fname: string
        """
        debugMsg("Called TiffImage.load()")

        # check to see if the file exists
        fileCheck(fname)

        self.renderer.runString("# TiffImage.load()")
        evalString = "_tiffReader.SetFileName(\"%s\")" % fname
        self.renderer.runString(evalString)
        return

    def render(self):
        """
        Does TiffImage object specific (pre)rendering stuff
        """
        debugMsg("Called TiffImage.render()")

        self.renderer.runString("# TiffImage.render()")
        self.renderer.runString("_imgActor = vtk.vtkImageActor()")
        self.renderer.runString(\
                "_imgActor.SetInput(_tiffReader.GetOutput())")
        self.renderer.runString("_renderer.AddActor(_imgActor)")
        return

class PnmImage(Image):
    """
    Subclass of Image class to explicitly handle pnm (ppm, pgm, pbm) images
    """
    def __init__(self, scene=None):
        """
        Initialises the PnmImage class object

        @param scene: The Scene object to add to
        @type scene: Scene object
        """
        debugMsg("Called PnmImage.__init__()")
        Image.__init__(self)

        if scene is not None:
            self.renderer = scene.renderer
            self.renderer.runString("# PnmImage.__init__()")
            self.renderer.runString("_pnmReader = vtk.vtkPNMReader()")
            self.readerName = "_pnmReader"
            # add the image to the scene
            scene.add(self)

        self.format = "pnm"
        
    def load(self, fname):
        """
        Loads pnm (ppm, pgm, pbm) image data from file.

        @param fname: The filename from which to load pnm image data
        @type fname: string
        """
        debugMsg("Called PnmImage.load()")

        # check to see if the file exists
        fileCheck(fname)

        self.renderer.runString("# PnmImage.load()")
        evalString = "_pnmReader.SetFileName(\"%s\")" % fname
        self.renderer.runString(evalString)
        return

    def render(self):
        """
        Does PnmImage object specific (pre)rendering stuff
        """
        debugMsg("Called PnmImage.render()")

        self.renderer.runString("# PnmImage.render()")
        self.renderer.runString("_imgActor = vtk.vtkImageActor()")
        self.renderer.runString(\
                "_imgActor.SetInput(_pnmReader.GetOutput())")
        self.renderer.runString("_renderer.AddActor(_imgActor)")
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
        debugMsg("Called PsImage.__init__()")
        Image.__init__(self)

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

        # need to check if the file exists
        fileCheck(fname)

        # this ability not handled by this renderer module
        unsupportedError()
        
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
        debugMsg("Called PdfImage.__init__()")
        Image.__init__(self)

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

        # need to check that the file exists
        fileCheck(fname)

        # this ability not handled by this renderer module
        unsupportedError()
        
        return

    def render(self):
        """
        Does PdfImage object specific (pre)rendering stuff
        """
        debugMsg("Called PdfImage.render()")

        return

# vim: expandtab shiftwidth=4:
