"""
Example of loading and viewing an image using pyvisi

Will hopefully help me write a decent interface.

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


method = 'pyvisi'
format = 'pnm'

if method == 'pyvisi':
    ## this is the pyvisi code
    import sys
    
    # import the python visualisation interface
    from esys.pyvisi import *
    # this is now where the renderer is specified
    from esys.pyvisi.renderers.vtk import *
    
    # set up a scene
    scene = Scene()
    
    if format == 'jpeg':
        # add a jpeg image to the scene, and then load the file
        jpegImage = JpegImage(scene)
        jpegImage.load(fname="Flinders_eval.jpg")
        jpegImage.render()  # this should be done at the scene.render step
    
    elif format == 'png':
        # add png image to the scene, and then load the file
        pngImage = PngImage(scene)
        pngImage.load(fname="Flinders_eval.png")
        pngImage.render()

    elif format == 'bmp':
        # add bmp image to the scene, and then load the file
        bmpImage = BmpImage(scene)
        bmpImage.load(fname="Flinders_eval.bmp")
        bmpImage.render()

    elif format == 'tiff':
        # add tiff image to the scene, and then load the file
        tiffImage = TiffImage(scene)
        tiffImage.load(fname="Flinders_eval.tiff")
        tiffImage.render()

    elif format == 'pnm':
        # add pnm (ppm, pgm, pbm) image to the scene, and then load the file
        pnmImage = PnmImage(scene)
        pnmImage.load(fname="Flinders_eval.pnm")
        pnmImage.render()

    else:
        raise ValueError, "Unknown format: %s" % format

    # render the scene, pausing so that the opengl window doesn't disappear
    scene.render(pause=True,interactive=True)

elif method == 'vtk':
    ## this is the original vtk code
    
    import vtk 
    _ren = vtk.vtkRenderer()
    _renWin = vtk.vtkRenderWindow()
    _renWin.AddRenderer(_ren)
    
    _imgActor = vtk.vtkImageActor()

    if format == 'jpeg':
        _jpegReader = vtk.vtkJPEGReader()
        _jpegReader.SetFileName("Flinders_eval.jpg")
        _imgActor.SetInput(_jpegReader.GetOutput())

    elif format == 'png':
        _pngReader = vtk.vtkPNGReader()
        _pngReader.SetFileName("Flinders_eval.png")
        _imgActor.SetInput(_pngReader.GetOutput())

    elif format == 'bmp':
        _bmpReader = vtk.vtkBMPReader()
        _bmpReader.SetFileName("Flinders_eval.bmp")
        _imgActor.SetInput(_bmpReader.GetOutput())

    elif format == 'tiff':
        _tiffReader = vtk.vtkTIFFReader()
        _tiffReader.SetFileName("Flinders_eval.tiff")
        _imgActor.SetInput(_tiffReader.GetOutput())

    elif format == 'pnm':
        _pnmReader = vtk.vtkPNMReader()
        _pnmReader.SetFileName("Flinders_eval.pnm")
        _imgActor.SetInput(_pnmReader.GetOutput())

    else:
        raise ValueError, "Unknown format: %s" % format
    
    _ren.AddActor(_imgActor)
    _renWin.SetSize(400,400)
    _ren.SetBackground(0.1,0.2,0.4)
    _renWin.Render()
    raw_input("Press any key to continue")

else:
    print "Eeek!  What plotting method am I supposed to use???"
    
# vim: expandtab shiftwidth=4:
