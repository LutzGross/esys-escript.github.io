# $Id: clipBox.py,v 1.1 2005/09/26 08:50:29 paultcochrane Exp $

"""
Example of plotting spheres with pyvisi and clipping them with a box

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


# what plotting method are we using?
#method = 'vtk'
method = 'pyvisi'

import os, sys
import random

# set up some data to plot
from Numeric import *

# the three axes in space
# this will give us 10 particles (_not_ 1000)
x = arange(10, typecode=Float)
y = arange(10, typecode=Float)
z = arange(10, typecode=Float)

# 3D position information
posArray = []
for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            posArray.append( (x[i], y[j], z[k]) )

# radius information
random.seed()
radiiArray = zeros(len(x)*len(y)*len(z), typecode=Float)
for i in range(len(x)*len(y)*len(z)):
    radiiArray[i] = random.random()*0.8

# tag information
random.seed()
tagsArray = zeros(len(x)*len(y)*len(z), typecode=Int)
for i in range(len(x)*len(y)*len(z)):
    tagsArray[i] = int(random.random()*10)

# plot it using one of the three methods
if method == 'pyvisi':

    # example code for how a user would write a script in pyvisi
    from esys.pyvisi import *          # base level visualisation stuff
    # import the objects to render the scene using the specific renderer
    from esys.pyvisi.renderers.vtk import *       # vtk
    #from esys.pyvisi.renderers.povray import *       # povray
    
    # define the scene object
    scene = Scene()
    
    # create a BallPlot object
    plot = BallPlot(scene)
    
    # add some helpful info to the plot
    plot.title = 'Example ball plot, clipped by a box'

    # assign some data to the plot
    plot.setData(fname="cp_test_0.xml", 
            format="vtk-xml", 
            radii="radius", 
            tags="particleTag")

    # clip the plot by a box defined by its bounds
    clipBox = ClipBox(plot)
    clipBox.setBounds(-1, 10, 10, 25, 2, 15)

    # render the scene
    scene.render(pause=True, interactive=True)

    # save the scene out to file
    ## as png
    plot.setData(fname="cp_test_0.xml", 
            format="vtk-xml", 
            radii="radius", 
            tags="particleTag")
    scene.save(fname="clipBox.png", format="png")

elif method == 'povray':
    ### povray code

    raise NotImplementedError, "Sorry, povray version not implemented yet."

    # load the data from the vtk file (yes, I know this is a bit dodgy)
    import vtk

    # create the reader of the file
    _reader = vtk.vtkXMLUnstructuredGridReader()
    _reader.SetFileName("cp_test_0.xml")
    #_reader.SetFileName("/home/cochrane/raid2/vis4people/steffen/frame_0.xml")
    _reader.Update()

    # read the output into an unstructured grid
    _grid = _reader.GetOutput()

    _modelCentre = _grid.GetCenter()
    _xMin, _xMax, _yMin, _yMax, _zMin, _zMax = _grid.GetBounds()

    # grab the points where the data sit
    _vtkPoints = _grid.GetPoints()
    
    # grab the radius data for the radii of the balls
    _vtkRadii = _grid.GetPointData().GetScalars("radius")

    # grab the tag data and use for colouring the balls
    _vtkTags = _grid.GetPointData().GetScalars("particleTag")

    # work out dynamically the number of different tags so that can use this
    # information to automatically set the scalar range for colouring
    _numPoints = _vtkTags.GetNumberOfTuples()
    _valueDict = {}
    for i in range(_numPoints):
        _tagValue = _vtkTags.GetValue(i)
        _valueDict[_tagValue] = 1

    _numTags = len(_valueDict.keys())

    _tagValues = _valueDict.keys()
    _tagValues.sort()

    # count the number of tags, and make an evenly spaced array of points
    # between zero and one, then use these as the scalars to colour by
    _vtkScaledTags = vtk.vtkFloatArray()
    _vtkScaledTags.SetNumberOfTuples(_numPoints)
    _vtkScaledTags.SetNumberOfComponents(1)
    _vtkScaledTags.SetName("tags")
    for i in range(_numPoints):
        _tagValue = _vtkTags.GetValue(i)
        for j in range(_numTags):
            if _tagValues[j] == _tagValue:
                _vtkScaledTags.InsertTuple1(i, float(j)/float(_numTags-1))

    # use vtk to generate the colour map, will have to do this myself at
    # some point
    _lut = vtk.vtkLookupTable()
    _lut.Build()

    _red = zeros(_numPoints, typecode=Float)
    _green = zeros(_numPoints, typecode=Float)
    _blue = zeros(_numPoints, typecode=Float)
    for i in range(_numPoints):
        _red[i], _green[i], _blue[i] = _lut.GetColor(_vtkScaledTags.GetValue(i))

    # now convert the information we want (radii, colours, positions) into
    # array objects so that I can play with them as per normal in python

    # note:  this is an inefficient way to do this, I can do it in one loop,
    # but this way makes the meaning of the code a lot clearer

    ### the points
    _xData = zeros(_numPoints, typecode=Float)
    _yData = zeros(_numPoints, typecode=Float)
    _zData = zeros(_numPoints, typecode=Float)
    for i in range(_numPoints):
        _xData[i], _yData[i], _zData[i] = _vtkPoints.GetPoint(i)

    ### the radii
    _radii = zeros(_numPoints, typecode=Float)
    for i in range(_numPoints):
        _radii[i] = _vtkRadii.GetValue(i)

    ### the tags
    _scaledTags = zeros(_numPoints, typecode=Float)
    _tags = zeros(_numPoints, typecode=Int)
    for i in range(_numPoints):
        _scaledTags[i] = _vtkScaledTags.GetValue(i)
        _tags[i] = _vtkTags.GetValue(i)

    ### generate the pov file

    # open the pov file to write to
    pov = open("ballPlot.pov", "w")

    # the include files to add
    pov.write("#include \"shapes.inc\"\n")
    pov.write("#include \"colors.inc\"\n")

    # the camera
    pov.write("camera {\n")
    pov.write("  location <%f, %f, -100>\n" % 
            (_modelCentre[0], _modelCentre[1]))
    pov.write("  direction <0, 0, 2>\n")
    pov.write("  up <0, 1, 0>\n")
    pov.write("  right <4/3, 0, 0>\n")
    pov.write("  look_at <%f, %f, %f>\n" % 
            (_modelCentre[0], _modelCentre[1], _modelCentre[2]))
    pov.write("}\n")

    # the light source
    pov.write("light_source {\n")
    pov.write("  <%f, %f, -300>\n" % (_modelCentre[0], _modelCentre[1]))
    pov.write("  colour White\n")
    pov.write("}\n")

    # the spheres
    for i in range(_numPoints):
        pov.write("sphere { \n")
        pov.write("  <%f, %f, %f>, %f\n" %
                (_xData[i], _yData[i], _zData[i], _radii[i]))
        pov.write("  pigment {\n")
        if _tags[i] != 20:
            pov.write("    colour red %f green %f blue %f\n" %
                    (_red[i], _green[i], _blue[i]))
        else:
            pov.write("    rgbt <%f, %f, %f, 0.90>\n" %
                    (_red[i], _green[i], _blue[i]))
        pov.write("  }\n")
        pov.write("}\n")

    # put a title on it
    pov.write("text {\n")
    pov.write("  ttf \"timrom.ttf\" \"Example ball plot\" 0.1, 0\n")
    pov.write("  pigment {\n")
    pov.write("    colour White\n")
    pov.write("  }\n")
    pov.write("  scale <3, 3, 1>\n")
    pov.write("  translate <%f, %f, 0>\n" % 
            (_modelCentre[0]-10, 1.2*_yMax))
    pov.write("}\n")

    # close the file
    pov.close()

    ### generate the ini file

    # open the ini file to write to
    ini = open("ballPlot.ini", "w")

    # the output resolution
    ini.write("Width=640\n")
    ini.write("Height=480\n")

    # anti-aliasing settings
    ini.write("Antialias=on\n")

    # generate png files
    ini.write("Output_File_Type=N\n")

    # the name of the input pov file
    ini.write("Input_File_Name=ballPlot.pov\n")

    # pause when done
    ini.write("Pause_When_Done=on\n")

    # close the file
    ini.close()

    # run povray on the file
    os.system("povray ballPlot.ini")
    
elif method == "vtk":
    import vtk

        # we're loading a file; make sure it exists first
    if os.path.exists("cp_test_0.xml"):
        pass
    else:
        raise IOError, "File not found"

    # create the reader of the file
    _reader = vtk.vtkXMLUnstructuredGridReader()
    _reader.SetFileName("cp_test_0.xml")
    _reader.Update()

    # read the output into an unstructured grid
    _grid = _reader.GetOutput()

    # grab the radius data for the radii of the balls
    _radii = _grid.GetPointData().GetScalars("radius")

    # grab the tag data and use for colouring the balls
    _tags = _grid.GetPointData().GetScalars("particleTag")

    # work out dynamically the number of different tags so that can use this
    # information to automatically set the scalar range for colouring
    _numPoints = _tags.GetNumberOfTuples()
    _valueDict = {}
    for i in range(_numPoints):
        _tagValue = _tags.GetValue(i)
        _valueDict[_tagValue] = 1

    _numTags = len(_valueDict.keys())

    _tagValues = _valueDict.keys()
    _tagValues.sort()

    # count the number of tags, and make an evenly spaced array of points
    # between zero and one, then use these as the scalars to colour by
    _scaledTags = vtk.vtkFloatArray()
    _scaledTags.SetNumberOfTuples(_numPoints)
    _scaledTags.SetNumberOfComponents(1)
    _scaledTags.SetName("tags")
    for i in range(_numPoints):
        _tagValue = _tags.GetValue(i)
        for j in range(_numTags):
            if _tagValues[j] == _tagValue:
                _scaledTags.InsertTuple1(i, float(j)/float(_numTags-1))

    # in vtk 4.2 have to set up an array of two components to get
    # the data through the glyph object to the mapper so create this now
    _data = vtk.vtkFloatArray()
    _data.SetNumberOfComponents(3)
    _data.SetNumberOfTuples(_radii.GetNumberOfTuples())
    _data.CopyComponent(0, _radii, 0)
    _data.CopyComponent(1, _tags, 0)
    _data.CopyComponent(2, _scaledTags, 0)
    _data.SetName("data")

    # add the data array to the grid (again, only specific to vtk 4.2 and 4.4)
    _grid.GetPointData().AddArray(_data)

    # make the data the active scalars
    _grid.GetPointData().SetActiveScalars("data")

    # to make sphere glyphs need a sphere source
    _sphere = vtk.vtkSphereSource()
    _sphere.SetRadius(1.0)  # set to 0.5 by default in vtk, we need it to be 1.0
    _sphere.SetThetaResolution(5)  # chunky, but will be good in large models
    _sphere.SetPhiResolution(5)

    # the spheres are 3D glyphs so set that up now
    _glyph = vtk.vtkGlyph3D()
    _glyph.ScalingOn()
    _glyph.SetScaleModeToScaleByScalar() # scale by the radii scalars
    _glyph.SetColorModeToColorByScalar() # colour by the tag scalars
    _glyph.SetScaleFactor(1.0)  # just in case
    _glyph.SetInput(_grid)
    _glyph.SetSource(_sphere.GetOutput())  # grab the 3D glyph to use
    _glyph.ClampingOff()  # so we even see the tiny spheres

    # set up a stripper (this will speed up rendering a lot)
    _stripper = vtk.vtkStripper()
    _stripper.SetInput(_glyph.GetOutput())

    # set up the planes for the box to clip by
    _planes = vtk.vtkPlanes()
    _planes.SetBounds(-1, 10, 10, 25, 2, 15)

    # now make the clipper
    _clipper = vtk.vtkClipPolyData()
    _clipper.SetClipFunction(_planes)
    _clipper.SetInput(_stripper.GetOutput())
    _clipper.InsideOutOff()

    # set up the mapper
    _mapper = vtk.vtkPolyDataMapper()
    _mapper.SetInput(_clipper.GetOutput())
    _mapper.ScalarVisibilityOn()
    _mapper.ColorByArrayComponent("data", 2)
    _mapper.SetScalarRange(0, 1)  

    # set up the actor
    _actor = vtk.vtkActor()
    _actor.SetMapper(_mapper)

    # text properties
    _font_size = 14
    _textProp = vtk.vtkTextProperty()
    _textProp.SetFontSize(_font_size)
    _textProp.SetFontFamilyToArial()
    _textProp.BoldOff()
    _textProp.ItalicOff()
    _textProp.ShadowOff()

    # add a title
    _titleMapper = vtk.vtkTextMapper()
    _title = "Example ball plot, clipped by a box"
    _titleMapper.SetInput(_title)

    _titleProp = _titleMapper.GetTextProperty()
    _titleProp.ShallowCopy(_textProp)
    _titleProp.SetJustificationToCentered()
    _titleProp.SetVerticalJustificationToTop()
    _titleProp.SetColor(0,0,0)
    _titleProp.BoldOn()
    _titleProp.SetFontSize(20)

    # set up the text actor
    _titleActor = vtk.vtkTextActor()
    _titleActor.SetMapper(_titleMapper)
    _titleActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    _titleActor.GetPositionCoordinate().SetValue(0.5, 0.95)

    # set up the renderer and the render window
    _ren = vtk.vtkRenderer()
    _renWin = vtk.vtkRenderWindow()
    _renWin.SetSize(640, 480)
    _renWin.AddRenderer(_ren)
    _ren.SetBackground(1,1,1)

    # add the actor
    _ren.AddActor(_actor)
    _ren.AddActor(_titleActor)

    # render the scene
    _iren = vtk.vtkRenderWindowInteractor()
    _iren.SetRenderWindow(_renWin)
    _iren.Initialize()
    _renWin.Render()
    _iren.Start()

    # convert the render window to an image
    _win2imgFilter = vtk.vtkWindowToImageFilter()
    _win2imgFilter.SetInput(_renWin)

    # save the image to file
    _outWriter = vtk.vtkPNGWriter()
    _outWriter.SetInput(_win2imgFilter.GetOutput())
    _outWriter.SetFileName("clipBox.png")
    _outWriter.Write()

else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:

