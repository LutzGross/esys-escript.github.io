
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
__url__="https://launchpad.net/escript-finley"


import vtk
import sys, os, re
import getopt

(opts, args) = getopt.getopt(sys.argv[1:],
        "d:f:ph",
        ["dirname=", "fnamestem=", "preproc", "help"],
        )

def usage():
    print "Usage:"
    print "  velVis.py -d <dirname> -f <fnamestem> -p -i\n"
    print "  Arguments:"
    print "   -d/--dirname: directory of vtk files to process (req)"
    print "   -f/--fnamestem: filename stem of vtk files to process (req)"
    print "   -p/--preproc: preprocess vtk files to find global vector max?  (opt)"
    print "   -i/--interactive: interactive use of first frame? (opt)"

#dirName = '/data/raid2/matt/results/050516/m3dpData/vis'
dirName = None
fnameStem = None
interactive = False
preproc = False   # preprocess to find the global max of the data

for option, arg in opts:
    if option in ('-d', '--dirname'):
        dirName = arg
    elif option in ('-f', '--fnamestem'):
        fnameStem = arg
    elif option in ('-p', '--preproc'):
        preproc = True
    elif option in ('-i', '--interactive'):
        interactive = True
    elif option in ('-h', '--help'):
        usage()
        sys.exit(0)

if dirName is None:
    print "You must supply a directory of vtk files to process\n"
    usage()
    sys.exit(1)

if fnameStem is None:
    print "You must supply the filename stem of the files to process"
    print "The filename stem is the part of the filename without the"
    print "frame number or extension\n"
    usage()
    sys.exit(1)

# now need to determine the number of files to process
dirList = os.listdir(dirName)
r = re.compile( "%s\\d+\\.\\w+" % fnameStem )
count = 0  # counter for the number of files found
fnames = []
for fname in dirList:
    if r.match(fname):
        fnames.append(fname)
        count += 1

# determine the first frame number and number of frames
fnames.sort()
firstFile = fnames[0]

r = re.compile("([a-zA-Z-])(\\d+)(\\.)(\\w+)")
firstNum = r.findall(firstFile)
firstNum = int(firstNum[0][1])

# grab the filename extension as well, just in case
extn = r.findall(firstFile)
extn = extn[0][3]

# the number of frames should be equal to the count
if interactive:
    numFrames = 1
    mesa = False
else:
    numFrames = count
    mesa = True

# try and work out the format of the zero-padded frame number
if count < 10:
    formatString = '"%01d"'
elif count < 100:
    formatString = '"%02d"'
elif count < 1000:
    formatString = '"%03d"'
elif count < 10000:
    formatString = '"%04d"'
elif count < 100000:
    formatString = '"%05d"'
else:
    print "Sorry, can't handle greater than 100000 input frames"
    sys.exit(1)

# use mesa if required
if mesa:
    factGraphics = vtk.vtkGraphicsFactory()
    factGraphics.SetUseMesaClasses(1)

    factImage = vtk.vtkImagingFactory()
    factImage.SetUseMesaClasses(1)

# perform some preprocessing to find the global max of the vector data norm
maxNorm = 0
if not interactive and preproc:
    print "Preprocessing vector data"
    sys.stdout.write("Completed: ")
    for i in range(numFrames):
        frameNum = eval(formatString) % i
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName("%s/%s%s.%s" % (dirName,fnameStem,frameNum,extn))
        reader.Update()

        grid = reader.GetOutput()
        norm = grid.GetPointData().GetVectors().GetMaxNorm()
        if norm > maxNorm:
            maxNorm = norm

	sys.stdout.write("%3.0f%%\b\b\b\b" % (float(i)/float(numFrames), ))
    print ""

for i in range(firstNum, numFrames+firstNum):

    frameNum = eval(formatString) % i
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("%s/%s%s.%s" % (dirName,fnameStem,frameNum,extn))
    reader.Update()

    print "Processing %s" % fnames[i]

    grid = reader.GetOutput()

    # grab the model centre and bounds
    centre = grid.GetCenter()
    bounds = grid.GetBounds()

    # grab the norm of the vectors
    norm = vtk.vtkVectorNorm()
    norm.SetInput(grid)

    if not preproc:
        maxNorm = grid.GetPointData().GetVectors().GetMaxNorm()

    # to make arrow glyphs need an arrow source
    arrow = vtk.vtkArrowSource()
    
    # the arrows are 3D glyphs so set that up now
    glyph = vtk.vtkGlyph3D()
    glyph.ScalingOn()
    glyph.SetScaleModeToScaleByScalar()
    glyph.SetColorModeToColorByScalar()
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleFactor(0.1/maxNorm)
    glyph.SetInput(norm.GetOutput())
    glyph.SetSource(arrow.GetOutput())
    glyph.ClampingOff()

    # set up a stripper to speed up rendening
    stripper = vtk.vtkStripper()
    stripper.SetInput(glyph.GetOutput())
    
    # make a lookup table for the colour map and invert it (colours look
    # better when it's inverted)
    lut = vtk.vtkLookupTable()
    refLut = vtk.vtkLookupTable()
    lut.Build()
    refLut.Build()
    for j in range(256):
        lut.SetTableValue(j, refLut.GetTableValue(255-j))
    
    # set up the mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(stripper.GetOutput())
    mapper.SetScalarRange(0,maxNorm)
    mapper.SetLookupTable(lut)
    
    # set up the actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
   
    # set up the text properties for nice text
    font_size = 20
    textProp = vtk.vtkTextProperty()
    textProp.SetFontSize(font_size)
    textProp.SetFontFamilyToArial()
    textProp.BoldOff()
    textProp.ItalicOff()
    textProp.ShadowOff()
    textProp.SetColor(0.0, 0.0, 0.0)
 
    # make a title
    title = vtk.vtkTextMapper()
    title.SetInput("Velocity flow vectors in mantle convection")

    # make the title text use the text properties
    titleProp = title.GetTextProperty()
    titleProp.ShallowCopy(textProp)
    titleProp.SetJustificationToCentered()
    titleProp.SetVerticalJustificationToTop()
    titleProp.BoldOn()

    # make the actor for the title
    titleActor = vtk.vtkTextActor()
    titleActor.SetMapper(title)
    titleActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    titleActor.GetPositionCoordinate().SetValue(0.5, 0.95)

    # make a frame counter
    counter = vtk.vtkTextMapper()
    counter.SetInput("frame: %d" % i)

    # make the counter use the text properties
    counterProp = counter.GetTextProperty()
    counterProp.ShallowCopy(textProp)
    counterProp.SetJustificationToLeft()
    counterProp.SetVerticalJustificationToTop()
    counterProp.SetFontSize(14)

    # make the actor for the frame counter
    counterActor = vtk.vtkTextActor()
    counterActor.SetMapper(counter)
    counterActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    counterActor.GetPositionCoordinate().SetValue(0.05, 0.8)

    # use a scalar bar
    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lut)
    scalarBar.SetWidth(0.1)
    scalarBar.SetHeight(0.7)
    scalarBar.SetPosition(0.9, 0.2)
    scalarBar.SetTitle("|v|")

    # set up its title text properties
    scalarBarProp = scalarBar.GetTitleTextProperty()
    scalarBarProp.ShallowCopy(textProp)
    scalarBarProp.SetFontSize(10)

    # set up the label text properties 
    scalarBarTextProp = scalarBar.GetLabelTextProperty()
    scalarBarTextProp.ShallowCopy(textProp)
    scalarBarTextProp.SetFontSize(10)

    # put an outline around the data
    outline = vtk.vtkOutlineSource()
    outline.SetBounds(bounds)

    # make its mapper
    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInput(outline.GetOutput())

    # make its actor
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(0,0,0)

    # set up the renderer and render window
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    
    renWin.SetSize(800,600)
    renWin.AddRenderer(ren)
    ren.SetBackground(1,1,1)
    
    # add the relevant actors
    ren.AddActor(actor)
    ren.AddActor(titleActor)
    ren.AddActor(counterActor)
    ren.AddActor(scalarBar)
    ren.AddActor(outlineActor)
   
    cam = ren.GetActiveCamera()
    cam.Azimuth(0)
    cam.Elevation(-90)
    cam.Zoom(1.2)
    ren.SetActiveCamera(cam)
    ren.ResetCameraClippingRange()

    # add some axes
    axes = vtk.vtkCubeAxesActor2D()
    axes.SetInput(grid)
    axes.SetCamera(ren.GetActiveCamera())
    axes.SetLabelFormat("%6.4g")
    axes.SetFlyModeToOuterEdges()
    axes.SetFontFactor(0.8)
    axes.SetAxisTitleTextProperty(textProp)
    axes.SetAxisLabelTextProperty(textProp)
    axes.SetXLabel("x")
    axes.SetYLabel("y")
    axes.SetZLabel("z")
    axes.SetNumberOfLabels(5)
    axes.GetProperty().SetColor(0,0,0)
    ren.AddProp(axes)

    if interactive:
        # set up stuff for interactive viewing
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        
        iren.Initialize()
        renWin.Render()
        iren.Start()
    else:
        renWin.OffScreenRenderingOn()
        renWin.Render()

        # the WindowToImageFilter is what one uses to save the window to an 
        # image file
        win2img = vtk.vtkWindowToImageFilter()
        win2img.SetInput(renWin)
        
        # set up the PNMWriter as we're saving to pnm
        writer = vtk.vtkPNMWriter()
        writer.SetFileName("%s%04d.pnm" % (fnameStem,i))
        writer.SetInput(win2img.GetOutput())
        writer.Write()
    
        print "Wrote %s%04d.pnm" % (fnameStem,i)

# vim: expandtab shiftwidth=4:
