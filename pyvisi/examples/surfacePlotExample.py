# $Id: surfacePlotExample.py,v 1.6 2005/05/05 00:46:44 paultcochrane Exp $

"""
Example of plotting surfaces with pyvisi 
"""

# what plotting method are we using?
method = 'pyvisi'

# set up some data to plot
from Numeric import *

# the x and y axes
x = arange(-2,2,0.2, typecode=Float)
y = arange(-2,3,0.2, typecode=Float)

# pick some interesting function to generate the data in the third dimension
# this is the one used in the matlab docs: z = x*exp(-x^2-y^2)
z = zeros((len(x),len(y)), typecode=Float)

# boy do *I* feel old fashioned writing it this way
# surely there's another way to do it: - something to do later
for i in range(len(x)):
    for j in range(len(y)):
	z[i,j] = x[i]*exp(-x[i]*x[i] - y[j]*y[j])

# plot it with either gnuplot, vtk or pyvisi
if method == 'pyvisi':
    #### pyvisi version of code

    # import the general pyvisi stuff
    from esys.pyvisi import *
    # import the gnuplot overrides of the interface
    #from esys.pyvisi.renderers.gnuplot import *
    from esys.pyvisi.renderers.plplot import *

    # define a scene object
    # a Scene is a container for all of the kinds of things you want to put
    # into your plot, for instance, images, meshes, arrow/vector/quiver
    # plots, contour plots, spheres etc.
    scene = Scene()

    # create a SurfacePlot object
    plot = SurfacePlot(scene)

    # add some helpful info to the plot
    plot.title = 'Example surface plot'
    plot.xlabel = 'x'
    plot.ylabel = 'y'
    plot.zlabel = 'z'

    # assign the data to the plot
    # this version assumes that we have x, then y, then z and that z is 2D
    # and that x and y are 1D arrays
    plot.setData(x,y,z)
    # alternative syntax
    #plot.setData(xData=x, yData=y, zData=z)
    # or (but more confusing depending upon one's naming conventions)
    #plot.setData(x=x, y=y, z=z)

    # render the scene to screen
    scene.render(pause=True, interactive=True)

    # save the scene to file
    plot.setData(x,y,z)  # have to do this now because we've already
                         # render()ed the scene.  This requirement will be
                         # removed in the future
    scene.save(fname="surfacePlot.png", format=PngImage())

elif method == 'vtk':
    print "vtk surface plotting not yet implemented"

elif method == 'plplot':
    import plplot

    # determine the min and max of x, y, and z in world coordinates
    xMin = min(x)
    xMax = max(x)
    
    yMin = min(y)
    yMax = max(y)

    zMin = min(z.flat)
    zMax = max(z.flat)

    # min and max of x and y variables in normalised coordinates
    # (these are values recommended by plplot in an example)
    xMin2D = -2.5
    xMax2D = 2.5

    yMin2D = -2.5
    yMax2D = 4.0

    # sides of box in normalised coordinates
    # (these are values recommended by plplot in an example)
    basex = 2.0
    basey = 4.0
    height = 3.0

    # angle to view box
    alt = 45.0
    az = 30.0

    side = 1
    opt = 3  # plots a net of lines

    plplot.plsdev("xwin")
    plplot.plinit()
    plplot.plenv(xMin2D, xMax2D, yMin2D, yMax2D, 0, -2)
    plplot.plw3d(basex, basey, height, 
            xMin, xMax, yMin, yMax, zMin, zMax, 
            alt, az)
    plplot.plmtex("t", 1.0, 0.5, 0.5, "Example surface plot")
    plplot.plbox3("bnstu", "x axis", 0.0, 0, 
            "bnstu", "y axis", 0.0, 0, 
            "bcdmnstuv", "z axis", 0.0, 0)
    plplot.plsurf3d(x, y, z, 0, ())
    plplot.plend()
    
    # to save as well, have to set everything up again, and replot
    # save as png
    plplot.plsdev("png")
    plplot.plsfnam("surfacePlot.png")
    plplot.plinit()
    plplot.plenv(xMin2D, xMax2D, yMin2D, yMax2D, 0, -2)
    plplot.plw3d(basex, basey, height, 
            xMin, xMax, yMin, yMax, zMin, zMax, 
            alt, az)
    plplot.plmtex("t", 1.0, 0.5, 0.5, "Example surface plot")
    plplot.plbox3("bnstu", "x axis", 0.0, 0, 
            "bnstu", "y axis", 0.0, 0, 
            "bcdmnstuv", "z axis", 0.0, 0)
    plplot.plsurf3d(x, y, z, 0, ())
    plplot.plend()

else:
    print "Eeek!  What plotting method am I supposed to use???"

# vim: expandtab shiftwidth=4:
