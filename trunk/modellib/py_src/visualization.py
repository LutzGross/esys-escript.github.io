# $Id$

from esys.escript.modelframe import Model
from esys.escript import saveVTK
import os

class Visualization(Model):
    """
    Generic visualization of scalar, vector and tensorial data 
    (not implemeted yet)

    @ivar t: current time
    @ivar scalar: scalar data set
    @ivar vector: vector data set
    @ivar tensor: tensor data set
    @ivar stride: visulaization is done every strides time step
    @ivar filename: name of the movie file
    """

    def __init__(self, debug=False):
        """
        Initialisation of the visualisation model object

        @param debug: Debugging flag
        """
        Model.__init__(self, debug=debug)
        self.declareParameter(t=0.,
                scalar=None, vector=None, tensor=None,
                stride=1, movie="movie.mpg", counter=0)

    def doInitialization(self):
        """
        Does some kind of initialisation
        """
        self.__n = 0
        self.__scene = None

    def doStepPostprocessing(self, dt):
        """
        Does any necessary postprocessing after each step

        @param dt:
        """
        self.__n += 1
        if self.__n % self.stride:
            data = self.scalar
            if data != None:
                pass
            data = self.vector
            if data != None:
                pass
            data = self.tensor
            if data != None:
                pass

    def doFinalization(self):
        """
        Finalises the visualisation.  For instance, makes a movie of the
        image files.
        """
        # make the movie into self.filename
        pass

class ShadePlot(Visualization):
    """
    Shaded contour plots
    """
    
    def __init(self, debug=False):
        """
        Initialisation
        """
        Visualization.__init__(self, debug)
        self.declareParameter(filename="shadePlot.%s.png")

    def doStepPostprocessing(self, dt):
        """
        Do any necessary postprocessing operations after a timestep.

        @param dt:
        """
        self.counter += 1
        if self.counter % self.stride == 0:
            n = self.counter/self.stride

            # look for a vtk xml file of the right name and plot it
            dataFname = "data.s.%d.xml" % n
            if not os.path.exists(dataFname):
                print "Data file doesn't exist!  Skipping frame generation."

            else:
                import pyvisi
                import pyvisi.renderers.vtk
    
                scene = pyvisi.renderers.vtk.Scene()
                plot = pyvisi.renderers.vtk.ContourPlot(scene)
                plot.setData(fname=dataFname, 
                        format='vtk-xml', 
                        scalars='escript_scalar_data')
                scene.save(fname="shadePlot.%05d.png" % n, format="png")


class ArrowPlot(Visualization):
    """
    Arrow/vector/quiver plots
    """
    
    def __init(self, debug=False):
        """
        Initialisation
        """
        Visualization.__init__(self, debug)
        self.declareParameter(filename="arrowPlot.%s.png")

    def doStepPostprocessing(self, dt):
        """
        Do any necessary postprocessing operations after a timestep.

        @param dt:
        """
        self.counter += 1
        if self.counter % self.stride == 0:
            n = self.counter/self.stride

            # look for a vtk xml file of the right name and plot it
            dataFname = "data.v.%d.xml" % n
            if not os.path.exists(dataFname):
                print "Data file doesn't exist!  Skipping frame generation."

            else:
                import pyvisi
                import pyvisi.renderers.vtk
    
                scene = pyvisi.renderers.vtk.Scene()
                plot = pyvisi.renderers.vtk.ArrowPlot3D(scene)
                plot.setData(fname=dataFname, 
                        format='vtk-xml')
                scene.save(fname="arrowPlot.%05d.png" % n, format="png")



class EllipsoidPlot(Visualization):
    """
    Ellipsoid plots
    """
    
    def __init(self, debug=False):
        """
        Initialisation
        """
        Visualization.__init__(self, debug)


class WriteVTK(Visualization):
    """
    Writes data into a VTK file for further processing. Currently data 
    are written in several files for each data type. This may change 
    in the future.

    scalar: scalar data set
    vector: vector data set
    tensor: tensor data set
    stride: file is written every stride-th time step
    filename: name of the data files. use %s for indication of data type 
    (s,v,t) and time step id.
    """

    def __init__(self, debug=False):
        """
        Initialisation of the WriteVTK object

        @param debug: Debugging flag
        """
        Visualization.__init__(self, debug=debug)
        self.declareParameter(filename="data.%s.xml")

    def doStepPostprocessing(self, dt):
        """
        Do any necessary postprocessing operations after a timestep.

        @param dt:
        """
        self.counter += 1
        if self.counter % self.stride == 0:
            n = self.counter/self.stride
            data = self.scalar
            if hasattr(data, "saveVTK"): 
                data.saveVTK(self.filename % ("s.%d" % n))
            data = self.vector
            if hasattr(data, "saveVTK"): 
                data.saveVTK(self.filename % ("v.%d" % n))
            data = self.tensor
            if hasattr(data, "saveVTK"): 
                data.saveVTK(self.filename % ("t.%d" % n))

# vim: expandtab shiftwidth=4:
