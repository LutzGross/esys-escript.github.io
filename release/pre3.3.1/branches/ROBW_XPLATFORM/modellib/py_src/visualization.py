# $Id$

from esys.escript.modelframe import Model
from esys.escript import saveVTK
import os

class Visualization(Model):
    """
    Generic visualization of scalar, vector and tensorial data 
    (not implemeted yet)

    @ivar t: current time
    @ivar n: frame counter
    @ivar scalar: scalar data set
    @ivar vector: vector data set
    @ivar tensor: tensor data set
    @ivar dt: increment for output 
    @ivar filename: name of the output file
    """

    def __init__(self, debug=False):
        """
        Initialisation of the visualisation model object

        @param debug: Debugging flag
        """
        super(Visualization,self).__init__(debug=debug)
        self.declareParameter(t=0.,
                              n=0,
                              scalar=None,
                              vector=None,
                              tensor=None,
                              dt=1,
                              filename="movie.mpg")

    def doInitialization(self):
        """
        Does some kind of initialisation
        """
        self.__last_t=self.t

    def writeFrame(self):
       out=self.t>=self.__last_t+self.dt
       if out: 
            self.__last_t+=self.dt
            self.n+=1
       return out

    def getFrameCounter(self):
        return self.n-1

    def getSafeTimeStepSize(self,dt):
           """
           returns new step size
                                                                                                                                                                                                    
           @param dt: last time step size used
           @type dt: C{float}
           @return: time step size that can savely be used
           @rtype: C{float}
           """
           return self.__last_t+self.dt-self.t

    def doStepPostprocessing(self, dt):
        """
        renders the scene

        @note: to be overwritten
        """
        if self.writeFrame():
            self.trace("%s-th frame at time %s"%(self.getFrameCounter(),self.t))
            if not self.scalar==None:
               self.trace("scalar data: (min,max) =(%s,%s)"%(inf(self.scalar),sup(self.scalar)))
            if not self.vector==None:
               self.trace("vector data: (min,max) =(%s,%s)"%(inf(self.vector),sup(self.vector)))
            if not self.tensor==None:
               self.trace("tensor data: (min,max) =(%s,%s)"%(inf(self.tensor),sup(self.tensor)))
           
    def doFinalization(self):
        """
        Finalises the visualisation.  For instance, makes a movie of the image files.

        @note: to be overwritten
        """
        pass

class WriteVTK(Visualization):
    """
    Writes data into VTK files for further processing. 
    """

    def __init__(self, debug=False):
        """
        Initialisation of the WriteVTK object

        @param debug: Debugging flag
        """
        super(WriteVTK,self).__init__(debug=debug)

    def doInitialization(self):
        """
        Does some kind of initialisation
        """
        super(WriteVTK,self).doInitialization()
        fnc=self.filename.split('.')
        if len(fnc)==0:
           self.__filename="data.%s.xml"
        else:
           n=fnc[0]
           for i in range(1,len(fnc)-1):
              n+="."+fnc[i]
           if len(fnc)==1:
              self.__filename=n+".%s"
           else:
              self.__filename=n+".%s."+fnc[-1]
        self.trace("output filename is %s."%self.__filename)
       

    def doStepPostprocessing(self, dt):
        """
        Do any necessary postprocessing operations after a timestep.

        @param dt:
        """
        if self.writeFrame():
            kwargs={}
            if not self.scalar==None: kwargs["scalar"] = self.scalar
            if not self.vector==None: kwargs["vector"] = self.vector
            if not self.tensor==None: kwargs["tensor"] = self.tensor
            saveVTK(self.__filename%self.getFrameCounter(),**kwargs)
            self.trace("%s-th frame at time %s is writen to %s"%(self.getFrameCounter(),self.t,self.__filename%self.getFrameCounter()))

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


# vim: expandtab shiftwidth=4:
