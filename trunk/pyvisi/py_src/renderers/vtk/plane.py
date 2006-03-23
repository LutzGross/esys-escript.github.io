"""
The classes associated with Planes

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
from common import debugMsg

# module specific imports
from item import Item

class Plane(Item):
    """
    Generic class for Plane objects
    """

    def __init__(self, scene):
        """
        Initialisation of the Plane object
        """
        debugMsg("Called Plane.__init__()")
        Item.__init__(self)

        self.renderer = scene.renderer

    def setOrigin(self, x, y, z):
        """
        Set the origin of the plane
        """
        self.origin = (x, y, z)
        return

    def getOrigin(self):
        """
        Get the current origin of the plane
        """
        return self.origin

    def setNormal(self, vx, vy, vz):
        """
        Set the normal vector to the plane
        """
        self.normal = (vx, vy, vz)
        return

    def getNormal(self):
        """
        Get the current normal vector to the plane
        """
        return self.normal

    def mapImageToPlane(self, image):
        # this really needs to go somewhere else!!!
        """
        Maps an Image object onto a Plane object
        """
        debugMsg("Called Plane.mapImageToPlane()")

        # need to work out the name of the internal image object name
        imgObjectName = image.readerName

        self.renderer.runString("# Plane.mapImageToPlane()")
        self.renderer.runString("_tex = vtk.vtkTexture()")
        evalString = "_tex.SetInput(%s.GetOutput())" % imgObjectName
        self.renderer.runString(evalString)
        self.renderer.runString("_plane = vtk.vtkPlaneSource()\n")
        self.renderer.runString(\
                "_planeMapper = vtk.vtkPolyDataMapper()\n")
        self.renderer.runString(\
                "_planeMapper.SetInput(_plane.GetOutput())\n")
        return

    def render(self):
        """
        Perform Plane object specific (pre)rendering tasks
        """
        debugMsg("Called Plane.mapImageToPlane()")

        self.renderer.runString("# Plane.render()\n")
        self.renderer.runString("_planeActor = vtk.vtkActor()\n")
        self.renderer.runString("_planeActor.SetMapper(_planeMapper)\n")
        self.renderer.runString("_planeActor.SetTexture(_tex)\n")
        self.renderer.runString("_renderer.AddActor(_planeActor)\n")

        return


class CutPlane(Plane):
    """
    Cut plane class: used to cut data sets with a plane

    Cut plane objects define a plane to cut a data set or plot by and return
    the data along the intersection between the data set or plot with the
    defined plane.
    """

    def __init__(self):
        """
        Intialisation of the CutPlane object
        """
        debugMsg("Called CutPlane.__init__()")
        Plane.__init__(self)


class ClipPlane(Plane):
    """
    Class for planes used to clip datasets
    """

    def __init__(self):
        """
        Intialisation of the ClipPlane object
        """
        debugMsg("Called ClipPlane.__init__()")
        Plane.__init__(self)

        # set the default inside out flag value
        self.insideOut = False

    def setInsideOut(self, insideOut):
        """
        Set the inside out flag
        """
        self.insideOut = insideOut
        return

    def getInsideOut(self):
        """
        Get the current value of the inside out flag
        """
        return self.insideOut

# vim: expandtab shiftwidth=4:
