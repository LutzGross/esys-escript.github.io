"""
defines a scene in which items are shown

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
__author__="Paul Cochrane, L. Gross"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision:$"
__date__="$Date:$"

from camera import IsometricCamera
from colors import White



def Scene(PyvisiObject):
   """
   a scene for visualization
   """
   def __init__(self,x_size,y_size,camera,background):
      self.setSize(x_size,y_size)
      self.setBackgroundColor(backgound)
      self.setCamera(backgound)

   def render(self,filename=None,wait=False):
      """
      renders the scene. If filename is present the file is rendered in the file filename
      """
      pass
   def setSize(self,x_size,y_size):
      """
      sets the size of the output window
      """
      pass

   def setBackgroundColor(self,color):
      """
      Sets the background color of the Scene

      @param color: The color to set the background to.  Can be RGB or CMYK
      @type color: L{Color}
      """
      pass

   def setCamera(self,camera):
      """
      sets the camera
`     """
      pass

def openScene(renderer="pyvisi",x_size=640,y_size=480,camera=IsometricCamera(),background=White()):
     """
     opens a scene 

     @param renderer: renderer to be used
     @param x_size: window width (in pixel)
     @param y_size: window height (in pixel)
     @param camera: camera 
     @param background: background color
     @returns: L{Scene} 
     """
     if renderer == "vtk_online":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"online", offscreen=False)
     elif renderer == "vtk_ps":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"ps", offscreen=True)
     elif renderer == "vtk_png":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"png", offscreen=True)
     elif renderer == "vtk_jpeg":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"jpeg", offscreen=True)
     elif renderer == "vtk_tiff":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"tiff", offscreen=True)
     elif renderer == "vtk_bmp":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"bmp", offscreen=True)
     elif renderer == "vtk_pnm":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,camera,background,"pnm", offscreen=True)
     elif renderer == "pyvisi" or renderer==None:
         import Scene
         return Scene(x_size,y_size,camera,background)
     else:
         raise ValueError, "Unknown renderer %s."%renderer

