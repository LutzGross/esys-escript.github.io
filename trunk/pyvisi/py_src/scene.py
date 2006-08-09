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



def Scene(Component):
   """
   a scene for visualization
   """
   def __init__(self,x_size,y_size)
      super(Scene, self).__init__()
      self.x_size=x_size
      self.y_size=y_size
      self.setBackgroundColor()
      self.setCamera()

   def render(self,filename=None,wait=False):
      """
      renders the scene. If filename is present the file is rendered in the file filename
      """
      pass

   def _addComponent(self,comp):
       pass

   def setBackgroundColor(self,color=White()):
      """
      Sets the background color of the Scene

      @param color: The color to set the background to.  Can be RGB or CMYK
      @type color: L{Color}
      """
      pass

def openScene(renderer="pyvisi",x_size=640,y_size=480):
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
         return VTKScene(x_size,y_size,"online", offscreen=False)
     elif renderer == "vtk_ps":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,"ps", offscreen=True)
     elif renderer == "vtk_png":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,"png", offscreen=True)
     elif renderer == "vtk_jpeg":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,"jpeg", offscreen=True)
     elif renderer == "vtk_tiff":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,"tiff", offscreen=True)
     elif renderer == "vtk_bmp":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,"bmp", offscreen=True)
     elif renderer == "vtk_pnm":
         import renderers.vtk.VTKScene as VTKScene
         return VTKScene(x_size,y_size,"pnm", offscreen=True)
     elif renderer == "pyvisi" or renderer==None:
         import Scene
         return Scene(x_size,y_size)
     else:
         raise ValueError, "Unknown renderer %s."%renderer
##########################################################################################
## vtk specific stuff (should be in renderers.vtk)
##########################################################################################
def renderArrow(comp):
      pass # 
class VTKScene(Scene):
    def __init__(x_size,y_size,format,offscreen):
         super(VTKScene, self).__init__()
         self.__comps=[]
         pass
    def _addComponent(self,comp):
        if isinstance(Arrow,comp):
           # ... some vtk set up here"
           comp._render=renderArrow
        self.__comps.append(comp)
 
   def render(self,filename=None,wait=False):
      """
      renders the scene. If filename is present the file is rendered in the file filename
      """
      for i in self.__comps: i.render()
      pass # more rendering work
      for i in self.__comps: i.markFeaturesAsUsed()
  
