"""
class that shows a vector field by arrows

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

class Property(object):
    def __init__(self,val=None):
       self.setValue(val)

    def setValue(self,val):
       self.__val=val
       self.__altered=True

    def getValue(self):
       return self.__val

    def isAltered(self):
       return self.__altered

    def markAsUsed(self):
       self.__altered=False
 
   
class Component(object):
    """
    shows a vector field by arrows
    """
    def __init__(self):
       self.features={} # item must be a Component or Property

    def render(self):
       for i in self.features:
          if isinstance(self.features[i],Component):
             self.features[i].render()
       self._render()

    def markFeaturesAsUsed(self):
       for i in self.features:
          if isinstance(self.features[i],Component):
             self.features[i].markAsUsed()
          else:
             self.features[i].markFeaturesAsUsed()
