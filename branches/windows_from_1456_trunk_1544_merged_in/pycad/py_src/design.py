#
# $Id$
#
#######################################################
#
#           Copyright 2003-2007 by ACceSS MNRF
#       Copyright 2007 by University of Queensland
#
#                http://esscc.uq.edu.au
#        Primary Business: Queensland, Australia
#  Licensed under the Open Software License version 3.0
#     http://www.opensource.org/licenses/osl-3.0.php
#
#######################################################
#

"""
template for the Design which defines a regions and features
for a mesh generator.

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""


__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2007 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"


from primitives import Primitive, ReversePrimitive, PropertySet, Point, Manifold1D, Manifold2D, Manifold3D
from xml.dom import minidom

class TagMap(object):
    """
    a class that allows to map tags to names
   
    tm=TagMap({5 : x })
    tm.setMap(a=1,x=4)
    assert tm.getTags("a") == [ 1 ]
    assert tm.getTags("x") == [ 5, 4 ]
    assert tm.map(x=10., a=20.) == { 5 : 10, 4: 10, 1 : 20 }

    """
    def __init__(self, mapping={}):
      """
      initizlizes the mapping. map defines an initial mapping from tag to a name.
      """
      self.__mapping={}
      for tag, name in mapping.items():
          if not isinstance(tag, int):
              raise TypeError("tag needs to be int")
          if not isinstance(name, str):
              raise TypeError("name needs to be a str.")
          self.__mapping[tag]=name
    def setMap(self,**kwargs):
      """
      set a new map where <name>=<tag> assigns the tag <tag> to name <name>. <tag> has to be integer.
      If <tag> has been assigned to a name before the mapping will be overwritten. Otherwise a new 
      mapping <tag> -> <name> is set. Notice that a single name can be assigned to different tags.
      """
      for  name, tag in kwargs.items():
          if not isinstance(tag, int):
             raise TypeError("tag needs to be int")
          self.__mapping[tag]=name
    def getTags(self,name=None):
        """
        returns a list of the tags assigned to name. If name is not present a list of tags is returned.
        """
        if name == None:
           out=self.__mapping.keys()
        else:
           out=[]
           for tag, arg in self.__mapping.items():
             if arg == name: out.append(tag)
        return out
    def getName(self,tag=None):
        """
        returns the name of a tag. If tag is not present a list of names is returned.
        """
        if tag == None:
           return list(set(self.__mapping.values()))
        else:
            return self.__mapping[tag]

    def getMapping(self):
        """
        returns a dictionary where the tags define the keys and the values the corresposnding names.
        """
        return self.__mapping

    def map(self,default=0,**kwargs):
        """
        returns a dictionary where the tags define the keys and the values give the values assigned to the tag via name 
        and kwargs:

        tm=TagMap(x=5)
        tm.setMap(a=1,x=4)
        print tm.map(x=10., a=20.) 
        { 5 : 10, 4: 10, 1 : 20 }
   
        the default is used for tags which map onto name with unspecified values
        """
        out={}
        for tag in self.__mapping:
           if kwargs.has_key(self.__mapping[tag]):
              out[tag]=kwargs[self.__mapping[tag]]
           else:
              out[tag]=default
        return out

    def insert(self,data,default=0,**kwargs):
        """
        inserts values into the L{esys.escript.Data} object according to the given values assigned to the keywords.
        the default is used for tags which map onto name with unspecified values
        """
        d=self.map(default=default,**kwargs)
        for t,v in d.items():
             data.setTaggedValue(t,v)
    def passToDomain(self,domain):
        """
        passes the tag map to  L{esys.escript.Domain} domain.
        """
        for tag, name in self.__mapping.items():
          domain.setTagMap(name,tag)
         
    def toDOM(self,dom):
         """
         adds object to dom
         """
         tm=dom.createElement("TagMap")
         dom.appendChild(tm)
         for tag,name in self.getMapping().items():
             item_dom=dom.createElement("map")
             tag_dom=dom.createElement("tag")
             name_dom=dom.createElement("name")
             tag_dom.appendChild(dom.createTextNode(str(tag)))
             name_dom.appendChild(dom.createTextNode(str(name)))
             item_dom.appendChild(tag_dom)
             item_dom.appendChild(name_dom)
             tm.appendChild(item_dom)
         return tm
    def fromDom(self,node):
        """
        fills from dom node
        """
        for node in node.childNodes:
           if isinstance(node, minidom.Element):
             if node.tagName == 'map':
               tag=int(node.getElementsByTagName("tag")[0].firstChild.nodeValue.strip())
               name=str(node.getElementsByTagName("name")[0].firstChild.nodeValue.strip())
               self.setMap(**{ name : tag })
        return

    def fillFromXML(self,iostream):
       """
       uses the xml file or string to set the mapping
       """
       if isinstance(iostream,str):
             dom=minidom.parseString(iostream)
       else:
           dom=minidom.parse(iostream)
       root=dom.getElementsByTagName('ESys')[0]
       for node in root.childNodes:
           if isinstance(node, minidom.Element):
              if node.tagName == 'TagMap':
                 self.fromDom(node)
                 return
             
          

    def writeXML(self,iostream=None):
         """
         writes XML serialization into the iostream or if not present returns the XML as string
         """
         dom=minidom.Document()
         esys=dom.createElement('ESys')
         esys.appendChild(self.toDOM(dom))
         dom.appendChild(esys)
         if iostream == None:
            return dom.toprettyxml()
         else:
            iostream.write(dom.toprettyxml())
        
class Design(object):
    """
    template for a design which defines the input for a mesh generator
    """
    def __init__(self,dim=3,element_size=1.,order=1,keep_files=False):
       """
       initializes a design 

       @param dim: patial dimension
       @param element_size: global element size
       @param order: element order
       @param keep_files: flag to keep work files.
       """ 
       self.clearItems()
       self.setElementSize(element_size)
       self.setDim(dim)
       self.setElementOrder(order)
       if keep_files:
          self.setKeepFilesOn()
       else:
          self.setKeepFilesOff()
    def setDim(self,dim=3):
        """
        sets the spatial dimension
        """
        if not dim  in [1,2,3]: 
           raise ValueError("only dimension 1, 2, 3 are supported.")
        self.__dim=dim
    def getDim(self,dim=3):
        """
        returns the spatial dimension
        """
        return self.__dim
    def setElementOrder(self,order=1):
        """
        sets the element order
        """
        if not order in [1,2]:
           raise ValueError("only element orser 1 or 2 is supported.")
        self.__order=order
        
    def getElementOrder(self):
        """
        returns the element order
        """
        return self.__order
        
    def setElementSize(self,element_size=1.):
        """
        set the global element size.
        """
        if element_size<=0.:
           raise ValueError("element size needs to be positive.")
        self.__element_size=element_size
        
    def getElementSize(self):
        """
        returns the global element size.
        """
        return self.__element_size
        
    def setKeepFilesOn(self):
        """
        work files are kept at the end of the generation
        """
        self.__keep_files=True
    def setKeepFilesOff(self):
        """
        work files are deleted at the end of the generation
        """
        self.__keep_files=False
    def keepFiles(self):
        """
        returns True if work files are kept
        """
        return self.__keep_files
    def addItems(self,*items):
       """
       adds items to the design
       """
       for i in range(len(items)):
          if not isinstance(items[i],(Primitive, ReversePrimitive)):
             raise TypeError("%s-th argument is not a Primitive object"%i)
       for i in items:
          self.__items.append(i)
    def getItems(self):
        """
        returns a list of the items used in the design
        """
        return self.__items
    def clearItems(self):
        """
        resets the items in design
        """
        self.__items=[]
    def getAllPrimitives(self):
        """
        returns a list of all primitives used to create the design.
        each primitve appears once. The primitives are ordered by their
        order of generation
        """
        prims=[]
        for i in self.getItems(): 
            for p in i.getPrimitives():
                if not p in prims: prims.append(p)
        prims.sort()
        return prims

    def setOptions(self,**kwargs):
        """
        sets options of the mesh generator

        @note: this method is typically overwritten by a particular Design implementation
        """
        pass
    def getMeshHandler(self):
        """
        returns a handle to a mesh meshing the design

        @note: this method has to be overwritten by a particular Design implementation
        """
        raise NotImplementedError()

    def getTagMap(self):
        """
        returns a L{TagMap} to map the names of L{PropertySet}s to tags
        """
        m={}
        for p in self.getAllPrimitives():
           if isinstance(p, PropertySet): m[ p.getTag() ] = p.getName()
        return TagMap(m)
