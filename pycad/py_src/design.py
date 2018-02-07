
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

"""
Template for the Design which defines regions and features
for a mesh generator.

:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"

from .primitives import Primitive, ReversePrimitive, PropertySet, Point, Manifold1D, Manifold2D, Manifold3D
from xml.dom import minidom
import tempfile, os

class TagMap(object):
    """
    A class that allows to map tags to names.

    Example::

        tm=TagMap({5 : x })
        tm.setMap(a=1,x=4)
        assert tm.getTags("a") == [ 1 ]
        assert tm.getTags("x") == [ 5, 4 ]
        assert tm.map(x=10., a=20.) == { 5 : 10, 4: 10, 1 : 20 }

    """
    def __init__(self, mapping={}):
      """
      Initializes the mapping. ``mapping`` defines an initial mapping from tag
      to a name.
      """
      self.__mapping={}
      for tag, name in sorted(mapping.items(), key=lambda x: x[1]):
          if not isinstance(tag, int):
              raise TypeError("tag needs to be an int")
          if not isinstance(name, str):
              raise TypeError("name needs to be a str.")
          self.__mapping[tag]=name

    def setMap(self,**kwargs):
      """
      Sets a new map where <name>=<tag> assigns the tag <tag> to name <name>.
      <tag> has to be an integer. If <tag> has been assigned to a name before
      the mapping will be overwritten. Otherwise a new mapping <tag> -> <name>
      is set. Notice that a single name can be assigned to different tags.
      """
      for name, tag in sorted(kwargs.items(), key=lambda x: x[0]):
          if not isinstance(tag, int):
             raise TypeError("tag needs to be an int")
          self.__mapping[tag]=name

    def getTags(self,name=None):
        """
        Returns a list of the tags assigned to ``name``. If name is not present
        a list of all tags is returned.
        """
        if name == None:
           out=sorted(self.__mapping.keys())
        else:
           out=[]
           for tag, arg in sorted(self.__mapping.items(), key=lambda x: x[0]):
             if arg == name: out.append(tag)
        return out

    def getName(self,tag=None):
        """
        Returns the name of a tag. If ``tag`` is not present a list of all names
        is returned.
        """
        if tag == None:
           return sorted(list(set(self.__mapping.values())))
        else:
            return self.__mapping[tag]

    def getMapping(self):
        """
        Returns a dictionary where the tags define the keys and the values the
        corresponding names.
        """
        return self.__mapping

    def map(self,default=0,**kwargs):
        """
        Returns a dictionary where the tags define the keys and the values give
        the values assigned to the tag via name and kwargs::

            tm=TagMap(x=5)
            tm.setMap(a=1,x=4)
            print tm.map(x=10., a=20.)
            { 5 : 10, 4: 10, 1 : 20 }

        The default is used for tags which map onto name with unspecified
        values.
        """
        out={}
        for tag in self.__mapping:
           if self.__mapping[tag] in kwargs:
              out[tag]=kwargs[self.__mapping[tag]]
           else:
              out[tag]=default
        return out

    def insert(self,data,default=0,**kwargs):
        """
        Inserts values into the `esys.escript.Data` object according to the
        given values assigned to the keywords. The default is used for tags
        which map onto name with unspecified values.
        """
        d=self.map(default=default,**kwargs)
        for t,v in sorted(d.items(), key=lambda x: x[0]):
             data.setTaggedValue(t,v)

    def passToDomain(self,domain):
        """
        Passes the tag map to the `esys.escript.Domain` ``domain``.
        """
        for tag, name in sorted(self.__mapping.items(), key=lambda x: x[1]):
          print("Tag",name, "is mapped to id ", tag)
          domain.setTagMap(name,tag)

    def toDOM(self,dom):
         """
         Adds object to ``dom``.
         """
         tm=dom.createElement("TagMap")
         dom.appendChild(tm)
         for tag,name in sorted(self.getMapping().items(), key=lambda x: x[1]):
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
        Fills names and tags from dom ``node``.
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
       Uses the XML file or string to set the mapping.
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
         Serializes self as XML into ``iostream`` or if not present returns the
         XML as string.
         """
         dom=minidom.Document()
         esys=dom.createElement('ESys')
         esys.appendChild(self.toDOM(dom))
         dom.appendChild(esys)
         if iostream == None:
            return dom.toprettyxml()
         else:
            iostream.write(dom.toprettyxml())

class AbstractDesign(object):
    """
    Template for a design which defines the input for a mesh generator.

    :note: class variable GMSH - `gmsh <http://www.geuz.org/gmsh/>`_ file format
    :note: class variable IDEAS - `I_DEAS <http://www.plm.automation.siemens.com/en_us/products/nx/>`_ universal file format
    :note: class variable VRML - `VRML <http://www.w3.org/MarkUp/VRML/>`_ file format
    :note: class variable STL - `STL <http://en.wikipedia.org/wiki/STL_(file_format)>`_ file format
    :note: class variable NASTRAN - `Nastran <http://simcompanion.mscsoftware.com/infocenter/index?page=content&channel=DOCUMENTATION>`_ bulk data format
    :note: class variable MEDIT - `Medit <http://www-rocq.inria.fr/OpenFEM/Doc/>`_ file format
    :note: class variable CGNS - `CGNS <http://cgns.sourceforge.net/>`_ file format
    :note: class variable PLOT3D -  `Plot3D <http://www.plot3d.net/>`_ file format
    :note: class variable DIFFPACK- `Diffpack <http://www.diffpack.com/>`_ file format
    """
    GMSH="msh"
    IDEAS="unv"
    VRML="vrml"
    STL="stl"
    NASTRAN="bdf"
    MEDIT="mesh"
    CGNS="cgns"
    PLOT3D="p3d"
    DIFFPACK="diff"
    def __init__(self,dim=3,element_size=1.,order=1,keep_files=False):
       """
       Initializes a design.

       :param dim: spatial dimension
       :param element_size: global element size
       :param order: element order
       :param keep_files: flag to keep work files
       """
       self.clearItems()
       self.setElementSize(element_size)
       self.setDim(dim)
       self.setElementOrder(order)
       self.setFileFormat()
       if keep_files:
          self.setKeepFilesOn()
       else:
          self.setKeepFilesOff()
       self.__mshname=""
       self.setMeshFileName()

    def setDim(self,dim=3):
        """
        Sets the spatial dimension.
        """
        if not dim  in [1,2,3]:
           raise ValueError("only dimension 1, 2, 3 are supported.")
        self.__dim=dim

    def getDim(self,dim=3):
        """
        Returns the spatial dimension.
        """
        return self.__dim

    def setElementOrder(self,order=1):
        """
        Sets the element order.
        """
        if not order in [1,2]:
           raise ValueError("only element order 1 or 2 is supported.")
        self.__order=order

    def getElementOrder(self):
        """
        Returns the element order.
        """
        return self.__order

    def setElementSize(self,element_size=1.):
        """
        Sets the global element size.
        """
        if element_size<=0.:
           raise ValueError("element size needs to be positive.")
        self.__element_size=element_size

    def getElementSize(self):
        """
        Returns the global element size.
        """
        return self.__element_size

    def setKeepFilesOn(self):
        """
        Work files are kept at the end of the generation.
        """
        self.__keep_files=True

    def setKeepFilesOff(self):
        """
        Work files are deleted at the end of the generation
        """
        self.__keep_files=False

    def keepFiles(self):
        """
        Returns True if work files are kept, False otherwise.
        """
        return self.__keep_files

    def addItems(self,*items):
       """
       Adds items to the design.
       """
       new_items=[]
       for i in range(len(items)):
          if not isinstance(items[i],(Primitive, ReversePrimitive)):
             raise TypeError("%s-th argument is not a Primitive object"%i)
          if isinstance(items[i],PropertySet):
             q=items[i]
          else:
             q=PropertySet("__%s__"%(items[i].getID()), items[i])
          for p in self.getAllPrimitives():
              if isinstance(p, PropertySet): 
                if q.getName() == p.getName():
                   raise ValueError("Property set name %s is allready in use."%q.getName())
          new_items.append(q)
       for q in new_items: self.__items.append(q)

    def getItems(self):
        """
        Returns a list of the items used in the design.
        """
        return self.__items

    def clearItems(self):
        """
        Removes all items from the design.
        """
        self.__items=[]

    def getAllPrimitives(self):
        """
        Returns a list of all primitives used to create the design.
        Each primitive appears once. The primitives are ordered by their
        order of generation.
        """
        prims=[]
        for i in self.getItems():
            for p in i.getPrimitives():
                if not p in prims: prims.append(p)
        prims.sort()
        return prims

    def setOptions(self,**kwargs):
        """
        Sets options of the mesh generator.

        :note: this method is typically overwritten by a particular design
               implementation.
        """
        pass
    def generate(self):
        """
        generate output file 
        
        :note: this method may be overwritten by a particular design
               implementation.
        """
        self.getMeshHandler()

    def getMeshHandler(self):
        """
        Returns a handle to a mesh meshing the design.

        :note: this method has to be overwritten by a particular design
               implementation.
        """
        raise NotImplementedError()

    def getTagMap(self):
        """
        Returns a `TagMap` to map the names of `PropertySet` s to tags.
        """
        m={}
        for p in self.getAllPrimitives():
           if isinstance(p, PropertySet): m[ p.getTag() ] = p.getName()
        return TagMap(m)

    def setFileFormat(self,format='msh'):
       """
       Sets the file format to be used.

       :param format: format to be used. needs to be one of 

       """
       if not format in [ self.GMSH, self.IDEAS, self.VRML, self.STL, self.NASTRAN, self.MEDIT, self.CGNS, self.PLOT3D, self.DIFFPACK] :
           raise ValueError("unknown file format %s."%format)
       self.__fileformat=format
           
    def getFileFormat(self):
       """
       Returns the file format
       """
       return self.__fileformat

    def setMeshFileName(self, name=None):
       """
       Sets the name for the mesh file. If no name is given a name is generated.
       """
       if self.__mshname:
           os.unlink(self.__mshname)
       if name == None:
           self.__mshname_set=False
           tmp_f_id=tempfile.mkstemp(suffix="."+self.getFileFormat())
           self.__mshname=tmp_f_id[1]
           os.close(tmp_f_id[0])
       else:
           self.__mshname=name
           self.__mshname_set=True

    def getMeshFileName(self):
       """
       Returns the name of the mesh file.
       """
       return self.__mshname



