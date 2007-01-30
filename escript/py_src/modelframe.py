# $Id$

"""
Environment for implementing models in escript 

@var __author__: name of author
@var __copyright__: copyrights
@var __license__: licence agreement
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


from types import StringType,IntType,FloatType,BooleanType,ListType,DictType
from sys import stdout
import numarray
import operator
import itertools

# import the 'set' module if it's not defined (python2.3/2.4 difference)
try:
    set
except NameError:
    from sets import Set as set

from xml.dom import minidom


def all(seq):
    for x in seq:
        if not x:
            return False
    return True

def any(seq):
    for x in seq:
        if x:
            return True
    return False

def importName(modulename, name):
    """ Import a named object from a module in the context of this function,
        which means you should use fully qualified module paths.
        Return None on failure.

        This function from: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52241
    """
    module = __import__(modulename, globals(), locals(), [name])
        
    try:
        return vars(module)[name]
    except KeyError:
        raise ImportError("Could not import %s from %s" % (name, modulename))

class ESySXMLParser(object):
    """
    parser for ESysXML file
    """
    def __init__(self,xml, debug=False):
      self.__dom = minidom.parseString(xml)
      self.__linkable_object_registry= {}
      self.__link_registry=  []
      self.__esys=self.__dom.getElementsByTagName('ESys')[0]
      self.debug=debug
  
    def getClassPath(self, node): 
        type = node.getAttribute("type")
        if (node.getAttribute("module")):
            module = node.getAttribute("module")
            return importName(module, type)
        else:
            return importName("__main__", type)

    def setLinks(self):
        for obj_id, link in self.__link_registry: 
            link.target = self.__linkable_object_registry[obj_id]

    def parse(self):
       """
       parser method for EsysXML and returns the list of generating ParameterSets 
       """
       found=[]
       for node in self.__esys.childNodes:
           if isinstance(node, minidom.Element):
               if node.tagName == 'Simulation':
                        found.append(Simulation.fromDom(self, node))
               elif node.tagName == 'Model':
                        found.append(self.getClassPath(node).fromDom(self, node))
               elif node.tagName == 'ParameterSet':
                        found.append(self.getClassPath(node).fromDom(self, node))
               else:
                  raise "Invalid type, %r" % node.getAttribute("type")
       self.setLinks()
       return found

    def registerLink(self,obj_id, link):
        self.__link_registry.append((int(obj_id),link))

    def registerLinkableObject(self,obj, node):
        id_str=node.getAttribute('id').strip()
        if len(id_str)>0:
           id=int(id_str)
           if self.__linkable_object_registry.has_key(id):
               raise ValueError("Object id %s already exists."%id)
           else:
               self.__linkable_object_registry[id]=obj

    def getComponent(self, node):
       """
       returns a single component + rank from a simulation 
       parser method for EsysXML and returns the list of generating ParameterSets 
       """
       rank  = int(node.getAttribute("rank"))
       for n in node.childNodes:
           if isinstance(n, minidom.Element):
               if n.tagName == 'Simulation':
                        return (rank, Simulation.fromDom(self, n))
               elif n.tagName == 'Model':
                        return (rank, self.getClassPath(n).fromDom(self, n))
               elif n.tagName == 'ParameterSet':
                        return (rank, self.getClassPath(n).fromDom(self, n))
               else:
                 raise ValueError("illegal component type %s"%n.tagName)
       raise ValueError("cannot resolve Component")

class ESySXMLCreator(object):
    """
    creates an XML Dom representation
    """
    def __init__(self):
       self.__dom=minidom.Document()
       self.__esys =self.__dom.createElement('ESys')
       self.__dom.appendChild(self.__esys)
       self.__linkable_object_registry={}
       self.__number_sequence = itertools.count(100)
    def getRoot(self):
       return self.__esys
    def createElement(self,name):
      return self.__dom.createElement(name)
    def createTextNode(self,name):
      return self.__dom.createTextNode(name) 
    def getElementById(self,name):
      return self.__dom.getElementById(name)
    def createDataNode(self, tagName, data):
          """
          C{createDataNode}s are the building blocks of the xml documents constructed in
          this module.  
    
          @param tagName: the associated xml tag
          @param data: the values in the tag
          """
          n = self.createElement(tagName)
          n.appendChild(self.createTextNode(str(data)))
          return n
    def getLinkableObjectId(self, obj):
        for id, o in self.__linkable_object_registry.items():
            if o == obj: return id
        id =self.__number_sequence.next()
        self.__linkable_object_registry[id]=obj
        return id
         
    def registerLinkableObject(self, obj, node):
        """
        returns a unique object id for object obj
        """
        id=self.getLinkableObjectId(obj)
        node.setAttribute('id',str(id))
        node.setIdAttribute("id")

    def includeTargets(self):
        target_written=True
        while target_written:
            targetsList =self.__dom.getElementsByTagName('Target')
            target_written=False
            for element in targetsList:
               targetId = int(element.firstChild.nodeValue.strip())
               if self.getElementById(str(targetId)): continue
               targetObj = self.__linkable_object_registry[targetId]
               targetObj.toDom(self, self.__esys)
               target_written=True

    def toprettyxml(self):
        self.includeTargets()
        return self.__dom.toprettyxml()
 
class Link:
    """
    A Link makes an attribute of an object callable:: 

          o.object()
          o.a=8
          l=Link(o,"a")
          assert l()==8
     """
    
    def __init__(self,target,attribute=None):
        """
        Creates a link to the object target. If attribute is given, the link is
        establised to this attribute of the target.  Otherwise the attribute is
        undefined.
        """
        self.target = target
        self.attribute = None
        self.setAttributeName(attribute)

    def getTarget(self):
        """
        returns the target
        """ 
        return self.target 
    def getAttributeName(self):
        """
        returns the name of the attribute the link is pointing to
        """
        return self.attribute
    
    def setAttributeName(self,attribute):
        """
        Set a new attribute name to be collected from the target object. The
        target object must have the attribute with name attribute.
        """
        if attribute and self.target:
            if isinstance(self.target,LinkableObject):
               if not self.target.hasAttribute(attribute):
                  raise AttributeError("%s: target %s has no attribute %s."%(self, self.target, attribute))
            else:
               if not hasattr(self.target,attribute):
                  raise AttributeError("%s: target %s has no attribute %s."%(self, self.target, attribute))
        self.attribute = attribute
    
    def hasDefinedAttributeName(self):
        """
        Returns true if an attribute name is set.
        """
        return self.attribute != None
    
    def __repr__(self):
        """
        Returns a string representation of the link.
        """
        if self.hasDefinedAttributeName():
            return "<Link to attribute %s of %s>" % (self.attribute, self.target)
        else:
            return "<Link to target %s>" % self.target
    
    def __call__(self,name=None):
        """
        Returns the value of the attribute of the target object. If the
        atrribute is callable then the return value of the call is returned.
        """
        if name:
            out=getattr(self.target, name)
        else:
            out=getattr(self.target, self.attribute)

        if callable(out):
            return out()
        else:
            return out

    def toDom(self, esysxml, node):
        """
        C{toDom} method of Link. Creates a Link node and appends it to the
	current XML esysxml.
        """
        link = esysxml.createElement('Link')
        assert (self.target != None), ("Target was none, name was %r" % self.attribute)
        link.appendChild(esysxml.createDataNode('Target', esysxml.getLinkableObjectId(self.target)))
        # this use of id will not work for purposes of being able to retrieve the intended
        # target from the xml later. I need a better unique identifier.
        assert self.attribute, "You can't xmlify a Link without a target attribute"
        link.appendChild(esysxml.createDataNode('Attribute', self.attribute))
        node.appendChild(link)

    def fromDom(cls, esysxml, node):
        targetid = int(node.getElementsByTagName("Target")[0].firstChild.nodeValue.strip())
        attribute =str(node.getElementsByTagName("Attribute")[0].firstChild.nodeValue.strip())
        l = cls(None, attribute)
        esysxml.registerLink(targetid, l)
        return l

    fromDom = classmethod(fromDom)
    
class LinkableObject(object):
    """
    An object that allows to link its attributes to attributes of other objects
    via a Link object. For instance::
           
           p = LinkableObject()
           p.x = Link(o,"name")
           print p.x
    
    links attribute C{x} of C{p} to the attribute name of object C{o}. 

    C{p.x} will contain the current value of attribute C{name} of object
    C{o}.  

    If the value of C{getattr(o, "name")} is callable, C{p.x} will return 
    the return value of the call. 
    """
   
    
    def __init__(self, id = None, debug=False):
        """
	Initializes LinkableObject so that we can operate on Links 
	"""
        self.debug = debug
        self.__linked_attributes={}
        
    def trace(self, msg):
        """
	If debugging is on, print the message, otherwise do nothing
        """
        if self.debug:
            print "%s: %s"%(str(self),msg)
    
    def __getattr__(self,name):
        """
	Returns the value of attribute name. If the value is a Link object the
        object is called and the return value is returned.
	"""
        out = self.getAttributeObject(name)
        if isinstance(out,Link):
            return out()
        else:
            return out
    
    def getAttributeObject(self,name):
        """
	Return the object stored for attribute name.
	"""

        if self.__dict__.has_key(name):
            return self.__dict__[name]

        if self.__linked_attributes.has_key(name):
            return self.__linked_attributes[name]

        if self.__class__.__dict__.has_key(name):
            return self.__class.__dict__[name]

        raise AttributeError,"No attribute %s."%name
    
    def hasAttribute(self,name):
        """
	Returns True if self as attribute name.
	"""
        return self.__dict__.has_key(name) or self.__linked_attributes.has_key(name) or  self.__class__.__dict__.has_key(name)

    def __setattr__(self,name,value):
        """
	Sets the value for attribute name. If value is a Link the target
        attribute is set to name if no attribute has been specified.
	"""

        if self.__dict__.has_key(name): 
            del self.__dict__[name]

        if isinstance(value,Link):
            if not value.hasDefinedAttributeName(): 
                value.setAttributeName(name)
            self.__linked_attributes[name] = value

            self.trace("attribute %s is now linked by %s."%(name,value))
        else:
            self.__dict__[name] = value
    
    def __delattr__(self,name):
        """
	Removes the attribute name.
	"""

        if self.__linked_attributes.has_key[name]:
            del self.__linked_attributes[name]
        elif self.__dict__.has_key(name):
            del self.__dict__[name]
        else:
            raise AttributeError,"No attribute %s."%name

class _ParameterIterator:
    def __init__(self,parameterset):

        self.__set=parameterset
        self.__iter=iter(parameterset.parameters)

    def next(self):
        o=self.__iter.next()
        return (o,self.__set.getAttributeObject(o))

    def __iter__(self):
        return self

class ParameterSet(LinkableObject):
    """
    A class which allows to emphazise attributes to be written and read to XML
       
    Leaves of an ESySParameters object can be:
    
	 - a real number
	 - a integer number
	 - a string
	 - a boolean value
	 - a ParameterSet object 
	 - a Simulation object
	 - a Model object 
	 - a numarray object
         - a list of booleans
        - any other object (not considered by writeESySXML and writeXML)
    
    Example how to create an ESySParameters object::
    
        p11=ParameterSet(gamma1=1.,gamma2=2.,gamma3=3.)
        p1=ParameterSet(dim=2,tol_v=0.001,output_file="/tmp/u.%3.3d.dx",runFlag=True,parm11=p11) 
        parm=ParameterSet(parm1=p1,parm2=ParameterSet(alpha=Link(p11,"gamma1")))
    
    This can be accessed as::
    
	parm.parm1.gamma=0.
        parm.parm1.dim=2
        parm.parm1.tol_v=0.001
        parm.parm1.output_file="/tmp/u.%3.3d.dx"
        parm.parm1.runFlag=True
        parm.parm1.parm11.gamma1=1.
        parm.parm1.parm11.gamma2=2.
        parm.parm1.parm11.gamma3=3.
        parm.parm2.alpha=1. (value of parm.parm1.parm11.gamma1)
    """
    def __init__(self, parameters=[], **kwargs):
        """
	Creates a ParameterSet with parameters parameters.
	"""
        LinkableObject.__init__(self, **kwargs)
        self.parameters = set()
        self.declareParameters(parameters)

    def __repr__(self):
        return "<%s %d>"%(self.__class__.__name__,id(self))
    
    def declareParameter(self,**parameters):
        """
	Declares a new parameter(s) and its (their) initial value.
	"""
        self.declareParameters(parameters)
    
    def declareParameters(self,parameters):
        """
	Declares a set of parameters. parameters can be a list, a dictionary 
	or a ParameterSet.
	"""
        if isinstance(parameters,ListType):
            parameters = zip(parameters, itertools.repeat(None))
        if isinstance(parameters,DictType):
            parameters = parameters.iteritems()

        for prm, value in parameters:
            setattr(self,prm,value)
            self.parameters.add(prm)

            self.trace("parameter %s has been declared."%prm)

    def releaseParameters(self,name):
        """
	Removes parameter name from the paramameters.
	"""
        if self.isParameter(name): 
            self.parameters.remove(name)
            self.trace("parameter %s has been removed."%name)
    
    def __iter__(self):
        """
	Creates an iterator over the parameter and their values.
	"""
        return _ParameterIterator(self)
    
    def showParameters(self):
        """
	Returns a descrition of the parameters.
	"""        
        out="{"
        notfirst=False
        for i,v in self:
            if notfirst: out=out+","
            notfirst=True
            if isinstance(v,ParameterSet):
                out="%s\"%s\" : %s"%(out,i,v.showParameters())
            else:
                out="%s\"%s\" : %s"%(out,i,v)
        return out+"}"
    
    def __delattr__(self,name):
        """
	Removes the attribute name.
	"""
        LinkableObject.__delattr__(self,name)
        try:
            self.releaseParameter(name) 
        except:
            pass

    def toDom(self, esysxml, node):
        """
	C{toDom} method of Model class
	"""
        pset = esysxml.createElement('ParameterSet')
        pset.setAttribute('type', self.__class__.__name__)
        pset.setAttribute('module', self.__class__.__module__)
        esysxml.registerLinkableObject(self, pset)
        self._parametersToDom(esysxml, pset)
        node.appendChild(pset)

    def _parametersToDom(self, esysxml, node):
        for name,value in self:
            # convert list to numarray when possible:
            if isinstance (value, list):
                elem_type=-1
                for i in value:
                    if isinstance(i, bool): 
                        elem_type = max(elem_type,0)
                    if isinstance(i, int) and not isinstance(i, bool): 
                        elem_type = max(elem_type,1)
                    if isinstance(i, float):
                        elem_type = max(elem_type,2)
                if elem_type == 0: value = numarray.array(value,numarray.Bool)
                if elem_type == 1: value = numarray.array(value,numarray.Int)
                if elem_type == 2: value = numarray.array(value,numarray.Float)

            param = esysxml.createElement('Parameter')
            param.setAttribute('type', value.__class__.__name__)

            param.appendChild(esysxml.createDataNode('Name', name))

            val = esysxml.createElement('Value')
            if isinstance(value,(ParameterSet,Link,DataSource)):
                value.toDom(esysxml, val)
                param.appendChild(val)
            elif isinstance(value, numarray.NumArray):
                shape = value.getshape()
                if isinstance(shape, tuple):
                    size = reduce(operator.mul, shape)
                    shape = ' '.join(map(str, shape))
                else:
                    size = shape
                    shape = str(shape)

                arraytype = value.type()
                if isinstance(arraytype, numarray.BooleanType):
                      arraytype_str="Bool"
                elif isinstance(arraytype, numarray.IntegralType):
                      arraytype_str="Int"
                elif isinstance(arraytype, numarray.FloatingType):
                      arraytype_str="Float"
                elif isinstance(arraytype, numarray.ComplexType):
                      arraytype_str="Complex"
                else:
                      arraytype_str=str(arraytype)
                numarrayElement = esysxml.createElement('NumArray')
                numarrayElement.appendChild(esysxml.createDataNode('ArrayType', arraytype_str))
                numarrayElement.appendChild(esysxml.createDataNode('Shape', shape))
                numarrayElement.appendChild(esysxml.createDataNode('Data', ' '.join(
                    [str(x) for x in numarray.reshape(value, size)])))
                val.appendChild(numarrayElement)
                param.appendChild(val)
            elif isinstance(value, list):
                param.appendChild(esysxml.createDataNode('Value', ' '.join([str(x) for x in value]) ))
            elif isinstance(value, (str, bool, int, float, type(None))):
                param.appendChild(esysxml.createDataNode('Value', str(value)))
            elif isinstance(value, dict):
                 dic = esysxml.createElement('dictionary')
                 if len(value.keys())>0:
                     dic.setAttribute('key_type', value.keys()[0].__class__.__name__)
                     dic.setAttribute('value_type', value[value.keys()[0]].__class__.__name__)
                 for k,v in value.items():
                    i=esysxml.createElement('item')
                    i.appendChild(esysxml.createDataNode('key', k))
                    i.appendChild(esysxml.createDataNode('value', v))
                    dic.appendChild(i)
                 param.appendChild(dic)
            else:
                print value
                raise ValueError("cannot serialize %s type to XML."%str(value.__class__))

            node.appendChild(param)

    def fromDom(cls, esysxml, node):
        # Define a host of helper functions to assist us.
        def _children(node):
            """
            Remove the empty nodes from the children of this node.
            """
            ret = []
            for x in node.childNodes:
                if isinstance(x, minidom.Text):
                    if x.nodeValue.strip():
                        ret.append(x)
                else:
                    ret.append(x)
            return ret

        def _floatfromValue(esysxml, node):
            return float(node.nodeValue.strip())

        def _stringfromValue(esysxml, node):
            return str(node.nodeValue.strip())
       
        def _intfromValue(esysxml, node):
            return int(node.nodeValue.strip())

        def _boolfromValue(esysxml, node):
            return _boolfromstring(node.nodeValue.strip())

        def _nonefromValue(esysxml, node):
            return None

        def _numarrayfromValue(esysxml, node):
            for node in _children(node):
                if node.tagName == 'ArrayType':
                    arraytype = node.firstChild.nodeValue.strip()
                if node.tagName == 'Shape':
                    shape = node.firstChild.nodeValue.strip()
                    shape = [int(x) for x in shape.split()]
                if node.tagName == 'Data':
                    data = node.firstChild.nodeValue.strip()
                    data = [float(x) for x in data.split()]
            return numarray.reshape(numarray.array(data, type=getattr(numarray, arraytype)),
                                    shape)
      
        def _listfromValue(esysxml, node):
            return [x for x in node.nodeValue.split()]

        def _boolfromstring(s):
            if s == 'True':
                return True
            else:
                return False
        # Mapping from text types in the xml to methods used to process trees of that type
        ptypemap = {"Simulation": Simulation.fromDom,
                    "Model":Model.fromDom,
                    "ParameterSet":ParameterSet.fromDom,
                    "Link":Link.fromDom,
                    "DataSource":DataSource.fromDom,
                    "float":_floatfromValue,
                    "int":_intfromValue,
                    "str":_stringfromValue,
                    "bool":_boolfromValue,
                    "list":_listfromValue,
                    "NumArray":_numarrayfromValue,
                    "NoneType":_nonefromValue,
                    }

        parameters = {}
        for n in _children(node):
            ptype = n.getAttribute("type")
            if not ptypemap.has_key(ptype):
               raise KeyError("cannot handle parameter type %s."%ptype)

            pname = pvalue = None
            for childnode in _children(n):
                if childnode.tagName == "Name":
                    pname = childnode.firstChild.nodeValue.strip()

                if childnode.tagName == "Value":
                    nodes = _children(childnode)
                    pvalue = ptypemap[ptype](esysxml, nodes[0])

            parameters[pname] = pvalue

        # Create the instance of ParameterSet
        o = cls(debug=esysxml.debug)
        o.declareParameters(parameters)
        esysxml.registerLinkableObject(o, node)
        return o
    
    fromDom = classmethod(fromDom)

    def writeXML(self,ostream=stdout):
        """
	Writes the object as an XML object into an output stream.
	"""
        esysxml=ESySXMLCreator()
        self.toDom(esysxml, esysxml.getRoot())
        ostream.write(esysxml.toprettyxml())
    
class Model(ParameterSet):
    """
    A Model object represents a processess marching over time until a
    finalizing condition is fullfilled. At each time step an iterative
    process can be performed and the time step size can be controlled. A
    Model has the following work flow::
          
          doInitialization()
          while not terminateInitialIteration(): doInitializationiStep()
          doInitialPostprocessing()
          while not finalize():
               dt=getSafeTimeStepSize(dt)
               doStepPreprocessing(dt)
               while not terminateIteration(): doStep(dt)
               doStepPostprocessing(dt)
          doFinalization()

    where C{doInitialization}, C{finalize}, C{getSafeTimeStepSize},
    C{doStepPreprocessing}, C{terminateIteration}, C{doStepPostprocessing},
    C{doFinalization} are methods of the particular instance of a Model. The
    default implementations of these methods have to be overwritten by the
    subclass implementing a Model.
    """

    UNDEF_DT=1.e300

    def __init__(self,parameters=[],**kwargs):
        """
	Creates a model.

        Just calls the parent constructor.
        """
        ParameterSet.__init__(self, parameters=parameters,**kwargs)

    def __str__(self):
       return "<%s %d>"%(self.__class__.__name__,id(self))


    def doInitialization(self):
        """
	Initializes the time stepping scheme.  
	
	This function may be overwritten.
	"""
        pass
    def doInitialStep(self):
        """
	performs an iteration step in the initialization phase

	This function may be overwritten.
	"""
        pass

    def terminateInitialIteration(self):
        """
	Returns True if iteration at the inital phase is terminated.
	"""
        return True

    def doInitialPostprocessing(self):
        """
	finalises the initialization iteration process 

	This function may be overwritten.
	"""
        pass
    
    def getSafeTimeStepSize(self,dt):
        """
	Returns a time step size which can safely be used. 

        C{dt} gives the previously used step size.

        This function may be overwritten.
	"""
        return self.UNDEF_DT
    
    def finalize(self):
        """
	Returns False if the time stepping is finalized. 
	
	This function may be overwritten.
	"""
        return False
       
    def doFinalization(self):
        """
	Finalizes the time stepping. 
	
	This function may be overwritten.
	"""
        pass
    
    def doStepPreprocessing(self,dt):
        """
	Sets up a time step of step size dt. 
	
	This function may be overwritten.
	"""
        pass
    
    def doStep(self,dt):
        """
	Executes an iteration step at a time step. 

        C{dt} is the currently used time step size.

        This function may be overwritten.
	"""
        pass
    
    def terminateIteration(self):
        """
	Returns True if iteration on a time step is terminated.
	"""
        return True

       
    def doStepPostprocessing(self,dt):
        """
	finalises the time step.

        dt is the currently used time step size.

        This function may be overwritten.
	"""
        pass

    def toDom(self, esysxml, node):
        """
	C{toDom} method of Model class
	"""
        pset = esysxml.createElement('Model')
        pset.setAttribute('type', self.__class__.__name__)
        pset.setAttribute('module', self.__class__.__module__)
        esysxml.registerLinkableObject(self, pset)
        node.appendChild(pset)
        self._parametersToDom(esysxml, pset)
    
class Simulation(Model): 
    """
    A Simulation object is special Model which runs a sequence of Models. 

    The methods C{doInitialization}, C{finalize}, C{getSafeTimeStepSize},
    C{doStepPreprocessing}, C{terminateIteration}, C{doStepPostprocessing},
    C{doFinalization} are executing the corresponding methods of the models in
    the simulation.
    """
    
    FAILED_TIME_STEPS_MAX=20
    MAX_ITER_STEPS=50
    MAX_CHANGE_OF_DT=2.
    
    def __init__(self, models=[], **kwargs):
        """
	Initiates a simulation from a list of models.
	"""
        super(Simulation, self).__init__(**kwargs)
        for m in models:
            if not isinstance(m, Model):
                 raise TypeError("%s is not a subclass of Model."%m)
        self.__models=[]
        for i in range(len(models)): 
            self[i] = models[i]
            

    def __repr__(self):
        """
        Returns a string representation of the Simulation.
        """
        return "<Simulation %r>" % self.__models

    def __str__(self):
        """
        Returning Simulation as a string.
        """
        return "<Simulation %d>"%id(self)
    
    def iterModels(self):
        """
	Returns an iterator over the models.
	"""
        return self.__models
    
    def __getitem__(self,i):
        """
	Returns the i-th model.
	"""
        return self.__models[i]
    
    def __setitem__(self,i,value):
        """
	Sets the i-th model.
	"""
        if not isinstance(value,Model):
            raise ValueError,"assigned value is not a Model but instance of %s"%(value.__class__.__name__,)
        for j in range(max(i-len(self.__models)+1,0)): 
            self.__models.append(None)
        self.__models[i]=value
    
    def __len__(self):
        """
	Returns the number of models.
	"""
        return len(self.__models)

    def getAllModels(self):
        """
        returns a list of all models used in the Simulation including subsimulations
        """
        out=[]
        for m in self.iterModels():
            if isinstance(m, Simulation):
               out+=m.getAllModels()
            else:
               out.append(m)
        return list(set(out))

    def checkModelLinks(self, models):
        """
        returns a list of (model,parameter, target model ) if the the parameter of model
        is linking to the target_model which is not in list of models.
        """
        out=[]
        for m in self.iterModels():
            if isinstance(m, Simulation):
               out+=[ (m,) + f for f in  m.checkModelLinks(models) ]
            else:
              for p in m:
                 if isinstance(p[1], Link):
                    l_m=p[1].getTarget()
                    if isinstance(l_m, Model) and not l_m in models: out.append( (m,p[0],l_m) )
        return out

    
    def getSafeTimeStepSize(self,dt):
        """
	Returns a time step size which can safely be used by all models.

        This is the minimum over the time step sizes of all models.
	"""
        out=min([o.getSafeTimeStepSize(dt) for o in self.iterModels()])
        return out
    
    def doInitialization(self):
        """
	Initializes all models.
	"""
        self.n=0
        self.tn=0.
        for o in self.iterModels(): 
             o.doInitialization()
    def doInitialStep(self):
        """
	performs an iteration step in the initialization step for all models
	"""
        iter=0
        while not self.terminateInitialIteration():
            if iter==0: self.trace("iteration for initialization starts")
            iter+=1
            self.trace("iteration step %d"%(iter))
            for o in self.iterModels(): 
                 o.doInitialStep()
            if iter>self.MAX_ITER_STEPS:
                 raise IterationDivergenceError("initial iteration did not converge after %s steps."%iter)
        self.trace("Initialization finalized after %s iteration steps."%iter)

    def doInitialPostprocessing(self):
        """
	finalises the initialization iteration process for all models.
	"""
        for o in self.iterModels(): 
            o.doInitialPostprocessing()
    def finalize(self):
        """
	Returns True if any of the models is to be finalized.
	"""
        return any([o.finalize() for o in self.iterModels()])
       
    def doFinalization(self):
        """
	finalises the time stepping for all models.
	"""
        for i in self.iterModels(): i.doFinalization()
        self.trace("end of time integation.")
    
    def doStepPreprocessing(self,dt):
        """
	Initializes the time step for all models.
	"""
        for o in self.iterModels(): 
            o.doStepPreprocessing(dt)
    
    def terminateIteration(self):
        """
	Returns True if all iterations for all models are terminated.
	"""
        out=all([o.terminateIteration() for o in self.iterModels()])
        return out

    def terminateInitialIteration(self):
        """
	Returns True if all initial iterations for all models are terminated.
	"""
        out=all([o.terminateInitialIteration() for o in self.iterModels()])
        return out
       
    def doStepPostprocessing(self,dt):
        """
	finalises the iteration process for all models.
	"""
        for o in self.iterModels(): 
            o.doStepPostprocessing(dt)
        self.n+=1
        self.tn+=dt
    
    def doStep(self,dt):
        """
	Executes the iteration step at a time step for all model::
 
            self.doStepPreprocessing(dt)
            while not self.terminateIteration(): 
	        for all models: 
		    self.doStep(dt)
                self.doStepPostprocessing(dt)
        """
        self.iter=0
        while not self.terminateIteration(): 
            if self.iter==0: self.trace("iteration at %d-th time step %e starts"%(self.n+1,self.tn+dt))
            self.iter+=1
            self.trace("iteration step %d"%(self.iter))
            for o in self.iterModels(): 
                  o.doStep(dt)
        if self.iter>0: self.trace("iteration at %d-th time step %e finalized."%(self.n+1,self.tn+dt))

    def run(self,check_point=None):
        """
	Run the simulation by performing essentially::
    
	    self.doInitialization()
            while not self.terminateInitialIteration(): self.doInitialStep()
            self.doInitialPostprocessing()
	    while not self.finalize():
	        dt=self.getSafeTimeStepSize()
	        self.doStep(dt)
	        if n%check_point==0: 
		    self.writeXML() 
	    self.doFinalization()

        If one of the models in throws a C{FailedTimeStepError} exception a 
	new time step size is computed through getSafeTimeStepSize() and the 
	time step is repeated.
   
        If one of the models in throws a C{IterationDivergenceError} 
	exception the time step size is halved and the time step is repeated.

        In both cases the time integration is given up after
	C{Simulation.FAILED_TIME_STEPS_MAX} attempts.
        """
        # check the completness of the models:
        # first a list of all the models involved in the simulation including subsimulations:
        # 
        missing=self.checkModelLinks(self.getAllModels())
        if len(missing)>0:
            msg=""
            for l in missing:
                 msg+="\n\t"+str(l[-1])+" at "+str(self)
                 for i in xrange(len(l)-1): msg+="."+str(l[i])
            raise MissingLink("link targets missing in the Simulation: %s"%msg)
        #==============================
        self.doInitialization()
        self.doInitialStep()
        self.doInitialPostprocessing()
        dt=self.UNDEF_DT
        while not self.finalize():
            step_fail_counter=0
            iteration_fail_counter=0
            if self.n==0:
                dt_new=self.getSafeTimeStepSize(dt)
            else:
                dt_new=min(max(self.getSafeTimeStepSize(dt),dt/self.MAX_CHANGE_OF_DT),dt*self.MAX_CHANGE_OF_DT)
            self.trace("%d. time step %e (step size %e.)" % (self.n+1,self.tn+dt_new,dt_new))
            end_of_step=False
            while not end_of_step:
               end_of_step=True
               if not dt_new>0:
                  raise NonPositiveStepSizeError("non-positive step size in step %d"%(self.n+1))
               try:
                  self.doStepPreprocessing(dt_new)
                  self.doStep(dt_new)
                  self.doStepPostprocessing(dt_new)
               except IterationDivergenceError:
                  dt_new*=0.5
                  end_of_step=False
                  iteration_fail_counter+=1
                  if iteration_fail_counter>self.FAILED_TIME_STEPS_MAX:
                           raise SimulationBreakDownError("reduction of time step to achieve convergence failed after %s steps."%self.FAILED_TIME_STEPS_MAX)
                  self.trace("Iteration failed. Time step is repeated with new step size %s."%dt_new)
               except FailedTimeStepError:
                  dt_new=self.getSafeTimeStepSize(dt)
                  end_of_step=False
                  step_fail_counter+=1
                  self.trace("Time step is repeated with new time step size %s."%dt_new)
                  if step_fail_counter>self.FAILED_TIME_STEPS_MAX:
                        raise SimulationBreakDownError("Time integration is given up after %d attempts."%step_fail_counter)
            dt=dt_new
            if not check_point==None:
                if n%check_point==0: 
                    self.trace("check point is created.")
                    self.writeXML()
        self.doFinalization()


    def toDom(self, esysxml, node):
        """
	C{toDom} method of Simulation class.
	"""
        simulation = esysxml.createElement('Simulation')
        esysxml.registerLinkableObject(self, simulation)
        for rank, sim in enumerate(self.iterModels()):
            component = esysxml.createElement('Component')
            component.setAttribute('rank', str(rank))
            sim.toDom(esysxml, component)
            simulation.appendChild(component)
        node.appendChild(simulation)


    def fromDom(cls, esysxml, node):
        sims = []
        for n in node.childNodes:
            if isinstance(n, minidom.Text):
                continue
            sims.append(esysxml.getComponent(n))
        sims.sort(cmp=_comp)
        sim=cls([s[1] for s in sims], debug=esysxml.debug)
        esysxml.registerLinkableObject(sim, node)
        return sim

    fromDom = classmethod(fromDom)

def _comp(a,b):
    if a[0]<a[1]:
      return 1
    elif a[0]>a[1]:
      return -1
    else:
      return 0

class IterationDivergenceError(Exception):
    """
    Exception which is thrown if there is no convergence of the iteration 
    process at a time step.

    But there is a chance that a smaller step could help to reach convergence.
    """
    pass

class FailedTimeStepError(Exception):
    """
    Exception which is thrown if the time step fails because of a step 
    size that have been choosen to be too large.
    """
    pass

class SimulationBreakDownError(Exception):
    """
    Exception which is thrown if the simulation does not manage to 
    progress in time.
    """
    pass

class NonPositiveStepSizeError(Exception):
    """
    Exception which is thrown if the step size is not positive.
    """
    pass

class MissingLink(Exception):
    """
    Exception thrown when a link is missing
    """
    pass

class DataSource(object):
    """
    Class for handling data sources, including local and remote files. This class is under development.
    """

    def __init__(self, uri="file.ext", fileformat="unknown"):
        self.uri = uri
        self.fileformat = fileformat

    def toDom(self, esysxml, node):
        """
        C{toDom} method of DataSource. Creates a DataSource node and appends it to the
	current XML esysxml.
        """
        ds = esysxml.createElement('DataSource')
        ds.appendChild(esysxml.createDataNode('URI', self.uri))
        ds.appendChild(esysxml.createDataNode('FileFormat', self.fileformat))
        node.appendChild(ds)

    def fromDom(cls, esysxml, node):
        uri= str(node.getElementsByTagName("URI")[0].firstChild.nodeValue.strip())
        fileformat= str(node.getElementsByTagName("FileFormat")[0].firstChild.nodeValue.strip())
        ds = cls(uri, fileformat)
        return ds

    def getLocalFileName(self):
        return self.uri

    fromDom = classmethod(fromDom)
    
# vim: expandtab shiftwidth=4:
