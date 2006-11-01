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
# import modellib  temporarily removed!!!

# import the 'set' module if it's not defined (python2.3/2.4 difference)
try:
    set
except NameError:
    from sets import Set as set

from xml.dom import minidom

def dataNode(document, tagName, data):
    """
    C{dataNode}s are the building blocks of the xml documents constructed in
    this module.  
    
    @param document: the current xml document
    @param tagName: the associated xml tag
    @param data: the values in the tag
    """
    t = document.createTextNode(str(data))
    n = document.createElement(tagName)
    n.appendChild(t)
    return n

def esysDoc():
    """
    Global method for creating an instance of an EsysXML document.
    """
    doc = minidom.Document()
    esys = doc.createElement('ESys')
    doc.appendChild(esys)
    return doc, esys

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

LinkableObjectRegistry = {}

def registerLinkableObject(obj_id, o):
    LinkableObjectRegistry[obj_id] = o

LinkRegistry = []

def registerLink(obj_id, l):
    LinkRegistry.append((obj_id,l))

def parse(xml):
    """
    Generic parse method for EsysXML.  Without this, Links don't work.
    """
    global LinkRegistry, LinkableObjectRegistry
    LinkRegistry = []
    LinkableObjectRegistry = {}

    doc = minidom.parseString(xml)
    sim = getComponent(doc.firstChild)
    for obj_id, link in LinkRegistry:
        link.target = LinkableObjectRegistry[obj_id]

    return sim

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

def getComponent(doc):
    """ 
    Used to get components of Simualtions, Models.
    """
    for node in doc.childNodes:
        
        if isinstance(node, minidom.Element):
            if node.tagName == 'Simulation':
                if node.getAttribute("type") == 'Simulation':
                    return Simulation.fromDom(node)
            if node.tagName == 'Model':
                if (node.getAttribute("module")):
                    model_module = node.getAttribute("module")
                    model_type = node.getAttribute("type")
                    return importName(model_module, model_type).fromDom(node)
                else:
                    model_type = node.getAttribute("type")
                    model_subclasses = Model.__subclasses__()
                    for model in model_subclasses:
                        if model_type == model.__name__:
                            return Model.fromDom(node)
            if node.tagName == 'ParameterSet':
                parameter_type = node.getAttribute("type")
                return ParameterSet.fromDom(node)
            raise "Invalid simulation type, %r" % node.getAttribute("type")
        

    raise ValueError("No Simulation Found")
            

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

    def toDom(self, document, node):
        """
        C{toDom} method of Link. Creates a Link node and appends it to the
	current XML document.
        """
        link = document.createElement('Link')
        assert (self.target != None), ("Target was none, name was %r" % self.attribute)
        link.appendChild(dataNode(document, 'Target', self.target.id))
        # this use of id will not work for purposes of being able to retrieve the intended
        # target from the xml later. I need a better unique identifier.
        assert self.attribute, "You can't xmlify a Link without a target attribute"
        link.appendChild(dataNode(document, 'Attribute', self.attribute))
        node.appendChild(link)

    def fromDom(cls, doc):
        targetid = doc.getElementsByTagName("Target")[0].firstChild.nodeValue.strip()
        attribute = doc.getElementsByTagName("Attribute")[0].firstChild.nodeValue.strip()
        l = cls(None, attribute)
        registerLink(targetid, l)
        return l

    fromDom = classmethod(fromDom)
    
    def writeXML(self,ostream=stdout):
        """
        Writes an XML representation of self to the output stream ostream.
        If ostream is nor present the standart output stream is used.  If
        esysheader==True the esys XML header is written
        """
        print 'I got to the Link writeXML method'
        document, rootnode = esysDoc()
        self.toDom(document, rootnode)

        ostream.write(document.toprettyxml())

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
   
    number_sequence = itertools.count(100)
    
    def __init__(self, debug=False):
        """
	Initializes LinkableObject so that we can operate on Links 
	"""
        self.debug = debug
        self.__linked_attributes={}
        self.id = self.number_sequence.next()
        registerLinkableObject(self.id, self)

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
        return "<%s %r>" % (self.__class__.__name__, 
                            [(p, getattr(self, p, None)) for p in self.parameters])
    
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

    def toDom(self, document, node):
        """
	C{toDom} method of ParameterSet class.
	"""
        pset = document.createElement('ParameterSet')
        node.appendChild(pset)
        self._parametersToDom(document, pset)

    def _parametersToDom(self, document, node):
        node.setAttribute('id', str(self.id))
        node.setIdAttribute("id")
        for name,value in self:
            param = document.createElement('Parameter')
            param.setAttribute('type', value.__class__.__name__)

            param.appendChild(dataNode(document, 'Name', name))

            val = document.createElement('Value')

            if isinstance(value,(ParameterSet,Link,DataSource)):
                value.toDom(document, val)
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
                numarrayElement = document.createElement('NumArray')
                numarrayElement.appendChild(dataNode(document, 'ArrayType', str(arraytype)))
                numarrayElement.appendChild(dataNode(document, 'Shape', shape))
                numarrayElement.appendChild(dataNode(document, 'Data', ' '.join(
                    [str(x) for x in numarray.reshape(value, size)])))
                val.appendChild(numarrayElement)
                param.appendChild(val)
            elif isinstance (value, list):
                param.appendChild(dataNode(document, 'Value', ' '.join(
                    [str(x) for x in value]) 
                ))
            else:
                param.appendChild(dataNode(document, 'Value', str(value)))

            node.appendChild(param)

    def fromDom(cls, doc):

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

        def _floatfromValue(doc):
            return float(doc.nodeValue.strip())

        def _stringfromValue(doc):
            return str(doc.nodeValue.strip())
       
        def _intfromValue(doc):
            return int(doc.nodeValue.strip())

        def _boolfromValue(doc):
            return _boolfromstring(doc.nodeValue.strip())

        def _nonefromValue(doc):
            return None

        def _numarrayfromValue(doc):
            for node in _children(doc):
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
      
        def _listfromValue(doc):
            return [_boolfromstring(x) for x in doc.nodeValue.split()]


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

#        print doc.toxml()

        parameters = {}
        for node in _children(doc):
            ptype = node.getAttribute("type")

            pname = pvalue = None
            for childnode in _children(node):

                if childnode.tagName == "Name":
                    pname = childnode.firstChild.nodeValue.strip()

                if childnode.tagName == "Value":
                    nodes = _children(childnode)
                #    if ptype == 'NumArray':
                 #       pvalue = _numarrayfromValue(nodes)
                 #   else:
                    pvalue = ptypemap[ptype](nodes[0])

            parameters[pname] = pvalue

        # Create the instance of ParameterSet
        o = cls()
        o.declareParameters(parameters)
        registerLinkableObject(doc.getAttribute("id"), o)
        return o
    
    fromDom = classmethod(fromDom)
    
    def writeXML(self,ostream=stdout):
        """
	Writes the object as an XML object into an output stream.
	"""
        # ParameterSet(d) with d[Name]=Value
        document, node = esysDoc()
        self.toDom(document, node)
        ostream.write(document.toprettyxml())

class Model(ParameterSet):
    """
    A Model object represents a processess marching over time until a
    finalizing condition is fullfilled. At each time step an iterative
    process can be performed and the time step size can be controlled. A
    Model has the following work flow::

          doInitialization()
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

    def __init__(self,parameters=[],**kwarg):
        """
	Creates a model.

        Just calls the parent constructor.
        """
        ParameterSet.__init__(self, parameters=parameters,**kwarg)

    def __str__(self):
       return "<%s %d>"%(self.__class__,id(self))

    def toDom(self, document, node):
        """
	C{toDom} method of Model class
	"""
        pset = document.createElement('Model')
        pset.setAttribute('type', self.__class__.__name__)
        if not self.__class__.__module__.startswith('esys.escript'):
            pset.setAttribute('module', self.__class__.__module__)
        node.appendChild(pset)
        self._parametersToDom(document, pset)

    def doInitialization(self):
        """
	Initializes the time stepping scheme.  
	
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
	Finalalizes the time step.

        dt is the currently used time step size.

        This function may be overwritten.
	"""
        pass
    
    def writeXML(self, ostream=stdout):
        document, node = esysDoc()
        self.toDom(document, node)
        ostream.write(document.toprettyxml())
    

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
        Model.__init__(self, **kwargs)
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

    def toDom(self, document, node):
        """
	C{toDom} method of Simulation class.
	"""
        simulation = document.createElement('Simulation')
        simulation.setAttribute('type', self.__class__.__name__)

        for rank, sim in enumerate(self.iterModels()):
            component = document.createElement('Component')
            component.setAttribute('rank', str(rank))

            sim.toDom(document, component)

            simulation.appendChild(component)

        node.appendChild(simulation)

    def writeXML(self,ostream=stdout):
        """
	Writes the object as an XML object into an output stream.
	"""
        document, rootnode = esysDoc()
        self.toDom(document, rootnode)
        targetsList = document.getElementsByTagName('Target')
        
        for element in targetsList:
            targetId = int(element.firstChild.nodeValue.strip())
            if document.getElementById(str(targetId)):
                continue
            targetObj = LinkableObjectRegistry[targetId]
            targetObj.toDom(document, rootnode)
        ostream.write(document.toprettyxml())
    
    def getSafeTimeStepSize(self,dt):
        """
	Returns a time step size which can safely be used by all models.

        This is the minimum over the time step sizes of all models.
	"""
        out=min([o.getSafeTimeStepSize(dt) for o in self.iterModels()])
        #print "%s: safe step size is %e."%(str(self),out)
        return out
    
    def doInitialization(self):
        """
	Initializes all models.
	"""
        self.n=0
        self.tn=0.
        for o in self.iterModels(): 
            o.doInitialization()
    
    def finalize(self):
        """
	Returns True if any of the models is to be finalized.
	"""
        return any([o.finalize() for o in self.iterModels()])
       
    def doFinalization(self):
        """
	Finalalizes the time stepping for all models.
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
       
    def doStepPostprocessing(self,dt):
        """
	Finalalizes the iteration process for all models.
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
        dt=self.UNDEF_DT
        self.doInitialization()
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

    def fromDom(cls, doc):
        sims = []
        for node in doc.childNodes:
            if isinstance(node, minidom.Text):
                continue

            sims.append(getComponent(node))

        return cls(sims)

    fromDom = classmethod(fromDom)


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

class DataSource(object):
    """
    Class for handling data sources, including local and remote files. This class is under development.
    """

    def __init__(self, uri="file.ext", fileformat="unknown"):
        self.uri = uri
        self.fileformat = fileformat

    def toDom(self, document, node):
        """
        C{toDom} method of DataSource. Creates a DataSource node and appends it to the
	current XML document.
        """
        ds = document.createElement('DataSource')
        ds.appendChild(dataNode(document, 'URI', self.uri))
        ds.appendChild(dataNode(document, 'FileFormat', self.fileformat))
        node.appendChild(ds)

    def fromDom(cls, doc):
        uri= doc.getElementsByTagName("URI")[0].firstChild.nodeValue.strip()
        fileformat= doc.getElementsByTagName("FileFormat")[0].firstChild.nodeValue.strip()
        ds = cls(uri, fileformat)
        return ds

    def getLocalFileName(self):
        return self.uri

    fromDom = classmethod(fromDom)
    
# vim: expandtab shiftwidth=4:
