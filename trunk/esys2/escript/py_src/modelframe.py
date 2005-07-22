# $Id$

from types import StringType,IntType,FloatType,BooleanType,ListType,DictType
from sys import stdout
import itertools

import modellib

# import the 'set' module if it's not defined (python2.3/2.4 difference)
try:
    set
except NameError:
    from sets import Set as set

from xml.dom import minidom

def dataNode(document, tagName, data):
    """
    dataNodes are the building blocks of the xml documents constructed in
    this module. document is the current xml document, tagName is the
    associated xml tag, and data is the values in the tag.
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

LinkableObjectRegistry = {}

def registerLinkableObject(obj_id, o):
    LinkableObjectRegistry[obj_id] = o

LinkRegistry = []

def registerLink(obj_id, l):
    LinkRegistry.append((obj_id,l))

def parse(xml):
    global LinkRegistry, LinkableObjectRegistry
    LinkRegistry = []
    LinkableObjectRegistry = {}

    doc = minidom.parseString(xml)

    sim = getComponent(doc.firstChild)
    for obj_id, link in LinkRegistry:
        link.target = LinkableObjectRegistry[obj_id]

    return sim

def getComponent(doc):
    for node in doc.childNodes:
        if isinstance(node, minidom.Element):
            if node.tagName == 'Simulation':
                if node.getAttribute("type") == 'Simulation':
                    return Simulation.fromDom(node)
                elif node.getAttribute("type") == 'ExplicitSimulation':
                    return ExplicitSimulation.fromDom(node)
            if node.tagName == 'Model':
                model_type = node.getAttribute("type")
                model_subclasses = Model.__subclasses__()
                for model in model_subclasses:
                    if model_type == model.__name__:
                        return model.fromDom(node)
                    
            raise "Invalid simulation type, %r" % node.getAttribute("type")

    raise ValueError("No Simulation Found")
            

class Link:
    """
    a Link makes an attribute of an object callable: 
          o.object()
          o.a=8
          l=Link(o,"a")
          assert l()==8
     """
    
    def __init__(self,target,attribute=None):
        """
        creates a link to the object target. If attribute is given, the link is
        establised to this attribute of the target.  Otherwise the attribute is
        undefined.
        """
        self.target = target
        self.attribute = None
        self.setAttributeName(attribute)
    
    def setAttributeName(self,attribute):
        """
        set a new attribute name to be collected from the target object. The
        target object must have the attribute with name attribute.
        """
        if attribute and self.target and not hasattr(self.target, attribute):
            raise AttributeError("%s: target %s has no attribute %s."%(self, self.target, attribute))
        self.attribute = attribute
    
    def hasDefinedAttributeName(self):
        """
        returns true if an attribute name is set
        """
        return self.attribute != None
    
    def __repr__(self):
        """
        returns a string representation of the link
        """
        if self.hasDefinedAttributeName():
            return "<Link to attribute %s of %s>" % (self.attribute, self.target)
        else:
            return "<Link to target %s>" % self.target
    
    def __call__(self,name=None):
        """
        returns the value of the attribute of the target object. If the
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
        toDom method of Link. Creates a Link node and appends it to the current XML 
        document 
        """
        link = document.createElement('Link')
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
        writes an XML representation of self to the output stream ostream.
        If ostream is nor present the standart output stream is used.  If
        esysheader==True the esys XML header is written
        """

        document, rootnode = esysDoc()
        self.toDom(document, rootnode)

        ostream.write(document.toprettyxml())

class LinkableObject(object):
    """
    An object that allows to link its attributes to attributes of other objects
    via a Link object. For instance
           
           p = LinkableObject()
           p.x = Link(o,"name")
           print p.x
    
    links attribute x of p to the attribute name of object o. 

    p.x will contain the current value of attribute name of object o.  

    If the value of getattr(o, "name") is callable, p.x will rturn the return
    value of the call. 
    """
   
    number_sequence = itertools.count(100)
    
    def __init__(self, debug=False):
        """ initializes LinkableObject so that we can operate on Links """
        self.debug = debug
        self.__linked_attributes={}
        self.id = self.number_sequence.next()

    def trace(self, msg):
        """ If debugging is on, print the message, otherwise do nothing
        """
        if self.debug:
            print msg
    
    def __getattr__(self,name):
        """returns the value of attribute name. If the value is a Link object the
        object is called and the return value is returned."""
        out = self.getAttributeObject(name)
        if isinstance(out,Link):
            return out()
        else:
            return out
    
    def getAttributeObject(self,name):
        """return the object stored for attribute name."""

        if self.__dict__.has_key(name):
            return self.__dict__[name]

        if self.__linked_attributes.has_key(name):
            return self.__linked_attributes[name]

        raise AttributeError,"No attribute %s."%name
    
    def __setattr__(self,name,value):
        """sets the value for attribute name. If value is a Link the target
        attribute is set to name if no attribute has been specified."""


        if self.__dict__.has_key(name): 
            del self.__dict__[name]

        if isinstance(value,Link):
            if not value.hasDefinedAttributeName(): 
                value.setAttributeName(name)
            self.__linked_attributes[name] = value

            self.trace("DEBUG: %s: attribute %s is now linked by %s."%(self,name,value))
        else:
            self.__dict__[name] = value
    
    def __delattr__(self,name):
        """removes the attribute name."""

        if self.__linked_attributes.has_key[name]:
            del self.__linked_attributes[name]
        elif self.__dict__.has_key(name):
            del self.__dict__[name]
        else:
            raise AttributeError,"No attribute %s."%name

class SimulationFrame(LinkableObject): 
    """A SimulationFrame objects represents a processess marching over time
    until a finalizing condition is fullfilled.  At each time step an iterative
    process can be performed and the time step size can be controlled
    """
    UNDEF_DT=1.e300
    MAX_TIME_STEP_REDUCTION=20
    MAX_ITER_STEPS=50
    
    def __init__(self,**kwargs):
        """
        Initialises a simulation
        
        Just calls the parent constructor.
        """
        LinkableObject.__init__(self,**kwargs)
    
    def doInitialization(self,t):
        """initializes the time stepping scheme. This function may be
        overwritten."""
        pass
    
    def getSafeTimeStepSize(self,dt):
        """returns a time step size which can savely be used. This function may
        be overwritten."""
        return self.UNDEF_DT
    
    def finalize(self):
        """returns True if the time stepping is finalized. This function may be
        overwritten."""
        return True
       
    def doFinalization(self):
        """finalizes the time stepping.  This function may be overwritten."""
        pass
    
    def doIterationInitialization(self,dt):
        """initializes the iteration at time step t. This function may be
        overwritten. (only called if doStep is not overwritten)"""
        pass
    
    def doIterationStep(self,dt):
        """executes the iteration step. This function may be overwritten. (only
        called if doStep is not overwritten)"""
        pass
    
    def terminate(self):
        """returns True if iteration on a time step is terminated. (only called
        if doStep is not overwritten)"""
        return True
       
    def doIterationFinalization(self,dt):
        """finalalizes the iteration process. (only called if doStep is not
        overwritten)"""
        pass
    
    def run(self,check_point=None):
        """run the simulation by performing essentially 
    
            self.doInitialization()
            while not self.finalize():
               dt=self.getSafeTimeStepSize()
               self.doStep(dt)
               if n%check_point==0: self.writeXML() 
            self.doFinalization()
    
        """
        self.__tn=0.
        self.__n=0
        self.__dt=None
        self.doInitialization(self.__tn)
        while not self.finalize():
            self.__n+=1
            self.__dt=self.getSafeTimeStepSize(self.__dt)
            if self.__dt==None: self.__dt=self.UNDEF_DT
            if not self.__dt>0:
                 raise NonPositiveStepSizeError("non-positive step size in step %d",self.__n)
            self.trace("%s: %d. time step %e (step size %e.)" % 
                          (self,self.__n,self.__tn+self.__dt,self.__dt))
            endoftimestep=False
            failcounter=0
            while not endoftimestep:
                endoftimestep=True 
                try:
                    self.doStep(self.__dt)
                except FailedTimeStepError:
                    self.__dt=self.getSafeTimeStepSize(self.__dt)
                    if self.__dt==None: self.__dt=self.UNDEF_DT
                    endoftimestep=False
                    self.trace("%s: time step is repeated with new step size %e."%(self,self.__dt))
                except IterationDivergenceError:
                    self.__dt*=0.5
                    endoftimestep=False
                    failcounter+=1
                    if failcounter>self.MAX_TIME_STEP_REDUCTION:
                        raise IterationBreakDownError("reduction of time step to achieve convergence failed.")

                    self.trace("%s: iteration failes. time step is repeated with new step size %e."
                        % (self,self.__dt))
            self.__tn+=self.__dt
            if not check_point==None:
                if self.__n%check_point==0: 
                    self.trace("%s: check point is created."%self)
                    self.writeXML()
        self.doFinalization()

    def writeXML(self):
        raise RuntimeError, "Not implemented"
    
    def doStep(self,dt):
        """executes a time step by iteration. This function may be overwritten.
           
           basicly it performs :
    
            self.doIterationInitialization(dt)
            while not self.terminate(): self.doIterationStep(dt)
            self.doIterationFinalization(dt)
    
        """
        self.doIterationInitialization(dt)
        self.__iter=0
        while not self.terminate(): 
            self.trace("%s: iteration step %d"%(self,self.__iter))
            self.doIterationStep(dt)
            self.__iter+=1
            if self.__iter>self.MAX_ITER_STEPS:
                raise IterationDivergenceError("%s: divergence in step %d"%(self,self.__iter))
        self.doIterationFinalization(dt)

class Simulation(SimulationFrame): 
    """A Simulation object is comprised by SimulationFrame(s) called subsimulations."""
    
    def __init__(self, subsimulations=[], **kwargs):
        """initiates a simulation from a list of subsimulations. """
        SimulationFrame.__init__(self, **kwargs)
        self.__subsimulations=[]

        for i in range(len(subsimulations)): 
            self[i] = subsimulations[i]
    
    def iterSubsimulations(self):
        """returns an iterator over the subsimulations"""
        return self.__subsimulations
    
    def __getitem__(self,i):
        """returns the i-th subsimulation"""
        return self.__subsimulations[i]
    
    def __setitem__(self,i,value):
        """sets the i-th subsimulation"""
        if not isinstance(value,SimulationFrame):
            raise ValueError("assigned value is not a Simulation")
        for j in range(max(i-len(self.__subsimulations)+1,0)): self.__subsimulations.append(None)
        self.__subsimulations[i]=value
    
    def __len__(self):
        """returns the number of subsimulations"""
        return len(self.__subsimulations)

    def toDom(self, document, node):
        """ toDom method of Simulation class """
        simulation = document.createElement('Simulation')
        simulation.setAttribute('type', self.__class__.__name__)

        for rank, sim in enumerate(self.iterSubsimulations()):
            component = document.createElement('Component')
            component.setAttribute('rank', str(rank))

            sim.toDom(document, component)

            simulation.appendChild(component)

        node.appendChild(simulation)

    def writeXML(self,ostream=stdout):
        """writes the object as an XML object into an output stream"""
        document, rootnode = esysDoc()
        self.toDom(document, rootnode)
        ostream.write(document.toprettyxml())
    
    def getSafeTimeStepSize(self,dt):
        """returns a time step size which can safely be used by all subsimulations"""
        out=self.UNDEF_DT
        for o in self.iterSubsimulations(): 
            dt_new = o.getSafeTimeStepSize(dt)
            if dt_new != None: 
                out = min(out,dt_new)

        return out
    
    def doInitialization(self,dt):
        """initializes all subsimulations """
        for o in self.iterSubsimulations(): 
            o.doInitialization(dt)
    
    def finalize(self):
        """returns True if all subsimulations are finalized"""
        return all([o.finalize() for o in self.iterSubsimulations()])
       
    def doFinalization(self):
        """finalalizes the time stepping for all subsimulations."""
        for i in self.iterSubsimulations(): i.doFinalization()
    
    def doIterationInitialization(self,dt):
        """initializes the iteration at time t for all subsimulations."""
        self.__iter=0
        self.trace("%s: iteration starts"%self)

        for o in self.iterSubsimulations(): 
            o.doIterationInitialization(dt)
    
    def terminate(self):
        """returns True if all iterations for all subsimulations are terminated."""
        return all([o.terminate() for o in self.iterSubsimulations()])
       
    def doIterationFinalization(self,dt):
        """finalalizes the iteration process for each of the subsimulations."""
        for o in self.iterSubsimulations(): 
            o.doIterationFinalization(dt)
        self.trace("%s: iteration finalized after %s steps"%(self,self.__iter+1))
    
    def doIterationStep(self,dt):
        """executes the iteration step at time step for each subsimulation"""
        self.__iter+=1
        self.trace("%s: iteration step %d"%(self,self.__iter))
        for o in self.iterSubsimulations(): 
            o.doIterationStep(dt)

    def fromDom(cls, doc):
        """
        Needs to be filled out.
        """
        sims = []
        for node in doc.childNodes:
            if isinstance(node, minidom.Text):
                continue

            sims.append(getComponent(node))


        return cls(sims)


            


    fromDom = classmethod(fromDom)

class ExplicitSimulation(Simulation): 
    """This is a modified form of the Simulation class. In fact it overwrites
    the doStep method by executing the doStep method of all subsimulations
    rather then iterating over all subsimulations."""

    def doStep(self,dt):
        """executes the time step for all subsimulation"""
        for i in self.iterSubsimulations(): 
            i.doStep(dt)

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
    """a class which allows to emphazise attributes to be written and read to XML
       
       Leaves of  an ESySParameters objects can be 
    
            a real number
            a integer number
            a string
            a boolean value
            a ParameterSet object 
            a Simulation object
            a Model object 
            any other object (not considered by writeESySXML and writeXML)
    
           Example how to create an ESySParameters object:
    
                 p11=ParameterSet(gamma1=1.,gamma2=2.,gamma3=3.)
                 p1=ParameterSet(dim=2,tol_v=0.001,output_file="/tmp/u.%3.3d.dx",runFlag=True,parm11=p11) 
                 parm=ParameterSet(parm1=p1,parm2=ParameterSet(alpha=Link(p11,"gamma1")))
    
           This can be accessed as 
    
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
        """creates a ParameterSet with parameters parameters"""
        LinkableObject.__init__(self, **kwargs)
        self.parameters = set()
        self.declareParameters(parameters)
    
    def declareParameter(self,**parameters):
        """declares a new parameter(s) and its (their) inital value."""
        self.declareParameters(parameters)
    
    def declareParameters(self,parameters):
        """declares a set of parameters. parameters can be a list, a dictonary or a ParameterSet."""
        if isinstance(parameters,ListType):
            parameters = zip(parameters, itertools.repeat(None))
        if isinstance(parameters,DictType):
            parameters = parameters.iteritems()

        for prm, value in parameters:
            setattr(self,prm,value)
            self.parameters.add(prm)

            self.trace("%s: parameter %s has been declared."%(self,prm))

    def releaseParameters(self,name):
        """removes parameter name from the paramameters"""
        if self.isParameter(name): 
            self.parameters.remove(name)
            self.debug("%s: parameter %s has been removed."%(self, name))
    
    def __iter__(self):
        """creates an iterator over the parameter and their values"""
        return _ParameterIterator(self)
    
    def showParameters(self):
        """returns a descrition of the parameters"""        
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
        """removes the attribute name."""
        LinkableObject.__delattr__(self,name)
        try:
            self.releaseParameter(name) 
        except:
            pass

    def toDom(self, document, node):
        """ toDom method of ParameterSet class """
        pset = document.createElement('ParameterSet')
        node.appendChild(pset)
        self._parametersToDom(document, pset)

    def _parametersToDom(self, document, node):
        node.setAttribute ('id', str(self.id))
        for name,value in self:
            param = document.createElement('Parameter')
            param.setAttribute('type', value.__class__.__name__)

            param.appendChild(dataNode(document, 'Name', name))

            val = document.createElement('Value')

            if isinstance(value,ParameterSet):
                value.toDom(document, val)
                param.appendChild(val)
            elif isinstance(value, Link):
                value.toDom(document, val)
                param.appendChild(val)
            elif isinstance(value,StringType):
                param.appendChild(dataNode(document, 'Value', value))
            else:
                param.appendChild(dataNode(document, 'Value', str(value)))

            node.appendChild(param)

   
    def fromDom(cls, doc):

        # Define a host of helper functions to assist us.
        def _children(node):
            """
            Remove the empty nodes from the children of this node
            """
            return [x for x in node.childNodes 
                    if not isinstance(x, minidom.Text) or x.nodeValue.strip()]

        def _floatfromValue(doc):
            return float(doc.nodeValue.strip())

        def _stringfromValue(doc):
            return str(doc.nodeValue.strip())
       
        def _intfromValue(doc):
            return int(doc.nodeValue.strip())

        def _boolfromValue(doc):
            return bool(doc.nodeValue.strip())
       
        # Mapping from text types in the xml to methods used to process trees of that type
        ptypemap = {"Simulation": Simulation.fromDom,
                    "Model":Model.fromDom,
                    "ParameterSet":ParameterSet.fromDom,
                    "Link":Link.fromDom,
                    "float":_floatfromValue,
                    "int":_intfromValue,
                    "str":_stringfromValue,
                    "bool":_boolfromValue
                    }

        parameters = {}
        for node in _children(doc):
            ptype = node.getAttribute("type")

            pname = pvalue = None
            for childnode in _children(node):

                if childnode.tagName == "Name":
                    pname = childnode.firstChild.nodeValue.strip()

                if childnode.tagName == "Value":
                    nodes = _children(childnode)
                    pvalue = ptypemap[ptype](nodes[0])

            parameters[pname] = pvalue

        # Create the instance of ParameterSet
        o = cls()
        o.declareParameters(parameters)
        registerLinkableObject(doc.getAttribute("id"), o)
        return o
    
    fromDom = classmethod(fromDom)
    
    def writeXML(self,ostream=stdout):
        """writes the object as an XML object into an output stream"""
        # ParameterSet(d) with d[Name]=Value
        document, node = esysDoc()
        self.toDom(document, node)
        ostream.write(document.toprettyxml())

class Model(ParameterSet,SimulationFrame):
    """a Model is a SimulationFrame which is also a ParameterSet.""" 

    def __init__(self,parameters=[],**kwargs):
        """creates a model"""
        ParameterSet.__init__(self, parameters=parameters)
        SimulationFrame.__init__(self,**kwargs)

    def toDom(self, document, node):
        """ toDom method of Model class """
        pset = document.createElement('Model')
        pset.setAttribute('type', self.__class__.__name__)
        node.appendChild(pset)
        self._parametersToDom(document, pset)


class IterationDivergenceError(Exception):
    """excpetion which should be thrown if there is no convergence of the iteration process at a time step but there is a chance taht a smaller step could help
       to reach convergence."""
    pass

class IterationBreakDownError(Exception):
    """excpetion which should be thrown if there is no conevregence and there is no chance that a time step reduction would help"""
    pass

class FailedTimeStepError(Exception):
    """excpetion which should be thrown if the time step fails because of a step size that have been choosen to be to large"""
    pass

class NonPositiveStepSizeError(Exception):
    """excpetion which is thrown if the step size is not positive"""
    pass

#
#   ignore this text:
#
""" 
the Model class provides a framework to run a time-dependent simulation. A
Model has a set of parameter which may be fixed or altered by the Model itself
or other Models over time.   

       The parameters of a models are declared at instantion, e.g.

           m=Model({"message" : "none" })

       creates a Model with parameters p1 and p2 with inital values 1 and 2.
       Typically a particular model is defined as a subclass of Model:

        class Messenger(Model):
            def __init__(self):
               Model.__init__(self,parameters={"message" : "none" })

        m=MyModel()

       There are various ways how model parameters can be changed:

       1) use object attributes:

          m.message="Hello World!"

       2) use setParamter method
        
          
          m.setParameters(message="Hello World!")

       3) or dictonaries
  
           d={ message : "Hello World!" }
           m.setParameters(**d)


       A model executed buy staring the run method of the model:
 
          m=Messenger()
          m.run()

       The run methods marches through time. It first calls the
       doInitialization() method of the Model to set up the process. In each
       time step the doStep() method is called to get from the current to the
       next time step. The step size is defined by calling the
       getSafeTimeStepSize() method.  the time integration process is
       terminated when the finalize() methods return true. Final the
       doFinalization() method is called to finalize the process. To implement
       a particular model a subclass of the Model class is defined. The
       subclass overwrites the default methods of Model. 

       The following class defines a messenger printing in the doStep method
       what ever the current value of its parameter message is:

       class Messenger(Model):
            def __init__(self):
               Model.__init__(self,parameters={"message" : "none" })

            def doInitialization(self):
               print "I start talking now!"

            def doStep(self,t):
               print "Message (time %e) : %s "%(t,self.message)

            def doFinalization(self):
               print "I have no more to say!"
       
       If a instance of the Messenger class is run, it will print the
       initialization and finalization message only.  This is because the
       default method for finalize() does always returns True and therefore the
       transition is terminated startcht away. 
        
       Following example for solving the ODE using a forward euler scheme:

                u(t=0)=u0
                u_t=a*u**2       for all 0<t<=ten

      exact solution is given by u(t)=1/(1/u0-a*t)

      class  Ode1(Model):
         def __init__(self,**args):
            Model.__init__(self,parameters={"tend" : 1., "dt" : 0.0001 ,"a" : 0.1 ,"u" : 1. },name="test",debug=True)

         def doInitialization(self):
             self._tn=0

         def doStep(self,t):
             self.u=self.u+(t-self._tn)*self.a*self.u**2
             self._tn=t

         def doFinalization(self):
             print "all done final error = ",abs(self.u-1./(1./3.-self.a*self._tn))

         def getSafeTimeStepSize(self):
             return self.dt

         def finalize(self):
             return self._tn>=self.tend

       In some cases at a given time step an iteration process has to be
       performed to get the state of the Model for the next time step. ` In
       this case the doStep() method is replaced by a sequance of methods which
       implements this iterative process.  The method then will control the
       iteration process by initializing the iteration through calling the
       doIterationInitialization() method. The iteration is preformed by
       calling the doIterationStep() method until the terminate() method
       returns True. The doIterationFinalization() method is called to end the
       iteration. 
       For a particular model these methods have to overwritten by a suitable
       subclass without touching the doStep() method.

       following example is a modification of the example above. Here an
       implicit euler scheme is used. in each time step the problem
           
           0= u_{n+1}-u_{n}+a*dt*u_{n+1}**2

       has to be solved for u_{n+1}. The Newton scheme is used to solve this non-linear problem.

 
      class  Ode2(Model):

       def __init__(self,**args):
           Model.__init__(self,{"tend" : 1., "dt" : 0.1 ,"a" : 10. ,"u" : 1. , "tol " : 1.e-8},"test","bla",None,True)

       def doInitialization(self):
           self.__tn=0

       def doIterationInitialization(self,t):
            self.__iter=0
            self.u_last=self.u            
            self.current_dt=t-self.tn
            self.__tn=t

       def doIterationStep(self):
          self.__iter+=1
          self.u_old=self.u
          self.u=(self.current_dt*self.a*self.u**2-self.u_last)/(2*self.current_dt*self.a*self.u-1.)

       def terminate(self):
           return abs(self.u_old-self.u)<self.tol*abs(self.u)

       def doIterationFinalization(self)
           print "all done"

       def getSafeTimeStepSize(self):
           return self.dt

       def finalize(self):
            return self.__tn>self.tend

       A model can be composed from subsimulations. Subsimulations are treated
       as model parameters. If a model parameter is set or a value of a model
       parameter is requested, the model will search for this parameter its
       subsimulations in the case the model does not have this parameter
       itself. The order in which the subsimulations are searched is critical.
       By default a Model initializes all its subsimulations, is finalized when
       all its subsimulations are finalized and finalizes all its
       subsimulations. In the case an iterative process is applied on a
       particular time step the iteration is initialized for all
       subsimulations, then the iteration step is performed for each
       subsimulation until all subsimulations indicate termination. Then the
       iteration is finalized for all subsimulations. Finally teh doStop()
       method for all submethods is called. 

       Here we are creating a model which groups ab instantiation of the Ode2 and the Messenger Model

       o=Ode2()
       m=Messenger()
       om=Model(subsimulations=[o,m],debug=True)
       om.dt=0.01
       om.u=1.
       m.message="it's me!"
       om.run()

       Notice that dt and u are parameters of class Ode2 and message is a
       parameter of the Messenger class. The Model formed from these models
       automatically hand the assignment of new values down to the
       subsimulation. om.run() starts this combined model where now the
       soStep() method of the Messenger object printing the value of its
       parameter message together with a time stamp is executed in each time
       step introduced by the Ode2 model. 

       A parameter of a Model can be linked to an attribute of onother object,
       typically an parameter of another Model object.     
       
       
       which is comprised by a set of subsimulations. 
       The simulation is run through its run method which in the simplest case has the form:

          s=Model()
          s.run()

       The run has an initializion and finalization phase. The latter is called
       if all subsimulations are to be finalized. The simulation is processing
       in time through calling the stepForward methods which updates the
       observables of each subsimulation.  A time steps size which is save for
       all subsimulation is choosen. 

       At given time step an iterative process may be performed to make sure
       that all observables are consistent across all subsimulations.  In this
       case, similar the time dependence, an initialization and finalization of
       the iteration is performed. 

       A Model has input and output parameters where each input parameter can
       be constant, time dependent or may depend on an output parameter of
       another model or the model itself. To create a parameter name of a model
       and to assign a value to it one can use the statement

           model.name=object


       At any time the current value of the parameter name can be obtained by

               value=model.name

       If the object that has been assigned to the paramter/attribute name has
       the attribute/parameter name isself the current value of this attribute
       of the object is returned (e.g. for model.name=object where object has
       an attribute name, the statement value=model.name whould assign the
       value object.name to value.). If the name of the parameters of a model
       and an object don't match the setParameter method of model can be used.
       So 

           model.setParameter(name,object,name_for_object)

       links the parameter name of model with the parameter name_for_object of
       object.

       The run method initiates checkpointing (it is not clear how to do this
       yet)
   =====
            
   """



if __name__=="__main__":
    import math
    #
    # test for parameter set
    #
    p11=ParameterSet()
    p11.declareParameter(gamma1=1.,gamma2=2.,gamma3=3.)
    p1=ParameterSet()
    p1.declareParameter(dim=2,tol_v=0.001,output_file="/tmp/u.%3.3d.dx",runFlag=True,parm11=p11) 
    parm=ParameterSet({ "parm1" : p1 , "parm2" : ParameterSet(["alpha"])})
    parm.parm2.alpha=Link(p11,"gamma1")
    parm.x="that should not be here!"
    print parm.showParameters()
    # should be something like: {"parm2" : {"alpha" : reference to attribute
    # gamma1 of <__main__.ParameterSet instance at 0xf6db51cc>},"parm1" : {"dim"
    # : 2,"runFlag" : True,"tol_v": 0.001,"parm11" : {"gamma3" : 3.0,"gamma2" :
    # 2.0,"gamma1" : 1.0},"output_file" : /tmp/u.%3.3d.dx}}
    assert parm.parm2.alpha==1.
    parm.writeXML()
    
    #=======================
    class Messenger(Model):
        def __init__(self):
            Model.__init__(self)
            self.declareParameter(message="none")

        def doInitialization(self,t):
            self.__t=t
            print "I start talking now!"

        def doStep(self,dt):
            self.__t+=dt
            print "Message (time %e) : %s "%(self.__t,self.message)

        def doFinalization(self):
            print "I have no more to say!"
   
    class ODETEST(Model):
        """ implements a solver for the ODE 

              du/dt=a*u+f(t)

           we use a implicit euler scheme :

              u_n-u_{n-1}= dt*a u_n + st*f(t_n)

           to get u_n we run an iterative process
 
               u_{n.k}=u_{n-1}+dt*(a u_{n.i-1} + f(t_n))


           input for this model are step size dt, end time tend and a value for
           a, f and initial value for u. we need also a tolerance tol for a
           stopping criterion.

        """

        def __init__(self):
            Model.__init__(self,debug=True)
            self.declareParameter(tend=1.,dt=0.1,a=0.9,u=10.,f=0.,message="",tol=1.e-8)

        def doInitialization(self,t):
            self.__tn=t
            self.__iter=0
           
        def doIterationInitialization(self,dt):
            self.__iter=0
            self.__u_last=self.u            

        def doIterationStep(self,dt):
            self.__iter+=1
            self.__u_old=self.u
            self.u=self.__u_last+dt*(self.a*self.__u_old+self.f)
         
        def terminate(self):
            if self.__iter<1:
                return False
            else:
                return abs(self.__u_old-self.u)<self.tol*abs(self.u)

        def doIterationFinalization(self,dt):
            self.__tn+=dt
            self.message="current error = %e"%abs(self.u-10.*math.exp((self.a-1.)*self.__tn))

        def getSafeTimeStepSize(self,dt):
            return min(self.dt,1./(abs(self.a)+1.))

        def finalize(self):
            return self.__tn>=self.tend

    #
    #   s solves the coupled ODE:
    #
    #     du/dt=a*u+  v 
    #     dv/dt=  u+a*v
    #  
    #   each equation is treated through the ODETEST class. The equations are
    #   linked and iteration over each time step is performed.  the current
    #   error of v is reported by the Messenger class.
    #
    o1=ODETEST()
    o1.u=10
    o2=ODETEST()
    o2.u=-10.
    o1.f=Link(o2,"u")
    o2.f=Link(o1,"u")
    m=Messenger()
    o1.dt=0.01
    m.message=Link(o1)
    s=ExplicitSimulation([Simulation([o1,o2],debug=True),m],debug=True)
    s.run() 
    s.writeXML()
