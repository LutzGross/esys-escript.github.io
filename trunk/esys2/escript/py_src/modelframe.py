# $Id$
from types import StringType

class Link:
  """ """
  def __init__(self,object,attribute=None):
     self.__object=object
     self.setAttributeName(attribute)

  def setAttributeName(self,name):
     if not name==None:
        if not hasattr(self.__object,name):
           raise AttributeError("Link: object %s has no attribute %s."%(self.__object,name))
     self.__attribute=name

  def hasAttributeName(self):
      if self.__attribute==None:
         return False
      else:
         return True

  def __str__(self):
      if self.hasAttributeName():
          return "reference to %s of %s"%(self.__attribute,self.__object)
      else:
          return "reference to object %s"%self.__object

  def getValue(self,name=None):
      if not self.hasAttributeName():
         out=getattr(self.__object,name)
      else:
         out=getattr(self.__object,self.__attribute)
      if callable(out):
          return out()
      else:
          return out
   
class Model: 
   """ the Model class provides a framework to run a time-dependent simulation. A Model has a set of parameter which 
       may be fixed or altered by the Model itself or other Models over time.   

       The parameters of a models are declared at instantion, e.g.

           m=Model({"message" : "none" })

       creates a Model with parameters p1 and p2 with inital values 1 and 2. Typically a particular model is defined as a subclass of Model:

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
       doInitialization() method of the Model to set up the process. In each time step the doStep() method is called 
       to get from the current to the next time step. The step size is defined by calling the getSafeTimeStepSize() method. 
       The time integration process is terminated when the finalize() methods return true. Final the doFinalization() method 
       is called to finalize the process. To implement a particular model a subclass
       of the Model class is defined. The subclass overwrites the default methods of Model. 

       The following class defines a messenger printing in the doStep method what ever the current value of its parameter message is:

       class Messenger(Model):
            def __init__(self):
               Model.__init__(self,parameters={"message" : "none" })

            def doInitialization(self):
               print "I start talking now!"

            def doStep(self,t):
               print "Message (time %e) : %s "%(t,self.message)

            def doFinalization(self):
               print "I have no more to say!"
       
       If a instance of the Messenger class is run, it will print the initialization and finalization message only. 
       This is because the default method for finalize() does always returns True and therefore the transition is 
       terminated startcht away. 
        
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

       In some cases at a given time step an iteration process has to be performed to get the state of the Model for the next time step. `
       In this case the doStep() method is replaced by a sequance of methods which implements this iterative process.
       The method then will control the iteration process by initializing the iteration through calling the 
       doIterationInitialization() method. The iteration is preformed by calling the doIterationStep() method until 
       the terminate() method returns True. The doIterationFinalization() method is called to end the iteration. 
       For a particular model these methods have to overwritten by a suitable subclass without touching the doStep() method.

       following example is a modification of the example above. Here an implicit euler scheme is used. in each time step the problem
           
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

       A model can be composed from submodels. Submodels are treated as model parameters. If a model parameter is set or a value of 
       a model parameter is requested, the model will search for this parameter its submodels in the case the model does not have this
       parameter itself. The order in which the submodels are searched is critical. By default a Model initializes all its submodels, 
       is finalized when all its submodels are finalized and finalizes all its submodels. In the case an iterative process is applied
       on a particular time step the iteration is initialized for all submodels, then the iteration step is performed for each submodel
       until all submodels indicate termination. Then the iteration is finalized for all submodels. Finally teh doStop() method for all 
       submethods is called. 

       Here we are creating a model which groups ab instantiation of the Ode2 and the Messenger Model

       o=Ode2()
       m=Messenger()
       om=Model(submodels=[o,m],debug=True)
       om.dt=0.01
       om.u=1.
       m.message="it's me!"
       om.run()

       Notice that dt and u are parameters of class Ode2 and message is a parameter of the Messenger class. The Model formed from these models
       automatically hand the assignment of new values down to the submodel. om.run() starts this combined model where now the soStep() method
       of the Messenger object printing the value of its parameter message together with a time stamp is executed in each time step introduced 
       by the Ode2 model. 

       A parameter of a Model can be linked to an attribute of onother object, typically an parameter of another Model object.     
       
       
       which is comprised by a set of submodels. 
       The simulation is run through its run method which in the simplest case has the form:

          s=Model()
          s.run()

       The run has an initializion and finalization phase. The latter is called if all submodels are to be finalized. The 
       simulation is processing in time through calling the stepForward methods which updates the observables of each submodel. 
       A time steps size which is save for all submodel is choosen. 

       At given time step an iterative process may be performed to make sure that all observables are consistent across all submodels.
       In this case, similar the time dependence, an initialization and finalization of the iteration is performed. 

       A Model has input and output parameters where each input parameter can be constant, time dependent or may depend on an 
       output parameter of another model or the model itself. To create a parameter name of a model and to 
       assign a value to it one can use the statement

           model.name=object


       At any time the current value of the parameter name can be obtained by

               value=model.name

       If the object that has been assigned to the paramter/attribute name has the attribute/parameter name isself the current value of this 
       attribute of the object is returned (e.g. for model.name=object where object has an attribute name, the statement value=model.name whould assign 
       the value object.name to value.). If the name of the parameters of a model and an object don't match the setParameter method of model can be used. So 

           model.setParameter(name,object,name_for_object)

       links the parameter name of model with the parameter name_for_object of object.

       The run method initiates checkpointing (it is not clear how to do this yet)
   =====
            
   """
   # step size used in case of an undefined value for the step size
   UNDEF_DT=1.e300

   def __init__(self,submodels=[],parameters={},name="model",description="none",check_pointing=None,debug=False):
      """initiates a model from a list of submodels. """
      self.setDebug(debug)
      self.__check_pointing=check_pointing
      self.__parameters={}
      self.setName(name)
      self.setDescription(description)
      self.declareParameter(**parameters)
      # get the models defined in parameters:
      self.__submodels=[]
      # submodels==None means no submodels used:
      if submodels==None:
         pass 
      # no submodel list given means all submodels are used as defined by the parameters dictionary:
      elif len(submodels)==0:
            for i in parameters.keys():
                if isinstance(parameters[i],Model): self.__submodels.append(i)
      # submodel list of strings and Models is given, submodels defines the order in which the 
      # submodels are processed. if new models are found in the list they are added to the parameter dictionary.
      else:
         c=0
         for i in submodels:
            if isinstance(i,StringType):
              m=self.getParameter(i)
              if not isinstance(m,Model):
                 raise ValueError,"submodel %s is not a model."%i
            else:
               if not isinstance(i,Model):
                 raise ValueError,"submodel list does contain item which is not a Model class object."
               m=i
               i="__submodel%d__"%c
               self.declareParameter(**{i : m})
               c+=1
            self.__submodels.append(i)
            if self.debug(): print "%s: model %s is added as parameter %s."%(self,m,i)
      if len(self.__submodels)>0 and self.debug(): print "%s: model ordering is %s"%(self,self.__submodels) 
   def setSubmodelOrder(submodels=[]):
      """sets a new ordering for submodels"""
      
     
   #
   # some basic fuctions:
   #
   def debugOn(self):
      """sets debugging to on"""
      self.__debug=True
   def debugOff(self):
      """sets debugging to off"""
      self.__debug=False
   def debug(self):
      """returns True if debug mode is set to on"""
      return self.__debug
   def setDebug(self,flag=False):
      """sets debugging to flag"""
      if flag:
         self.debugOn()
      else:
         self.debugOff()
   def setDebug(self,flag=False):
      """sets debugging to flag"""
      if flag:
         self.debugOn()
      else:
         self.debugOff()
   # name and description handling
   def __str__(self):
       """returns the name of the model"""
       return self.getName()

   def getName(self):
       """returns the name of the model"""
       return self.__name

   def getFullName(self):
       """returns the full name of the model including all the names of the submodels"""
       out=str(self)+"("
       notfirst=False
       for i in self.__submodels:
            if notfirst: out=out+"," 
            out=out+i.getFullName()
            notfirst=True
       return out+")"

   def setName(self,name):
       """sets the name of the model"""
       self.__name=name

   def setDescription(self,description="none"):
       """sets new description"""
       self.__description=description
       if self.debug(): print "%s: description is set to %s."%(self,description)
   def getDescription(self):
       """returns the description of the model"""
       return self.__description
   #
   #    parameter/attribute handling:
   #
   def declareParameter(self,**parameters):
      """declares a new parameter and its inital value."""
      for prm in parameters.keys():
         if prm in self.__dict__.keys():
             raise ValueError,"object attribute %s of %s cannot be used as a model parameter."%(prm,self)
         self.__parameters[prm]=parameters[prm]
         if self.debug(): print "%s: parameter %s has been declared."%(self,prm)



   def showParameters(self):
      """returns a descrition of the parameters"""
      out=""
      notfirst=False
      for i in self.__parameters:
          if notfirst: out=out+","
          notfirst=True
          out="%s%s=%s"%(out,i,self.__parameters[i])
      return out


   def deleteParameter(self,name):
      """removes parameter name from the model"""
      raise IllegalParameterError("Cannot delete parameter %s."%name)

   def getParameter(self,name):
      """returns the value of parameter name. If the parameter is not declared in self, the submodels are searched.
         if the parameter is a Link, the current value of the obejective is returned."""
      if self.__parameters.has_key(name):
          if isinstance(self.__parameters[name],Link):
             out=self.__parameters[name].getValue(name)
          else:
             out=self.__parameters[name]
      else:
          out=None
          for i in self.__submodels:
             try:
                out=self.__parameters[i].getParameter(name)
             except IllegalParameterError:
                pass
          if out==None: raise IllegalParameterError("Cannot find parameter %s."%name)
      return out

   def setParameter(self,**parameters):
      """sets parameter name to value. If the initial value for the parameter is a Model, the new value has to be a Model."""
      for name in parameters.keys():
         if self.__parameters.has_key(name):
            if not isinstance(parameters[name],Model) and isinstance(self.__parameters[name],Model):
                raise ValueError,"%s: parameter %s can assigned to a Model object only."%(self,name)
            if isinstance(parameters[name],Model) and not isinstance(self.__parameters[name],Model):
                raise ValueError,"%s: parameter %s is not declared as a Model."%(self,name)
            self.__parameters[name]=parameters[name]
            if isinstance(self.__parameters[name],Link):
                 if not self.__parameters[name].hasAttributeName(): self.__parameters[name].setAttributeName(name) 
            if self.debug(): print "%s: parameter %s has now value %s"%(self,name,self.__parameters[name])
         else:
            set=False
            for i in self.__submodels:
                try: 
                   self.__parameters[i].setParameter(**{name : parameters[name]})
                   set=True
                except IllegalParameterError:
                    pass
            if not set: raise IllegalParameterError("%s: Attempt to set undeclared parameter %s."%(self,name))

   def hasParameter(self,name):
      """returns True if self or one of the submodels has parameter name"""
      if self.__parameters.has_key(name):
         out=True
      else:
         out=False
         for i in self.__submodels: out= out or self.__parameters[i].hasParameter(name)
      return out

   def checkParameter(self,name):
      """checks if self has the parameter name. Otherewise ParameterError is thrown."""
      if not self.hasParameter(name):
           raise ParameterError("%s has no parameter %s."%(str(self),name))
  
   def __getattr__(self,name):
      """returns the value for attribute name. If name is in the Link list, the corresponding attribute is returned.""" 
      if self.__dict__.has_key(name):
         return self.__dict__[name]
      elif self.__dict__.has_key("_Model__parameters") and self.__dict__.has_key("_Model__submodels"): 
         return self.getParameter(name)
      else:
         raise AttributeError,"No attribute %s."%name

   def __setattr__(self,name,value):
      """returns the value for attribute name."""
      if self.__dict__.has_key("_Model__parameters") and self.__dict__.has_key("_Model__submodels"):
         if self.hasParameter(name): 
            self.setParameter(**{ name : value })
         else:
            self.__dict__[name]=value
      else: 
         self.__dict__[name]=value

   def __delattr__(self,name):
      """removes the attribute name."""
      if self.__dict__.has_key(name):
         del self.__dict__[name]
      elif self.__dict__.has_key("_Model__parameters"): 
         self.deleteParameter(name)
      else:
         raise AttributeError,"No attribute %s."%name

   # 
   #    submodel handeling:
   #
   def doInitializationOfSubmodels(self):
      """initializes the time stepping for all submodels."""
      for i in self.__submodels: self.getParameter(i).doInitialization()

   def getSafeTimeStepSizeFromSubmodels(self):
      """returns a time step size which can savely be used by all submodels. To avoid a big increase in the step size, 
         the new step size is restricted to the double of the precious step size."""
      out=None
      for i in self.__submodels: 
          dt=self.getParameter(i).getSafeTimeStepSize()
          if not dt==None: 
              if out==None:
                 out=dt
              else:
                 out=min(out,dt)
      return out

   def doStepOfSubmodels(self,t):
      """executes the time step for each submodel"""
      for i in self.__submodels: self.getParameter(i).doStep(t)

   def finalizeAllSubmodels(self):
      """returns True if all submodels can be finalized"""
      out=True
      for i in self.__submodels: out = out and self.getParameter(i).finalize()
      return out
      
   def doFinalizationOfSubmodels(self):
      """finalalizes the time stepping for each of the submodels."""
      for i in self.__submodels: self.getParameter(i).doFinalization()

   def doIterationInitializationOfSubmodels(self,t):
      """initializes the iteration for each of the submodels."""
      for i in self.__submodels: self.getParameter(i).doIterationInitialization(t)

   def doIterationStepOfSubmodels(self):
      """executes the iteration step at time step for each submodel"""
      for i in self.__submodels: self.getParameter(i).doIterationStep()

   def terminateAllSubmodels(self):
      """returns True if all iterations for all submodels are terminated."""
      out=True
      for i in self.__submodels: out = out and self.getParameter(i).terminate()
      return out
      
   def doIterationFinalizationOfSubmodels(self):
      """finalalizes the iteration process for each of the submodels."""
      for i in self.__submodels: self.getParameter(i).doIterationFinalization()

   def checkPointSubmodels(self):
      """performs check pointing for each submodel"""
      for i in self.__submodels: self.getParameter(i).checkPoint()

   #
   #   these methods control the time stepping
   #  
   def doInitialization(self):
      """initializes the time stepping"""
      self.doInitializationOfSubmodels()

   def getSafeTimeStepSize(self):
      """returns a time step size which can savely be used"""
      return self.getSafeTimeStepSizeFromSubmodels()

   def doStep(self,t):
      """executes the time step by first iterating over time step t and then step forward"""
      # run iteration on simulation until terminated:
      self.doIterationInitialization(t)
      while not self.terminate(): self.doIterationStep()
      self.doIterationFinalization()
      self.doStepOfSubmodels(t)

   def finalize(self):
      """returns True if all submodels are to be finalized"""
      return self.finalizeAllSubmodels()
      
   def doFinalization(self):
      """finalizes the time stepping."""
      self.doFinalizationOfSubmodels()
   #
   #   methods deal with iterations:
   #
   def doIterationInitialization(self,t):
      """initializes the iteration on a time step"""
      self.__iter=0
      if self.debug(): print "%s: iteration starts"%self
      self.doIterationInitializationOfSubmodels(t)

   def doIterationStep(self):
      """executes the iteration step"""
      self.__iter+=1
      if self.debug(): print "%s: iteration step %d"%(self,self.__iter)
      try:
         self.doIterationStepOfSubmodels()
      except IterationDivergenceError,e:
         raise IterationDivergenceError("divergence at time step %s in iteration step %s by reason: \n%s."%(self.__n,self.__iter,e.value))

   def terminate(self):
      """returns True if time steping is terminated"""
      return self.terminateAllSubmodels() 
      
   def doIterationFinalization(self):
      """finalalizes the iteration process."""
      self.doIterationFinalizationOfSubmodels()
      if self.debug(): print "%s: iteration finalized after %s step"%(self,self.__iter)
   #
   #   sum other method:
   #
   def checkPoint(self):
      """performs check pointing for each submodel"""
      if not self.__check_pointing==None:
         if self.__n%self.__check_pointing==0: self.checkPointsSubmodels()

   def run(self):
      """After check_pointing time steps the model will start to create checkpoint files for each of the submodels"""
      self.__tn=0.
      self.__n=0
      self.__dt=None
      self.doInitialization()
      while not self.finalize():
         self.__n+=1
         self.__dt=self.getSafeTimeStepSize()
         if self.__dt==None: self.__dt=self.UNDEF_DT
         if self.debug(): print "%s: %d. time step %e (step size %e.)"%(self,self.__n,self.__tn+self.__dt,self.__dt)
         endoftimestep=False
         while not endoftimestep:
              endoftimestep=True 
              try:
                 self.doStep(self.__tn+self.__dt)
              except FailedTimeStepError:
                 self.__dt=self.getSafeTimeStepSize()
                 if self.__dt==None: self.__dt=self.UNDEF_DT
                 endoftimestep=False
                 if self.debug(): print "%s: time step is repeated with new step size %e."%(self,self.__dt) 
              except IterationDivergenceError:
                 self.__dt*=0.5
                 endoftimestep=False
                 if self.debug(): print "%s: iteration failes. time step is repeated with new step size %e."%(self,self.__dt) 
         self.checkPoint()
         self.__tn+=self.__dt
      self.doFinalization()

class IterationDivergenceError(Exception):
    """excpetion which should be thrown if an iteration at a time step fails"""
    pass

class FailedTimeStepError(Exception):
    """excpetion which should be thrown if the time step fails because of a step size that have been choosen to be to large"""
    pass

class IllegalParameterError(Exception):
    """excpetion which is thrown if model has not the desired parameter"""
    pass


if __name__=="__main__":
   class Messenger(Model):
      def __init__(self):
         Model.__init__(self,parameters={"message" : "none" },name="messenger")

      def doInitialization(self):
         print "I start talking now!"

      def doStep(self,t):
         print "Message (time %e) : %s "%(t,self.message)

      def doFinalization(self):
         print "I have no more to say!"
   
   # explicit scheme
   class  Ode1(Model):
      def __init__(self,**args):
           Model.__init__(self,parameters={"tend" : 1., "dt" : 0.0001 ,"a" : 0.1 ,"u" : 1. , "message" : "none" },name="Ode1",debug=True)

      def doInitialization(self):
           self._tn=0

      def doStep(self,t):
           self.u=self.u+(t-self._tn)*self.a*self.u**2
           self._tn=t

      def doFinalization(self):
           self.message="current error = %e"%abs(self.u-1./(1./3.-self.a*self._tn))
           print self.message

      def getSafeTimeStepSize(self):
           return self.dt

      def finalize(self):
           return self._tn>=self.tend
   # explicit scheme
   class  Ode2(Model):

       def __init__(self,**args):
           Model.__init__(self,parameters={"tend" : 1., "dt" : 0.0001 ,"a" : 0.1 ,"u" : 10000. },name="Ode2",debug=True)
           self.declareParameter(tol=1.e-8,message="none")
           

       def doInitialization(self):
           self._tn=0
           self._iter=0

       def doIterationInitialization(self,t):
            self._iter=0
            self._u_last=self.u            
            self._dt=t-self._tn
            self._tn=t

       def doIterationStep(self):
          self._iter+=1
          self._u_old=self.u
          self.u=(self._dt*self.a*self.u**2-self._u_last)/(2*self._dt*self.a*self.u-1.)

       def terminate(self):
          if self._iter<1:
              return False
          else:
             return abs(self._u_old-self.u)<self.tol*abs(self.u)

       def doIterationFinalization(self):
           self.message="current error = %e"%abs(self.u-1./(1-self.a*self._tn))
           print self.message

       def getSafeTimeStepSize(self):
           return self.dt

       def finalize(self):
            return self._tn>=self.tend

   # a simple model with paramemter tend, dt, p1, p2, and p3
   class Test1(Model):

       def __init__(self,**args):
           Model.__init__(self,{"tend" : 1., "dt" : 0.1 ,"p1" : 0 ,"p2" : 0 ,"p3" : 0 },"test","bla",None,True)
           self.setParameters(args)

       def doInitialization(self):
           self.__tn=0
           self.__n=0

       def doStep(self,t):
           self.p3=self.p1+t*self.p2
           self.__tn=t
           print "test1 set the value out1 to ",self.p3

       def doFinalization(self):
           pass

       def getSafeTimeStepSize(self):
           return self.dt

       def finalize(self):
            return self._tn>self.tend


   class Test2(Model):

       def __init__(self):
           Model.__init__(self,{"q1": None},"test2","",None,True)


       def doInitialization(self):
           print "the whole thing starts"

       def doStep(self,t):
           print "test2 things that out1 is now ",self.out1

       def doFinalization(self):
           print "all done"

       def finalize(self):
            return True

   class Test12(Model):
     """model build from two models in a transperent way"""
     def __init__(self):
         Model.__init__(self,{"sm1": None, a : 0, "sm2": None},"test2","",None,True)
         self.setExecutionOrder(["sm2","sm1"])

   # test messenger
   m=Messenger()
   m.run()
   # ode1      
   o=Ode1()
   o.dt=0.001
   o.u=3.
   o.run()
   # ode1
   o=Ode2()
   o.dt=0.01
   o.a=0.1
   o.u=1.
   o.run()
   # and they are linked together:
   o=Ode2()
   m=Messenger()
   om=Model(submodels=[o,m],debug=True)
   om.dt=0.01
   om.u=1.
   m.message=Link(o)
   om.run()
   print om.showParameters()
   1/0

   t=Test1()
   t.tend=1.
   t.dt=0.25
   t.in1=1.
   t.in2=3.
   t.run()
   # and a coupled problem:
   t2=Test2()
   t2.out1=Link(t)
   Model([t,t2],debug=True).run()
