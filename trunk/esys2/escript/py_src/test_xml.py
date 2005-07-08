#!/usr/bin/python


import unittest
from modelframe import *
import math
from cStringIO import StringIO
from xml.dom import minidom

class XMLDocumentTestCase(unittest.TestCase):

    def setUp(self):

        o1=ODETEST(debug=False)
        o1.u=10
        o2=ODETEST(debug=False)
        o2.u=-10.
        o1.f=Link(o2,"u")
        o2.f=Link(o1,"u")
        m=Messenger()
        o1.dt=0.01
        m.message=Link(o1)
        s=ExplicitSimulation([Simulation([o1,o2],debug=False),m],debug=False)
        s.run() 
        output = StringIO()
        s.writeXML(output)
        output.reset()
        outputList = output.readlines()
        self.xmlList = outputList

    def testFirstLine(self):
        firstLine = self.xmlList[0]
        self.assertEqual('<?xml version="1.0" ?>\n', firstLine)
        
    def testEsysHeader(self):
        header = self.xmlList[1]
        self.assertEqual('<ESys>\n', header)

    def testEsysFooter(self):
        footer = self.xmlList[-1]
        self.assertEqual('</ESys>\n', footer)

    def testSimulationHeader(self):
        pass

    def testSimulationFooter(self):
        pass

class SimulationTestCase(unittest.TestCase):
    def setUp(self):
        o1=ODETEST(debug=False)
        o1.u=10
        o2=ODETEST(debug=False)
        o2.u=-10.
        o1.f=Link(o2,"u")
        o2.f=Link(o1,"u")
        m=Messenger()
        o1.dt=0.01
        m.message=Link(o1)
        self.s=s=ExplicitSimulation([Simulation([o1,o2],debug=False),m],debug=False)
        s.run() 
        output = StringIO()
        s.writeXML(output)
        output.reset()
        self.xml = output.read()

    def testSimulation(self):
        assert "<Link" in self.xml, "I should see a link"

        o = parse(self.xml)
        output = StringIO()
        o.writeXML(output)
        output.reset()
        xml = output.read()
        assert '<Link' in xml, "Missing a link! It should be in this!"




class LinkTestCase(unittest.TestCase):
   

    def setUp(self):

        self.o1=ODETEST(debug=False)
        #self.o1.u=10
        self.o2=ODETEST(debug=False)
        self.o2.u=-10.
        self.o1.f=Link(self.o2,"u")
        self.o2.f=Link(self.o1,"u")
        self.o2.declareParameter(child=self.o1)

    def testLinkCreation(self):
       self.o1.f=Link(self.o2,"u")
       assert self.o1.f


    def testLinkValue(self):
       self.assertEqual(self.o1.f, -10)

    def testLinkTarget(self):
       pass

    def testLinkDefaultAttribute(self):
        Link(self.o2)

    def testLinkXML(self):
        s = StringIO()
        self.o2.writeXML(s)
        s.reset()
        xmlout = s.read()
        assert '<Link' in xmlout

    def testLinkTargetXML(self):
       pass
        
class ParamaterSetTestCase(unittest.TestCase):


    def setUp(self):
        self.p = ParameterSet()
        self.p.declareParameter(gamma1=1.,gamma2=2.,gamma3=3.)

    def testParameterSetCreation(self):
        self.assertEqual(self.p.gamma1, 1.)

    def testParameterSetXMLCreation(self):
        s = StringIO()
        self.p.writeXML(s)
        s.reset()
        xmlout = s.read()
        assert ("gamma1" in xmlout)
        assert ("gamma2" in xmlout)
        assert ("gamma3" in xmlout)
        parsable = minidom.parseString(xmlout)
        assert parsable.getElementsByTagName("ParameterSet")

    def testParameterSetFromXML(self):
        doc = self._dom()
        pset = ParameterSet.fromDom(doc.getElementsByTagName("ParameterSet")[0])
        assert (isinstance(pset, ParameterSet))
        self.assertEqual(self.p.gamma1,pset.gamma1)


    def testParameterSetWithChildrenFromXML(self):
        p2 = ParameterSet()
        p2.declareParameter(s="abc", f=3.)
        self.p.declareParameter(child=p2)
        doc = self._dom()
        pset = ParameterSet.fromDom(doc.getElementsByTagName("ParameterSet")[0])
        self.assertEqual(self.p.child.f, pset.child.f) 

    def testParameterSetChild(self):
        p2 = ParameterSet()
        p2.declareParameter(s="abc", f=3.)
        self.p.declareParameter(child=p2)
        self.assertEqual(self.p.child.s, "abc")
        self.assertEqual(self.p.child.f, 3.)

    def _dom(self):
        s = StringIO()
        self.p.writeXML(s)
        s.reset()
        xmlout = s.read()
        doc = minidom.parseString(xmlout)
        return doc

class ModeltoDomTestCase(unittest.TestCase):
    
    def _dom(self):
        s = StringIO()
        self.o1.writeXML(s)
        s.reset()
        xmlout = s.read()
        doc = minidom.parseString(xmlout)
        return doc

    def setUp(self):
        self.o1=ODETEST(debug=False)

    def testModelExists(self):
        assert self._dom().getElementsByTagName("Model")
        #print self._dom().toxml()

    def testModelhasID(self):
        assert int(self._dom().getElementsByTagName("Model")[0].getAttribute("id"))>99
        
class Messenger(Model):
    def __init__(self, *args, **kwargs):
        Model.__init__(self, *args, **kwargs)
        self.declareParameter(message="none")

    def doInitialization(self,t):
        self.__t=t
        #print "I start talking now!"

    def doStep(self,dt):
        self.__t+=dt
        #print "Message (time %e) : %s "%(self.__t,self.message)

    def doFinalization(self):
        #print "I have no more to say!"
        pass



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

    def __init__(self, *args, **kwargs):
        Model.__init__(self, *args, **kwargs)
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

    
if __name__ == "__main__":
            unittest.main()

