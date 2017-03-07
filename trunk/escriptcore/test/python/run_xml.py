# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2017 by The University of Queensland
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

__copyright__="""Copyright (c) 2003-2017 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript.modelframe import Model,Link,Simulation,ParameterSet,ESySXMLParser,DataSource
import math
from io import StringIO
from xml.dom import minidom
import numpy
import sys

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
        s=Simulation([o1,o2,m],debug=False)
        s.run() 
        output = StringIO()
        s.writeXML(output)
        output.seek(0)
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
        self.s=Simulation([o1,o2,m],debug=False)
        self.s.run() 
        output = StringIO()
        self.s.writeXML(output)
        output.seek(0)
        self.xml = output.read()

    def testSimulation(self):
        assert "<Simulation" in self.xml, "I should see a Simulation"
        
    def testParseAndInstanceOfSimulation(self):
        
        xml = ESySXMLParser(self.xml)
        newSim = xml.parse()[0]
        assert isinstance(newSim, Simulation)
        newout = StringIO()
        newSim.writeXML(newout)
        newout.seek(0)
        xml = newout.read()
        assert '<Simulation' in xml, "Missing a Simulation! It should be in this!"




class Test_LinkTestCase(unittest.TestCase):
   

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
        s.seek(0)
        xmlout = s.read()
        assert '<Link' in xmlout

    def testLinkTargetXML(self):
       pass
        
class Test_ParameterSetTestCase(unittest.TestCase):


    def setUp(self):
        self.p = ParameterSet()
        self.p.declareParameter(gamma1=1.,gamma2=2.,gamma3=3.)

    def testParameterSetCreation(self):
        self.assertEqual(self.p.gamma1, 1.)

    def testParameterSetXMLCreation(self):
        s = StringIO()
        self.p.writeXML(s)
        s.seek(0)
        xmlout = s.read()
        assert ("gamma1" in xmlout)
        assert ("gamma2" in xmlout)
        assert ("gamma3" in xmlout)
        esysxml=ESySXMLParser(xmlout)
        parsable = esysxml.parse()[0]
        assert (isinstance (parsable, ParameterSet))
        assert (self._dom(self.p).getElementsByTagName("ParameterSet"))

    def testParameterSetFromXML(self):
        s = StringIO()
        self.p.writeXML(s)
        s.seek(0)
        xmlout = s.read()
        esysxml=ESySXMLParser(xmlout)
        doc = self._class(self.p)
        pset = ParameterSet.fromDom(esysxml, self._dom(self.p).getElementsByTagName("ParameterSet")[0])
        assert (isinstance(pset, ParameterSet))
        assert (isinstance(doc, ParameterSet))
        self.assertEqual(self.p.gamma1,doc.gamma1)


    def testParameterSetWithChildrenFromXML(self):
        p2 = ParameterSet()
        p2.declareParameter(s="abc", f=3.)
        self.p.declareParameter(child=p2)
        doc = self._class(self.p)
        #pset = ParameterSet.fromDom(doc.getElementsByTagName("ParameterSet")[0])
        self.assertEqual(self.p.child.f, doc.child.f) 

    def testParameterSetChild(self):
        p2 = ParameterSet()
        p2.declareParameter(s="abc", f=3.)
        self.p.declareParameter(child=p2)
        self.assertEqual(self.p.child.s, "abc")
        self.assertEqual(self.p.child.f, 3.)

    def _dom(self, input):
        s = StringIO()
        input.writeXML(s)
        s.seek(0)
        xmlout = s.read()
        doc = minidom.parseString(xmlout)
        return doc

    def _class(self, input):
        s = StringIO()
        input.writeXML(s)
        s.seek(0)
        xmlout = s.read()
        esysxml=ESySXMLParser(xmlout)
        doc =esysxml.parse()[0]
        return doc

    def testFromDomInt(self):
        p3 = ParameterSet()
        p3.declareParameter(inttest=1)
        doc = self._class(p3)
        assert type(doc.inttest)==int

    def testFromDomNumarrayVector(self):
        p3 = ParameterSet()
        mynumpy = numpy.array([3., 4., 5.], dtype=numpy.float64)
        p3.declareParameter(numtest=mynumpy)
        doc = self._class(p3)
        assert doc.numtest.dtype == numpy.float64
        assert type(doc.numtest) == numpy.ndarray

    def testFromDomNumarrayMulti(self):
        p3 = ParameterSet()
        mynumpy = numpy.array([[1., 2., 3.], [3., 4., 5.]], dtype=numpy.float64)
        p3.declareParameter(numtest=mynumpy)
        doc = self._class(p3)
        assert doc.numtest.dtype == numpy.float64
        assert type(doc.numtest) == numpy.ndarray

    def testBoolLists(self):
        p4 = ParameterSet()
        mylist = [True, False, False, True]
        p4.declareParameter(listest=mylist)
        doc = self._class(p4)
        assert type(doc.listest) == numpy.ndarray
        assert doc.listest.dtype == numpy.bool_
        assert len(doc.listest) == len(mylist)
        assert min([ mylist[i] == doc.listest[i] for i in range(len( doc.listest)) ])

    def testIntLists(self):
        p4 = ParameterSet()
        mylist = [1,2,4]
        p4.declareParameter(listest=mylist)
        doc = self._class(p4)
        assert type(doc.listest) == numpy.ndarray
        assert doc.listest.dtype == numpy.int_
        assert len(doc.listest) == len(mylist)
        assert min([ mylist[i] == doc.listest[i] for i in range(len( doc.listest)) ])

    def testFloatLists(self):
        p4 = ParameterSet()
        mylist = [1.,2.,4.]
        p4.declareParameter(listest=mylist)
        doc = self._class(p4)
        assert type(doc.listest) == numpy.ndarray
        assert doc.listest.dtype == numpy.float_
        assert len(doc.listest) == len(mylist)
        assert min([ mylist[i] == doc.listest[i] for i in range(len( doc.listest)) ])

    def testStringLists(self):
        p4 = ParameterSet()
        mylist = ["a", "b", "c"]
        p4.declareParameter(listest=mylist)
        doc = self._class(p4)
        assert type(doc.listest) == list
        assert len(doc.listest) == len(mylist)
        assert min([ mylist[i] == doc.listest[i] for i in range(len( doc.listest)) ])
        
    def testDatasource(self):
        p5 = ParameterSet()
        myURI = DataSource("somelocalfile.txt", "text")
        p5.declareParameter(uritest=myURI)
        doc = self._class(p5)
        self.assertEqual(myURI.uri, doc.uritest.uri)
        self.assertEqual(myURI.fileformat, doc.uritest.fileformat)
        assert type(doc.uritest) == DataSource
        

       
class Test_ModeltoDomTestCase(unittest.TestCase):
    
    def _class(self):
        # returns a modelframe class, generated from the xml
        s = StringIO()
        self.o1.writeXML(s)
        s.seek(0)
        self.xmlout = s.read()
        esysxml=ESySXMLParser(self.xmlout)
        doc =esysxml.parse()[0]
        return doc

    def _dom(self):
        # returns a minidom dom element, generated from the xml
        s = StringIO()
        self.o1.writeXML(s)
        s.seek(0)
        self.xmlout = s.read()
        doc = minidom.parseString(self.xmlout)
        return doc
    
    def setUp(self):
        self.o1=ODETEST(debug=False)
        self.o1.message='blah'

    def testModelExists(self):
        modeldoc = self._class()
        assert isinstance(modeldoc, Model)
        assert self._dom().getElementsByTagName("Model")
    
    def testModelhasID(self):
        assert int(self._dom().getElementsByTagName("Model")[0].getAttribute("id"))>99

class ModeltoDomTestCase(unittest.TestCase):
    def _xml(self, modulename, modelname):
        # returns a modelframe class, generated from the xml
        return '''<?xml version="1.0" ?>
<ESys> <Simulation type="Simulation"> <Component rank="0"> 

    <Model id="127" module="%s" type="%s"> 

<Parameter type="float"> <Name> a </Name> <Value> 0.9 </Value> </Parameter>
<Parameter type="Link"> <Name> f </Name> <Value> <Link> <Target> 128 </Target>
<Attribute> u </Attribute> </Link> </Value> </Parameter> <Parameter
type="float"> <Name> tend </Name> <Value>
1.0 </Value> </Parameter> <Parameter type="int"> <Name> u </Name> <Value> 10
  </Value> </Parameter> <Parameter type="float"> <Name> tol </Name> <Value>
  1e-08 </Value> </Parameter> <Parameter type="float"> <Name> dt </Name>
  <Value>
0.01 </Value> </Parameter> <Parameter type="str"> <Name> message </Name>
  <Value> current error = 9.516258e-01 </Value> </Parameter> </Model>
  </Component> <Component rank="1"> <Model id="128" type="ODETEST"> <Parameter
  type="float"> <Name> a </Name> <Value>
0.9 </Value> </Parameter> <Parameter type="Link"> <Name> f </Name> <Value>
  <Link> <Target> 127 </Target> <Attribute> u </Attribute> </Link> </Value>
  </Parameter> <Parameter type="float"> <Name> tend </Name> <Value>
1.0 </Value> </Parameter> <Parameter type="float"> <Name> u </Name> <Value>
  -10.0 </Value> </Parameter> <Parameter type="float"> <Name> tol </Name>
  <Value> 1e-08 </Value> </Parameter> <Parameter type="float"> <Name> dt
  </Name> <Value>
0.1 </Value> </Parameter> <Parameter type="str"> <Name> message </Name> <Value>
  current error = 1.904837e+01 </Value> </Parameter> </Model> </Component>
  <Component rank="2"> <Model id="129" type="Messenger"> <Parameter
  type="Link"> <Name> message </Name> <Value> <Link> <Target> 127 </Target>
  <Attribute> message </Attribute> </Link> </Value> </Parameter> </Model>
  </Component> </Simulation> <Model id="150" type="ODETEST"> <Parameter
  type="float"> <Name> a </Name> <Value>
0.9 </Value> </Parameter> <Parameter type="Link"> <Name> f </Name> <Value>
  <Link> <Target> 127 </Target> <Attribute> u </Attribute> </Link> </Value>
  </Parameter> <Parameter type="float"> <Name> tend </Name> <Value>
1.0 </Value> </Parameter> <Parameter type="float"> <Name> u </Name> <Value>
  -10.0 </Value> </Parameter> <Parameter type="float"> <Name> tol </Name>
  <Value> 1e-08 </Value> </Parameter> <Parameter type="float"> <Name> dt
  </Name> <Value>
0.1 </Value> </Parameter> <Parameter type="str"> <Name> message </Name> <Value>
  current error = 1.904837e+01 </Value> </Parameter> </Model> <Model id="130"
  type="ODETEST"> <Parameter type="float"> <Name> a </Name> <Value>
0.9 </Value> </Parameter> <Parameter type="Link"> <Name> f </Name> <Value>
  <Link> <Target> 128 </Target> <Attribute> u </Attribute> </Link> </Value>
  </Parameter> <Parameter type="float"> <Name> tend </Name> <Value>
1.0 </Value> </Parameter> <Parameter type="int"> <Name> u </Name> <Value> 10
  </Value> </Parameter> <Parameter type="float"> <Name> tol </Name> <Value>
  1e-08 </Value> </Parameter> <Parameter type="float"> <Name> dt </Name>
  <Value>
0.01 </Value> </Parameter> <Parameter type="str"> <Name> message </Name>
  <Value> current error = 9.516258e-01 </Value> </Parameter> </Model> <Model
  id="170" type="ODETEST"> <Parameter type="float"> <Name> a </Name> <Value>
0.9 </Value> </Parameter> <Parameter type="Link"> <Name> f </Name> <Value>
  <Link> <Target> 128 </Target> <Attribute> u </Attribute> </Link> </Value>
  </Parameter> <Parameter type="float"> <Name> tend </Name> <Value>
1.0 </Value> </Parameter> <Parameter type="int"> <Name> u </Name> <Value> 10
  </Value> </Parameter> <Parameter type="float"> <Name> tol </Name> <Value>
  1e-08 </Value> </Parameter> <Parameter type="float"> <Name> dt </Name>
  <Value>
0.01 </Value> </Parameter> <Parameter type="str"> <Name> message </Name>
  <Value> current error = 9.516258e-01 </Value> </Parameter> </Model> </ESys>
''' % (modulename, modelname)

    def testModuleAttribute(self):
        esysxml=ESySXMLParser(self._xml('run_xml', 'ODETEST'))
        modeldoc=esysxml.parse()[0]

    def testModuleAttributeFails(self):
        try:
            esysxml=ESySXMLParser(self._xml('a', 'b'))
            modeldoc=esysxml.parse()[0]
        except ImportError:
            return # correct

        assert False, "This test should have resulted in an ImportError"
        
class Messenger(Model):
    def __init__(self, *args, **kwargs):
        Model.__init__(self, *args, **kwargs)
        self.declareParameter(message="none")

    def doInitialization(self):
        self.__t=0
        #print "I start talking now!"

    def doStepPostprocessing(self,dt):
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

    def doInitialization(self):
        self.__tn=0
        self.__iter=0
       
    def doStepPreprocessing(self,dt):
        self.__iter=0
        self.__u_last=self.u            

    def doStep(self,dt):
        self.__iter+=1
        self.__u_old=self.u
        self.u=self.__u_last+dt*(self.a*self.__u_old+self.f)
     
    def terminate(self):
        if self.__iter<1:
            return False
        else:
            return abs(self.__u_old-self.u)<self.tol*abs(self.u)

    def doStepPostprocessing(self,dt):
        self.__tn+=dt
        self.message="current error = %e"%abs(self.u-10.*math.exp((self.a-1.)*self.__tn))

    def getSafeTimeStepSize(self,dt):
        return min(self.dt,1./(abs(self.a)+1.))

    def finalize(self):
        return self.__tn>=self.tend


    
if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
