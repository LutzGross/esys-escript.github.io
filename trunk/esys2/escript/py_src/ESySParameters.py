# $Id$
""" A simple tool to handel parameters for a simulation in easy way
    The idea is that all parameters are stored in a single object in a hierachical form
    and can accessed using python's attribute notation.  For instance:

        parm.parm2.alpha=814.
        parm.parm1.gamma=0.
        parm.parm1.dim=2
        parm.parm1.tol_v=0.001
        parm.parm1.output_file="/tmp/u.%3.3d.dx"
        parm.parm1.runFlag=True
        parm.parm1.T=1.
        parm.parm1.x1=[-1.,2]
        parm.parm1.x2=(10.,11)
        parm.parm1.x3=(-10.,)
        parm.parm1.parm11.gamma1=1.
        parm.parm1.parm11.gamma2=2.
        parm.parm1.parm11.gamma3=3.

This structure can be stored/defined through an XML file parm.xml:

<?xml version="1.0"?>
<ESyS>
   <Component type="Geodynamics">
   <Name>parm1</Name>
   <Description>
     a few examples of parameters
   </Description>
   <Parameter><Item>gamma</Item><Value>0.</Value></Parameter>
   <Parameter type="int"><Item>dim</Item><Value>2</Value></Parameter>
   <Parameter type="real"><Item>tol_v</Item><Value>0.001</Value></Parameter>
   <Parameter type="string"><Item>output_file</Item><Value>/tmp/u.%3.3d.dx</Value></Parameter>
   <Parameter type="bool"><Item>runFlag</Item><Value>true</Value></Parameter>
   <Parameter type="real" sequence="single"><Item>T</Item><Value>1.</Value><Value>2</Value></Parameter>
   <Parameter type="real" sequence="list"><Item>x1</Item><Value>-1.</Value><Value>2</Value></Parameter>
   <Parameter type="real" sequence="tuple"><Item>x2</Item><Value>10</Value><Value>11</Value></Parameter>
   <Parameter sequence="tuple"><Item>x3</Item><Value>-10.</Value></Parameter>
   <Component>
           <Name>parm11</Name>
           <Description>
          a sub compoment
           </Description>
            <Parameter><Item>gamma1</Item><Value>1.</Value></Parameter>
            <Parameter><Item>gamma2</Item><Value>2.</Value></Parameter>
            <Parameter><Item>gamma3</Item><Value>3.</Value></Parameter>
           </Component>   
   </Component>
    <Component type="Geodynamics">
   <Name>parm2</Name>
   <Description>
     another component
   </Description>
   <Parameter><Item>alpha</Item><Value>0814</Value></Parameter>
   </Component>
</ESyS>

"""

import types
from xml.dom import minidom
from string import atoi,atof

class ESySParameters:
   """is an object to store simulation parameters in the form of a tree and
      access their values in an easy form
      
      Leaves of  an ESySParameters objects can be 

           a single real number or a list or tuple of real numbers
           a single integer number or a list or tuple of integer  numbers
           a single strings or a list or tuple of strings
           a single boolean value or a list or tuple of boolean values
           a ESySParameters object
           any other object (not considered by writeESySXML and writeXML)

          Example how to create an ESySParameters object:

                parm=ESySParameters()
                parm.parm1=ESySParameters()
                parm.parm1.gamma=0.
                parm.parm1.dim=2
                parm.parm1.tol_v=0.001
                parm.parm1.output_file="/tmp/u.%3.3d.dx"
                parm.parm1.runFlag=True
                parm.parm1.T=1.
                parm.parm1.x1=[-1.,2]
                parm.parm1.x2=(10.,11)
                parm.parm1.x3=(-10.,)
                parm.parm1.parm11=ESySParameters()
                parm.parm1.parm11.gamma1=1.
                parm.parm1.parm11.gamma2=2.
                parm.parm1.parm11.gamma3=3.
                parm.parm2=ESySParameters()
                parm.parm2.alpha=814.

                print parm

            Output is

            (parm1=(dim=2,output_file=/tmp/u.%3.3d.dx,parm11=(gamma3=3.0,gamma2=2.0,gamma1=1.0),
                    tol_v=0.001,T=1.0,x2=(10,),x3=(-10,),runFlag=True,x1=[-1.0, 2.0],gamma=0.0),parm2=(alpha=814.0))

            Notice that parm.parm1.x1 is now a list of two floats although it is defined by a list of a float and an integer.
            ESySParameter is trying to use the same type for all items in a list or a tuple.
            
   """

   def __init__(self,description="none",type=None):
      self.setDescription(description)
      self.setType(type)
      
   def getTypeName(self):
       if self.__type==None:
          return None
       else:
          return self.__type.__str__()

   def getDescription(self):
       return self.__description

   def setType(self,type=None):
       self.__type=type
         
   def setDescription(self,description="none"):
       self.__description=description

   def __str__(self):
       """returns a string representation"""
       out=""
       for name,value in self.__dict__.iteritems():
           if name[0]!="_":
               if out=="":
                   out=name+"="+str(value)
               else:
                   out=out+","+name+"="+str(value)
       return "("+out+")"
        
   def __setattr__(self,name,value):
     """defines attribute name and assigns value. if name does not start
        with an underscore value has to be a valid Parameter."""
     name=name.replace(" ","_")
     if name[0]!="_":
       if value==None:
          self.__dict__[name]=value
       elif isinstance(value,ESySParameters):
          self.__dict__[name]=value
       elif isinstance(value,types.BooleanType):
          self.__dict__[name]=value   
       elif isinstance(value,types.ListType):
          self.__dict__[name]=_mkSameType(value)
       elif isinstance(value,types.TupleType):
          self.__dict__[name]=tuple(_mkSameType(value))
       elif isinstance(value,types.BooleanType):
          self.__dict__[name]=value
       elif isinstance(value,types.IntType):
          self.__dict__[name]=value
       elif isinstance(value,types.FloatType):
          self.__dict__[name]=value
       elif isinstance(value,types.StringType) or isinstance(value,types.UnicodeType):
          self.__dict__[name]=str(value)
       else:
          self.__dict__[name]=value
     else:
       self.__dict__[name]=value
           
   def writeXML(self,iostream):
     """writes the object as an XML object into an IO stream"""
     for name,value in self.__dict__.iteritems():
        if name[0]!="_":
           if isinstance(value,ESySParameters):
              sequence=_PARAMETER_SEQUENCE_UNKNOWN
              type=value.getTypeName()
              iostream.write("<%s"%_COMPONENT)
              if type!=None: iostream.write("%s=\"%s\""%(_COMPONENT_TYPE_ATTRIBUTE,type))
              iostream.write(">\n<%s>%s</%s>\n<%s>%s</%s>\n"%(_NAME,name,_NAME,_DESCRIPTION,value.getDescription(),_DESCRIPTION))
              value.writeXML(iostream)
              iostream.write("</%s>"%_COMPONENT)
           else:
               if isinstance(value,types.ListType):
                  sequence=_PARAMETER_SEQUENCE_LIST
                  type=_getTypeNameOfList(value)
               elif isinstance(value,types.TupleType):
                  sequence=_PARAMETER_SEQUENCE_TUPLE
                  type=_getTypeNameOfList(value)
               else:
                  sequence=_PARAMETER_SEQUENCE_SINGLE
                  type=_getTypeName(value)
               iostream.write("<%s %s=\"%s\" %s=\"%s\"><%s>%s</%s>\n"% \
                           (_PARAMETER,_PARAMETER_TYPE_ATTRIBUTE,type, \
                                      _PARAMETER_SEQUENCE_ATTRIBUTE,sequence, \
                                      _PARAMETER_ITEM,name,_PARAMETER_ITEM))
               if type!=_PARAMETER_TYPE_UNKNOWN:
                   if sequence==_PARAMETER_SEQUENCE_LIST or sequence==_PARAMETER_SEQUENCE_TUPLE:
                       for i in value:
                           iostream.write("<%s>%s</%s>"%(_PARAMETER_VALUE,i.__str__(),_PARAMETER_VALUE))
                   elif sequence==_PARAMETER_SEQUENCE_SINGLE:
                       iostream.write("<%s>%s</%s>\n"%(_PARAMETER_VALUE,value.__str__(),_PARAMETER_VALUE))
               iostream.write("</%s>\n"%_PARAMETER)


   def writeProperties(self, iostream, nameSpace = None):
     """writes the object as a property list to an IO stream"""
     if nameSpace != None and nameSpace != "":
        nameSpace += ".";
     for name,value in self.__dict__.iteritems():
        if name[0]!="_":
           if isinstance(value,ESySParameters):
              value.writeProperties(iostream, name)
           else:
               if isinstance(value,types.ListType):
                  sequence=_PARAMETER_SEQUENCE_LIST
                  type=_getTypeNameOfList(value)
               elif isinstance(value,types.TupleType):
                  sequence=_PARAMETER_SEQUENCE_TUPLE
                  type=_getTypeNameOfList(value)
               else:
                  sequence=_PARAMETER_SEQUENCE_SINGLE
                  type=_getTypeName(value)
               iostream.write("%s = %s\n" % (nameSpace + name, value.__str__()))

   def writeESySXML(self,iostream):
        """writes an ESyS XML file"""
        iostream.write("<?xml version=\"1.0\"?><ESyS>")
        self.writeXML(iostream)
        iostream.write("</ESyS>")

def readESySXMLFile(filename):
       """reads an ESyS XML file and returns it as a ESySParameter object"""
       return _readParametersFromDOM(minidom.parse(filename).getElementsByTagName(_ESYS)[0])


_ESYS="ESyS"
_COMPONENT="Component"
_COMPONENT_TYPE_ATTRIBUTE="type"
_NAME="Name"
_DESCRIPTION="Description"
_PARAMETER="Parameter"
_PARAMETER_ITEM="Item"
_PARAMETER_VALUE="Value"
_PARAMETER_TYPE_ATTRIBUTE="type"
_PARAMETER_TYPE_REAL="real"
_PARAMETER_TYPE_INT="int"
_PARAMETER_TYPE_STRING="string"
_PARAMETER_TYPE_BOOL="bool"
_PARAMETER_TYPE_UNKNOWN="unknown"
_PARAMETER_SEQUENCE_ATTRIBUTE="sequence"
_PARAMETER_SEQUENCE_UNKNOWN="unknown"
_PARAMETER_SEQUENCE_SINGLE="single"
_PARAMETER_SEQUENCE_LIST="list"
_PARAMETER_SEQUENCE_TUPLE="tuple"


def _mkSameType(list):
    """returns list where all items in the list have the same type"""
    out=[]
    if len(list)>0:
        type=0
        for i in list:
            if isinstance(i,types.BooleanType):
                type=max(type,0)
            elif isinstance(i,types.IntType):
                type=max(type,1)
            elif isinstance(i,types.FloatType):
                type=max(type,2)
            elif isinstance(i,types.StringType):
                type=max(type,3)
            else:
                raise TypeError,"illegal item type"
            
        for i in list:
            if isinstance(i,types.BooleanType):
                if type==0:
                    out.append(i)
                elif type==1:
                    out.append(int(i))
                elif type==2:
                    out.append(float(i))
                else:
                    out.append(i.__str__()) 
            elif isinstance(i,types.IntType):
                if type==1:
                    out.append(i)
                elif type==2:
                    out.append(float(i))
                else:
                    out.append(i.__str__()) 
            elif isinstance(i,types.FloatType):
                if type==2:
                    out.append(i)
                else:
                    out.append(i.__str__()) 
            else:
                out.append(i)
    return out

def _getTypeNameOfList(values):
    """returns the type of the parameters in list values"""
    if len(values)==0:
        type=_PARAMETER_TYPE_UNKNOWN
    else:
        type=_getTypeName(values[0])
    return type

def _getTypeName(value):
    """returns the type of the parameter value"""
    if isinstance(value,types.FloatType):
           type=_PARAMETER_TYPE_REAL
    elif isinstance(value,types.BooleanType):
           type=_PARAMETER_TYPE_BOOL
    elif isinstance(value,types.IntType):
           type=_PARAMETER_TYPE_INT
    elif isinstance(value,types.StringType):
           type=_PARAMETER_TYPE_STRING
    else:
           type=_PARAMETER_TYPE_UNKNOWN
    return type

    
def _extractStrippedValue(dom):
    """exracts a string from a DOM node"""
    out=""
    for i in dom.childNodes:
        s=i.nodeValue.strip()
        if s!="\n": out+=s
    return str(out)
    
def _readParametersFromDOM(dom):
    out=ESySParameters()
    for node in dom.childNodes:
        if node.nodeType==node.ELEMENT_NODE:
            if node.nodeName==_COMPONENT:
                name=None
                description="none"
                type=None
                if node.hasAttribute(_COMPONENT_TYPE_ATTRIBUTE): type=node.getAttribute(_COMPONENT_TYPE_ATTRIBUTE)
                # find description and name:
                for c_dom in node.childNodes:
                    if c_dom.nodeType==c_dom.ELEMENT_NODE:
                        if c_dom.tagName==_NAME: name=_extractStrippedValue(c_dom)
                        if c_dom.tagName==_DESCRIPTION: description=_extractStrippedValue(c_dom)
                if name==None:
                    raise IOError,"name of component missing"
                p=_readParametersFromDOM(node)
                p.setDescription(description)
                p.setType(type)
                out.__setattr__(name,p)
            elif node.nodeName==_PARAMETER:
                if node.hasAttribute(_PARAMETER_TYPE_ATTRIBUTE):
                    type=node.getAttribute(_PARAMETER_TYPE_ATTRIBUTE)
                    if type==_PARAMETER_TYPE_UNKNOWN: type=_PARAMETER_TYPE_REAL
                else:
                    type=_PARAMETER_TYPE_REAL
                if node.hasAttribute(_PARAMETER_SEQUENCE_ATTRIBUTE):
                    sequence=node.getAttribute(_PARAMETER_SEQUENCE_ATTRIBUTE)
                    if sequence==_PARAMETER_SEQUENCE_UNKNOWN: sequence=_PARAMETER_SEQUENCE_SINGLE
                else:
                    sequence=_PARAMETER_SEQUENCE_SINGLE
                # get the name and values as list:
                name=None
                p=[]
                for c_dom in node.childNodes:
                    if c_dom.nodeType==c_dom.ELEMENT_NODE:
                        if c_dom.nodeName==_PARAMETER_ITEM: name=_extractStrippedValue(c_dom)
                        if c_dom.nodeName==_PARAMETER_VALUE:
                            value=_extractStrippedValue(c_dom)
                            if type==_PARAMETER_TYPE_REAL:
                                p.append(atof(value))
                            elif type==_PARAMETER_TYPE_INT:
                                p.append(atoi(value))
                            elif type==_PARAMETER_TYPE_BOOL:
                                if value=="true" or value=="True" or value=="TRUE":
                                    p.append(True)
                                elif value=="false" or value=="FALSE" or value=="False":
                                    p.append(False)
                                else:
                                    raise IOError,"cannot convert %s to bool"%value
                            elif type==_PARAMETER_TYPE_STRING: 
                               p.append(value)
                            else:
                                raise IOError,"unknown parameter type %s"%type
                if name==None: raise IOError,"Item tag missing"
                if sequence==_PARAMETER_SEQUENCE_SINGLE:
                   if len(p)==0:
                      p=None
                   else:
                      p=p[0]
                elif sequence==_PARAMETER_SEQUENCE_TUPLE:
                   if len(p)==0:
                      p=tuple()
                   else:
                      p=tuple(p)
                elif sequence==_PARAMETER_SEQUENCE_LIST:
                     pass
                else:
                     raise IOError,"unknown sequence attribute %s"%sequence
                out.__setattr__(name,p)  
    return out

# test section:
if  (__name__=="__main__"):
    def test(parm):
        if parm.parm1.gamma!=0. : raise IOError,"unexpected value for parm.parm1.gamma"
        if parm.parm1.dim!=2: raise IOError,"unexpected value for parm.parm1.dim"
        if parm.parm1.tol_v!=0.001: raise IOError,"unexpected value for parm.parm1.tol_v"
        if parm.parm1.output_file!="/tmp/u.%3.3d.dx": raise IOError,"unexpected value for parm.parm1.output_file"
        if parm.parm1.runFlag!=True: raise IOError,"unexpected value for parm.parm1.runFlag"
        if parm.parm1.T!=1.: raise IOError,"unexpected value for parm.parm1.T"
        if parm.parm1.x1[0]!=-1.: raise IOError,"unexpected value for parm.parm1.x1[0]"
        if parm.parm1.x1[1]!=2.: raise IOError,"unexpected value for parm.parm1.x1[1]"
        if parm.parm1.x2[0]!=10.: raise IOError,"unexpected value for parm.parm1.x2[0]"
        if parm.parm1.x2[1]!=11.: raise IOError,"unexpected value for parm.parm1.x2[1]"
        if parm.parm1.x3[0]!=-10: raise IOError,"unexpected value for parm.parm1.x3[0]"
        if parm.parm1.parm11.gamma1!=1.: raise IOError,"unexpected value for parm.parm1.parm11.gamma1"
        if parm.parm1.parm11.gamma2!=2.: raise IOError,"unexpected value for parm.parm1.parm11.gamma2"
        if parm.parm1.parm11.gamma3!=3.: raise IOError,"unexpected value for parm.parm1.parm11.gamma3"
        if parm.parm2.alpha!=814.: raise IOError,"unexpected value for parm.parm2.alpha"

    print "@@@ explicit construction"  
    parm=ESySParameters()
    parm.parm1=ESySParameters()
    parm.parm1.gamma=0.
    parm.parm1.dim=2
    parm.parm1.tol_v=0.001
    parm.parm1.output_file="/tmp/u.%3.3d.dx"
    parm.parm1.runFlag=True
    parm.parm1.T=1.
    parm.parm1.x1=[-1.,2]
    parm.parm1.x2=(10,11.)
    parm.parm1.x3=(-10.,)
    parm.parm1.parm11=ESySParameters()
    parm.parm1.parm11.gamma1=1.
    parm.parm1.parm11.gamma2=2.
    parm.parm1.parm11.gamma3=3.
    parm.parm2=ESySParameters()
    parm.parm2.alpha=814.
    print parm
    test(parm)
    print "@@@ read and write:"
    parm.writeESySXML(file("/tmp/test.xml",mode="w"))
    parm.writeProperties(file("/tmp/test.dat",mode="w"))
    parm2=readESySXMLFile("/tmp/test.xml")
    print parm2
    test(parm2)
    print "@@@ file"
    file("/tmp/test2.xml","w").write("""<?xml version="1.0"?>
<ESyS>
   <Component type="Geodynamics">
   <Name>parm1</Name>
   <Description>
     a few examples of parameters
   </Description>
   <Parameter><Item>gamma</Item><Value>0.</Value></Parameter>
   <Parameter type="int"><Item>dim</Item><Value>2</Value></Parameter>
   <Parameter type="real"><Item>tol_v</Item><Value>0.001</Value></Parameter>
   <Parameter type="string"><Item>output_file</Item><Value>/tmp/u.%3.3d.dx</Value></Parameter>
   <Parameter type="bool"><Item>runFlag</Item><Value>true</Value></Parameter>
   <Parameter type="real" sequence="single"><Item>T</Item><Value>1.</Value><Value>2</Value></Parameter>
   <Parameter type="real" sequence="list"><Item>x1</Item><Value>-1.</Value><Value>2</Value></Parameter>
   <Parameter type="real" sequence="tuple"><Item>x2</Item><Value>10</Value><Value>11</Value></Parameter>
   <Parameter sequence="tuple"><Item>x3</Item><Value>-10</Value></Parameter>
   <Component>
           <Name>parm11</Name>
           <Description>
          a sub compoment
           </Description>
            <Parameter><Item>gamma1</Item><Value>1.</Value></Parameter>
            <Parameter><Item>gamma2</Item><Value>2.</Value></Parameter>
            <Parameter><Item>gamma3</Item><Value>3.</Value></Parameter>
           </Component>   
   </Component>
    <Component type="Geodynamics">
   <Name>parm2</Name>
   <Description>
     another component
   </Description>
   <Parameter><Item>alpha</Item><Value>0814</Value></Parameter>
   </Component>
</ESyS>
""")
    parm3=readESySXMLFile("/tmp/test2.xml")
    print parm3
    test(parm3)
