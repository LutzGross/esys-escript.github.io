
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

from __future__ import print_function, division

__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://launchpad.net/escript-finley"

from esys.escript.modelframe import Model, Link
from esys.escript import Data
from esys.weipa import saveVTK
import os

class Visualization(Model):
    """
    Generic visualization Model

    :ivar t: current time (in)
    :ivar n: frame counter (in)
    :ivar dt: increment for output  (in)
    :ivar filename: name of the output file (in)
    """

    def __init__(self, **kwargs):
        """
        Initialisation of the visualisation model object

        :keyword debug: Debugging flag
        """
        super(Visualization,self).__init__(**kwargs)
        self.declareParameter(t=0.,
                              n=0,
                              dt=1,
                              filename="movie.mpg")

    def doInitialization(self):
       """
       does some kind of initialisation
       """
       self.__last_t=self.t
       fnc=self.filename.split('.')
       if len(fnc)==0:
           self.__frame_filename="data.%s.xml"
       else:
           n=fnc[0]
           for i in range(1,len(fnc)-1):
              n+="."+fnc[i]
           if len(fnc)==1:
              self.__frame_filename=n+".%s"
           else:
              self.__frame_filename=n+".%s."+fnc[-1]
       self.trace("output filename is %s."%self.__frame_filename)

    def getFrameFileName(self):
         return self.__frame_filename%self.getFrameCounter()

    def writeFrame(self):
       """
       returns True if the time stamp for writing frame is reached.
       """
       out=self.t>=self.__last_t+self.dt
       if out: 
            self.__last_t+=self.dt
            self.n+=1
       return out

    def getFrameCounter(self):
       """
       returns a frame counter
       """
       return self.n

    def getSafeTimeStepSize(self,dt):
       """
       returns new step size

       :param dt: last time step size used
       :type dt: ``float``
       :return: time step size that can savely be used
       :rtype: ``float``
       """
       return self.__last_t+self.dt-self.t

class WriteVTK(Visualization):
    """
    Writes data into a VTK file.

    The Model can handel up to 20 data sets that are written into a single file tagged with the given name. If no name is given and 
    the data are defined by a `Link` the name of the target attribute is used as a tag. 

    :ivar data0: data set 0 to be written
    :type data0: `escript.Data`
    :ivar name0: name tag for data set 0
    :type name0: ``str`` or ``None``
    :ivar data1: data set 1 to be written
    :type data1: `escript.Data`
    :ivar name1: name tag for data set 1
    :type name1: ``str`` or ``None``
    :ivar data2: data set 2 to be written
    :type data2: `escript.Data`
    :ivar name2: name tag for data set 2
    :type name2: ``str`` or ``None``
    :ivar data3: data set 3 to be written
    :type data3: `escript.Data`
    :ivar name3: name tag for data set 3
    :type name3: ``str`` or ``None``
    :ivar data4: data set 4 to be written
    :type data4: `escript.Data`
    :ivar name4: name tag for data set 4
    :type name4: ``str`` or ``None``
    :ivar data5: data set 5 to be written
    :type data5: `escript.Data`
    :ivar name5: name tag for data set 5
    :type name5: ``str`` or ``None``
    :ivar data6: data set 6 to be written
    :type data6: `escript.Data`
    :ivar name6: name tag for data set 6
    :type name6: ``str`` or ``None``
    :ivar data7: data set 7 to be written
    :type data7: `escript.Data`
    :ivar name7: name tag for data set 7
    :type name7: ``str`` or ``None``
    :ivar data8: data set 8 to be written
    :type data8: `escript.Data`
    :ivar name8: name tag for data set 8
    :type name8: ``str`` or ``None``
    :ivar data9: data set 9 to be written
    :type data9: `escript.Data`
    :ivar name9: name tag for data set 9
    :type name9: ``str`` or ``None``
    :ivar data10: data set 10 to be written
    :type data10: `escript.Data`
    :ivar name10: name tag for data set 10
    :type name10: ``str`` or ``None``
    :ivar data11: data set 11 to be written
    :type data11: `escript.Data`
    :ivar name11: name tag for data set 11
    :type name11: ``str`` or ``None``
    :ivar data12: data set 12 to be written
    :type data12: `escript.Data`
    :ivar name12: name tag for data set 12
    :type name12: ``str`` or ``None``
    :ivar data13: data set 13 to be written
    :type data13: `escript.Data`
    :ivar name13: name tag for data set 13
    :type name13: ``str`` or ``None``
    :ivar data14: data set 14 to be written
    :type data14: `escript.Data`
    :ivar name14: name tag for data set 14
    :type name14: ``str`` or ``None``
    :ivar data15: data set 15 to be written
    :type data15: `escript.Data`
    :ivar name15: name tag for data set 15
    :type name15: ``str`` or ``None``
    :ivar data16: data set 16 to be written
    :type data16: `escript.Data`
    :ivar name16: name tag for data set 16
    :type name16: ``str`` or ``None``
    :ivar data17: data set 17 to be written
    :type data17: `escript.Data`
    :ivar name17: name tag for data set 17
    :type name17: ``str`` or ``None``
    :ivar data18: data set 18 to be written
    :type data18: `escript.Data`
    :ivar name18: name tag for data set 18
    :type name18: ``str`` or ``None``
    :ivar data19: data set 19 to be written
    :type data19: `escript.Data`
    :ivar name19: name tag for data set 19
    :type name19: ``str`` or ``None``
    """
    def __init__(self, **kwargs):
        """
        Initialisation of the WriteVTK object

        :keyword debug: debugging flag
        :type debug: ``bool``
        """
        super(WriteVTK,self).__init__(**kwargs)
        self.declareParameter(data0=None,name0=None,
                              data1=None,name1=None,
                              data2=None,name2=None,
                              data3=None,name3=None,
                              data4=None,name4=None,
                              data5=None,name5=None,
                              data6=None,name6=None,
                              data7=None,name7=None,
                              data8=None,name8=None,
                              data9=None,name9=None,
                              data10=None,name10=None,
                              data11=None,name11=None,
                              data12=None,name12=None,
                              data13=None,name13=None,
                              data14=None,name14=None,
                              data15=None,name15=None,
                              data16=None,name16=None,
                              data17=None,name17=None,
                              data18=None,name18=None,
                              data19=None,name19=None)
    def collectData(self):
        kwargs={}
        if isinstance(self.data0, Data):
            if self.name0 == None:
               if isinstance(self.getAttributeObject("data0"),Link):
                  kwargs[self.getAttributeObject("data0").getAttributeName()]=self.data0
               else:
                  kwargs["data0"]=self.data0
            else:
               kwargs[str(self.name0)]=self.data0


        if isinstance(self.data1, Data):
            if self.name1 == None:
               if isinstance(self.getAttributeObject("data1"),Link):
                  kwargs[self.getAttributeObject("data1").getAttributeName()]=self.data1
               else:
                  kwargs["data1"]=self.data1
            else:
               kwargs[str(self.name1)]=self.data1


        if isinstance(self.data2, Data):
            if self.name2 == None:
               if isinstance(self.getAttributeObject("data2"),Link):
                  kwargs[self.getAttributeObject("data2").getAttributeName()]=self.data2
               else:
                  kwargs["data2"]=self.data2
            else:
               kwargs[str(self.name2)]=self.data2


        if isinstance(self.data3, Data):
            if self.name3 == None:
               if isinstance(self.getAttributeObject("data3"),Link):
                  kwargs[self.getAttributeObject("data3").getAttributeName()]=self.data3
               else:
                  kwargs["data3"]=self.data3
            else:
               kwargs[str(self.name3)]=self.data3


        if isinstance(self.data4, Data):
            if self.name4 == None:
               if isinstance(self.getAttributeObject("data4"),Link):
                  kwargs[self.getAttributeObject("data4").getAttributeName()]=self.data4
               else:
                  kwargs["data4"]=self.data4
            else:
               kwargs[str(self.name4)]=self.data4


        if isinstance(self.data5, Data):
            if self.name5 == None:
               if isinstance(self.getAttributeObject("data5"),Link):
                  kwargs[self.getAttributeObject("data5").getAttributeName()]=self.data5
               else:
                  kwargs["data5"]=self.data5
            else:
               kwargs[str(self.name5)]=self.data5


        if isinstance(self.data6, Data):
            if self.name6 == None:
               if isinstance(self.getAttributeObject("data6"),Link):
                  kwargs[self.getAttributeObject("data6").getAttributeName()]=self.data6
               else:
                  kwargs["data6"]=self.data6
            else:
               kwargs[str(self.name6)]=self.data6


        if isinstance(self.data7, Data):
            if self.name7 == None:
               if isinstance(self.getAttributeObject("data7"),Link):
                  kwargs[self.getAttributeObject("data7").getAttributeName()]=self.data7
               else:
                  kwargs["data7"]=self.data7
            else:
               kwargs[str(self.name7)]=self.data7


        if isinstance(self.data8, Data):
            if self.name8 == None:
               if isinstance(self.getAttributeObject("data8"),Link):
                  kwargs[self.getAttributeObject("data8").getAttributeName()]=self.data8
               else:
                  kwargs["data8"]=self.data8
            else:
               kwargs[str(self.name8)]=self.data8


        if isinstance(self.data9, Data):
            if self.name9 == None:
               if isinstance(self.getAttributeObject("data9"),Link):
                  kwargs[self.getAttributeObject("data9").getAttributeName()]=self.data9
               else:
                  kwargs["data9"]=self.data9
            else:
               kwargs[str(self.name9)]=self.data9


        if isinstance(self.data10, Data):
            if self.name10 == None:
               if isinstance(self.getAttributeObject("data10"),Link):
                  kwargs[self.getAttributeObject("data10").getAttributeName()]=self.data10
               else:
                  kwargs["data10"]=self.data10
            else:
               kwargs[str(self.name10)]=self.data10


        if isinstance(self.data11, Data):
            if self.name11 == None:
               if isinstance(self.getAttributeObject("data11"),Link):
                  kwargs[self.getAttributeObject("data11").getAttributeName()]=self.data11
               else:
                  kwargs["data11"]=self.data11
            else:
               kwargs[str(self.name11)]=self.data11


        if isinstance(self.data12, Data):
            if self.name12 == None:
               if isinstance(self.getAttributeObject("data12"),Link):
                  kwargs[self.getAttributeObject("data12").getAttributeName()]=self.data12
               else:
                  kwargs["data12"]=self.data12
            else:
               kwargs[str(self.name12)]=self.data12


        if isinstance(self.data13, Data):
            if self.name13 == None:
               if isinstance(self.getAttributeObject("data13"),Link):
                  kwargs[self.getAttributeObject("data13").getAttributeName()]=self.data13
               else:
                  kwargs["data13"]=self.data13
            else:
               kwargs[str(self.name13)]=self.data13


        if isinstance(self.data14, Data):
            if self.name14 == None:
               if isinstance(self.getAttributeObject("data14"),Link):
                  kwargs[self.getAttributeObject("data14").getAttributeName()]=self.data14
               else:
                  kwargs["data14"]=self.data14
            else:
               kwargs[str(self.name14)]=self.data14


        if isinstance(self.data15, Data):
            if self.name15 == None:
               if isinstance(self.getAttributeObject("data15"),Link):
                  kwargs[self.getAttributeObject("data15").getAttributeName()]=self.data15
               else:
                  kwargs["data15"]=self.data15
            else:
               kwargs[str(self.name15)]=self.data15


        if isinstance(self.data16, Data):
            if self.name16 == None:
               if isinstance(self.getAttributeObject("data16"),Link):
                  kwargs[self.getAttributeObject("data16").getAttributeName()]=self.data16
               else:
                  kwargs["data16"]=self.data16
            else:
               kwargs[str(self.name16)]=self.data16


        if isinstance(self.data17, Data):
            if self.name17 == None:
               if isinstance(self.getAttributeObject("data17"),Link):
                  kwargs[self.getAttributeObject("data17").getAttributeName()]=self.data17
               else:
                  kwargs["data17"]=self.data17
            else:
               kwargs[str(self.name17)]=self.data17


        if isinstance(self.data18, Data):
            if self.name18 == None:
               if isinstance(self.getAttributeObject("data18"),Link):
                  kwargs[self.getAttributeObject("data18").getAttributeName()]=self.data18
               else:
                  kwargs["data18"]=self.data18
            else:
               kwargs[str(self.name18)]=self.data18


        if isinstance(self.data19, Data):
            if self.name19 == None:
               if isinstance(self.getAttributeObject("data19"),Link):
                  kwargs[self.getAttributeObject("data19").getAttributeName()]=self.data19
               else:
                  kwargs["data19"]=self.data19
            else:
               kwargs[str(self.name19)]=self.data19
        return kwargs

    def doInitialPostprocessing(self):
        """
        writes vtk file at the end of initial iteration
        """
        super(WriteVTK,self).doInitialPostprocessing()
        kwargs=self.collectData()
        if len(kwargs)>0:
           saveVTK(self.getFrameFileName(),**kwargs)
           self.trace("%s-th frame at time %s is writen to %s"%(self.getFrameCounter(),self.t,self.getFrameFileName()))

    def doStepPostprocessing(self, dt):
        """
        writes vtk file at the end of time iteration
        """
        super(WriteVTK,self).doStepPostprocessing(dt)
        if self.writeFrame():
            kwargs=self.collectData()
            if len(kwargs)>0:
               saveVTK(self.getFrameFileName(),**kwargs)
               self.trace("%s-th frame at time %s is writen to %s"%(self.getFrameCounter(),self.t,self.getFrameFileName()))

# vim: expandtab shiftwidth=4:
