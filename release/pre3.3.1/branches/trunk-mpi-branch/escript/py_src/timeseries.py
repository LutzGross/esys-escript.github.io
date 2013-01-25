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
Time serieas analysis

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
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision$"
__date__="$Date$"


import numarray
from types import SliceType
DEFAULT_BUFFER_SIZE=1000
DEFAULT_FLOAT_TYPE=numarray.Float64

class TimeSeriesBase:
   """The TimeSeriesBase class is the base class for all class of the TimeSeries module."""

   def __init__(self,debug=False,description="TimeSeriesBase"):
       self.__debug=debug
       self.setDescription(description)

   def __str__(self):
       return self.__description
   
   def setDescription(self,text):
       self.__description=text

   def setDebugOn(self):
      """switch on degugging mode"""
      self.__debug=True

   def setDebugOff(self):
      """switch off degugging mode"""
      self.__debug=False

   def setDebug(self,flag=False):
      """sets debug mode to flag"""
      if flag:
         self.setDebugOn()
      else:
         self.setDebugOff()
         
   def debug(self):
      """returns true if debug mode is on"""
      return self.__debug

#============================================================================================================
class TimeSeriesBaseDataset(TimeSeriesBase):
   """provides an interface for accessing a set of linearly ordered data."""
   def __init__(self,buffer,offset=0,debug=False,description="TimeSeriesDataset"):
       TimeSeriesBase.__init__(self,debug,description)
       self.__buffer=buffer
       self.__offset=offset
       if self.debug(): print "Debug: %s: offset %d to buffer"%(self,self.getOffset())

   def __len__(self):
       """needed to handle negative indexing in slicing"""
       return 0

   def getNumComponents(self):
       """returns the number of components of the data (may be overwritten by subclass)"""
       return self.getBaseBuffer().getNumComponents()

   def getIdOfLastDatum(self):
      """returns the identification number of the last datum in the data set (may be overwritten by subclass)"""
      return self.getBaseBuffer().getIdOfLastDatum()-self.getOffset()

   def getIdOfFirstDatum(self):
      """returns the identification number of the first datum (may be overwritten by subclass)"""
      return self.getBaseBuffer().getIdOfFirstDatum()-self.getOffset()

   def getIdOfFirstAvailableDatum(self):
      """returns the identification number of the first avaiable datum (may be overwritten by subclass)"""
      return self.getBaseBuffer().getIdOfFirstAvailableDatum()-self.getOffset()

   def getOffsetInBaseBuffer(self):
      """returns the offset to access elements in getBaseBuffer() (may be overwritten by subclass)"""
      return  self.getOffset()

   def getIdOfLastUnreferencedDatum(self):
       """returns the identification number of the last datum which has been unused by all TimeSeries refering to the TimeSeriesBaseDataset (may be overwritten by subclass)"""
       return self.getBaseBuffer().getIdOfLastUnreferencedDatum()-self.getOffset()

   def updateIdOfLastUnreferencedDatum(self,last_unreferenced_datum):
       """updates the identification number of the last unused datum (to be overwritten by subclass)"""
       self.getBaseBuffer().updateIdOfLastUnreferencedDatum(last_unreferenced_datum+self.getOffset())

   def append(self,values):
       """appends data to the buffer. If the buffer would be full the buffer is rearranged before the data are appended  (to be overwritten by subclass)"""
       self.getBaseBuffer().append(values)

   def getBaseBufferSize(self):
       """returns the size of the buffer (to be overwritten by subclass)"""
       return self.getBaseBuffer().getBaseBufferSize()
   
   def needsRearrangement(self,num_new_data=0):
       """returns True if the buffer will be full after num_new_data have been appended (to be overwritten by subclass)"""
       return self.getBaseBuffer().needsRearrangement(num_new_data)

   def isEmpty(self):
      """returns true if no data are appeneded to buffer"""
      return self.getNumData()<=0
   
   def getNumData(self):
      """returns the number of data (not all of them are accessible)"""
      return self.getIdOfLastDatum()-self.getIdOfFirstDatum()+1

   def getBaseBuffer(self):
      """return the buffer referenced by the TimeSeriesBaseDataset"""
      return self.__buffer

   def getOffset(self):
      """return the offset when referring to dataset elements"""
      return self.__offset

   def __getitem__(self,index):
      """returns the datum index"""
      if type(index)==SliceType:
         start=index.start
         end=index.stop
         if start==end: 
            return self[start]
         else:
             if start<self.getIdOfFirstDatum() or start>self.getIdOfLastDatum() or \
                 end-1<self.getIdOfFirstDatum() or end-1>self.getIdOfLastDatum(): raise IndexError,"%s: Index [%d:%d] out of range"%(self,start,end)
             return self.getBaseBuffer()[start+self.getOffsetInBaseBuffer():end+self.getOffsetInBaseBuffer()]
      else:
         if index<self.getIdOfFirstDatum() or index>self.getIdOfLastDatum(): raise IndexError,"%s: Index %d out of range"%(self,index)
         return self.getBaseBuffer()[index+self.getOffsetInBaseBuffer()]

class TimeSeriesBaseBuffer(TimeSeriesBaseDataset):
   """An inplementation of TimeSeriesBaseDataset which actually is storing data into a numarray buffer"""
   def __init__(self,buffer_size=DEFAULT_BUFFER_SIZE,numComponents=1,type=DEFAULT_FLOAT_TYPE,id_of_first_datum=0,debug=False,description="TimeSeriesBaseBuffer"):
       if numComponents<2:
          buffer=numarray.zeros((buffer_size,),type)
       else:
          buffer=numarray.zeros((buffer_size,numComponents),type)
       TimeSeriesBaseDataset.__init__(self,buffer,id_of_first_datum-1,debug,description)
       self.__num_data_in_buffer=0
       self.__id_last_unreferenced_datum=id_of_first_datum-1
       self.__id_last_datum=id_of_first_datum-1
       self.__id_first_datum=id_of_first_datum
       if self.debug(): print "Debug: %s : buffer of size %d with %d components allocated (first datum is %d)."% \
                       (self,self.getBaseBufferSize(),self.getNumComponents(),id_of_first_datum)


   def getBaseBufferSize(self):
       """returns the size of the buffer"""
       return self.getBaseBuffer().shape[0]
   
   def getNumComponents(self):
       """returns the number of components of the data (overwrites TimeSeriesBaseDataset method)"""
       if self.getBaseBuffer().rank==1:
          return 1
       else: 
          self.getBaseBuffer().shape[1]

   def getNumDataInBaseBuffer(self):
       """returns the number of data currently in the buffer"""
       return self.__num_data_in_buffer

   def getIdOfLastDatum(self):
      """returns the identification number of the last datum in the data set (overwrites method from TimeSeriesBaseDataset)"""
      return self.__id_last_datum

   def getIdOfFirstDatum(self):
      """returns the identification number of the first datum (overwrites method from TimeSeriesBaseDataset)"""
      return self.__id_first_datum

   def getOffsetInBaseBuffer(self):
      """returns the offset to access elements in the buffer (overwrites method from TimeSeriesBaseDataset)"""  
      return -self.getIdOfLastDatum()+self.getNumDataInBaseBuffer()-1  

   def getIdOfLastUnreferencedDatum(self):
       """returns the identification number of the last datum which has been unused by all TimeSeries refering to the TimeSeriesBaseDataset (overwrites method from TimeSeriesBaseDataset)"""
       return self.__id_last_unreferenced_datum

   def updateIdOfLastUnreferencedDatum(self,last_unreferenced_datum):
       """updates the identification number of the last unused datum (to be overwritten by subclass)"""
       self.getBaseBuffer().updateIdOfLastUnreferencedDatum(last_unreferenced_datum-self.getOffset())

   def updateIdOfLastUnreferencedDatum(self,last_unreferenced_datum):
       """updates the identification number of the last unused datum (overwrites TimeSeriesBaseDataset method)"""
       if self.__id_last_unreferenced_datum>last_unreferenced_datum:
           self.__id_last_unreferenced_datum=last_unreferenced_datum
           if self.debug(): print "Debug: %s: last unused datum is now %s"%(self,last_unreferenced_datum)

   def needsRearrangement(self,num_new_data=0):
       """returns True if the buffer will be full after num_new_data have been appended"""
       return self.getNumDataInBaseBuffer()+num_new_data>self.getBaseBufferSize()
         
   def getIdOfFirstAvailableDatum(self):
      """returns the identification number of the first avaiable datum (overwrites TimeSeriesBaseDataset method)"""
      return self.getIdOfLastDatum()-self.__num_data_in_buffer+1

   def append(self,data):
      """appends data to the buffer. If the buffer would be full the buffer is rearranged before the data are appended (overwrites TimeSeriesBaseDataset method)"""
      data=numarray.array(data)
      nc=self.getNumComponents()
      if data.rank==0:
        if nc==1:
           num_new_data=1
        else:
           raise ValueError,"%s: illegal data shape"%self
      elif data.rank==1:
        if nc==1:
             num_new_data=data.shape[0]
        else:
             num_new_data=1  
      elif data.rank==2: 
        if not nc==data.shape[1]: raise ValueError,"%s: illegal data shape"%self
        num_new_data=data.shape[0]
      else:
         raise ValueError,"%s: illegal rank"%self

      # check is buffer will be overflown when data are appended:
      if self.needsRearrangement(num_new_data):
        nn=self.getNumDataInBaseBuffer()
        num_protected_data=self.getIdOfLastDatum()-self.getIdOfLastUnreferencedDatum()
        if num_protected_data+num_new_data>self.getBaseBufferSize():
              raise ValueError,"%s: buffer overflow: buffer size has to be bigger than %d"%(self,num_protected_data+num_new_data)
        if num_protected_data>0: self.getBaseBuffer()[0:num_protected_data]=self.getBaseBuffer()[nn-num_protected_data:nn]
        self.__num_data_in_buffer=num_protected_data
        self.__id_last_unreferenced_datum=self.__id_last_datum
        if self.debug(): 
             print "Debug: %s: rearrangement: first data in buffer is %d."%(self,self.getIdOfLastDatum()-self.getNumDataInBaseBuffer()+1)
      # copy data over:
      nn=self.getNumDataInBaseBuffer()
      self.getBaseBuffer()[nn:nn+num_new_data]=data
      self.__num_data_in_buffer+=num_new_data
      self.__id_last_datum+=num_new_data
      self.__id_last_unreferenced_datum+=num_new_data
      if self.debug(): print "Debug: %s: %d data appended. Last unreferenced datum is now %d."%(self,num_new_data,self.__id_last_unreferenced_datum)

# ======================================
class TimeSeriesControlerView(TimeSeriesBase):
      """A TimeSeriesControlerView is attached to a Controler and moves forward in time by increasing the id of the last processed datum.
         Any implementation of a TimeSeriesControlerView must provide the getControler method which returns the controler"""
      def __init__(self,id_first_datum=0,debug=False,description="TimeSeries"):
        TimeSeriesBase.__init__(self,debug,description)
        self.__id_last_processed_datum=id_first_datum-1
        if self.debug(): print "Debug: %s  created with first datum %d"%(str(self),id_first_datum)

      def getIdOfLastProcessedDatum(self):
          return self.__id_last_processed_datum

      def updateIdOfLastProcessedDatum(self,id_last_processed_datum):
          self.__id_last_processed_datum=id_last_processed_datum

      # def getControler(self):
      #      """returns the Controler of the time series (to be overwritten by subclass)"""
      #      pass

class TimeSeries(TimeSeriesBaseDataset,TimeSeriesControlerView):
      """makes TimeSeriesBaseDataset look like a TimeSeries and introduces operations
         Any implementation of a TimeSeriesControlerView must provide the getControler method which returns the controler"""
      def __init__(self,dataset,debug=False,description="TimeSeries"):
        TimeSeriesControlerView.__init__(self,dataset.getIdOfFirstDatum(),debug,description)
        TimeSeriesBaseDataset.__init__(self,dataset,0,debug,description)
      
      def getDataset(self):
          """returns the TimeSeriesBaseDataset of the time series"""
          return self.getBaseBuffer()

      # def getControler(self):
      #      """returns the Controler of the time series (to be overwritten by subclass)"""
      #      pass

      def __add__(self,arg):
         if isinstance(arg,TimeSeriesBaseDataset):
            return TimeSeriesAdd(self,arg)
         else:
            return TimeSeriesAddScalar(self,arg)

      def __sub__(self,arg):
         return self+(-1.)*arg

      def __mul__(self,arg):
         if isinstance(arg,TimeSeriesBaseDataset):
            return TimeSeriesMult(self,arg)
         else:
            return TimeSeriesMultScalar(self,arg)

      def __div__(self,arg):
         if isinstance(arg,TimeSeriesBaseDataset):
            return TimeSeriesDiv(self,arg)
         else:
            return TimeSeriesMultScalar(self,1./arg)

      def __pow__(self,arg):
         if isinstance(arg,TimeSeriesBaseDataset):
            return TimeSeriesPower(self,arg)
         else:
            return TimeSeriesPowerScalar(self,arg)
      
      def __radd__(self,arg):
         return self.__add__(arg)

      def __rsub__(self,arg):
         return arg+(-1.)*self

      def __rmul__(self,arg):
         return self.__mul__(arg)

      def __rdiv__(self,arg):
         if isinstance(arg,TimeSeriesBaseDataset):
            return TimeSeriesDiv(arg,self)
         else:
            return TimeSeriesDivScalar(self,arg)

      def __rpow__(self,arg):
         if isinstance(arg,TimeSeriesBaseDataset):
            return TimeSeriesPower(arg,self)
         else:
            return Exp(numarray.log(arg)*self)

      def __lshift__(self,arg):
         return TimeSeriesShift(self,-arg)

      def __rshift__(self,arg):
         return TimeSeriesShift(self,arg)

      def __neg__(self):
         return (-1.0)*self

      def __pos__(self):
         return (1.0)*self

class TimeSeriesOperator(TimeSeriesControlerView):
      """a TimeSeriesOperator decribes an operation acting on list of TimeSeries time_series_args. It allows to update its output (if there is any)
         through the update method which is overwritten by a particular implementation of the class. The update method is called to process the data [start:end] using
         [start-left_wing_size:end+right_wing_size] of its arguments"""
      def __init__(self,controler,time_series_args=[],left_wing_size=0,right_wing_size=0,debug=False,description="TimeSeriesOperator"):
          id_first_datum=controler.getIdOfFirstDatum()
          for i in time_series_args: id_first_datum=max(id_first_datum,i.getIdOfFirstDatum())
          TimeSeriesControlerView.__init__(self,id_first_datum+left_wing_size,debug,description)
          self.__left_wing_size=left_wing_size
          self.__right_wing_size=right_wing_size
          self.__time_series_args=time_series_args
          self.__controler=controler
          controler.appendOperatorToUpdateList(self)
          if self.debug(): print "Debug: %s: with left/right wing size %d/%d and %d arguments."%(str(self),left_wing_size,right_wing_size,len(time_series_args))

      def __del__(self):
          self.getControler().removeOperatorFromUpdateList(self)

      def getControler(self):
          """returns the Controler updating the TimeSeriesOperator"""
          return self.__controler

      def getLeftWingSize(self):
          """returns the left wing size"""  
          return self.__left_wing_size

      def getRightWingSize(self):
          """returns the right wing size""" 
          return self.__right_wing_size

      def getArguments(self,index=None):
          """returns the list of arguments or, index is present, the argument with index index. In the latter case None is returned if no arguments are present"""
          if index==None:
             return self.__time_series_args
          else:
             if len(self.__time_series_args)>0:
                return self.__time_series_args[index]
             else:
                return None

      def getArgumentDataset(self,index):
          """returns the dataset of in the argument with index index"""
          arg=self.getArguments(index)
          if arg==None: 
             return None
          else:
              return self.getArguments(index).getDataset()

      def flush(self):
          """calls the update method with all the maximum processable range. It also updates the id of unused datum for all arguments"""
          start=self.getIdOfLastProcessedDatum()+1
          end=self.getControler().getIdOfLastDatum()
          for i in self.getArguments(): end=min(end,i.getIdOfLastDatum())
          if start<=end-self.getRightWingSize():
             if self.debug(): print "Debug: %s: range [%d:%d] is updated."%(self,start,end-self.getRightWingSize())
             self.update(start,end-self.getRightWingSize()+1)      
             for i in self.getArguments(): i.updateIdOfLastUnreferencedDatum(end-self.getLeftWingSize())
             self.updateIdOfLastProcessedDatum(end)

      def update(self,start,end):
          """updates the the data [start:end] using [start-left_wing_size:end+right_wing_size] of its arguments (is overwritten by a particular TimeSeriesOperator)"""
          pass


class TimeSeriesFilter(TimeSeries,TimeSeriesOperator):
      """a TimeSeriesFilter is a TimeSeries taht is created trough a TimeSeriesOperator"""
      def __init__(self,controler,dataset,time_series_args=[],left_wing_size=0,right_wing_size=0,debug=False,description="TimeSeriesFilter"):
         TimeSeriesOperator.__init__(self,controler,time_series_args,left_wing_size,right_wing_size,debug,description)
         TimeSeries.__init__(self,dataset,debug,description)

      def update(self,start,end):
          """appends zeros to the dataset. This method should be overwritten by a particular TimeSeriesFilter"""
          nc=self.getNumComponents()
          if nc>1:
             self.getDataset().append(numarray.zeros([nc,end-start]))
          else:
             self.getDataset().append(numarray.zeros(end-start))

class Controler(TimeSeries):
   """controls a set of TimeSeries"""
   def __init__(self,buffer_size=DEFAULT_BUFFER_SIZE,debug=False,description="TimeSeriesControler"): 
        TimeSeries.__init__(self,TimeSeriesBaseBuffer(buffer_size,1,DEFAULT_FLOAT_TYPE,0,debug,"node buffer of "+description),debug,"nodes of "+description)
        self.setFlushRate()   
        self.__update_time_series=list()
      
   def getControler(self):
       """returns the Controler of the time series (overwrites method of by TimeSeries)"""
       return self

   def setFlushRate(self,rate=50):
       """set the flush rate, i.e. after rate new time nodes have been checked in the flush method is called."""
       self.__flush_rate=rate
       if self.debug(): print "Debug: %s: flush rate is set to %d"%(self,rate)
 
   def needsFlushing(self):
      """returns true if the depending TimeSeriesFilters needs to be flushed becuase the time nodes buffer is full or because of the set flush rate"""
      return self.needsRearrangement(1) or (self.getNumData()+1)%self.__flush_rate==0

   def flush(self):
       """flushes all dependend TimeSeriesFilters by processing their flush method"""
       if self.debug(): print "Debug: %s: start flushing"%self
       for time_serie in self.__update_time_series: time_serie.flush()

   def appendOperatorToUpdateList(self,time_serie):
       if not time_serie.getControler()==self: raise ValueError,"%s: TimeSeries %s is not defined on this controler."%(self,time_serie)
       if not self.isEmpty(): raise ValueError,"%s: you can only check in a time series time_serie is controler is empty."%self
       self.__update_time_series.append(time_serie)
       if self.debug(): print "Debug: %s: %s has been added to update list."%(self,time_serie)

   def removeOperatorFromUpdateList(self,time_serie):
       self.__update_time_series.remove(time_serie)
       if self.debug(): print "Debug: %s: %s has been removed from update list."%(self,time_serie)

   def nextTime(self,value):
       if self.needsFlushing(): self.flush()
       self.getDataset().append(value)
       if self.debug(): print "Debug: %s: new time node %e has been added."%(self,value)

class TimeSeriesShift(TimeSeries):
      """creates a shift of the time series, i.e. if d[n] is the datum at time t[n], the value at t[n] becomes v[n+shift] on the output"""
      def __init__(self,time_serie,shift=1):
          if shift<0:
              dsc="(%s)<<%d"%(time_serie,-shift)
          else:
              dsc="(%s)>>%d"%(time_serie,shift)
          self.__controler=time_serie.getControler()
          TimeSeries.__init__(self,TimeSeriesBaseDataset(time_serie.getDataset(),-shift,time_serie.debug(),"buffer view to "+dsc),time_serie.debug(),dsc)

      def getControler(self):
          return self.__controler

class TimeSeriesAdd(TimeSeriesFilter):
      """adds two TimeSeries"""
      def __init__(self,time_serie_1,time_serie_2):
          dsc="(%s)+(%s)"%(time_serie_1,time_serie_2)
          dbg=time_serie_1.debug() or time_serie_2.debug()
          cntrl=time_serie_1.getControler()
          if not cntrl==time_serie_2.getControler(): 
                  raise ValueError("TimeSeriesAdd: %s and %s have different controler."%(time_serie_1,time_serie_2))
          id_first_datum=max(time_serie_1.getIdOfFirstDatum(),time_serie_2.getIdOfFirstDatum())
          TimeSeriesFilter.__init__(self,cntrl, \
                              TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                              [time_serie_1,time_serie_2],0,0,dbg,dsc)

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]+self.getArgumentDataset(1)[start:end])

class TimeSeriesAddScalar(TimeSeriesFilter):
      """adds a single value to a TimeSeries"""
      def __init__(self,time_serie,scalar):
          dsc="(%s)+(%s)"%(time_serie,scalar)
          dbg=time_serie.debug()
          cntrl=time_serie.getControler()
          id_first_datum=time_serie.getIdOfFirstDatum()
          TimeSeriesFilter.__init__(self,cntrl, \
                       TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                       [time_serie],0,0,dbg,dsc)
          self.__scalar=scalar

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]+self.__scalar)

class TimeSeriesMult(TimeSeriesFilter):
      """multiplies two TimeSeries"""
      def __init__(self,time_serie_1,time_serie_2):
          dsc="(%s)*(%s)"%(time_serie_1,time_serie_2)
          dbg=time_serie_1.debug() or time_serie_2.debug()
          cntrl=time_serie_1.getControler()
          if not cntrl==time_serie_2.getControler(): 
                  raise ValueError("TimeSeriesMult: %s and %s have different controler."%(time_serie_1,time_serie_2))
          id_first_datum=max(time_serie_1.getIdOfFirstDatum(),time_serie_2.getIdOfFirstDatum())
          TimeSeriesFilter.__init__(self,cntrl, \
                   TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                   [time_serie_1,time_serie_2],0,0,dbg,dsc)

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]*self.getArgumentDataset(1)[start:end])

class TimeSeriesMultScalar(TimeSeriesFilter):
      """multiplies a TimeSeries with a single value"""
      def __init__(self,time_serie,scalar):
          dsc="(%s)*%s"%(time_serie,scalar)
          dbg=time_serie.debug()
          cntrl=time_serie.getControler()
          id_first_datum=time_serie.getIdOfFirstDatum()
          TimeSeriesFilter.__init__(self,cntrl, \
                       TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                       [time_serie],0,0,dbg,dsc)
          self.__scalar=scalar

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]*self.__scalar)

class TimeSeriesDiv(TimeSeriesFilter):
      """divides two TimeSeries"""
      def __init__(self,time_serie_1,time_serie_2):
          dsc="(%s)/(%s)"%(time_serie_1,time_serie_2)
          dbg=time_serie_1.debug() or time_serie_2.debug()
          cntrl=time_serie_1.getControler()
          if not cntrl==time_serie_2.getControler(): 
                  raise ValueError("TimeSeriesDiv: %s and %s have different controler."%(time_serie_1,time_serie_2))
          id_first_datum=max(time_serie_1.getIdOfFirstDatum(),time_serie_2.getIdOfFirstDatum())
          TimeSeriesFilter.__init__(self,cntrl, \
                     TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                     [time_serie_1,time_serie_2],0,0,dbg,dsc)

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]/self.getArgumentDataset(1)[start:end])

class TimeSeriesDivScalar(TimeSeriesFilter):
      """divides a scalar be a TimeSerie"""
      def __init__(self,time_serie,scalar):
          dsc="(%s)/(%s)"%(scalar,time_serie)
          dbg=time_serie.debug()
          cntrl=time_serie.getControler()
          id_first_datum=time_serie.getIdOfFirstDatum()
          TimeSeriesFilter.__init__(self,cntrl, \
                       TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                       [time_serie],0,0,dbg,dsc)
          self.__scalar=scalar

      def update(self,start,end):
          self.append(self.__scalar/self.getArgumentDataset(0)[start:end])

class TimeSeriesPower(TimeSeriesFilter):
      """raise one TimeSeries to the power of an other TimeSeries"""
      def __init__(self,time_serie_1,time_serie_2):
          dsc="(%s)**(%s)"%(time_serie_1,time_serie_2)
          dbg=time_serie_1.debug() or time_serie_2.debug()
          cntrl=time_serie_1.getControler()
          if not cntrl==time_serie_2.getControler(): 
                  raise ValueError("TimeSeriesPower: %s and %s have different controler."%(time_serie_1,time_serie_2))
          id_first_datum=max(time_serie_1.getIdOfFirstDatum(),time_serie_2.getIdOfFirstDatum())
          TimeSeriesFilter.__init__(self,cntrl, \
                TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                [time_serie_1,time_serie_2],0,0,dbg,dsc)

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]**self.getArgumentDataset(1)[start:end])

class TimeSeriesPowerScalar(TimeSeriesFilter):
      """raises a TimeSerie to the power of a scalar"""
      def __init__(self,time_serie,scalar):
          dsc="(%s)**(%s)"%(time_serie,scalar)
          dbg=time_serie.debug()
          cntrl=time_serie.getControler()
          id_first_datum=time_serie.getIdOfFirstDatum()
          TimeSeriesFilter.__init__(self,cntrl, \
                       TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                       [time_serie],0,0,dbg,dsc)
          self.__scalar=scalar

      def update(self,start,end):
          self.append(self.getArgumentDataset(0)[start:end]**self.__scalar)

class Exp(TimeSeriesFilter):
      """"""
      def __init__(self,time_serie):
          dsc="exp(%s)"%(time_serie)
          dbg=time_serie.debug()
          cntrl=time_serie.getControler()
          id_first_datum=time_serie.getIdOfFirstDatum()
          TimeSeriesFilter.__init__(self,cntrl, \
                     TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                     [time_serie],0,0,dbg,dsc)

      def update(self,start,end):
          self.append(numarray.exp(self.getArgumentDataset(0)[start:end]))

class Writer(TimeSeriesOperator):
      """writes the time series into an output strim ostream which mast have the writeline method. The values are seperated by the string seperator."""
      def __init__(self,time_serie,ostream,seperator=",",commend_tag="#"):
         dsc="write %s to %s"%(time_serie,ostream)
         dbg=time_serie.debug()
         cntrl=time_serie.getControler()
         self.__ostream=ostream
         self.__seperator=seperator
         TimeSeriesOperator.__init__(self,cntrl,[time_serie],0,0,dbg,dsc)
         ostream.writelines("%s time series %s\n"%(commend_tag,str(self)))

      def update(self,start,end):
         cntrl=self.getControler()
         arg=self.getArguments(0)
         n=arg.getNumComponents()
         if n<2:
            for i in range(start,end): self.__ostream.writelines("%s%s%s\n"%(cntrl[i],self.__seperator,arg[i]))
         else:
            for i in range(start,end): 
               l="%s"%cntrl[i]
               for j in range(n): l=l+"%s%s"(self.__seperator,arg[i][j])
               self.__ostream.writelines("%s\n"%l)

class DataCatcher(TimeSeries):
      """collects data into a time series."""
      def __init__(self,controler,numComponents=1,description="DataCatcher"):
         self.__controler=controler
         dbg=controler.debug()
         TimeSeries.__init__(self,TimeSeriesBaseBuffer(controler.getBaseBufferSize(),numComponents,DEFAULT_FLOAT_TYPE,controler.getIdOfFirstDatum(),dbg,"buffer for "+description),dbg,description)

      def getControler(self):
          return self.__controler

      def nextValue(self,value):
          """append a value to the time series"""
          id_last=self.getIdOfLastDatum()
          id_current=self.getControler().getIdOfLastDatum()
          if id_last+1==id_current:
             self.getDataset().append(value)
          elif id_last+1<id_current:
               if self.isEmpty():
                   self.getDataset().append(value)
                   id_last+=1
               t_last=self.getControler()[id_last]
               t_current=self.getControler()[id_current]
               value_last=self[id_last]
               out=(value_last-value)/(t_last-t_current)*(self.getControler()[id_last+1:id_current+1]-t_current)+value
               self.getDataset().append(out)
          else :
             raise ValueError,"%s: a new time node must be introduced before a new value can be added."
          self.updateIdOfLastUnreferencedDatum(id_last)
          
  
class TimeSeriesCumulativeSum(TimeSeriesFilter):
      """cummulative sum of the time series values"""
      def __init__(self,time_serie):
         dsc="cumsum(%s)"%(time_serie)
         dbg=time_serie.debug()
         cntrl=time_serie.getControler()
         id_first_datum=time_serie.getIdOfFirstDatum()
         TimeSeriesFilter.__init__(self,cntrl, \
                     TimeSeriesBaseBuffer(cntrl.getBaseBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                     [time_serie],0,0,dbg,dsc)
         self.__last_value=0

      def update(self,start,end):
          out=numarray.cumsum(self.getArgumentDataset(0)[start:end])+self.__last_value
          self.__last_value=out[end-start-1]
          self.append(out)
         

class Reader(TimeSeriesBase):
      """reads a list of input streams and creates a time series for each input stream but on the same Controler where the first column
         is used to create the time nodes"""
      def __init__(self,list_of_istreams,buffer_size=DEFAULT_BUFFER_SIZE,seperator=",",commend_tag="#",debug=False):
         TimeSeriesBase.__init__(self,debug=debug,description="reader")
         if not isinstance(list_of_istreams,list): 
              self.__list_of_istreams=[list_of_istreams]
         else:
              self.__list_of_istreams=list_of_istreams
         self.__cntrl=Controler(buffer_size,debug,"reader controler")
         self.__seperator=seperator
         self.__commend_tag=commend_tag
         self.__time_series={}
         self.__t={}
         self.__v={}
         # set up the time series:
         for i in self.__list_of_istreams:
           line=self.__commend_tag
           while  not line=="" and line[0]==self.__commend_tag:
               line=i.readline().strip()
           if line=="": 
              list_of_istreams.remove(i)
           else:
              d=line.split(self.__seperator)
              self.__t[i]=float(d[0])
              tmp=[]
              for j in d[1:]: tmp.append(float(j))
              self.__v[i]=numarray.array(tmp)
              self.__time_series[i]=DataCatcher(self.__cntrl,len(d)-1,str(i))

         # 
      def run(self):
         while len(self.__list_of_istreams)>0:
            if len(self.__time_series)>0:
               # find list all times with minumum time node:
               tminargs=[]
               for i in self.__time_series:
                   if len(tminargs)==0:
                       tminargs.append(i)
                   elif abs(t[tminargs[0]]-self.__t[i])<1.e-8*abs(self.__t[i]):
                       tminargs.append(i)
                   elif self.__t[i]<t[tminargs[0]]:
                       tminargs=[i]
               # find list all times with minumum time node:
               self.__cntrl.nextTime(self.__t[tminargs[0]])
               for i in tminargs:
                   self.__time_series[i].nextValue(self.__v[i])
                   # find next line without leading "#"
                   line="#"
                   while not line=="" and line[0]==self.__commend_tag:
                       line=i.readline().strip()
                   # if eof reached iostream is removed for searching
                   if line=="": 
                      self.__list_of_istreams.remove(i)
                   else:
                      d=line.split(self.__seperator)
                      self.__t[i]=float(d[0])
                      tmp=[]
                      for j in d[1:]: tmp.append(float(j))
                      self.__v[i]=numarray.array(tmp)

      def getControler(self):
         """returns the controler shared by all time series created through the input streams"""
         return self.__cntrl

      def getTimeSeries(self,istream=None):
         """returns the time series as a tuple. If istream is present its time series is returned"""
         if istream==None:
            out=self.__time_series.values()
            if len(out)>1:
               return tuple(out)
            elif len(out)>0:
               return out[0]
            else:
               return None
         else:
            return self.__time_series[istream]


class Plotter(TimeSeriesOperator):
    def __init__(self,time_series,window_size=DEFAULT_BUFFER_SIZE/4,file_name=None,format=None):
         if isinstance(time_series,list):
             dbg=time_series[0].getControler().debug()
             text=""
             for i in time_series: 
               if len(text)==0:
                  text=str(i)
               else:
                  text=text+","+str(i)
             TimeSeriesOperator.__init__(self,time_series[0].getControler(),time_series,window_size,0,dbg,"plot(%s)"%text)
         else:
             dbg=time_series.getControler().debug()
             text=str(time_series)
             TimeSeriesOperator.__init__(self,time_series.getControler(),[time_series],window_size,0,dbg,"plot(%s)"%text)
         from pyvisi.renderers.gnuplot import LinePlot,Scene,PsImage
         self.__renderer=Scene()
         self.__line_plot=LinePlot(self.__renderer)
         self.__line_plot.setTitle(text)
         self.__line_plot.setLineStyle("lines")
         self.__line_plot.setXLabel("time")
         self.__line_plot.setYLabel("values")
         self.__file_name=file_name
         if format==None:
             self.__format=PsImage()
         else:
             self.__format=format
         self.__window_size=window_size

    def update(self,start,end):
         s=max(end-self.__window_size,self.getControler().getIdOfFirstAvailableDatum())
         args=[self.getControler()[s:end]]
         for arg in self.getArguments(): args.append(arg[s:end])
         self.__line_plot.setData(*args)
         self.__line_plot.render()
         if self.__file_name==None:
             raise SystemError,"Online viewing is not avilabel yet!"
         else:
             self.__renderer.save(fname=self.__file_name, format=self.__format)
             

def viewer(time_serie,seperator=","):
      """creates a viewer for a time series"""
      import sys
      return Writer(time_serie,sys.stdout,seperator)

def differential(time_serie):
      """calculates the derivative Dv of the time series v:
         
            Dv[n]=(v[n]-v[n-1])/(t[n]-t[n-1])

      """
      out=(((time_serie<<1)-time_serie)/((time_serie.getControler()<<1)-time_serie.getControler())+ \
           ((time_serie>>1)-time_serie)/((time_serie.getControler()>>1)-time_serie.getControler()))/2.
      out.setDescription("d(%s)/dt"%str(time_serie))
      out.setDebug(time_serie.debug())
      return out

def integral(time_serie):
      """calculates the intagral Iv of the time series v using the trapozidal rule:
         
            Iv[n]=int_{t_0}^{t_n} v ~ sum_{0<i<=n} n (v[i]+v[i-1])/2*(t[i]-t[i-1])

      """
      out=TimeSeriesCumulativeSum(((time_serie>>1)+time_serie)/2.*(time_serie.getControler()-(time_serie.getControler()>>1)))
      out.setDescription("I (%s) dt"%str(time_serie))
      out.setDebug(time_serie.debug())
      return out

def smooth(time_serie,range=5):
     """smoothes a time series using the at each time the previous and next range values"""
     i=integral(time_serie)
     out=((i>>range)-(i<<range))/((time_serie.getControler()>>range)-(time_serie.getControler()<<range))
     out.setDescription("smooth(%s,-%d:%d) dt"%(str(time_serie),range,range))
     out.setDebug(time_serie.debug())
     return out

def leakySmooth(time_serie,l=0.99):
     """leaky smoother: s(t)=int_{t_0}^{t} v(r) l^{t-r} dr/ int_{t_0}^{t} l^{t-r} dr """
     w=l**(-time_serie.getControler())
     out=integrate(time_serie*w)/integrate(w)
     out.setDescription("leaky smoother(%s)"%str(time_serie))
     return out

# test

if __name__=="__main__":
   # tests the interfaces to data sets:
   print "Test of Datasets:"
   print "================="
   bf=TimeSeriesBaseBuffer(buffer_size=5,numComponents=1,debug=True,description="TestBaseBuffer")
   bfv_l=TimeSeriesBaseDataset(bf,offset=1,debug=True,description="offset 1")
   bfv_r=TimeSeriesBaseDataset(bf,offset=-1,debug=True,description="offset -1")
   bf.append([1.,2.,3.,4.])
   print "should be all 2. :",bfv_l[0]
   print bf[1]
   print bfv_r[2]
   bf.append([5.,6.,7.])
   print "should be all 5. :",bfv_l[3],bf[4],bfv_r[5]
   print "should be all 6. :",bfv_l[4],bf[5],bfv_r[6]
   print "should be all 7. :",bfv_l[5],bf[6],bfv_r[7]
   print "should be all [6., 7.] :",bfv_l[4:6],bf[5:7],bfv_r[6:8]

   print "Test of Controler"
   print "================="
   b=Controler(buffer_size=15,debug=True)
   s3=b>>3
   s1=b>>1
   s_3=b<<3 
   print s_3
   print b
   print b+s3
   sum=(s_3+b)+(b+s3)
   
   for i in range(30):
       b.nextTime(i*1.)
   b.flush()
   print "should be all 28. :",s_3.getDataset()[25],b.getDataset()[28],s3.getDataset()[31]
   print "should be all 29. :",s_3.getDataset()[26],b.getDataset()[29],s3.getDataset()[32]
   print "should be all 96. :",sum.getDataset()[24]
  
   print "Test of operators"
   print "================="
   b=Controler(buffer_size=15,debug=True)
   b.setFlushRate(2)
   q=DataCatcher(b)
   b1=b<<1
   a=b+b1
   a_s=b1+1.
   s_a=1.+b1
   d=b-b1
   d_s=b1-1.
   s_d=1.-b1
   m=b*b1
   m_s=b1*2.
   s_m=2.*b1
   dv=b/b1
   dv_s=b1/2.
   s_dv=2./b1
   p=b**b1
   p_s=b1**2.
   s_p=2.**b1
   pb=+b
   mb=-b
   sum=TimeSeriesCumulativeSum(b)
   diff=differential(b)
   smt=smooth(b,2)
   int=integral(b*2)
   fl=file("/tmp/test.csv","w")
   w=Writer(q,fl)
   v=viewer(q)
   plo=Plotter([a,a_s],window_size=4,file_name="s.ps")
   for i in range(30):
       b.nextTime(i*1.)
       if i%2==1: q.nextValue(i*28.)
   b.flush()
   print "a[28] should be %e: %e"%(28.+29.,a[28])
   print "a_s[28] should be %e: %e"%(29.+1.,a_s[28])
   print "s_a[28] should be %e: %e"%(29.+1.,s_a[28])
   print "d[28] should be %e: %e"%(28.-29.,d[28])
   print "d_s[28] should %e: %e"%(29.-1.,d_s[28])
   print "s_d[28] should %e: %e"%(1.-29.,s_d[28])
   print "m[28] should be %e: %e"%(28.*29.,m[28])
   print "m_s[28] should be %e: %e"%(29.*2.,m_s[28])
   print "s_m[28] should be %e: %e"%(29.*2.,s_m[28])
   print "dv[28] should be %e: %e"%(28./29.,dv[28])
   print "dv_s[28] should be %e: %e"%(29./2.,dv_s[28])
   print "s_dv[28] should be %e: %e"%(2./29.,s_dv[28])
   print "p[28] should be %e: %e"%(28.**29.,p[28])
   print "p_s[28] should be %e: %e"%(29.**2,p_s[28])
   print "s_p[28] should be %e: %e"%(2.**29.,s_p[28])
   print "pb[28] should be %e: %e"%(28.,pb[28])
   print "mb[28] should be %e: %e"%(-28.,mb[28])
   print "sum[28] should be %e: %e"%(28*29./2,sum[28])
   print "diff[28] should be %e: %e"%(1.,diff[28])
   print "smt[27] should be %e: %e"%(27.,smt[27])
   print "int[28] should be %e: %e"%(28.**2,int[28])
   print "q[27] should be %e: %e"%(27*28.,q[27])
   print "q[28] should be %e: %e"%(28*28.,q[28])
   print "q[29] should be %e: %e"%(29*28.,q[29])
   fl.flush()
   
   rin=Reader(file("/tmp/test.csv","r+"),buffer_size=15,debug=True)
   rin.run()
   inp=rin.getTimeSeries()
   print "inp[27] should be %e: %e"%(27*28.,inp[27])
   print "inp[28] should be %e: %e"%(28*28.,inp[28])
   print "inp[29] should be %e: %e"%(29*28.,inp[29])

