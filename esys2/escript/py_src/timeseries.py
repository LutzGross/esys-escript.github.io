# $Id$

import numarray
from types import SliceType
DEFAULT_BUFFER_SIZE=9
DEFAULT_FLOAT_TYPE=numarray.Float64

class TimeSeriesBase:
   """The TimeSeriesBase class is the base class for all class of the TimeSeries module."""

   def __init__(self,debug=False,description="timeseries.Base"):
       self.__debug=debug
       self.__description=description

   def __str__(self):
       return self.__description

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
class TimeSeriesDataset(TimeSeriesBase):
   """provides an interface for accessing a set of linearly ordered data."""
   def __init__(self,buffer,offset=0,debug=False,description="timeseries.Dataset"):
       TimeSeriesBase.__init__(self,debug,description)
       self.__buffer=buffer
       self.__offset=offset
       if self.debug(): print "Debug: %s: offset %d to buffer"%(self,self.getOffset())

   def __len__(self):
       """needed to handle negative indexing in slicing"""
       return 0

   def getNumComponents(self):
       """returns the number of components of the data (may be overwritten by subclass)"""
       return self.getBuffer().getNumComponents()

   def getIdOfLastDatum(self):
      """returns the identification number of the last datum in the data set (may be overwritten by subclass)"""
      return self.getBuffer().getIdOfLastDatum()-self.getOffset()

   def getIdOfFirstDatum(self):
      """returns the identification number of the first datum (may be overwritten by subclass)"""
      return self.getBuffer().getIdOfFirstDatum()-self.getOffset()

   def getOffsetInBuffer(self):
      """returns the offset to access elements in getBuffer() (may be overwritten by subclass)"""
      return  self.getOffset()

   def getIdOfLastUnusedDatum(self):
       """returns the identification number of the last datum which has been unused by all TimeSeries refering to the TimeSeriesDataset (may be overwritten by subclass)"""
       return self.getBuffer().getIdOfLastUnusedDatum()-self.getOffset()

   def updateIdOfLastUnusedDatum(self,last_unused_datum):
       """updates the identification number of the last unused datum (to be overwritten by subclass)"""
       self.getBuffer().updateIdOfLastUnusedDatum(last_unused_datum+self.getOffset())

   def append(self,values):
       """appends data to the buffer. If the buffer would be full the buffer is rearranged before the data are appended  (to be overwritten by subclass)"""
       self.getBuffer().append(values)

   def getBufferSize(self):
       """returns the size of the buffer (to be overwritten by subclass)"""
       return self.getBuffer().getBufferSize()
   
   def needsRearrangement(self,num_new_data=0):
       """returns True if the buffer will be full after num_new_data have been appended (to be overwritten by subclass)"""
       return self.getBuffer().needsRearrangement(num_new_data)

   def isEmpty(self):
      """returns true if no data are appeneded to buffer"""
      return self.getNumData()<=0
   
   def getNumData(self):
      """returns the number of data (not all of them are accessible)"""
      return self.getIdOfLastDatum()-self.getIdOfFirstDatum()+1

   def getBuffer(self):
      """return the buffer referenced by the TimeSeriesDataset"""
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
             return self.getBuffer()[start+self.getOffsetInBuffer():end+self.getOffsetInBuffer()]
      else:
         if index<self.getIdOfFirstDatum() or index>self.getIdOfLastDatum(): raise IndexError,"%s: Index %d out of range"%(self,index)
         return self.getBuffer()[index+self.getOffsetInBuffer()]

class TimeSeriesBuffer(TimeSeriesDataset):
   """An inplementation of TimeSeriesDataset which actually is storing data into a numarray buffer"""
   def __init__(self,buffer_size=DEFAULT_BUFFER_SIZE,numComponents=1,type=DEFAULT_FLOAT_TYPE,id_of_first_datum=0,debug=False,description="timeseries.Buffer"):
       if numComponents<2:
          buffer=numarray.zeros((buffer_size,),type)
       else:
          buffer=numarray.zeros((buffer_size,numComponents),type)
       TimeSeriesDataset.__init__(self,buffer,id_of_first_datum-1,debug,description)
       self.__num_data_in_buffer=0
       self.__id_last_unused_datum=id_of_first_datum-1
       self.__id_last_datum=id_of_first_datum-1
       self.__id_first_datum=id_of_first_datum
       if self.debug(): print "Debug: %s : buffer of size %d with %d components allocated (first datum is %d)."% \
                       (self,self.getBufferSize(),self.getNumComponents(),id_of_first_datum)


   def getBufferSize(self):
       """returns the size of the buffer"""
       return self.getBuffer().shape[0]
   
   def getNumComponents(self):
       """returns the number of components of the data (overwrites TimeSeriesDataset method)"""
       if self.getBuffer().rank==1:
          return 1
       else: 
          self.getBuffer().shape[1]

   def getNumDataInBuffer(self):
       """returns the number of data currently in the buffer"""
       return self.__num_data_in_buffer

   def getIdOfLastDatum(self):
      """returns the identification number of the last datum in the data set (overwrites method from TimeSeriesDataset)"""
      return self.__id_last_datum

   def getIdOfFirstDatum(self):
      """returns the identification number of the first datum (overwrites method from TimeSeriesDataset)"""
      return self.__id_first_datum

   def getOffsetInBuffer(self):
      """returns the offset to access elements in the buffer (overwrites method from TimeSeriesDataset)"""  
      return -self.getIdOfLastDatum()+self.getNumDataInBuffer()-1  

   def getIdOfLastUnusedDatum(self):
       """returns the identification number of the last datum which has been unused by all TimeSeries refering to the TimeSeriesDataset (overwrites method from TimeSeriesDataset)"""
       return self.__id_last_unused_datum

   def updateIdOfLastUnusedDatum(self,last_unused_datum):
       """updates the identification number of the last unused datum (to be overwritten by subclass)"""
       self.getBuffer().updateIdOfLastUnusedDatum(last_unused_datum-self.getOffset())

   def updateIdOfLastUnusedDatum(self,last_unused_datum):
       """updates the identification number of the last unused datum (overwrites TimeSeriesDataset method)"""
       if self.__id_last_unused_datum>last_unused_datum:
           self.__id_last_unused_datum=last_unused_datum
           if self.debug(): print "Debug: %s: last unused datum is now %s"%(self,last_unused_datum)

   def needsRearrangement(self,num_new_data=0):
       """returns True if the buffer will be full after num_new_data have been appended"""
       return self.getNumDataInBuffer()+num_new_data>self.getBufferSize()
         
   def append(self,data):
      """appends data to the buffer. If the buffer would be full the buffer is rearranged before the data are appended (overwrites TimeSeriesDataset method)"""
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
        nn=self.getNumDataInBuffer()
        num_protected_data=self.getIdOfLastDatum()-self.getIdOfLastUnusedDatum()
        if num_protected_data+num_new_data>self.getBufferSize():
              raise ValueError,"%s: buffer overflow: buffer size has to be bigger than %d"%(self,num_protected_data+num_new_data)
        if num_protected_data>0: self.getBuffer()[0:num_protected_data]=self.getBuffer()[nn-num_protected_data:nn]
        self.__num_data_in_buffer=num_protected_data
        self.__id_last_unused_datum=self.__id_last_datum
        if self.debug(): 
             print "Debug: %s: rearrangement: first data in buffer is %d."%(self,self.getIdOfLastDatum()-self.getNumDataInBuffer()+1)
      # copy data over:
      nn=self.getNumDataInBuffer()
      self.getBuffer()[nn:nn+num_new_data]=data
      self.__num_data_in_buffer+=num_new_data
      self.__id_last_datum+=num_new_data
      self.__id_last_unused_datum+=num_new_data
      if self.debug(): print "Debug: %s: %d data appended. Last unused datum is now %d."%(self,num_new_data,self.__id_last_unused_datum)

# ======================================
class TimeSeries(TimeSeriesDataset):
      """a TimeSeries glues a Controler controler and a TimeSeriesDataset dataset together. It also provides a TimeSeriesDataset view to the datset"""
      def __init__(self,dataset,debug=False,description="timeseries."):
           TimeSeriesDataset.__init__(self,dataset,0,debug,description)
           self.__id_last_processed_datum=dataset.getIdOfFirstDatum()-1         
      
      def getDataset(self):
          """returns the TimeSeriesDataset of the time series"""
          return self.getBuffer()

      def getControler(self):
          """returns the Controler of the time series (to be overwritten by subclass)"""
          pass

      def getIdOfLastProcessedDatum(self):
          return self.__id_last_processed_datum

      def updateIdOfLastProcessedDatum(self,id_last_processed_datum):
          self.__id_last_processed_datum=id_last_processed_datum

      def __add__(self,arg):
         if isinstance(arg,TimeSeriesDataset):
            return TimeSeriesSum(self,arg)
         else:
            return TimeSeriesAddScalar(self,arg)

      def __sub__(self,arg):
         return self+(-1.)*arg

      def __mul__(self,arg):
         if isinstance(arg,TimeSeriesDataset):
            return TimeSeriesMult(self,arg)
         else:
            return TimeSeriesMultScalar(self,arg)

      def __div__(self,arg):
         if isinstance(arg,TimeSeriesDataset):
            return TimeSeriesDiv(self,arg)
         else:
            return TimeSeriesMultScalar(self,1./arg)

      def __pow__(self,arg):
         if isinstance(arg,TimeSeriesDataset):
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
         if isinstance(arg,TimeSeriesDataset):
            return TimeSeriesDiv(arg,self)
         else:
            return TimeSeriesDivScalar(self,arg)

      def __rpow__(self,arg):
         if isinstance(arg,TimeSeriesDataset):
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

class TimeSeriesFilter(TimeSeries):
      """a TimeSeriesFilter is a TimeSeriesDataset attached to a Controler where the TimeSeriesDataset provides data 
          at the time nodes defined by the Controler. Additional to a TimeSeries a TimeSeriesFilter allows to update
          the underlying TimeSeriesDataset through the update method which is overwritten by a particular implementation of the 
          class. The update method is called to append the data [start:end] to the attached dataset by the the attached TimeSerieControler"""
      def __init__(self,controler,dataset,args=[],left_wing_size=0,right_wing_size=0,debug=False,description="timeseries.Filter"):
          TimeSeries.__init__(self,dataset,debug,description)
          self.__left_wing_size=left_wing_size
          self.__right_wing_size=right_wing_size
          self.__args=args
          self.__controler=controler
          controler.appendFilterToUpdateList(self)

      def getControler(self):
          """returns the Controler of the time series (overwrites method of by TimeSeries)"""
          return self.__controler

      def update(self,start,end):
          """appends zeros to the dataset. This method should be overwritten by a particular TimeSeriesFilter"""
          nc=self.getNumComponents()
          if nc>1:
             self.getDataset().append(numarray.zeros([nc,end-start]))
          else:
             self.getDataset().append(numarray.zeros(end-start))
      def getLeftWingSize(self):
          """returns the left wing size"""  
          return self.__left_wing_size

      def getRightWingSize(self):
          """returns the right wing size""" 
          return self.__right_wing_size

      def getArguments(self,index=None):
          """returns the list of arguments or, index is present, the argument with index index. In the latter case None is returned if no arguments are present"""
          if index==None:
             return self.__args
          else:
             if len(self.__args)>0:
                return self.__args[index]
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
          end=None
          for i in self.getArguments(): 
             if end==None:
                end=i.getIdOfLastDatum()
             else:
                end=min(end,i.getIdOfLastDatum())
          if not end==None:
              if self.debug(): print "Debug: %s: range [%d:%d] is updated."%(self,start,end-self.getRightWingSize())
              self.update(start,end-self.getRightWingSize()+1)      
              for i in self.getArguments(): i.updateIdOfLastUnusedDatum(end-self.getLeftWingSize())
              self.updateIdOfLastProcessedDatum(end)

class Controler(TimeSeries):
   """controls a set of TimeSeries"""
   def __init__(self,buffer_size=DEFAULT_BUFFER_SIZE,debug=False,description="timeseries.Controler"): 
        TimeSeries.__init__(self,TimeSeriesBuffer(buffer_size,1,DEFAULT_FLOAT_TYPE,0,debug,"Time nodes buffer of "+description),\
                                                                                   debug,"Time nodes of "+description)
        self.setFlushRate()   
        self.__update_time_series=list()
      
   def __del__(self):
       self.flush()

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

   def appendFilterToUpdateList(self,time_serie):
       if not time_serie.getControler()==self: raise ValueError,"%s: time series time_serie %s is not linked defined on %s."%(self,time_serie,self)
       if not self.isEmpty(): raise ValueError,"%s: you can only check in a time series time_serie is controler is empty."%self
       self.__update_time_series.append(time_serie)
       if self.debug(): print "Debug: %s: %s has been added to update list."%(self,time_serie)

   def newTimeNode(self,value):
       if self.needsFlushing(): self.flush()
       self.getDataset().append(value)
       if self.debug(): print "Debug: %s: new time node %e has been added."%(self,value)

# ============================================
class TimeSeriesShift(TimeSeries):
      """creates a shift of the time series, i.e. if d[n] is the datum at time t[n], the value at t[n] becomes v[n+shift] on the output"""
      def __init__(self,time_serie,shift=1):
          if shift<0:
              dsc="(%s)<<%d"%(time_serie,-shift)
          else:
              dsc="(%s)>>%d"%(time_serie,shift)
          self.__controler=time_serie.getControler()
          TimeSeries.__init__(self,TimeSeriesDataset(time_serie.getDataset(),-shift,time_serie.debug(),"buffer view to "+dsc),\
                                          time_serie.debug(),dsc)
      def getControler(self):
          return self.__controler

class TimeSeriesSum(TimeSeriesFilter):
      """adds two TimeSeries"""
      def __init__(self,time_serie_1,time_serie_2):
          dsc="(%s)+(%s)"%(time_serie_1,time_serie_2)
          dbg=time_serie_1.debug() or time_serie_2.debug()
          cntrl=time_serie_1.getControler()
          if not cntrl==time_serie_2.getControler(): 
                  raise ValueError("TimeSeriesSum: %s and %s have different controler."%(time_serie_1,time_serie_2))
          id_first_datum=max(time_serie_1.getIdOfFirstDatum(),time_serie_2.getIdOfFirstDatum())
          TimeSeriesFilter.__init__(self,cntrl, \
                              TimeSeriesBuffer(cntrl.getBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                       TimeSeriesBuffer(cntrl.getBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                   TimeSeriesBuffer(cntrl.getBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                       TimeSeriesBuffer(cntrl.getBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                     TimeSeriesBuffer(cntrl.getBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                       TimeSeriesBuffer(cntrl.getBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                TimeSeriesBuffer(cntrl.getBufferSize(),time_serie_1.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                       TimeSeriesBuffer(cntrl.getBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
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
                     TimeSeriesBuffer(cntrl.getBufferSize(),time_serie.getNumComponents(),DEFAULT_FLOAT_TYPE,id_first_datum,dbg,"buffer for "+dsc), \
                     [time_serie],0,0,dbg,dsc)

      def update(self,start,end):
          self.append(numarray.exp(self.getArgumentDataset(0)[start:end]))

class TimeSeriesCumulativeSum(TimeSeriesFilter):
      """creates a shift of the time series, i.e. if d[n] is the datum at time t[n], the value at t[n] becomes v[n+shift] on the output"""
      def __init__(self,time_series):
         TimeSeriesFilter.__init__(self,1)
         TimeSeries.__init__(self,frame_size=time_series.getDatasetSize(),buffer_size=time_series.getBufferSize(), \
                                                                         numComponents=time_series.getNumComponents())
         self.setDebug(time_series.debug())
         time_series.checkInUpdate(self)
         self.__integral=0

      def __str__(self):
         return "timeseries.Integrator"

      def update(self,times,data):
          l=times.shape[0]
          self.append(times[1:l],(data[0:l-1]+data[1:l])/2.*(times[1:l]-times[0:l-1]))
         

class TimeSeriesCollector(TimeSeries):
      """timeseries.Collector collects data at time nodes"""
      def __init__(self):
         TimeSeries.__init__(self)

      def __str__(self):
         return "timeseries.Collector"

      def add(self,time_mark,value):
           """adds the value at time time_mark to the time series"""
           self.append(numarray.array([time_mark]),numarray.array([value]))

      def read(self,istream,seperator=","):
        """reads pairs from iostream istream"""
        for l in istream:
           d=l.strip().split(seperator)
           self.add(float(d[0]),float(d[1]))

def Differential(time_series):
      """calculates the derivative Dv of the time series v:
         
            Dv[n]=(v[n]-v[n-1])/(t[n]-t[n-1])

      """
      return (time_series<<1-time_series)/(time_series.getControler()<<1-time_series.getControler())

def Integral(time_series):
      """calculates the intagral Iv of the time series v using the trapozidal rule:
         
            Iv[n]=sum_i<n (v[n]+v[n-1])/2*(t[n]-t[n-1])

      """
      return TimeSeriesCumulativeSum((time_series<<1+time_series)/2.*(time_series.getControler()-(time_series.getControler<<1)),0.)


class TimeSeriesViewer(TimeSeriesFilter):
      def __init__(self,time_series):
         TimeSeriesFilter.__init__(self,0)
         time_series.checkInUpdate(self)

      def __str__(self):
         return "timeseries.Viewer"

      def update(self,times,data):
          for i in range(times.shape[0]): print "[%s: %s]"%(times[i],data[i])

class TimeSeriesWriter(TimeSeriesFilter):
      def __init__(self,time_series,ostream,seperator=","):
         TimeSeriesFilter.__init__(self,0)
         time_series.checkInUpdate(self)
         self.setDebug(time_series.debug())
         self.__ostream=ostream
         self.__seperator=seperator

      def __str__(self):
         return "timeseries.Writer"

      def update(self,times,data):
        for i in range(times.shape[0]): self.__ostream.writelines("%s,%s\n"%(times[i],data[i]))

# test

if __name__=="__main__":
   # tests the interfaces to data sets:
   print "Test of Datasets:"
   print "================="
   bf=TimeSeriesBuffer(buffer_size=5,numComponents=1,debug=True,description="TestBuffer")
   bfv_l=TimeSeriesDataset(bf,offset=1,debug=True,description="offset 1")
   bfv_r=TimeSeriesDataset(bf,offset=-1,debug=True,description="offset -1")
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
   sum=(s_3+b)+(b+s3)
   
   for i in range(30):
       b.newTimeNode(i*1.)
   b.flush()
   print "should be all 28. :",s_3.getDataset()[25],b.getDataset()[28],s3.getDataset()[31]
   print "should be all 29. :",s_3.getDataset()[26],b.getDataset()[29],s3.getDataset()[32]
   print "should be all 96. :",sum.getDataset()[24]
  
   print "Test of operators"
   print "================="
   b=Controler(buffer_size=15,debug=True)
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
   for i in range(30):
       b.newTimeNode(i*1.)
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

   1/0
   c=TimeSeriesCollector(b)
   c.setDebugOn()
   ii=TimeSeriesIntegrator(c)
   d=TimeSeriesDifferential(c)
   v=TimeSeriesViewer(ii)
   w=TimeSeriesWriter(d,file("test.csv","w"))

   for i in range(15):
      b.newTime(i*1.)
      c.add(i+1.)
   
