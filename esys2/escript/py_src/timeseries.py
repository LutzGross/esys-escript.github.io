# $Id$

import numarray

class TimeSeriesBase:
   """The TimeSeriesBase class is the base class for all class of the TimeSeries module.
      It takes care of the updating depending TimeSeriesBase objects and the debuging mechnism"""

   def __init__(self):
       self.__debug=False

   def __str__(self):
       return "TimeSeriesBase"

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
      
class TimeSeriesFilter(TimeSeriesBase):
   """TimeSeriesFilter objects are applied to TimeSeries objects to filer out information or to convert it.
      A TimeSeriesFilter objects is called by the TimeSeries object it is depending on to consider the values currently in the buffer for
      updating. Some TimeSeriesFilter may require values outside the buffer. The TimeSeries object maintains the last buffer_overlap values
      in the buffer so they can be used to process (not neccesarily all) value in the buffer."""

   def __init__(self,buffer_overlap=0):
       self.__left_required_extension=buffer_overlap

   def __str__(self):
       return "TimeSeriesFilter"

   def getBufferOverlapNeededForUpdate(self):
       return self.__left_required_extension

   def update(self,times,values):
       pass

_DEFAULT_CACHE_SIZE=9
_DEFAULT_BUFFER_SIZE=5
_FLOATING_TYPE=numarray.Float64

class TimeSeries(TimeSeriesBase):
   def __init__(self,buffer_overlap=0,buffer_size=_DEFAULT_BUFFER_SIZE,cache_size=_DEFAULT_CACHE_SIZE,numComponents=1):
       if buffer_size>cache_size: raise ValueError,"buffer size has to be less or equal cache size"
       TimeSeriesBase.__init__(self)
       self.__updates=list()
       self.__max_buffer_overlap=0
       self.__buffer_overlap=0
       self.__numNodes=0
       self.__numNodesInBuffer=0
       self.__numNodesInCache=0
       self.__firstNodeInBuffer=0
       self.__firstNodeInCache=0
       self.__buffer_size=buffer_size
       self.__node_cache=numarray.zeros((cache_size,),_FLOATING_TYPE)
       self.__attachment_cache=[]
       if numComponents<2:
          self.__value_cache=numarray.zeros((cache_size,),_FLOATING_TYPE)
       else:
          self.__value_cache=numarray.zeros((cache_size,numComponents),_FLOATING_TYPE)
       self.resizeMaxBufferOverlap(buffer_overlap)

   def __del__(self):
       self.flush()

   def __str__(self):
       return "TimeSeries"

   def getNumComponents(self):
       if self.__value_cache.rank==1:
          return 1
       else: 
          self.__value_cache.shape[1]

   def getNumNodes(self):
       """returns the number of time nodes in the time series"""
       return self.__numNodes

   def getCacheSize(self):
       """returns the cache size"""
       return self.__node_cache.shape[0]

   def getBufferSize(self):
       """returns the cache size"""
       return self.__buffer_size

   def getNumNodesInCache(self):
       """returns the number of nodes in cache"""
       return self.__numNodesInCache

   def getNumNodesInBuffer(self):
       """returns the number of nodes in cache"""
       return self.__numNodesInBuffer
    
   def getFirstNodeInCache(self):
       """returns the id number of the first node in the cache"""
       return self.__firstNodeInCache

   def getFirstNodeInBuffer(self):
       """returns the id number of the first node in the buffer"""
       return self.__firstNodeInBuffer

   def getFirstNodeOfBufferInCache(self):
       """returns the first location of the first node in the buffer relative to the cache"""
       return self.getFirstNodeInBuffer()-self.getFirstNodeInCache()

   def getBufferOverlap(self):
       """returns the current size of the left extension"""
       return self.__buffer_overlap

   def getMaxBufferOverlap(self):
       """returns the maximum size of the left extension"""
       return self.__max_buffer_overlap

   def resizeMaxBufferOverlap(self,new_buffer_overlap=0):
       if new_buffer_overlap>self.__max_buffer_overlap:
          if self.getNumNodes()>0: raise ValueError,"left extension can only be resized for empty time series"
          if self.getCacheSize()<self.getBufferSize()+new_buffer_overlap:
               raise ValueError,"Cache size is too small! required cache size is %s"%self.getBufferSize()+new_buffer_overlap
          self.__max_buffer_overlap=new_buffer_overlap
          if self.debug(): print "Debug: %s: left extension is increased to %d"%(self,new_buffer_overlap)

   def getLastNode(self):
       if self.getNumNodesInCache()>0:
          return self.__node_cache[self.getNumNodesInCache()-1]
       else: 
          return -1.e300

   def getLastValue(self):
       if self.getNumNodesInCache()>0:
          return self.__node_cache[self.getNumNodesInCache()-1]
       else: 
          raise ValueError,"No value available"

   def checkInUpdate(self,time_series_filter):
       """checks in a time_series_filter object to be updated when buffer is full"""
       if self.getNumNodes()>0:
          raise TypeError,"Check in of TimeSeries requires empty buffer."
       self.__updates.append(time_series_filter)
       self.resizeMaxBufferOverlap(time_series_filter.getBufferOverlapNeededForUpdate())
       if self.debug(): print "Debug: %s: %s checked in successfully."%(self,time_series_filter)

   def append(self,time_nodes,values,attachments=None):
       """appends the time_nodes and values into the buffer"""
       num_additional_nodes=time_nodes.shape[0]
       if num_additional_nodes<1: return
       if self.debug(): 
            if num_additional_nodes>1:
               print "Debug: %s: values %d to %d are added to time series."%(self,self.getNumNodes(),self.getNumNodes()+num_additional_nodes-1)
            else:
               print "Debug: %s: value %d is added to time series."%(self,self.getNumNodes())
       if not num_additional_nodes==values.shape[0]:
          raise ValueError,"Number time nodes and number of values don't match."
       if self.getLastNode()>=time_nodes[0]:
          raise ValueError,"first time node to be checked in is less than last previously checked in node"

       if num_additional_nodes>1:
            if min(time_nodes[1:num_additional_nodes]-time_nodes[0:num_additional_nodes-1])<=0:
              raise ValueError,"time nodes have to be strictly increasing"
      
       # full cache requires a shift:
       if self.getNumNodesInCache()+num_additional_nodes>self.getCacheSize():
           new_num_nodes_in_cache=self.getNumNodesInBuffer()+self.getBufferOverlap()
           if new_num_nodes_in_cache+num_additional_nodes>self.getCacheSize():
              raise ValueError,"Cache overflow: Expected size is bigger than %d"%(new_num_nodes_in_cache+num_additional_nodes)
           start=self.getNumNodesInCache()-new_num_nodes_in_cache
           end=start+new_num_nodes_in_cache
           self.__node_cache[0:new_num_nodes_in_cache]=self.__node_cache[start:end]
           self.__value_cache[0:new_num_nodes_in_cache]=self.__value_cache[start:end]
           self.__attachment_cache[0:new_num_nodes_in_cache]=self.__attachment_cache[start:end]
 
           self.__firstNodeInCache+=start
           self.__numNodesInCache=new_num_nodes_in_cache
           if self.debug(): print "Debug: %s: %d values from %d onwards are moved to the beginning of the cache (first node in cache is now %d)."% \
                                                                                    (self,new_num_nodes_in_cache,start,self.__firstNodeInCache)
           
       # copy values into cache:
       if self.getNumNodesInCache()+num_additional_nodes>self.getCacheSize():
           raise ValueError,"Cache overflow: Expected size is bigger than %d"%(self.getNumNodesInCache()+num_additional_nodes)
       if self.debug(): 
           if num_additional_nodes>1:
              print "Debug: %s: values %d to %d of cache are updated"%(self,self.getNumNodesInCache(),self.getNumNodesInCache()+num_additional_nodes-1)
           else:
              print "Debug: %s: value %d of cache is updated."%(self,self.getNumNodesInCache())
       self.__node_cache[self.getNumNodesInCache():self.getNumNodesInCache()+num_additional_nodes]=time_nodes
       self.__value_cache[self.getNumNodesInCache():self.getNumNodesInCache()+num_additional_nodes]=values
       self.__numNodes+=num_additional_nodes
       self.__numNodesInBuffer+=num_additional_nodes
       self.__numNodesInCache+=num_additional_nodes
       print self.__node_cache
       print self.__value_cache
       # copy values into cache:
       if self.getNumNodesInBuffer()>=self.getBufferSize():
              if self.debug() and len(self.__updates)>0: print "Debug: %s: buffer is full. Updating process is started"%self
              self.processBuffer()

   def flush(self):
      self.processBuffer()

   def processBuffer(self):
        if self.getNumNodesInBuffer()>0:
           for i in self.__updates: 
             if self.debug(): print "Debug: %s: update for %s started"%(self,i)
             if i.getBufferOverlapNeededForUpdate()>self.getBufferOverlap():
                s=self.getFirstNodeOfBufferInCache()
                l=self.getNumNodesInBuffer()
             else:
                s=self.getFirstNodeOfBufferInCache()-i.getBufferOverlapNeededForUpdate()
                l=self.getNumNodesInBuffer()+i.getBufferOverlapNeededForUpdate()
             i.update(self.__node_cache[s:s+l],self.__value_cache[s:s+l])
           self.__firstNodeInBuffer+=self.__numNodesInBuffer
           self.__numNodesInBuffer=0
        self.__buffer_overlap=self.getMaxBufferOverlap()
        if self.debug(): print "Debug: %s: first node in buffer is now %d"%(self,self.__firstNodeInBuffer)

   


class TimeSeriesCollector(TimeSeries):
      """TimeSeriesCollector collects values at time nodes"""
      def __init__(self):
         TimeSeries.__init__(self)

      def __str__(self):
         return "TimeSeriesCollector"

      def add(self,time_mark,value):
           """adds the value at time time_mark to the time series"""
           self.append(numarray.array([time_mark]),numarray.array([value]))

      def read(self,istream,seperator=","):
        """reads pairs from iostream istream"""
        for l in istream:
           d=l.strip().split(seperator)
           self.add(float(d[0]),float(d[1]))

class TimeSeriesIntegrator(TimeSeries,TimeSeriesFilter):
      def __init__(self,time_series):
         TimeSeriesFilter.__init__(self,1)
         TimeSeries.__init__(self,buffer_size=time_series.getBufferSize(),cache_size=time_series.getCacheSize(), \
                                                                         numComponents=time_series.getNumComponents())
         self.setDebug(time_series.debug())
         time_series.checkInUpdate(self)
         self.__integral=0

      def __str__(self):
         return "TimeSeriesIntegrator"

      def update(self,times,values):
          l=times.shape[0]
          self.append(times[1:l],(values[0:l-1]+values[1:l])/2.*(times[1:l]-times[0:l-1]))
         

class TimeSeriesDifferential(TimeSeries,TimeSeriesFilter):
      def __init__(self,time_series):
         TimeSeriesFilter.__init__(self,1)
         TimeSeries.__init__(self,buffer_size=time_series.getBufferSize(),cache_size=time_series.getCacheSize(), \
                                                                         numComponents=time_series.getNumComponents())
         self.setDebug(time_series.debug())
         time_series.checkInUpdate(self)

      def __str__(self):
         return "TimeSeriesDifferential"

      def update(self,times,values):
          l=times.shape[0]
          self.append((times[0:l-1]+times[1:l])/2,(values[0:l-1]-values[1:l])/(times[0:l-1]-times[1:l]))

class TimeSeriesViewer(TimeSeriesFilter):
      def __init__(self,time_series):
         TimeSeriesFilter.__init__(self,0)
         time_series.checkInUpdate(self)

      def __str__(self):
         return "TimeSeriesViewer"

      def update(self,times,values):
          for i in range(times.shape[0]): print "[%s: %s]"%(times[i],values[i])

class TimeSeriesWriter(TimeSeriesFilter):
      def __init__(self,time_series,ostream,seperator=","):
         TimeSeriesFilter.__init__(self,0)
         time_series.checkInUpdate(self)
         self.setDebug(time_series.debug())
         self.__ostream=ostream
         self.__seperator=seperator

      def __str__(self):
         return "TimeSeriesWriter"

      def update(self,times,values):
        for i in range(times.shape[0]): self.__ostream.writelines("%s,%s\n"%(times[i],values[i]))

# test

if __name__=="__main__":

   c=TimeSeriesCollector()
   c.setDebugOn()
   ii=TimeSeriesIntegrator(c)
   d=TimeSeriesDifferential(c)
   v=TimeSeriesViewer(ii)
   w=TimeSeriesWriter(d,file("test.csv","w"))

   for i in range(15):
      c.add(i*1.,i+1.)
