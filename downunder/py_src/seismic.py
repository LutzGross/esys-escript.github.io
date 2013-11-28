
##############################################################################
#
# Copyright (c) 2003-2013 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development since 2012 by School of Earth Sciences
#
##############################################################################

__copyright__="""Copyright (c) 2003-2013 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

__all__ = ['SimpleSEGYWriter', 'Ricker', 'WaveBase', 'SonicWave' ]


from math import pi
import numpy as np
import sys
import time
from esys.escript import *
import esys.escript.unitsSI as U
from esys.escript.linearPDEs import LinearPDE



             
class Wavelet(object):
	"""
	place holder for source wavelet
	"""
	pass
	
class Ricker(Wavelet):
	"""
	The Ricker Wavelet w=f(t)
	"""
	def __init__(self, f_dom=40, t_dom=None):
	        """
	        set up Ricker wavelet wih dominant frequence f_dom and center at time t_dom. If t_dom
	        is not given an estimate for suitable t_dom is calculated so f(0)~0.
	        """
		drop=18
		self.__f=f_dom
		self.__f_max=sqrt(7)*f_dom
		self.__s=pi*self.__f
		if t_dom == None:
			t_dom=sqrt(drop)/self.__s
		self.__t0=t_dom
		
	def getCenter(self):
	       """
	       return value of wavelet center 
	       """
               return self.__t0	
        
        def getTimeScale(self):
                """
                returns the time scale which is the inverse of the largest freqence with a significant 
                spectral component. 
                """
		return 1/self.__f_max

	def getValue(self, t):
	        """
	        get value of wavelet at time t
	        """
		e2=(self.__s*(t-self.__t0))**2
		return (1-2*e2)*exp(-e2)

	def getAcceleration(self, t):
	        """
	        get the acceleration f''(t) at time t
	        """
		e2=(self.__s*(t-self.__t0))**2
                return 2*self.__s**2*(-4*e2**2 + 12*e2 - 3)*exp(-e2)
                

class SimpleSEGYWriter(object):
	"""
	as simple writer for 2D and 3D seimic lines in particular for synthetic data

	Typical usage:

    	   from esys.escript import unitsSI as U
	   sw=SimpleSEGYWriter([0.,100*U.m,200*U,m,300.], source=200*U.m, sampling_interval=4*U.msec)
	   while n < 10:
	        sw.addRecord([i*2., i*0.67, i**2, -i*7])
	   sw.write('example.segy')

	:note: the writer uses `obsy`
	"""
	COORDINATE_SCALE = 1000.
	def __init__(self, receiver_group=None, source=0., sampling_interval=4*U.msec, text="some seimic data"):
		"""
		initalize writer

		:param receiver_group: list of receiver coordinates (in meters). For the 2D case a list of floats is given, for the 3D case a list of coordinate tuples
		                       are given
		:param source: coordinates of the source (in meters).  For the 2D case a single floats is given, for the 3D case a coordinate tuples
		:param sampling_interval: sample rate in seconds
		:param text: a text for the header file (e.g a description)
		"""
		if isinstance(source, float) or  isinstance(source, int) :
			DIM=1
			source = (source, 0.)
		elif hasattr(source, '__len__'):
			DIM=len(source)
			if DIM == 1:
			    source = (source[0], 0.)
		        elif DIM==2:
			    source = (source[0], source[1] )
			else:
		            raise ValueError("Only 1D or 2D arrays accepted.")
	        else:
		    raise TypeError("illegal source type.")
	        if receiver_group== None: 
			if DIM == 1:
				receiver_group=[0.]
			else:
				receiver_group=[(0.,0.)]

		if isinstance(receiver_group, float) or  isinstance(receiver_group, int) :
			if DIM == 1:
			    rg = [ (receiver_group, 0. ) ]
			else:
		            raise TypeError("illegal receiver_group type.")
		elif hasattr(receiver_group, '__len__'):
			if DIM == 1:
			    rg = [(c,0)  for c in receiver_group]
			else:
			    rg = [(c[0],c[1])  for c in receiver_group]
	        else:
		    raise TypeError("illegal receiver_group type.")
                            
		self.__source=source
		self.__receiver_group=rg
		self.__text=text
                self.__trace = [ [] for i in self.__receiver_group ]
                self.__sampling_interval = sampling_interval   

	def addRecord(self, record):
	       """
	       adds a record to the traces. a time difference of sample_interval between two records is assumed. 
	       The record mast be a list of as many values as given receivers or a float if a single receiver is used.

               :param record: list of tracks to be added to the record. 
	       """
	       if not len(self.__receiver_group) == len(record): 
	                 raise ValueError("expected number of records is %s but %s given."%(len(self.__receiver_group), len(record)))
	       if len(self.__receiver_group) == 1:
		       if isinstance(self.__receiver_group, float) or isinstance(self.__receiver_group, int):
			         self.__trace[0].append(record)
		       else:
			         self.__trace[0].append(record[0])
	       else:
		       for i in xrange(len(record)):
			       self.__trace[i].append(record[i])

	def getSamplingInterval(self):
	       """
	       returns the sampling interval in seconds. 
	       """
               return self.__sampling_interval

	def write(self, filename):
            """
            writes to segy file 
            
            :param filename: file name
            :note: the function uses the `obspy` module. 
            """
            from obspy import Trace, Stream, UTCDateTime
	    from obspy.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
            from obspy.core import AttribDict

	    stream=Stream()
	    
	    for i in xrange(len(self.__receiver_group)):
		    trace = Trace(data=np.array(self.__trace[i], dtype='float32'))
                    # Attributes in trace.stats will overwrite everything in
                    # trace.stats.segy.trace_header (in Hz)
                    trace.stats.sampling_rate = 1./self.getSamplingInterval()
                    #trace.stats.starttime = UTCDateTime(2011,11,11,0,0,0)
                    if not hasattr(trace.stats, 'segy.trace_header'):
                        trace.stats.segy = {}
                    trace.stats.segy.trace_header = SEGYTraceHeader()
                    trace.stats.segy.trace_header.trace_identification_code=1
                    trace.stats.segy.trace_header.trace_sequence_number_within_line = i + 1
                    trace.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = int(self.COORDINATE_SCALE)
                    trace.stats.segy.trace_header.coordinate_units=1
                    trace.stats.segy.trace_header.source_coordinate_x=int(self.__source[0] * self.COORDINATE_SCALE)
                    trace.stats.segy.trace_header.source_coordinate_y=int(self.__source[1] * self.COORDINATE_SCALE)
                    trace.stats.segy.trace_header.group_coordinate_x=int(self.__receiver_group[i][0] * self.COORDINATE_SCALE)
                    trace.stats.segy.trace_header.group_coordinate_y=int(self.__receiver_group[i][1] * self.COORDINATE_SCALE)

                    # Add trace to stream
                    stream.append(trace)
            # A SEGY file has file wide headers. This can be attached to the stream
            # object.  If these are not set, they will be autocreated with default
            # values.
            stream.stats = AttribDict()
            stream.stats.textual_file_header = 'C.. '+self.__text+'\nC.. with esys.escript.downunder r%s at %s'%(getVersion(),time.asctime())
            stream.stats.binary_file_header = SEGYBinaryFileHeader()

            if getMPIRankWorld()<1:
                stream.write(filename, format="SEGY", data_encoding=1,byteorder=sys.byteorder)

class WaveBase(object):
      """
      Base for wave propagation using the Verlet scheme. 
            
            u_tt = A(t,u), u(t=t0)=u0, u_t(t=t0)=v0

      with a given  acceleration force as function of time. 
      
      a_n=A(t_{n-1})
      v_n=v_{n-1} +  dt *  a_n
      u_n=u_{n-1} +  dt *  v_n
      
      
      """
      def __init__(self, dt, u0, v0, t0=0.):
         """
         set up the wave base
         
         :param dt: time step size (need to be sufficiently small)
         :param u0: initial value
         :param v0: initial velocity
         :param t0: initial time 
         """
         self.u=u0
         self.v=v0
         self.t=t0
         self.t_last=t0
         self.__dt=dt
         
      def getTimeStepSize(self):
    	   return self.__dt
    	   
      def _getAcceleration(self, t, u):
           """
           returns the acceleraton for time t and solution u at time t
           : note: this function needs to be overwritten by a particular wave problem
           """
           pass
           
      def update(self, t):
             """
             returns the solution for the next time marker t which needs to greater than the time marker from the 
             previous call.
             """
             if not self.t_last <= t:
                  raise ValueError("You can not move backward in time.")
                  
             dt = self.getTimeStepSize()
             # apply Verlet scheme until 
             while self.t < t:
                  a=self._getAcceleration(self.t, self.u)
		  self.v += dt * a
		  self.u += dt * self.v 
                  self.t += dt
             
             # now we work backwards
             self.t_last=t
             return t, self.u + self.v * (t-self.t)
             
def createAbsorbtionLayerFunction(x, absorption_zone=300*U.m, absorption_cut=1.e-2):
    """
    creating a distribution which is one in the interior of the domain of x
    and is falling down to the value 'absorption_cut' over a margain of thickness 'absorption_zone' 
    toward each boundary except the top of the domain.
    
    :param x: location of points in the domain
    :type x: `Data`
    :param absorption_zone: thickness of the aborption zone 
    :param absorption_cut: value of decay function on domain boundary
    :return: function on 'x' which is one in the iterior and decays to almost zero over a margin 
             toward the boundary.            
    """    
    dom=x.getDomain()  
    bb=boundingBox(dom)
    DIM=dom.getDim()
    decay=-log(absorption_cut)/absorption_zone**2
    f=1
    for i in xrange(DIM):
        x_i=x[i]
        x_l=x_i-(bb[i][0]+absorption_zone)
	m_l=whereNegative(x_l)
	f=f*( (exp(-decay*(x_l*m_l)**2)-1) * m_l+1 )
	if not DIM-1 == i:
	    x_r=(bb[i][1]-absorption_zone)-x_i
	    m_r=whereNegative(x_r)
	    f=f*( (exp(-decay*(x_r*m_r)**2)-1) * m_r+1 )
    return f
     
class SonicWave(WaveBase):
	"""
	Solving the sonic wave equation 
	
	p_tt = (v_p**2 * p_i)_i  + f(t) * delta_s   where (p-) velocity v_p.
	
	f(t) is wavelet acting at a point source term at positon s 
	"""
	def __init__(self, domain, v_p, wavelet, source_tag, dt=None, p0=None, p0_t=None, absorption_zone=300*U.m, absorption_cut=1e-2):
           """
           initialize the sonic wave solver
           
           :param domain: domain of the problem
           :type domain: `Doamin`        
           :param v_p: p-velocity field    
           :type v_p: `Scalar`       
           :param wavelet: wavelet to describe the time evolution of source term 
           :type wavelet: `Wavelet`          
           :param source_tag: tag of the source location
           :type source_tag: 'str' or 'int'           
           :param dt: time step size. If not present a suitable time step size is calculated.           
           :param p0: initial solution. If not present zero is used.           
           :param p0_t: initial solution change rate. If not present zero is used.           
           :param absorption_zone: thickness of absorption zone           
           :param absorption_cut: boundary value of absorption decay factor
           """
           f=createAbsorbtionLayerFunction(Function(domain).getX(), absorption_zone, absorption_cut)
           v_p=v_p*f

	   if p0 == None:
	      p0=Scalar(0.,Solution(domain))
	   else:
	      p0=interpolate(p0, Solution(domain ))
	      
	   if p0_t == None:
	      p0_t=Scalar(0.,Solution(domain))
	   else:
	      p0_t=interpolate(p0_t, Solution(domain ))
	   
	   if dt == None:
                  dt=min(inf((1./5.)*domain.getSize()/v_p), wavelet.getTimeScale())
            
           super(SonicWave, self).__init__( dt, u0=p0, v0=p0_t, t0=0.)
           
           self.__wavelet=wavelet
           self.__mypde=LinearPDE(domain)
           self.__mypde.getSolverOptions().setSolverMethod(self.__mypde.getSolverOptions().HRZ_LUMPING)
           self.__mypde.setValue(D=1.)
	   self.__source_tag=source_tag
           self.__r=Scalar(0., DiracDeltaFunctions(self.__mypde.getDomain()))
	   self.__vp2=v_p**2


        def  _getAcceleration(self, t, u):
             """
             returns the acceleraton for time t and solution u at time t
             """
             self.__r.setTaggedValue(self.__source_tag, self.__wavelet.getAcceleration(t))
             self.__mypde.setValue(X=-self.__vp2*grad(u), y_dirac= self.__r)
             return self.__mypde.getSolution()

            

