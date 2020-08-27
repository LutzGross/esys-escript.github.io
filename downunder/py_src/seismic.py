
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

__all__ = ['SimpleSEGYWriter', 'Ricker', 'WaveBase', 'SonicWave', 'VTIWave',
'HTIWave', 'createAbsorptionLayerFunction', 'createAbsorptionLayerFunction',
'SonicHTIWave' , "TTIWave"]

import math
import numpy as np
import sys
import time
import esys.escript as escript
import esys.escript.unitsSI as U
import esys.escript.linearPDEs as lpde


OBSPY_AVAILABLE = False
try:
    from obspy import Trace, Stream, UTCDateTime
    try:
        # new interface
        from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
    except:
        from obspy.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader
    from obspy.core import AttribDict
    OBSPY_AVAILABLE = True
except:
    pass

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
                Sets up a Ricker wavelet wih dominant frequence `f_dom` and
                center at time `t_dom`. If `t_dom` is not given an estimate
                for suitable `t_dom` is calculated so f(0)~0.

                :note: maximum frequence is about 2 x the dominant frequence.
                """
                drop=18
                self.__f=f_dom
                self.__f_max=escript.sqrt(7)*f_dom
                self.__s=math.pi*self.__f
                if t_dom == None:
                        t_dom=escript.sqrt(drop)/self.__s
                self.__t0=t_dom

        def getCenter(self):
               """
               Return value of wavelet center
               """
               return self.__t0

        def getTimeScale(self):
                """
                Returns the time scale which is the inverse of the largest
                frequence with a significant spectral component.
                """
                return 1/self.__f_max

        def getValue(self, t):
                """
                get value of wavelet at time `t`
                """
                e2=(self.__s*(t-self.__t0))**2
                return (1-2*e2)*escript.exp(-e2)

        def getVelocity(self, t):
                """
                get the velocity f'(t) at time `t`
                """
                e2=(self.__s*(t-self.__t0))**2
                return (-3+2*e2)*escript.exp(-e2)*2*self.__s**2*(t-self.__t0)

        def getAcceleration(self, t):
                """
                get the acceleration f''(t) at time `t`
                """
                e2=(self.__s*(t-self.__t0))**2
                return 2*self.__s**2*(-4*e2**2 + 12*e2 - 3)*escript.exp(-e2)

class SimpleSEGYWriter(object):
        """
        A simple writer for 2D and 3D seismic lines, in particular for synthetic data

        Typical usage::

           from esys.escript import unitsSI as U
           sw=SimpleSEGYWriter([0.,100*U.m,200*U,m,300.], source=200*U.m, sampling_interval=4*U.msec)
           while n < 10:
               sw.addRecord([i*2., i*0.67, i**2, -i*7])
           sw.write('example.segy')

        :note: the writer uses `obspy`
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
               Adds a record to the traces. A time difference of sample_interval between two records is assumed.
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
                       for i in range(len(record)):
                               self.__trace[i].append(record[i])

        def getSamplingInterval(self):
               """
               returns the sampling interval in seconds.
               """
               return self.__sampling_interval

        def obspy_available(self):
            """
            for checking if the obspy module is available
            """
            return OBSPY_AVAILABLE

        def write(self, filename):
            """
            writes to segy file

            :param filename: file name
            :note: the function uses the `obspy` module.
            """

            if not OBSPY_AVAILABLE:
                raise RuntimeError("This feature (SimpleSEGYWriter.write())"+\
                        " depends on obspy, which is not installed, see "+\
                        "https://github.com/obspy/obspy for install guide")

            if escript.getMPISizeWorld() > 1:
                raise RuntimeError("Writing segy files with multiple ranks is"+\
                        " not yet supported.")

            stream=Stream()

            for i in range(len(self.__receiver_group)):
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
                    trace.stats.segy.trace_header.scalar_to_be_applied_to_all_coordinates = -int(self.COORDINATE_SCALE)
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
            stream.stats.textual_file_header = 'C.. '+self.__text+'\nC.. with esys.escript.downunder r%s\nC.. %s'%(escript.getVersion(),time.asctime())
            stream.stats.binary_file_header = SEGYBinaryFileHeader()

            if escript.getMPIRankWorld()<1:
                stream.write(filename, format="SEGY", data_encoding=1,byteorder=sys.byteorder)

class WaveBase(object):
      """
      Base for wave propagation using the Verlet scheme.

            ``u_tt = A(t,u), u(t=t0)=u0, u_t(t=t0)=v0``

      with a given acceleration force as function of time.

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
         self.n=0
         self.t_last=t0
         self.__dt=dt

      def getTimeStepSize(self):
           return self.__dt

      def _getAcceleration(self, t, u):
           """
           returns the acceleraton for time `t` and solution `u` at time `t`
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
                  self.n += 1

             # now we work backwards
             self.t_last=t
             return t, self.u + self.v * (t-self.t)

def createAbsorptionLayerFunction(x, absorption_zone=300*U.m,
        absorption_cut=1.e-2, top_absorption=False):
    """
    Creates a distribution which is one in the interior of the domain of `x`
    and is falling down to the value 'absorption_cut' over a margin of thickness 'absorption_zone'
    toward each boundary except the top of the domain.

    :param x: location of points in the domain
    :type x: `escript.Data`
    :param absorption_zone: thickness of the absorption zone
    :param absorption_cut: value of decay function on domain boundary
    :return: function on 'x' which is one in the iterior and decays to almost zero over a margin
             toward the boundary.
    """

    if absorption_zone is None or absorption_zone == 0:
        return 1

    dom=x.getDomain()
    bb=escript.boundingBox(dom)
    DIM=dom.getDim()
    decay=-escript.log(absorption_cut)/absorption_zone**2
    f=1
    for i in range(DIM):
        x_i=x[i]
        x_l=x_i-(bb[i][0]+absorption_zone)
        m_l=escript.whereNegative(x_l)
        f=f*( (escript.exp(-decay*(x_l*m_l)**2)-1) * m_l+1 )
        if  top_absorption or not DIM-1 == i:
            x_r=(bb[i][1]-absorption_zone)-x_i
            m_r=escript.whereNegative(x_r)
            f=f*( (escript.exp(-decay*(x_r*m_r)**2)-1) * m_r+1 )
    return f

class SonicWave(WaveBase):
        """
        Solving the sonic wave equation

        `p_tt = (v_p**2 * p_i)_i  + f(t) * delta_s`   where (p-) velocity v_p.

        f(t) is wavelet acting at a point source term at positon s
        """
        def __init__(self, domain, v_p, wavelet, source_tag, dt=None, p0=None, p0_t=None, absorption_zone=300*U.m, absorption_cut=1e-2, lumping=True):
           """
           initialize the sonic wave solver

           :param domain: domain of the problem
           :type domain: `Domain`
           :param v_p: p-velocity field
           :type v_p: `escript.Scalar`
           :param wavelet: wavelet to describe the time evolution of source term
           :type wavelet: `Wavelet`
           :param source_tag: tag of the source location
           :type source_tag: 'str' or 'int'
           :param dt: time step size. If not present a suitable time step size is calculated.
           :param p0: initial solution. If not present zero is used.
           :param p0_t: initial solution change rate. If not present zero is used.
           :param absorption_zone: thickness of absorption zone
           :param absorption_cut: boundary value of absorption decay factor
           :param lumping: if True mass matrix lumping is being used. This is accelerates the computing but introduces some diffusion.
           """
           f=createAbsorptionLayerFunction(escript.Function(domain).getX(), absorption_zone, absorption_cut)
           v_p=v_p*f

           if p0 == None:
              p0=escript.Scalar(0.,escript.Solution(domain))
           else:
              p0=escript.interpolate(p0, escript.Solution(domain ))

           if p0_t == None:
              p0_t=escript.Scalar(0.,escript.Solution(domain))
           else:
              p0_t=escript.interpolate(p0_t, escript.Solution(domain ))

           if dt == None:
                  dt=min(escript.inf((1./5.)*domain.getSize()/v_p), wavelet.getTimeScale())

           super(SonicWave, self).__init__( dt, u0=p0, v0=p0_t, t0=0.)

           self.__wavelet=wavelet
           self.__mypde=lpde.LinearSinglePDE(domain)
           if lumping: self.__mypde.getSolverOptions().setSolverMethod(lpde.SolverOptions.HRZ_LUMPING)
           self.__mypde.setSymmetryOn()
           self.__mypde.setValue(D=1./v_p**2)
           self.__source_tag=source_tag
           self.__r=escript.Scalar(0., escript.DiracDeltaFunctions(self.__mypde.getDomain()))


        def  _getAcceleration(self, t, u):
             """
             returns the acceleraton for time t and solution u at time t
             """
             self.__r.setTaggedValue(self.__source_tag, self.__wavelet.getValue(t))
             self.__mypde.setValue(X=-escript.grad(u,escript.Function(self.__mypde.getDomain())), y_dirac= self.__r)
             return self.__mypde.getSolution()


class VTIWave(WaveBase):
    """
    Solving the VTI wave equation

    :note: In case of a two dimensional domain the second spatial dimenion is depth.
    """
    def __init__(self, domain, v_p, v_s, wavelet, source_tag,
            source_vector = [0.,0.,1.], eps=0., gamma=0., delta=0., rho=1.,
            dt=None, u0=None, v0=None, absorption_zone=None,
            absorption_cut=1e-2, lumping=True, disable_fast_assemblers=False):
        """
        initialize the VTI wave solver

        :param domain: domain of the problem
        :type domain: `Domain`
        :param v_p: vertical p-velocity field
        :type v_p: `escript.Scalar`
        :param v_s: vertical s-velocity field
        :type v_s: `escript.Scalar`
        :param wavelet: wavelet to describe the time evolution of source term
        :type wavelet: `Wavelet`
        :param source_tag: tag of the source location
        :type source_tag: 'str' or 'int'
        :param source_vector: source orientation vector
        :param eps: first Thompsen parameter
        :param delta: second Thompsen parameter
        :param gamma: third Thompsen parameter
        :param rho: density
        :param dt: time step size. If not present a suitable time step size is calculated.
        :param u0: initial solution. If not present zero is used.
        :param v0: initial solution change rate. If not present zero is used.
        :param absorption_zone: thickness of absorption zone
        :param absorption_cut: boundary value of absorption decay factor
        :param lumping: if True mass matrix lumping is being used. This is accelerates the computing but introduces some diffusion.
        :param disable_fast_assemblers: if True, forces use of slower and more general PDE assemblers
        :type disable_fast_assemblers: `boolean`
        """
        DIM=domain.getDim()
        self.fastAssembler = hasattr(domain, "createAssembler") and not disable_fast_assemblers
        f=createAbsorptionLayerFunction(escript.Function(domain).getX(), absorption_zone, absorption_cut)

        f = escript.interpolate(f, escript.Function(domain))

        v_p=v_p*f
        v_s=v_s*f

        if u0 == None:
          u0=escript.Vector(0.,escript.Solution(domain))
        else:
          u0=escript.interpolate(p0, escript.Solution(domain ))

        if v0 == None:
          v0=escript.Vector(0.,escript.Solution(domain))
        else:
          v0=escript.interpolate(v0, escript.Solution(domain ))

        if dt == None:
              dt=min((1./5.)*min(escript.inf(domain.getSize()/v_p), escript.inf(domain.getSize()/v_s)), wavelet.getTimeScale())

        super(VTIWave, self).__init__( dt, u0=u0, v0=v0, t0=0.)

        self.__wavelet=wavelet

        self.c33=v_p**2 * rho
        self.c44=v_s**2 * rho
        self.c11=(1+2*eps) * self.c33
        self.c66=(1+2*gamma) * self.c44
        self.c13=escript.sqrt(2*self.c33*(self.c33-self.c44) * delta + (self.c33-self.c44)**2)-self.c44
        self.c12=self.c11-2*self.c66

        if self.fastAssembler:
            C = [("c11", self.c11),
                    ("c12", self.c12), ("c13", self.c13), ("c33", self.c33),
                    ("c44", self.c44), ("c66", self.c66)]
            if "speckley" in domain.getDescription().lower():
                C = [(n, escript.interpolate(d, escript.ReducedFunction(domain))) for n,d in C]
            self.__mypde=lpde.WavePDE(domain, C)
        else:
            self.__mypde=lpde.LinearPDESystem(domain)
            self.__mypde.setValue(X=self.__mypde.createCoefficient('X'))

        if lumping:
            self.__mypde.getSolverOptions().setSolverMethod(lpde.SolverOptions.HRZ_LUMPING)
        self.__mypde.setSymmetryOn()
        self.__mypde.setValue(D=rho*escript.kronecker(DIM))
        self.__source_tag=source_tag

        if DIM ==2 :
          source_vector= [source_vector[0],source_vector[2]]

        self.__r=escript.Vector(0, escript.DiracDeltaFunctions(self.__mypde.getDomain()))
        self.__r.setTaggedValue(self.__source_tag, source_vector)


    def setQ(self,q):
        """
        sets the PDE q value

        :param q: the value to set
        """
        self.__mypde.setValue(q=q)

    def _getAcceleration(self, t, u):
        """
        returns the acceleraton for time `t` and solution `u` at time `t`
        """
        du = escript.grad(u)
        if not self.fastAssembler:
            sigma=self.__mypde.getCoefficient('X')
            if self.__mypde.getDim() == 3:
                e11=du[0,0]
                e22=du[1,1]
                e33=du[2,2]
                sigma[0,0]=self.c11*e11+self.c12*e22+self.c13*e33
                sigma[1,1]=self.c12*e11+self.c11*e22+self.c13*e33
                sigma[2,2]=self.c13*(e11+e22)+self.c33*e33

                s=self.c44*(du[2,1]+du[1,2])
                sigma[1,2]=s
                sigma[2,1]=s

                s=self.c44*(du[2,0]+du[0,2])
                sigma[0,2]=s
                sigma[2,0]=s

                s=self.c66*(du[0,1]+du[1,0])
                sigma[0,1]=s
                sigma[1,0]=s


            else:
                e11=du[0,0]
                e22=du[1,1]
                sigma[0,0]=self.c11*e11+self.c13*e22
                sigma[1,1]=self.c13*e11+self.c33*e22
                s=self.c44*(du[1,0]+du[0,1])
                sigma[0,1]=s
                sigma[1,0]=s
            self.__mypde.setValue(X=-sigma, y_dirac= self.__r * self.__wavelet.getValue(t))

        else:
            self.__mypde.setValue(du=du, y_dirac= self.__r * self.__wavelet.getValue(t))

        return self.__mypde.getSolution()


class HTIWave(WaveBase):
    """
    Solving the HTI wave equation (along the x_0 axis)

    :note: In case of a two dimensional domain a horizontal domain is considered, i.e. the depth component is dropped.
    """

    def __init__(self, domain, v_p, v_s,   wavelet, source_tag,
            source_vector = [1.,0.,0.], eps=0., gamma=0., delta=0., rho=1.,
            dt=None, u0=None, v0=None, absorption_zone=None,
            absorption_cut=1e-2, lumping=True, disable_fast_assemblers=False):
       """
       initialize the VTI wave solver

       :param domain: domain of the problem
       :type domain: `Domain`
       :param v_p: vertical p-velocity field
       :type v_p: `escript.Scalar`
       :param v_s: vertical s-velocity field
       :type v_s: `escript.Scalar`
       :param wavelet: wavelet to describe the time evolution of source term
       :type wavelet: `Wavelet`
       :param source_tag: tag of the source location
       :type source_tag: 'str' or 'int'
       :param source_vector: source orientation vector
       :param eps: first Thompsen parameter
       :param delta: second Thompsen parameter
       :param gamma: third Thompsen parameter
       :param rho: density
       :param dt: time step size. If not present a suitable time step size is calculated.
       :param u0: initial solution. If not present zero is used.
       :param v0: initial solution change rate. If not present zero is used.
       :param absorption_zone: thickness of absorption zone
       :param absorption_cut: boundary value of absorption decay factor
       :param lumping: if True mass matrix lumping is being used. This is accelerates the computing but introduces some diffusion.
       :param disable_fast_assemblers: if True, forces use of slower and more general PDE assemblers
       """
       DIM=domain.getDim()
       self.fastAssembler = hasattr(domain, "createAssembler") and not disable_fast_assemblers
       f=createAbsorptionLayerFunction(v_p.getFunctionSpace().getX(), absorption_zone, absorption_cut)

       v_p=v_p*f
       v_s=v_s*f

       if u0 == None:
          u0=escript.Vector(0.,escript.Solution(domain))
       else:
          u0=escript.interpolate(p0, escript.Solution(domain ))

       if v0 == None:
          v0=escript.Vector(0.,escript.Solution(domain))
       else:
          v0=escript.interpolate(v0, escript.Solution(domain ))

       if dt == None:
            dt=min((1./5.)*min(escript.inf(domain.getSize()/v_p), escript.inf(domain.getSize()/v_s)), wavelet.getTimeScale())

       super(HTIWave, self).__init__( dt, u0=u0, v0=v0, t0=0.)

       self.__wavelet=wavelet

       self.c33 = v_p**2 * rho
       self.c44 = v_s**2 * rho
       self.c11 = (1+2*eps) * self.c33
       self.c66 = (1+2*gamma) * self.c44
       self.c13 = escript.sqrt(2*self.c33*(self.c33-self.c44) * delta + (self.c33-self.c44)**2)-self.c44
       self.c23 = self.c33-2*self.c66

       if self.fastAssembler:
            C = [("c11", self.c11),
                ("c23", self.c23), ("c13", self.c13), ("c33", self.c33),
                ("c44", self.c44), ("c66", self.c66)]
            if "speckley" in domain.getDescription().lower():
                C = [(n, escript.interpolate(d, escript.ReducedFunction(domain))) for n,d in C]
            self.__mypde=lpde.WavePDE(domain, C)
       else:
            self.__mypde=lpde.LinearPDESystem(domain)
            self.__mypde.setValue(X=self.__mypde.createCoefficient('X'))

       if lumping:
            self.__mypde.getSolverOptions().setSolverMethod(lpde.SolverOptions.HRZ_LUMPING)
       self.__mypde.setSymmetryOn()
       self.__mypde.setValue(D=rho*escript.kronecker(DIM))
       self.__source_tag=source_tag

       if DIM == 2:
          source_vector= [source_vector[0],source_vector[2]]

       self.__r=escript.Vector(0, escript.DiracDeltaFunctions(self.__mypde.getDomain()))
       self.__r.setTaggedValue(self.__source_tag, source_vector)


    def setQ(self,q):
        """
        sets the PDE q value

        :param q: the value to set
        """
        self.__mypde.setValue(q=q)

    def  _getAcceleration(self, t, u):
         """
         returns the acceleraton for time `t` and solution `u` at time `t`
         """
         du = escript.grad(u)
         if self.fastAssembler:
            self.__mypde.setValue(du=du, y_dirac= self.__r * self.__wavelet.getValue(t))
         else:
             sigma=self.__mypde.getCoefficient('X')

             if self.__mypde.getDim() == 3:
                e11=du[0,0]
                e22=du[1,1]
                e33=du[2,2]

                sigma[0,0]=self.c11*e11+self.c13*(e22+e33)
                sigma[1,1]=self.c13*e11+self.c33*e22+self.c23*e33
                sigma[2,2]=self.c13*e11+self.c23*e22+self.c33*e33

                s=self.c44*(du[2,1]+du[1,2])
                sigma[1,2]=s
                sigma[2,1]=s

                s=self.c66*(du[2,0]+du[0,2])
                sigma[0,2]=s
                sigma[2,0]=s

                s=self.c66*(du[0,1]+du[1,0])
                sigma[0,1]=s
                sigma[1,0]=s

             else:
                e11=du[0,0]
                e22=du[1,1]
                sigma[0,0]=self.c11*e11+self.c13*e22
                sigma[1,1]=self.c13*e11+self.c33*e22

                s=self.c66*(du[1,0]+du[0,1])
                sigma[0,1]=s
                sigma[1,0]=s
             self.__mypde.setValue(X=-sigma, y_dirac= self.__r * self.__wavelet.getValue(t))

         return self.__mypde.getSolution()

class TTIWave(WaveBase):
        """
        Solving the 2D TTI wave equation with

        `sigma_xx= c11*e_xx + c13*e_zz + c15*e_xz`
        `sigma_zz= c13*e_xx + c33*e_zz + c35*e_xz`
        `sigma_xz= c15*e_xx + c35*e_zz + c55*e_xz`

        the coefficients `c11`, `c13`, etc are calculated from the tompsen parameters `eps`, `delta` and the tilt `theta`

        :note: currently only the 2D case is supported.
        """

        def __init__(self, domain, v_p, v_s,   wavelet, source_tag,
                source_vector = [0.,1.], eps=0., delta=0., theta=0., rho=1.,
                dt=None, u0=None, v0=None, absorption_zone=300*U.m,
                absorption_cut=1e-2, lumping=True):
           """
           initialize the TTI wave solver

           :param domain: domain of the problem
           :type domain: `Domain`
           :param v_p: vertical p-velocity field
           :type v_p: `escript.Scalar`
           :param v_s: vertical s-velocity field
           :type v_s: `escript.Scalar`
           :param wavelet: wavelet to describe the time evolution of source term
           :type wavelet: `Wavelet`
           :param source_tag: tag of the source location
           :type source_tag: 'str' or 'int'
           :param source_vector: source orientation vector
           :param eps: first Thompsen parameter
           :param delta: second Thompsen parameter
           :param theta: tilting (in Rad)
           :param rho: density
           :param dt: time step size. If not present a suitable time step size is calculated.
           :param u0: initial solution. If not present zero is used.
           :param v0: initial solution change rate. If not present zero is used.
           :param absorption_zone: thickness of absorption zone
           :param absorption_cut: boundary value of absorption decay factor
           :param lumping: if True mass matrix lumping is being used. This is accelerates the computing but introduces some diffusion.
           """
           cos=escript.cos
           sin=escript.sin
           DIM=domain.getDim()
           if not DIM == 2:
                raise ValueError("Only 2D is supported.")
           f=createAbsorptionLayerFunction(escript.Function(domain).getX(), absorption_zone, absorption_cut)

           v_p=v_p*f
           v_s=v_s*f

           if u0 == None:
              u0=escript.Vector(0.,escript.Solution(domain))
           else:
              u0=escript.interpolate(p0, escript.Solution(domain ))

           if v0 == None:
              v0=escript.Vector(0.,escript.Solution(domain))
           else:
              v0=escript.interpolate(v0, escript.Solution(domain ))

           if dt == None:
                  dt=min((1./5.)*min(escript.inf(domain.getSize()/v_p), escript.inf(domain.getSize()/v_s)), wavelet.getTimeScale())

           super(TTIWave, self).__init__( dt, u0=u0, v0=v0, t0=0.)

           self.__wavelet=wavelet

           self.__mypde=lpde.LinearPDESystem(domain)
           if lumping: self.__mypde.getSolverOptions().setSolverMethod(lpde.SolverOptions.HRZ_LUMPING)
           self.__mypde.setSymmetryOn()
           self.__mypde.setValue(D=rho*escript.kronecker(DIM), X=self.__mypde.createCoefficient('X'))
           self.__source_tag=source_tag

           self.__r=escript.Vector(0, escript.DiracDeltaFunctions(self.__mypde.getDomain()))
           self.__r.setTaggedValue(self.__source_tag, source_vector)

           c0_33=v_p**2 * rho
           c0_66=v_s**2 * rho
           c0_11=(1+2*eps) * c0_33
           c0_13=escript.sqrt(2*c0_33*(c0_33-c0_66) * delta + (c0_33-c0_66)**2)-c0_66

           self.c11= c0_11*cos(theta)**4 - 2*c0_13*cos(theta)**4 + 2*c0_13*cos(theta)**2 + c0_33*sin(theta)**4 - 4*c0_66*cos(theta)**4 + 4*c0_66*cos(theta)**2
           self.c13= -c0_11*cos(theta)**4 + c0_11*cos(theta)**2 + c0_13*sin(theta)**4 + c0_13*cos(theta)**4 - c0_33*cos(theta)**4 + c0_33*cos(theta)**2 + 4*c0_66*cos(theta)**4 - 4*c0_66*cos(theta)**2
           self.c16= (-2*c0_11*cos(theta)**2 - 4*c0_13*sin(theta)**2 + 2*c0_13 + 2*c0_33*sin(theta)**2 - 8*c0_66*sin(theta)**2 + 4*c0_66)*sin(theta)*cos(theta)/2
           self.c33= c0_11*sin(theta)**4 - 2*c0_13*cos(theta)**4 + 2*c0_13*cos(theta)**2 + c0_33*cos(theta)**4 - 4*c0_66*cos(theta)**4 + 4*c0_66*cos(theta)**2
           self.c36= (2*c0_11*cos(theta)**2 - 2*c0_11 + 4*c0_13*sin(theta)**2 - 2*c0_13 + 2*c0_33*cos(theta)**2 + 8*c0_66*sin(theta)**2 - 4*c0_66)*sin(theta)*cos(theta)/2
           self.c66= -c0_11*cos(theta)**4 + c0_11*cos(theta)**2 + 2*c0_13*cos(theta)**4 - 2*c0_13*cos(theta)**2 - c0_33*cos(theta)**4 + c0_33*cos(theta)**2 + c0_66*sin(theta)**4 + 3*c0_66*cos(theta)**4 - 2*c0_66*cos(theta)**2

        def  _getAcceleration(self, t, u):
             """
             returns the acceleraton for time `t` and solution `u` at time `t`
             """
             du = escript.grad(u)
             sigma=self.__mypde.getCoefficient('X')

             e_xx=du[0,0]
             e_zz=du[1,1]
             e_xz=du[0,1]+du[1,0]


             sigma[0,0]= self.c11 * e_xx + self.c13 * e_zz + self.c16 * e_xz
             sigma[1,1]= self.c13 * e_xx + self.c33 * e_zz + self.c36 * e_xz
             sigma_xz  = self.c16 * e_xx + self.c36 * e_zz + self.c66 * e_xz

             sigma[0,1]=sigma_xz
             sigma[1,0]=sigma_xz

             self.__mypde.setValue(X=-sigma, y_dirac= self.__r * self.__wavelet.getValue(t))
             return self.__mypde.getSolution()

class SonicHTIWave(WaveBase):
        """
        Solving the HTI wave equation (along the x_0 axis) with azimuth (rotation around verticle axis)
        under the assumption of zero shear wave velocities
        The unknowns are the transversal (along x_0) and vertial stress (Q, P)

        :note: In case of a two dimensional domain the second spatial dimenion is depth.
        """
        def __init__(self, domain, v_p, wavelet, source_tag, source_vector = [1.,0.], eps=0., delta=0., azimuth=0.,
                     dt=None, p0=None, v0=None, absorption_zone=300*U.m, absorption_cut=1e-2, lumping=True):
           """
           initialize the HTI wave solver

           :param domain: domain of the problem
           :type domain: `Doamin`
           :param v_p: vertical p-velocity field
           :type v_p: `escript.Scalar`
           :param v_s: vertical s-velocity field
           :type v_s: `escript.Scalar`
           :param wavelet: wavelet to describe the time evolution of source term
           :type wavelet: `Wavelet`
           :param source_tag: tag of the source location
           :type source_tag: 'str' or 'int'
           :param source_vector: source orientation vector
           :param eps: first Thompsen parameter
           :param azimuth: azimuth (rotation around verticle axis)
           :param gamma: third Thompsen parameter
           :param rho: density
           :param dt: time step size. If not present a suitable time step size is calculated.
           :param p0: initial solution (Q(t=0), P(t=0)). If not present zero is used.
           :param v0: initial solution change rate. If not present zero is used.
           :param absorption_zone: thickness of absorption zone
           :param absorption_cut: boundary value of absorption decay factor
           :param lumping: if True mass matrix lumping is being used. This is accelerates the computing but introduces some diffusion.
           """
           DIM=domain.getDim()
           f=createAbsorptionLayerFunction(v_p.getFunctionSpace().getX(), absorption_zone, absorption_cut)

           self.v2_p=v_p**2
           self.v2_t=self.v2_p*escript.sqrt(1+2*delta)
           self.v2_n=self.v2_p*(1+2*eps)

           if p0 == None:
              p0=escript.Data(0.,(2,),escript.Solution(domain))
           else:
              p0=escript.interpolate(p0, escript.Solution(domain ))

           if v0 == None:
              v0=escript.Data(0.,(2,),escript.Solution(domain))
           else:
              v0=escript.interpolate(v0, escript.Solution(domain ))

           if dt == None:
                  dt=min(min(escript.inf(domain.getSize()/escript.sqrt(self.v2_p)), escript.inf(domain.getSize()/escript.sqrt(self.v2_t)), escript.inf(domain.getSize()/escript.sqrt(self.v2_n))) , wavelet.getTimeScale())*0.2

           super(SonicHTIWave, self).__init__( dt, u0=p0, v0=v0, t0=0.)

           self.__wavelet=wavelet

           self.__mypde=lpde.LinearPDESystem(domain)
           if lumping: self.__mypde.getSolverOptions().setSolverMethod(lpde.SolverOptions.HRZ_LUMPING)
           self.__mypde.setSymmetryOn()
           self.__mypde.setValue(D=escript.kronecker(2), X=self.__mypde.createCoefficient('X'))
           self.__source_tag=source_tag


           self.__r=escript.Vector(0, escript.DiracDeltaFunctions(self.__mypde.getDomain()))
           self.__r.setTaggedValue(self.__source_tag, source_vector)

        def  _getAcceleration(self, t, u):
            """
            returns the acceleraton for time `t` and solution `u` at time `t`
            """
            dQ = escript.grad(u[0])[0]
            dP = escript.grad(u[1])[1:]
            sigma=self.__mypde.getCoefficient('X')

            sigma[0,0] = self.v2_n*dQ
            sigma[0,1:] = self.v2_t*dP
            sigma[1,0] = self.v2_t*dQ
            sigma[1,1:] = self.v2_p*dP

            self.__mypde.setValue(X=-sigma, y_dirac= self.__r * self.__wavelet.getValue(t))
            return self.__mypde.getSolution()            
