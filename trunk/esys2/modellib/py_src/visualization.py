# $Id$


from esys.modelframe import Model

class Visualization(Model):
       """

             generic visualization of scalar, vector and tensorial data (not implemeted yet)

             scalar: scalar data set
             vector: vector data set
             tensor: tensor data set
             stride: visulaization is done every strides time step
             filename: name of the movie file

       """
       def __init__(self,debug=False):
            Model.__init__(self,debug=debug)
            self.declareParameter(scalar=None,vector=None,tensor=None,stride=1,filename="movie.mpg")

       def doInitialization(self,t):
           self.__tn=t
           self.__n=0
           self.__scene=None

       def doStep(self,dt):
           self.__tn+=dt
           self.__n+=1
           if self.__n%self.stride:
               data=self.scalar
               if data!=None:
                   pass
               data=self.vector
               if data!=None:
                   pass
               data=self.tensor
               if data!=None:
                   pass

       def doFinalization(self):
            # make the movie into self.filename
            pass

class WriteVTK(Visualization):
       """ writes data into a VTK file for further processing. Currently data are written in several files for
           each data type. This may change in the future.

             scalar: scalar data set
             vector: vector data set
             tensor: tensor data set
             stride: file is written every stride-th time step
             filename: name of the data files. use %s for indication of data type (s,v,t) and time step id.

       """
       def __init__(self,debug=False):
            Visualization.__init__(self,debug=debug)
            self.declareParameter(scalar=None,vector=None,tensor=None,stride=1,filename="data.%s.vtk",counter=0)

       def doInitialization(self,t):
           self.__tn=t

       def doStep(self,dt):
           self.__tn+=dt
           self.counter+=1
           print self.counter,self.stride,self.counter%self.stride
           if self.counter%self.stride==0:
               n=self.counter/self.stride
               data=self.scalar
               if hasattr(data,"saveVTK"): data.saveVTK(self.filename%("s%d"%n))
               data=self.vector
               if hasattr(data,"saveVTK"): data.saveVTK(self.filename%("v%d"%n))
               data=self.tensor
               if hasattr(data,"saveVTK"): data.saveVTK(self.filename%("t%d"%n))
