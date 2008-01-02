from esys.escript import *

class TransportPDE(object):
     def __init__(self,domain,num_equations=1,theta=0.,dt_max=-1.,trace=True):
        self.__domain=domain
        self.__num_equations=num_equations
        self.__theta=theta
        self.__dt_max=dt_max
        self.__transport_problem=None
	self.__trace=trace
	self.__matrix_type=0
     def trace(self,text):
	     if self.__trace: print text

     def getDomain(self):
        return self.__domain
     def getTheta(self):
        return self.__theta
     def getDt_max(self):
        return self.__dt_max
     def getNumEquations(self):
        return self.__num_equations
     def reduced(self):
	     return False
     def getFunctionSpace(self):
        if self.reduced():
           return ReducedSolution(self.getDomain())
        else:
           return Solution(self.getDomain())


     def __getNewTransportProblem(self):
       """
       returns an instance of a new operator
       """
       self.trace("New Transport problem is allocated.")
       return self.getDomain().newTransportProblem( \
                               self.getTheta(),
                               self.getDt_max(),
                               self.getNumEquations(), \
                               self.getFunctionSpace(), \
                               self.__matrix_type)

     def setValue(self,M=Data(),A=Data(),B=Data(),C=Data(),D=Data(),X=Data(),Y=Data(),
	          d=Data(),y=Data(),d_contact=Data(),y_contact=Data()):
         self.__transport_problem=self.__getNewTransportProblem()
	 if self.getNumEquations() ==1 :
		self.__source=Data(0.0,(),self.getFunctionSpace()) 
         else:
	         self.__source=Data(0.0,(self.getNumEquations(),),self.getFunctionSpace())
	 self.getDomain().addPDEToTransportProblem(
	             self.__transport_problem,
		     self.__source,
		     M,A,B,C,D,X,Y,d,y,d_contact,y_contact)

     def setInitialSolution(self,u):
	     self.__transport_problem.setInitialValue(interpolate(u,self.getFunctionSpace()))
     def solve(self,dt):
	   return self.__transport_problem.solve(self.__source,dt,{})
from esys.finley import Rectangle

dom=Rectangle(10,10)
fc=TransportPDE(dom,num_equations=1,theta=0.0,dt_max=0.)
fc.setValue(M=Scalar(1.,Function(dom)),C=Scalar(1.,Function(dom))*[-1.,0.])
x=dom.getX()
u=whereNonPositive(abs(x[0]-0.3)-0.1)*whereNonPositive(abs(x[1]-0.3)-0.1)
c=0
saveVTK("u.%s.xml"%c,u=u)
fc.setInitialSolution(u)
dt=0.301
t=0.
while t<dt:
    print "time step t=",t+dt	
    u=fc.solve(dt)	
    print "range u",inf(u),sup(u)
    c+=1
    saveVTK("u.%s.xml"%c,u=u)
    t+=dt
