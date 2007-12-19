from esys.escript import *

class TransportPDE(object):
    def __init__(domain,num_equations=1,theta=0.,dt_max=-1.,trace=True):
        self.__domain=domain
        self.__num_equations=num_equations
        self.__theta=theta
        self.__dt_max=dt_max
        self.__transport_problem=None

     def getDomain(self):
        return self.__domain
     def getTheta(self):
        return self.__theta
     def getDt_max(self):
        return self.__dt_max
     def getNumEquations(self):
        return self.__getNumEquations
     def getFunctionSpace(self):
        if self.redueced():
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

     def setValue(M=Data(),A=Data(),B=Data(),C=Data(),D=Data(),X=Data(),Y=Data()):
         self.__transport_problem=self.__getNewTransportProblem()
         self.getDomain().

     def solve(X=Data(),Y=Data()):
