##############################################################################
#
# Copyright (c) 2003-2022 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
##############################################################################

__copyright__ = """Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__ = """Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__ = "https://github.com/LutzGross/esys-escript.github.io"

"""

This example illustrates the use of the esys.escript.minimizers. The example
is minimization of the Rosenbrock function which is commonly used test problem
for optimization.
 
"""


from esys.escript.minimizer import MinimizerLBFGS, CostFunction
import numpy as np
from scipy.optimize import rosen, rosen_der, rosen_hess
import logging
import matplotlib
# For interactive use, you can comment out this line
matplotlib.use('agg')

mylogger = logging.getLogger('esys')
mylogger.setLevel(logging.DEBUG)


# This defines the cost function to be minimized
class RosenFunc(CostFunction):
    def __init__(self, N=4):
        super().__init__()
        self.N = N

    def getDualProduct(self,  m, g):
        """
        returns the dot product of gradient g with solution vector m:
        """
        return np.dot(g, m)

    def getArguments(self, m):
        """
        return precalculated values for current approximation m. Here nothing is returned:
        """
        return None,

    def getValue(self, m, *args):
        """
        returns the value of the cost function:
        """
        return rosen(m)

    def getGradient(self, m, *args):
        """
        returns the value of the cost function:
        """
        return rosen_der(m)

    def getInverseHessianApproximation(self, r, m, *args, initializeHessian=False):
        """
        returns a vectore of an approximation of the Hessian inverse at m for vector r.
        """
        # here the Hessian is updated when at the first iteration step and when
        # the BFGS is restarted.
        if initializeHessian:
            self.invH = np.linalg.inv(rosen_hess(m))
        return self.invH.dot(r)

    def getNorm(self, m):
        """
        return the solution of m
        """
        return np.linalg.norm(m, np.inf)

# ... callback function is used to collect the table of iteration count vs. value of the cost function
TABLE = []
def myCallback(iterCount, m, dm, Fm, grad_Fm, norm_m, args_m, failed):
    if not failed:
        TABLE.append( [ iterCount, Fm ] )

#
# Create an instance of the cost function
#
F = RosenFunc(N=20)

#... Create an instance of the solver
#
solve = MinimizerLBFGS(F, iterMax=300, logger=mylogger)
solve.setTolerance(m_tol=1e-7)
solve.setCallback(myCallback)
#... run the solver:
m0 = np.zeros((F.N,))
solve.run(m0)
m = solve.getResult()
print("m = ", m)
print("close to true solution? ", np.allclose(m, np.ones( (F.N, )) ) )

# plotting the convergency history:

import matplotlib.pyplot as plt
plt.plot([t[0] for t in TABLE], [t[1] for t in TABLE], '-')
plt.ylabel('F(m)')
plt.xlabel('iteration count')
plt.yscale('log')
plt.savefig("rosenhistory.png")
