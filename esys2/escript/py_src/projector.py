# $Id$

"""
@brief A simple projector for rank 0 and rank 1 Data to a continuous
function space.
"""

from esys.escript import *
from esys.linearPDEs import LinearPDE
import numarray

class Projector:

  def __init__(self, domain, rank = 0, reduce = False):
    """
    @brief Create a continuous function space projector for a domain.

    @param domain Domain of the projection.
    @param rank   Rank of data to project.
    @param reduce Flag to reduce projection order.
    """
    dim = rank * domain.getDim() + (1 - rank)
    self.pde = LinearPDE(domain, numEquations = dim, numSolutions = dim)
    self.pde.setLumping(True)
    self.pde.setReducedOrderTo(reduce)
    if rank == 0:
      D = 1
    elif rank == 1:
      D = numarray.identity(domain.getDim())
    else:
      raise Exception("Projection restricted to rank 0 or rank 1 Data.")
    self.pde._setValue(D = D)
    return

  def __del__(self):
    return

  def __call__(self, value):
    """
    @brief Execute a projection to a continous function space on value.

    @param value  The value to be projected.
    """
    self.pde._setValue(Y = value)
    return self.pde.getSolution()
