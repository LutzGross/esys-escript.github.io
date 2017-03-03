
##############################################################################
#
# Copyright (c) 2003-2016 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

from __future__ import print_function, division

from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from numpy.random import permutation
M = N = 10

for nnz in [0, 1, 2, 5, 8, 10, 15, 20, 30, 50, 80, 100]:
    P = permutation(M * N)[:nnz]
    I = P / N
    J = P % N
    V = permutation(nnz) + 1

    A = coo_matrix( (V,(I,J)) , shape=(M,N))
    filename = '%03d_nonzeros.mtx' % (nnz,)
    mmwrite(filename, A)

