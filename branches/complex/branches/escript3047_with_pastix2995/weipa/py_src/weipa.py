
########################################################
#
# Copyright (c) 2003-2010 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2010 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

def __checkAndFilterData(domain, data):
    # (Reduced)Solution is not directly supported so interpolate to different
    # function space and ensure all data is defined on same domain
    from esys.escript import Solution, ReducedSolution
    from esys.escript import ContinuousFunction, ReducedContinuousFunction
    from esys.escript.util import interpolate
    new_data={}
    for n,d in data.items():
        if not d.isEmpty():
            fs=d.getFunctionSpace()
            domain2=fs.getDomain()
            if fs == Solution(domain2):
                new_data[n]=interpolate(d, ContinuousFunction(domain2))
            elif fs == ReducedSolution(domain2):
                new_data[n]=interpolate(d, ReducedContinuousFunction(domain2))
            else:
                new_data[n]=d
            if domain==None:
                domain=domain2
            elif not domain==domain2:
                raise ValueError, "save: Data must be on same domain!"
    return domain, new_data

def saveSilo(filename, cycle=0, time=0., domain=None, **data):
    """
    Writes `Data` objects and their mesh into a file using the Silo file
    format.

    Example::

        ## Within a loop using 'cycle' as the loop counter:
        time=...
        tmp=Scalar(..)
        v=Vector(..)
        saveSilo("timestep_%02d.silo"%(cycle), cycle, time, domain,
                 temperature=tmp, velocity=v)

    ``tmp`` and ``v`` are written into "timestep_XX.silo" where ``tmp`` is
    named "temperature" and ``v`` is named "velocity". The Silo file will
    contain the cycle number and time value as metadata.

    :param filename: file name of the output file
    :type filename: ``str``
    :param cycle: cycle number for this data set
    :type cycle: ``int``
    :param time: time value for this data set
    :type time: ``float``
    :param domain: domain of the `Data` objects. If not specified, the domain
                   of the given `Data` objects is used.
    :type domain: `escript.Domain`
    :note: The data objects have to be defined on the same domain although
           they don't have to share the same `FunctionSpace`.
    :note: If escript was compiled without Silo support then this function
           will throw an exception.
    """
    from weipacpp import _saveSilo
    domain, newData = __checkAndFilterData(domain, data)
    if domain==None:
        raise ValueError,"saveSilo: no domain detected."

    _saveSilo(filename, cycle, time, domain, newData)

