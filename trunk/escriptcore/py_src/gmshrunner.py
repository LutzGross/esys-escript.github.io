# -*- coding: utf-8 -*-

##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################

__copyright__="""Copyright (c) 2003-2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

"""
interface to gmsh
"""

__all__ = ['gmshGeo2Msh']

from .util import getMPIRankWorld, getMPIWorldMax
from .escriptcpp import getEscriptParamInt

try:
    import gmshpy
    HAVE_GMSHPY=True
except ImportError:
    HAVE_GMSHPY=False

HAVE_GMSH = getEscriptParamInt("GMSH_SUPPORT")
GMSH_MPI = HAVE_GMSH and getEscriptParamInt("GMSH_MPI")

def _runGmshPy(geoFile, mshFile, numDim, order, verbosity):
    if getMPIRankWorld() == 0:
        gmshpy.Msg_SetVerbosity(verbosity)
        gmshpy.GModel_readGEO(geoFile)
        model=gmshpy.GModel_current()
        linear=False
        incomplete=False
        model.setOrderN(order, linear, incomplete)
        model.mesh(numDim)
        ret = model.writeMSH(mshFile)==0 # 0 indicates error for gmshpy
        gmshpy.GmshClearProject()
    else:
        ret = 0
    ret = getMPIWorldMax(ret)
    return ret

def _runGmshSerial(geoFile, mshFile, numDim, order, verbosity):
    if getMPIRankWorld() == 0:
        import shlex, subprocess
        cmdline = "gmsh -format msh -%s -order %s -o '%s' '%s'"%(numDim, order, mshFile, geoFile)
        args = shlex.split(cmdline)
        try:
            ret = subprocess.call(args)
        except Exception as e:
            ret = 1
    else:
        ret = 0
   
    ret=getMPIWorldMax(ret)
    return ret

def _runGmshMPI(geoFile, mshFile, numDim, order, verbosity):
    import shlex
    from .escriptcpp import runMPIProgram
    from time import sleep

    cmdline = "gmsh -format msh -%s -order %s -v %s -o '%s' '%s'"%(numDim, order, verbosity, mshFile, geoFile)
    args = shlex.split(cmdline)
    ret = runMPIProgram(args)
    # runMPIProgram returns immediately so wait 'a bit' to let gmsh finish
    # (which is still not guaranteed, see bug #33)
    sleep(10)
    return ret

def gmshGeo2Msh(geoFile, mshFile, numDim, order=1, verbosity=0):
    """
    Runs gmsh to mesh input `geoFile`.
    Returns 0 on success.
    """
    if HAVE_GMSHPY:
        return _runGmshPy(geoFile, mshFile, numDim, order, verbosity)
    elif GMSH_MPI:
        return _runGmshMPI(geoFile, mshFile, numDim, order, verbosity)
    else: # we try our luck even if gmsh was not available at build time
        return _runGmshSerial(geoFile, mshFile, numDim, order, verbosity)

