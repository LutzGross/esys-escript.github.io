# $Id$

## @file util.py

"""
@brief Utility functions for escript
"""

import numarray
#
#   escript constants:
#
FALSE=0
TRUE=1
UNKNOWN=-1
EPSILON=1.e-15
Pi=3.1415926535897931
# matrix types
CSC=0   
CSR=1
LUMPED=10
# some solver options:
NO_REORDERING=0
MINIMUM_FILL_IN=1
NESTED_DISSECTION=2
DEFAULT_METHOD=0
PCG=1
CR=2
CGS=3
BICGSTAB=4
SSOR=5
ILU0=6
ILUT=7
JACOBI=8
# supported file formats:
VRML=1
PNG=2
JPEG=3
JPG=3
PS=4
OOGL=5
BMP=6
TIFF=7
OPENINVENTOR=8
RENDERMAN=9
PNM=10
#
# wrapper for various functions: if the argument has attribute the function name
# as an argument it calls the correspong methods. Otherwise the coresponsing numarray
# function is called.
#
def L2(arg):
    """
    @brief

    @param arg
    """
    return arg.L2()

def grad(arg,where=None):
    """
    @brief

    @param arg
    @param where
    """
    if where==None:
       return arg.grad()
    else:
       return arg.grad(where)

def integrate(arg):
    """
    @brief

    @param arg
    """
    return arg.integrate()

def interpolate(arg,where):
    """
    @brief

    @param arg
    @param where
    """
    return arg.interpolate(where)

def transpose(arg):
    """
    @brief

    @param arg
    """
    if hasattr(arg,"transpose"):
       return arg.transpose()
    else:
       return numarray.transpose(arg,axis=None)

def trace(arg):
    """
    @brief

    @param arg
    """
    if hasattr(arg,"trace"):
       return arg.trace()
    else:
       return numarray.trace(arg,k=0)

def exp(arg):
    """
    @brief

    @param arg
    """
    if hasattr(arg,"exp"):
       return arg.exp()
    else:
       return numarray.exp(arg)

def sqrt(arg):
    """
    @brief

    @param arg
    """
    if hasattr(arg,"sqrt"):
       return arg.sqrt()
    else:
       return numarray.sqrt(arg)

def sin(arg):
    """
    @brief

    @param arg
    """
    if hasattr(arg,"sin"):
       return arg.sin()
    else:
       return numarray.sin(arg)

def cos(arg):
    """
    @brief

    @param arg
    """
    if hasattr(arg,"cos"):
       return arg.cos()
    else:
       return numarray.cos(arg)

def maxval(arg):
    """
    @brief

    @param arg
    """
    return arg.maxval()

def minval(arg):
    """
    @brief

    @param arg
    """
    return arg.minval()

def sup(arg):
    """
    @brief

    @param arg
    """
    return arg.sup()

def inf(arg):
    """
    @brief

    @param arg
    """
    return arg.inf()

def Lsup(arg):
    """
    @brief

    @param arg
    """
    return arg.Lsup()

def length(arg):
    """
    @brief

    @param arg
    """
    return arg.length()

def sign(arg):
    """
    @brief

    @param arg
    """
    return arg.sign()
#
# $Log$
# Revision 1.2  2004/10/27 00:23:36  jgs
# fixed minor syntax error
#
# Revision 1.1.1.1  2004/10/26 06:53:56  jgs
# initial import of project esys2
#
# Revision 1.1.2.3  2004/10/26 06:43:48  jgs
# committing Lutz's and Paul's changes to brach jgs
#
# Revision 1.1.4.1  2004/10/20 05:32:51  cochrane
# Added incomplete Doxygen comments to files, or merely put the docstrings that already exist into Doxygen form.
#
# Revision 1.1  2004/08/05 03:58:27  gross
# Bug in Assemble_NodeCoordinates fixed
#
#
