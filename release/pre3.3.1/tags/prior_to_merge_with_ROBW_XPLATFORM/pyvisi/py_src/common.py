"""
Variables common to all classes and functions

@var __author__: name of author
@var __license__: licence agreement
@var __copyright__: copyrights
@var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__copyright__="""  Copyright (c) 2006 by ACcESS MNRF
                    http://www.access.edu.au
                Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
             http://www.opensource.org/licenses/osl-3.0.php"""
__author__="Paul Cochrane"
__url__="http://www.iservo.edu.au/esys"
__version__="$Revision$"
__date__="$Date$"


_debug = False
_pyvisiVersion = __version__
_pyvisiRevision = __version__
__revision__ = _pyvisiRevision

import os.path

def overrideWarning(methodName):
    """
    Print a warning message for functions that need to be overridden but are
    called.

    @param methodName: The method name as a string.
    """
    # print a warning message if get to here
    print "\nWarning!!  If you are reading this message, then the renderer"
    print "           you have chosen hasn't overridden this method as"
    print "           they should have.  Please contact the maintainer of"
    print "           the renderer module, quoting the method name given"
    print "           below, to get this problem fixed."
    print "\nMethod: %s()\n" % methodName

    # barf
    raise NotImplementedError, "Method not implemented by renderer"

def unsupportedError(rendererName):
    """
    Print an error message when a method is called that is defined in pyvisi
    but is not supported at the renderer module level.

    @param rendererName: the name of the renderer module
    @type rendererName: string
    """
    errorString = "Sorry, but %s doesn't support this method." % rendererName
    raise NotImplementedError, errorString

def fileCheck(fname):
    """
    Check to see if the specified file exists, if not, raise an exception

    @param fname: the name of the file to check for
    @type fname: string
    """
    if os.path.exists(fname) is None:
        raise ValueError, "File %s doesn't exist" % fname

    return

def debugMsg(message):
    """
    Convenience function for debugging messages.

    This function will print out a debugging message if the debug variable
    is set.

    @param message: the message to output if the debug flag is set
    @type message: string
    """
    if _debug:
        print "\tBASE: %s" % message

# vim: expandtab shiftwidth=4:
