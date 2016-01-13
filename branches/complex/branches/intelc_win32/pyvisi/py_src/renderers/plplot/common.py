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


from esys.pyvisi.common import _debug
_rendererName = 'PLPLOT'
_rendererVersion = __version__
_rendererRevision = __version__

def debugMsg(message):
    """
    Convenience function for debugging messages.

    This function will print out a debugging message if the debug variable
    is set.

    @param message: the message to output if the debug flag is set
    @type message: string
    """
    if _debug:
        print "\t%s: %s" % (_rendererName, message)

def unsupportedError():
    """
    Print an error message when a method is called that is defined in pyvisi
    but is not supported at the renderer module level.
    """
    errorString = "Sorry, but %s doesn't support this method." % _rendererName
    raise NotImplementedError, errorString

# vim: expandtab shiftwidth=4:
